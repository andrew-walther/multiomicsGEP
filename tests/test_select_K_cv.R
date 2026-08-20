# ============================================================
# Script:  test_select_K_cv.R
# Purpose: Tests for select_K_cv() — K selection via cross-validated C-index
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-06
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_select_K_cv.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival)
})
if (!exists("fit_supervised_mf_modular")) {
  source("code/update_beta.R"); source("code/update_L.R")
  source("code/update_F.R");    source("code/update_tau.R")
  source("code/compute_elbo.R")
  # fit_modular.R has a DATA_MODE runner block that stops when real_Y is NULL;
  # tryCatch ensures the function definition above the block is captured.
  suppressMessages(tryCatch(source("code/fit_modular.R"), error = function(e) invisible(NULL)))
}
if (!exists("predict_supervised_mf"))
  source("code/predict.R")
if (!exists("create_stratified_folds"))
  source("code/train_test_split.R")
if (!exists("select_K_cv"))
  source("code/select_K.R")

cat("\n--- test_select_K_cv.R ---\n")

# Small synthetic dataset — fast enough to run 3-fold CV over a short K grid
.kcv_dat <- local({
  set.seed(77)
  n <- 60; p <- 40; K_true <- 3
  L <- matrix(rexp(n * K_true), n, K_true)
  F <- matrix(rexp(p * K_true), p, K_true)
  Y <- L %*% t(F) + matrix(rnorm(n * p, sd = 0.5), n, p)
  lp <- L[, 1] - 0.5 * L[, 2]
  tv  <- rexp(n, rate = exp(scale(lp)[, 1]))
  sv  <- rep(c(1L, 0L), length.out = n)
  list(Y = Y, time = tv, status = sv)
})

# ============================================================
# T1: Output structure ----
# ============================================================

run_test("KCV-T1: returns K_opt, cv_table, fold_results, selection_rule", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 3L, 5L), n_folds = 3L,
                     alpha = 0.5, max_iter = 8L)
  assert_true(!is.null(res$K_opt),          "K_opt missing")
  assert_true(!is.null(res$cv_table),        "cv_table missing")
  assert_true(!is.null(res$fold_results),    "fold_results missing")
  assert_true(!is.null(res$selection_rule),  "selection_rule missing")
})

run_test("KCV-T2: cv_table has one row per K with required columns", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 3L, 5L), n_folds = 3L,
                     alpha = 0.5, max_iter = 8L)
  assert_equal(nrow(res$cv_table), 3L, "cv_table should have one row per K")
  required <- c("K", "mean_cindex", "se_cindex", "n_folds")
  for (col in required)
    assert_true(col %in% names(res$cv_table),
                sprintf("cv_table missing column '%s'", col))
})

run_test("KCV-T3: fold_results has one row per (K, fold) combination", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 4L), n_folds = 3L,
                     alpha = 0.5, max_iter = 8L)
  # 2 K values × 3 folds = 6 rows
  assert_equal(nrow(res$fold_results), 6L,
               sprintf("Expected 6 fold rows, got %d", nrow(res$fold_results)))
})

run_test("KCV-T4: K_opt is always a member of K_grid", {
  d <- .kcv_dat
  K_grid <- c(2L, 3L, 5L)
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = K_grid, n_folds = 3L,
                     alpha = 0.5, max_iter = 8L)
  assert_true(res$K_opt %in% K_grid,
              sprintf("K_opt=%d not in K_grid", res$K_opt))
})

# ============================================================
# T2: Selection rules ----
# ============================================================

run_test("KCV-T5: use_1se=TRUE reports selection_rule='1se'", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 3L, 5L), n_folds = 3L,
                     use_1se = TRUE, alpha = 0.5, max_iter = 8L)
  assert_equal(res$selection_rule, "1se", "Expected selection_rule='1se'")
})

run_test("KCV-T6: use_1se=FALSE reports selection_rule='max'", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 3L, 5L), n_folds = 3L,
                     use_1se = FALSE, alpha = 0.5, max_iter = 8L)
  assert_equal(res$selection_rule, "max", "Expected selection_rule='max'")
})

run_test("KCV-T7: 1-SE rule selects K <= K chosen by max rule", {
  d <- .kcv_dat
  res_1se <- select_K_cv(d$Y, d$time, d$status,
                          K_grid = c(2L, 3L, 5L), n_folds = 3L,
                          use_1se = TRUE,  alpha = 0.5, max_iter = 8L, seed = 7L)
  res_max <- select_K_cv(d$Y, d$time, d$status,
                          K_grid = c(2L, 3L, 5L), n_folds = 3L,
                          use_1se = FALSE, alpha = 0.5, max_iter = 8L, seed = 7L)
  assert_true(res_1se$K_opt <= res_max$K_opt,
              sprintf("1-SE K (%d) should be <= max K (%d)",
                      res_1se$K_opt, res_max$K_opt))
})

run_test("KCV-T8: 1-SE K_opt mean C-index >= (best_mean - se_at_best)", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 3L, 5L), n_folds = 3L,
                     use_1se = TRUE, alpha = 0.5, max_iter = 8L)
  best_idx  <- which.max(res$cv_table$mean_cindex)
  se_best   <- res$cv_table$se_cindex[best_idx]
  se_best   <- ifelse(is.na(se_best), 0, se_best)
  threshold <- res$cv_table$mean_cindex[best_idx] - se_best
  chosen_c  <- res$cv_table$mean_cindex[res$cv_table$K == res$K_opt]
  assert_true(chosen_c >= threshold,
              sprintf("K_opt mean C (%.4f) must be >= threshold (%.4f)",
                      chosen_c, threshold))
})

# ============================================================
# T3: Reproducibility ----
# ============================================================

run_test("KCV-T9: same seed gives identical cv_table and fold_results", {
  d <- .kcv_dat
  r1 <- select_K_cv(d$Y, d$time, d$status,
                    K_grid = c(2L, 3L), n_folds = 3L,
                    alpha = 0.5, max_iter = 8L, seed = 99L)
  r2 <- select_K_cv(d$Y, d$time, d$status,
                    K_grid = c(2L, 3L), n_folds = 3L,
                    alpha = 0.5, max_iter = 8L, seed = 99L)
  assert_true(identical(r1$cv_table, r2$cv_table),
              "cv_table must be identical for same seed")
  assert_true(identical(r1$fold_results, r2$fold_results),
              "fold_results must be identical for same seed")
})

# ============================================================
# T4: Input validation ----
# ============================================================

run_test("KCV-T10: negative K in K_grid raises error", {
  d <- .kcv_dat
  err <- tryCatch(
    select_K_cv(d$Y, d$time, d$status,
                K_grid = c(-1L, 3L), n_folds = 3L, alpha = 0.5, max_iter = 4L),
    error = function(e) conditionMessage(e)
  )
  assert_true(is.character(err), "Negative K should raise an error")
  assert_true(grepl("positive", err), "Error should mention 'positive'")
})

run_test("KCV-T11: n_folds < 2 raises error", {
  d <- .kcv_dat
  err <- tryCatch(
    select_K_cv(d$Y, d$time, d$status,
                K_grid = c(2L, 3L), n_folds = 1L, alpha = 0.5, max_iter = 4L),
    error = function(e) conditionMessage(e)
  )
  assert_true(is.character(err), "n_folds=1 should raise an error")
  assert_true(grepl("n_folds", err), "Error should mention n_folds")
})

run_test("KCV-T12: C-index values in fold_results lie in [0, 1]", {
  d <- .kcv_dat
  res <- select_K_cv(d$Y, d$time, d$status,
                     K_grid = c(2L, 3L), n_folds = 3L,
                     alpha = 0.5, max_iter = 8L)
  c_obs <- res$fold_results$cindex[!is.na(res$fold_results$cindex)]
  assert_true(all(c_obs >= 0 & c_obs <= 1),
              "All observed C-indices should lie in [0, 1]")
})

# ============================================================
# T5: YFB model path ----
# ============================================================

# Source YFB functions if not already available
if (!exists("fit_cox_on_yf", mode = "function")) {
  tryCatch(source("code/fit_cox_on_yf.R"),    error = function(e) invisible(NULL))
  tryCatch(source("code/predict_cox_on_yf.R"), error = function(e) invisible(NULL))
}

if (exists("fit_cox_on_yf", mode = "function") &&
    exists("predict_cox_on_yf", mode = "function")) {

  run_test("KCV-T13: model='YFB' returns valid structure with model field", {
    d <- .kcv_dat
    res <- select_K_cv(d$Y, d$time, d$status,
                       K_grid = c(2L, 3L), n_folds = 3L,
                       model = "YFB", max_iter = 8L)
    assert_true(!is.null(res$K_opt),         "K_opt missing")
    assert_true(!is.null(res$cv_table),       "cv_table missing")
    assert_true(!is.null(res$fold_results),   "fold_results missing")
    assert_equal(res$model, "YFB",            "model field should echo 'YFB'")
    assert_true(res$K_opt %in% c(2L, 3L),    "K_opt must be in K_grid")
  })

  run_test("KCV-T14: model='YFB' fold_results has correct row count", {
    d <- .kcv_dat
    res <- select_K_cv(d$Y, d$time, d$status,
                       K_grid = c(2L, 3L), n_folds = 3L,
                       model = "YFB", max_iter = 8L)
    # 2 K values × 3 folds = 6 rows
    assert_equal(nrow(res$fold_results), 6L,
                 sprintf("Expected 6 fold rows, got %d", nrow(res$fold_results)))
    c_obs <- res$fold_results$cindex[!is.na(res$fold_results$cindex)]
    assert_true(all(c_obs >= 0 & c_obs <= 1),
                "All YFB C-indices should lie in [0, 1]")
  })

  run_test("KCV-T15: model='YFB' mean C-indices are >= 0.5 (sign correction disabled in folds)", {
    # Before fix: fit_cox_on_yf always applied Phase C sign correction, orienting
    # EBeta so ZF·EBeta is concordant (C>0.5 on training). select_K_cv then evaluated
    # I(-risk_scores), double-flipping sign → C<0.5 in all folds.
    # After fix: sign_correction=FALSE inside folds. Raw SVD orientation keeps EBeta
    # in the anti-concordant direction, so I(-risk_scores) yields C>=0.5.
    # Use a high-SNR dataset (sd=0.1, beta=3, all events, n=100) so YFB reliably
    # recovers the factor in few iterations.
    d_t15 <- local({
      set.seed(31)
      n <- 100L; p <- 40L
      L  <- matrix(rexp(n * 2L), n, 2L)
      F2 <- matrix(rexp(p * 2L), p, 2L)
      Y  <- L %*% t(F2) + matrix(rnorm(n * p, sd = 0.1), n, p)
      Y  <- scale(Y, center = TRUE, scale = FALSE)
      lp <- 3.0 * scale(L[, 1L])[, 1L]   # strong signal: beta=3 on factor 1
      tv <- rexp(n, rate = exp(lp))
      ct <- rexp(n, rate = 0.5)           # random censoring (~30% censored)
      list(Y = Y, time = pmin(tv, ct), status = as.integer(tv <= ct))
    })
    res <- select_K_cv(d_t15$Y, d_t15$time, d_t15$status,
                       K_grid = c(2L, 3L), n_folds = 3L,
                       model = "YFB", max_iter = 40L)
    c_obs <- res$cv_table$mean_cindex[!is.na(res$cv_table$mean_cindex)]
    assert_true(length(c_obs) > 0, "No non-NA mean C-indices returned")
    assert_true(all(c_obs >= 0.5),
                sprintf("YFB mean C-indices should be >= 0.5 after sign fix; got: %s",
                        paste(round(c_obs, 3), collapse = ", ")))
  })

} else {
  cat("  [KCV-T13, T14, T15] Skipped: fit_cox_on_yf not available\n")
}

# ---------------------------------------------------------------------------
# T16: cohort_id is subsetted per fold — no dimension mismatch
# ---------------------------------------------------------------------------
# ============================================================
# T6: classify_factors() ----
# ============================================================

run_test("KCV-T17: classify_factors() labels survival_active/genomics_only/dead correctly", {
  # 4 factors, hand-built so PVE and EBeta land unambiguously in each category:
  #   factor 1: high |EBeta|, high PVE     -> survival_active
  #   factor 2: |EBeta| below thresh, high PVE -> genomics_only
  #   factor 3: |EBeta| below thresh, PVE below thresh -> dead
  #   factor 4: high |EBeta|, PVE below thresh -> survival_active (survival wins over dead genomics)
  n <- 20; p_ <- 10
  EL <- matrix(0, n, 4)
  EF <- matrix(0, p_, 4)
  EL[, 1] <- 1; EF[, 1] <- 1   # PVE_1 large
  EL[, 2] <- 1; EF[, 2] <- 1   # PVE_2 large
  EL[, 3] <- 0.001; EF[, 3] <- 0.001  # PVE_3 ~ 0
  EL[, 4] <- 0.001; EF[, 4] <- 0.001  # PVE_4 ~ 0
  res <- list(EL = EL, EF = EF, EBeta = c(0.5, 0.0001, 0.0001, 0.5))
  Y   <- matrix(rnorm(n * p_), n, p_)

  out <- classify_factors(res, Y, beta_thresh = 0.001, pve_thresh = 0.01)

  assert_equal(nrow(out), 4L, "classify_factors should return one row per factor")
  assert_equal(out$category[1], "survival_active", "factor 1 should be survival_active")
  assert_equal(out$category[2], "genomics_only",    "factor 2 should be genomics_only")
  assert_equal(out$category[3], "dead",             "factor 3 should be dead")
  assert_equal(out$category[4], "survival_active",  "factor 4 should be survival_active (beta wins over dead PVE)")
})

run_test("KCV-T18: classify_factors() borderline threshold edge case (values exactly at thresholds are NOT active)", {
  # Comparisons in classify_factors() use strict '>', so a value exactly equal
  # to the threshold is NOT classified as active — matches auto_prune_K()'s convention.
  n <- 10; p_ <- 5
  EL <- matrix(1, n, 1)
  EF <- matrix(1, p_, 1)
  res <- list(EL = EL, EF = EF, EBeta = 0.001)  # EBeta == beta_thresh exactly
  Y   <- matrix(1, n, p_)  # total_var = n*p_ = 50; PVE = (n*p_)^2 / 50^2... compute directly

  pve_val <- compute_pve(res, Y)
  out <- classify_factors(res, Y, beta_thresh = 0.001, pve_thresh = pve_val)  # PVE == pve_thresh exactly

  assert_true(!out$surv_active[1], "EBeta exactly at beta_thresh should not be surv_active (strict >)")
  assert_true(!out$geno_active[1], "PVE exactly at pve_thresh should not be geno_active (strict >)")
  assert_equal(out$category[1], "dead", "borderline factor (both exactly at threshold) should be classified dead")
})

run_test("KCV-T19: classify_factors() rel_thresh excludes a small-but-nonzero spurious factor relative to the max", {
  # Mirrors the pattern found in the 2026-08-19 K-recovery simulation: one
  # factor with a large, real EBeta, and a second whose EBeta clears
  # beta_thresh but is only ~29% of the max (matching the real-data K=7 fit's
  # own ratio, 0.0115/0.0404) -- rel_thresh=0.65 should exclude the second.
  n <- 10; p_ <- 5
  EL <- matrix(1, n, 2)
  EF <- matrix(1, p_, 2)
  res <- list(EL = EL, EF = EF, EBeta = c(0.0404, 0.0115))
  Y   <- matrix(rnorm(n * p_), n, p_)

  out_no_rel <- classify_factors(res, Y, beta_thresh = 0.001, pve_thresh = 0.01)
  assert_equal(out_no_rel$category[2], "survival_active",
               "without rel_thresh, the smaller factor clears beta_thresh and counts as active")

  out_rel <- classify_factors(res, Y, beta_thresh = 0.001, pve_thresh = 0.01, rel_thresh = 0.65)
  assert_equal(out_rel$category[1], "survival_active",
               "the largest factor should still be survival_active under rel_thresh")
  assert_true(out_rel$category[2] != "survival_active",
              "rel_thresh=0.65 should exclude a factor at only ~28% of the max")
})

run_test("KCV-T16: cohort_id in ... is correctly row-subsetted per fold (LB)", {
  set.seed(21)
  n <- 60L; p <- 20L
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, 0.1)
  status <- rbinom(n, 1, 0.7)
  cid    <- rep(c("A", "B"), each = n / 2)   # length n — would mismatch without fix

  # If cohort_id is NOT row-subsetted, fit_supervised_mf_modular sees n_fold rows in Y
  # but n rows in cohort_id → error. The fix extracts cohort_id and passes
  # cohort_id[train_idx] to each fold fit.
  res <- tryCatch(
    select_K_cv(Y, time, status,
                K_grid  = c(2L, 3L),
                n_folds = 3L,
                model   = "LB",
                max_iter = 5L,
                cohort_id = cid),
    error = function(e) e
  )
  assert_true(!inherits(res, "error"),
              msg = sprintf("select_K_cv with cohort_id should not error; got: %s",
                            if (inherits(res, "error")) conditionMessage(res) else "ok"))
  assert_true(is.data.frame(res$cv_table),
              msg = "cv_table should be a data.frame when cohort_id is supplied")
})
