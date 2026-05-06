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
  source("code/compute_elbo.R"); source("code/fit_modular.R")
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
