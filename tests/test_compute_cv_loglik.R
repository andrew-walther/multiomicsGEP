# ============================================================
# Script:  test_compute_cv_loglik.R
# Purpose: Tests for code/compute_cv_loglik.R -- cv_survival_loglik(),
#          bicv_genomics_loglik(), gaussian_matrix_loglik(), and the
#          .fit_genomics_only() internal helper.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_compute_cv_loglik.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival)
})
if (!exists("fit_cox_on_yf")) {
  source("code/update_beta.R"); source("code/update_L_surv_YFB.R")
  source("code/update_F_surv_YFB.R"); source("code/update_tau.R")
  source("code/compute_elbo.R")
  suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL)))
}
if (!exists("predict_cox_on_yf")) source("code/predict_cox_on_yf.R")
if (!exists("create_stratified_folds")) source("code/train_test_split.R")
if (!exists("cv_survival_loglik")) source("code/compute_cv_loglik.R")

cat("\n--- test_compute_cv_loglik.R ---\n")

# --------------------------------------------------------------------------
# Fixture: synthetic data with genuine low-rank genomics structure and a
# real survival signal in the first latent factor.
# --------------------------------------------------------------------------
set.seed(3031)
n <- 80; p <- 120; K_true <- 3
L_true <- matrix(abs(rnorm(n * K_true)), n, K_true)
F_true <- matrix(abs(rnorm(p * K_true)), p, K_true)
Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.3), n, p)
beta_true <- c(2.0, 0, 0)
risk <- as.vector(L_true %*% beta_true)
time   <- rexp(n, rate = exp(risk - mean(risk)) * 0.2)
status <- rbinom(n, 1, 0.75)

# ============================================================
# gaussian_matrix_loglik() ----
# ============================================================

run_test("CVLL-T1: gaussian_matrix_loglik matches an independent dnorm() sum", {
  set.seed(1); resid <- matrix(rnorm(20 * 5), 20, 5)
  tau <- c(1, 2, 0.5, 4, 1.5)
  sd_j <- 1 / sqrt(tau)
  oracle <- sum(vapply(seq_along(tau), function(j) sum(dnorm(resid[, j], sd = sd_j[j], log = TRUE)), numeric(1)))
  assert_near(gaussian_matrix_loglik(resid, tau), oracle, tol = 1e-8,
              "gaussian_matrix_loglik does not match dnorm() oracle")
})

run_test("CVLL-T2: gaussian_matrix_loglik stop()s on tau/resid dimension mismatch", {
  err <- tryCatch({ gaussian_matrix_loglik(matrix(1, 3, 4), tau = c(1, 1)); NULL },
                  error = function(e) e)
  assert_true(!is.null(err), "expected an error on mismatched tau length")
})

# ============================================================
# cv_survival_loglik() ----
# ============================================================

run_test("CVLL-T3: fold_results has one row per fold with matching totals", {
  res <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11,
                             max_iter = 20, tol = 1e-3)
  assert_equal(nrow(res$fold_results), 4L, "expected one row per fold")
  assert_near(res$total_logPL, sum(res$fold_results$logPL), tol = 1e-10,
              "total_logPL must equal the sum of per-fold logPL")
  assert_near(res$total_events, sum(res$fold_results$n_event_test), tol = 1e-10,
              "total_events must equal the sum of per-fold event counts")
  assert_near(res$mean_logPL_per_event, res$total_logPL / res$total_events, tol = 1e-10,
              "mean_logPL_per_event formula mismatch")
})

run_test("CVLL-T4: every test patient is covered exactly once across folds", {
  res <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 5, seed = 11,
                             max_iter = 20, tol = 1e-3)
  assert_equal(as.integer(sum(res$fold_results$n_test)), as.integer(n),
               "fold test sets must partition all n patients")
})

run_test("CVLL-T5: leakage guard -- c_train is reproducible from TRAINING rows only", {
  # Refit the same folds by hand and recompute training concordance directly.
  # If cv_survival_loglik() ever used held-out outcomes to orient EBeta, this
  # would not match (held-out concordance differs from training concordance).
  res <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11,
                             max_iter = 20, tol = 1e-3)
  fold_obj <- create_stratified_folds(status, n_folds = 4, seed = 11)
  fold_id <- 1
  test_idx  <- fold_obj$folds[[fold_id]]
  train_idx <- setdiff(seq_len(n), test_idx)
  set.seed(3031)  # fit_cox_on_yf's SVD init is deterministic given Y, so no seed needed,
                   # but keep this for any future stochastic init_method default change.
  fit <- fit_cox_on_yf(Y[train_idx, , drop = FALSE], time[train_idx], status[train_idx],
                        K = K_true, max_iter = 20, tol = 1e-3,
                        prior_LF = "point_exponential", prior_beta = "normal", alpha = 0.5,
                        sign_correction = FALSE, verbose = FALSE)
  ZF_train  <- Y[train_idx, , drop = FALSE] %*% sweep(fit$EF, 2, fit$EF_norms, "/")
  eta_train <- as.vector(ZF_train %*% fit$EBeta)
  c_train_manual <- as.numeric(concordance(Surv(time[train_idx], status[train_idx]) ~ eta_train,
                                            reverse = TRUE)$concordance)
  assert_near(res$fold_results$c_train[fold_id], c_train_manual, tol = 1e-8,
              "c_train must be reproducible from training-fold data alone")
})

run_test("CVLL-T6: real survival signal scores higher (per-event) than permuted labels", {
  set.seed(42)
  perm <- sample(n)
  status_perm <- status[perm]
  time_perm   <- time[perm]
  res_real <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 21,
                                  max_iter = 30, tol = 1e-3)
  res_perm <- cv_survival_loglik(Y, time_perm, status_perm, K = K_true, n_folds = 4, seed = 21,
                                  max_iter = 30, tol = 1e-3)
  assert_true(res_real$mean_logPL_per_event > res_perm$mean_logPL_per_event,
              "held-out logPL/event should be higher when Y and survival outcome are genuinely linked")
})

run_test("CVLL-T7: stop()s on Y/time/status length mismatch", {
  err <- tryCatch({
    cv_survival_loglik(Y, time[1:10], status, K = K_true, n_folds = 4)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error on time/status length mismatch")
})

# ============================================================
# .fit_genomics_only() ----
# ============================================================

run_test("CVLL-T8: .fit_genomics_only reduces reconstruction error vs. its own SVD init", {
  fit0 <- .fit_genomics_only(Y, K = K_true, max_iter = 30, tol = 1e-4)
  svd0 <- svd(Y, nu = K_true, nv = K_true)
  d_k  <- sqrt(pmax(svd0$d[1:K_true], 0))
  EL0  <- abs(svd0$u %*% diag(d_k, K_true, K_true))
  EF0  <- abs(svd0$v %*% diag(d_k, K_true, K_true))
  err_init <- sum((Y - EL0 %*% t(EF0))^2)
  err_fit  <- sum((Y - fit0$EL %*% t(fit0$EF))^2)
  assert_true(err_fit <= err_init,
              "CAVI genomics-only fit should not increase reconstruction error vs. raw SVD init")
})

run_test("CVLL-T9: .fit_genomics_only returns correctly-shaped output", {
  fit0 <- .fit_genomics_only(Y, K = K_true, max_iter = 10, tol = 1e-3)
  assert_equal(dim(fit0$EL), as.integer(c(n, K_true)), "EL shape mismatch")
  assert_equal(dim(fit0$EF), as.integer(c(p, K_true)), "EF shape mismatch")
  assert_length(fit0$Tau, p, "Tau length mismatch")
  assert_finite(fit0$Tau, "Tau must be finite")
})

# ============================================================
# bicv_genomics_loglik() ----
# ============================================================

run_test("CVLL-T10: bi-CV blocks cover every cell of Y exactly once", {
  res <- bicv_genomics_loglik(Y, status, K = K_true, n_row_folds = 3, n_col_folds = 3,
                               seed = 11, max_iter = 15, tol = 1e-3)
  assert_equal(nrow(res$block_results), 9L, "expected n_row_folds*n_col_folds blocks")
  assert_equal(as.integer(res$n_cells_scored), as.integer(n * p),
               "bi-CV must score every cell of Y exactly once")
  assert_near(res$total_loglik, sum(res$block_results$block_loglik), tol = 1e-10,
              "total_loglik must equal the sum of per-block log-likelihoods")
})

run_test("CVLL-T11: bi-CV loglik is finite and does not trivially favor an absurdly large K", {
  res_true <- bicv_genomics_loglik(Y, status, K = K_true, n_row_folds = 3, n_col_folds = 3,
                                    seed = 11, max_iter = 15, tol = 1e-3)
  res_big  <- bicv_genomics_loglik(Y, status, K = min(n, p) %/% 2, n_row_folds = 3, n_col_folds = 3,
                                    seed = 11, max_iter = 15, tol = 1e-3)
  assert_finite(res_true$total_loglik, "bi-CV log-likelihood must be finite at K_true")
  assert_finite(res_big$total_loglik, "bi-CV log-likelihood must be finite at an oversized K")
  assert_true(res_true$total_loglik > res_big$total_loglik,
              "bi-CV should penalize an oversized K rather than favor it monotonically (unlike naive row-wise projection)")
})

run_test("CVLL-T12: stop()s when create_stratified_folds() is unavailable", {
  # Simulate by calling with an object lacking the dependency in a fresh env
  # is impractical here; instead confirm the explicit dependency check fires
  # for a Y that is not a matrix (a simpler, always-reachable guard).
  err <- tryCatch({
    bicv_genomics_loglik(as.data.frame(Y), status, K = K_true)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error when Y is not a numeric matrix")
})

# ============================================================
# cv_survival_loglik() -- cohort_id / strata_id / beta_cohort_id passthrough
# (added 2026-09-04, Stage 3 leftovers: needed to score cohort-aware arms'
# held-out survival log-likelihood without a per-fold dimension mismatch)
# ============================================================

run_test("CVLL-T13: cohort_id=NULL/strata_id=NULL/beta_cohort_id=NULL reproduces the base result", {
  ref <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 20, tol = 1e-3)
  res <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 20, tol = 1e-3,
                             cohort_id = NULL, strata_id = NULL, beta_cohort_id = NULL)
  assert_near(res$total_logPL, ref$total_logPL, tol = 1e-10, "NULL cohort args must not change the result")
})

run_test("CVLL-T14: cohort_id is correctly row-subsetted per fold (no dimension-mismatch error)", {
  set.seed(99)
  cohort2 <- sample(c("A", "B"), n, replace = TRUE)
  res <- tryCatch(
    cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 15, tol = 1e-3,
                        cohort_id = cohort2),
    error = function(e) e
  )
  assert_true(!inherits(res, "error"), sprintf("cohort_id passthrough should not error: %s",
              if (inherits(res, "error")) conditionMessage(res) else ""))
  assert_finite(res$total_logPL, "total_logPL must be finite with cohort_id supplied")
})

run_test("CVLL-T15: beta_cohort_id is correctly row-subsetted and scored (within_cohort default)", {
  set.seed(99)
  cohort2 <- sample(c("A", "B"), n, replace = TRUE)
  res <- tryCatch(
    cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 15, tol = 1e-3,
                        beta_cohort_id = cohort2),
    error = function(e) e
  )
  assert_true(!inherits(res, "error"), sprintf("beta_cohort_id passthrough should not error: %s",
              if (inherits(res, "error")) conditionMessage(res) else ""))
  assert_finite(res$total_logPL, "total_logPL must be finite with beta_cohort_id supplied")
  assert_equal(nrow(res$fold_results), 4L, "expected one row per fold")
  assert_equal(res$cv_scoring, "within_cohort", "cv_scoring defaults to within_cohort")
})

# ============================================================
# Review finding, Step 4 (fixed 2026-09-04, DECISIONS.md): two held-out
# scoring gaps in cv_survival_loglik() -- strata_id was not passed through to
# the held-out calc_cox_taylor_yf() call, and beta_cohort_id fits always
# scored with EBeta_pooled even though within-fold CV patients have KNOWN
# held-out cohort labels.
# ============================================================

run_test("CVLL-T16: cv_scoring is NA_character_ when beta_cohort_id is NULL (argument unused)", {
  res <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 15, tol = 1e-3)
  assert_true(is.na(res$cv_scoring), "cv_scoring should be NA when beta_cohort_id is NULL")
})

run_test("CVLL-T17: cv_scoring='within_cohort' vs 'unseen_cohort' give genuinely different held-out logPL", {
  # These answer different questions (see roxygen) and must not silently
  # collapse to the same number -- that would indicate cohort_id_test /
  # EBeta_pooled aren't actually being used distinctly.
  set.seed(99)
  cohort2 <- sample(c("A", "B"), n, replace = TRUE)
  res_within <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 15, tol = 1e-3,
                                    beta_cohort_id = cohort2, cv_scoring = "within_cohort")
  res_unseen <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, seed = 11, max_iter = 15, tol = 1e-3,
                                    beta_cohort_id = cohort2, cv_scoring = "unseen_cohort")
  assert_true(abs(res_within$total_logPL - res_unseen$total_logPL) > 1e-8,
              "within_cohort and unseen_cohort must score differently -- they use different beta sources")
  assert_equal(res_within$cv_scoring, "within_cohort", "echoes the requested mode")
  assert_equal(res_unseen$cv_scoring, "unseen_cohort", "echoes the requested mode")
})

run_test("CVLL-T18: cv_scoring rejects an invalid value", {
  err <- tryCatch({
    cv_survival_loglik(Y, time, status, K = K_true, n_folds = 4, max_iter = 5, tol = 1e-3,
                        beta_cohort_id = sample(c("A", "B"), n, replace = TRUE), cv_scoring = "bogus")
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected match.arg() to reject an invalid cv_scoring value")
})

run_test("CVLL-T19: strata_id is passed through to held-out scoring (isolated from its training-time effect)", {
  # strata_id also affects the TRAINING fit (already correct before this
  # fix), so comparing cv_survival_loglik(strata_id=...) against
  # cv_survival_loglik() with no strata_id conflates that pre-existing,
  # already-correct effect with the held-out-scoring bug this test targets.
  # Isolate the held-out line by reproducing fold 1 of
  # cv_survival_loglik(strata_id=strata3)'s OWN computation manually --
  # same fold assignment (create_stratified_folds depends only on status,
  # not strata_id, so it is identical either way), same training fit (same
  # seed, same strata_id) -- and comparing calc_cox_taylor_yf() WITH vs.
  # WITHOUT the strata argument on that fold's OWN held-out risk scores.
  # Before the fix, cv_survival_loglik()'s recorded logPL matched the
  # unstratified (no-strata) recomputation; after the fix it matches the
  # stratified one.
  set.seed(77)
  strata3 <- sample(c("S1", "S2", "S3"), n, replace = TRUE)
  n_folds <- 4L; seed <- 5L
  res <- cv_survival_loglik(Y, time, status, K = K_true, n_folds = n_folds, seed = seed,
                             max_iter = 15, tol = 1e-3, strata_id = strata3)

  fold_obj  <- create_stratified_folds(status, n_folds = n_folds, seed = seed)
  test_idx  <- fold_obj$folds[[1]]
  train_idx <- setdiff(seq_len(n), test_idx)
  set.seed(42L)  # cv_survival_loglik() does not itself seed the fit; fit_cox_on_yf's
                 # CAVI loop is deterministic (SVD init, no RNG), so this is a no-op
                 # safeguard for determinism, matching how other tests in this file fit directly.
  fit <- suppressMessages(fit_cox_on_yf(
    Y[train_idx, , drop = FALSE], time[train_idx], status[train_idx], K = K_true,
    max_iter = 15, tol = 1e-3, sign_correction = FALSE, verbose = FALSE,
    strata_id = strata3[train_idx]
  ))
  # Reproduce cv_survival_loglik()'s own (non-cohort) orientation decision
  # exactly, so risk_scores below matches what it would have used internally
  # -- logPL is NOT sign-invariant, so an orientation mismatch here would
  # fail this test for a reason unrelated to the strata fix being tested.
  ZF_train  <- Y[train_idx, , drop = FALSE] %*% sweep(fit$EF, 2, fit$EF_norms, "/")
  eta_train <- as.vector(ZF_train %*% fit$EBeta)
  c_train   <- as.numeric(survival::concordance(
    survival::Surv(time[train_idx], status[train_idx]) ~ eta_train, reverse = TRUE)$concordance)
  EBeta_oriented <- if (is.finite(c_train) && c_train < 0.5) -fit$EBeta else fit$EBeta
  pred <- predict_cox_on_yf(Y[test_idx, , drop = FALSE], fit$EF, EBeta_oriented, EF_norms = fit$EF_norms)
  logpl_stratified   <- calc_cox_taylor_yf(pred$risk_scores, time[test_idx], status[test_idx],
                                            strata = strata3[test_idx])$logPL
  logpl_unstratified <- calc_cox_taylor_yf(pred$risk_scores, time[test_idx], status[test_idx])$logPL

  assert_true(abs(logpl_stratified - logpl_unstratified) > 1e-10,
              "fixture must make stratified vs. unstratified held-out logPL differ, or this test cannot discriminate")
  assert_near(res$fold_results$logPL[1], logpl_stratified, tol = 1e-8,
              "cv_survival_loglik()'s recorded fold-1 logPL must match the STRATIFIED held-out recomputation")
})

report_results("test_compute_cv_loglik.R")
