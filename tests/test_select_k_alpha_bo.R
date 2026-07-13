# ============================================================
# Script:  test_select_k_alpha_bo.R
# Purpose: Tests for select_k_alpha_bayesopt() (code/select_k_alpha_bo.R)
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_select_k_alpha_bo.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival); library(rBayesianOptimization)
})
if (!exists("fit_cox_on_yf")) {
  source("code/update_beta.R"); source("code/update_L_surv_YFB.R")
  source("code/update_F_surv_YFB.R"); source("code/update_tau.R")
  source("code/compute_elbo.R"); source("code/update_F_cohort.R")
  suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL)))
}
if (!exists("predict_cox_on_yf")) source("code/predict_cox_on_yf.R")
if (!exists("create_stratified_folds")) source("code/train_test_split.R")
if (!exists("select_K_cv")) source("code/select_K.R")
if (!exists("select_k_alpha_bayesopt")) source("code/select_k_alpha_bo.R")

cat("\n--- test_select_k_alpha_bo.R ---\n")

set.seed(77)
n <- 60; p <- 80
Y   <- matrix(abs(rnorm(n * p, mean = 3)), n, p)
tv  <- rexp(n, 0.1)
sv_ <- rbinom(n, 1, 0.6)

# T1: runs end to end and returns the expected structure
run_test("BO-T1: returns Best_Par with K and alpha, and a History table", {
  res <- select_k_alpha_bayesopt(Y, tv, sv_,
           K_bounds = c(2L, 4L), alpha_bounds = c(0.2, 0.8),
           n_folds = 2, init_points = 4, n_iter = 1,
           max_iter = 8, model = "YFB", verbose = FALSE)
  assert_true(all(c("K", "alpha") %in% names(res$Best_Par)), "Best_Par should contain K and alpha")
  assert_true(!is.null(res$History), "should return a History table")
  assert_true(nrow(res$History) >= 3, "History should have at least init_points + n_iter rows")
})

# T2: selected K is an integer within bounds
run_test("BO-T2: Best_Par K is an integer within K_bounds", {
  res <- select_k_alpha_bayesopt(Y, tv, sv_,
           K_bounds = c(2L, 4L), alpha_bounds = c(0.2, 0.8),
           n_folds = 2, init_points = 4, n_iter = 1,
           max_iter = 8, model = "YFB", verbose = FALSE)
  k_sel <- res$Best_Par[["K"]]
  assert_true(abs(k_sel - round(k_sel)) < 1e-8, "K should be (numerically) an integer")
  assert_true(round(k_sel) >= 2 && round(k_sel) <= 4, "K should be within K_bounds")
})

# T3: selected alpha is within bounds
run_test("BO-T3: Best_Par alpha is within alpha_bounds", {
  res <- select_k_alpha_bayesopt(Y, tv, sv_,
           K_bounds = c(2L, 4L), alpha_bounds = c(0.2, 0.8),
           n_folds = 2, init_points = 4, n_iter = 1,
           max_iter = 8, model = "YFB", verbose = FALSE)
  a_sel <- res$Best_Par[["alpha"]]
  assert_true(a_sel >= 0.2 - 1e-8 && a_sel <= 0.8 + 1e-8, "alpha should be within alpha_bounds")
})

# T4: the objective function is consistent with a direct select_K_cv() call --
#     regression guard that the BO wrapper isn't scoring something different
#     from what select_K_cv() itself would report for the same (K, alpha).
run_test("BO-T4: objective score matches a direct select_K_cv() call at the same point", {
  obj <- select_k_alpha_bo_objective(Y, tv, sv_, K = 3, alpha = 0.5,
           n_folds = 2, seed = 42, model = "YFB", max_iter = 8)
  direct <- select_K_cv(Y, tv, sv_, K_grid = 3, n_folds = 2, use_1se = FALSE,
             seed = 42, verbose = FALSE, model = "YFB", alpha = 0.5, max_iter = 8)
  assert_near(obj$Score, direct$cv_table$mean_cindex[1], tol = 1e-8,
              "BO objective score should match select_K_cv()'s mean_cindex exactly")
})

# T5: non-integer K passed to the objective is rounded before fitting (BO's
#     GP proposes continuous K values internally)
run_test("BO-T5: objective rounds a fractional K before calling select_K_cv", {
  obj <- select_k_alpha_bo_objective(Y, tv, sv_, K = 3.4, alpha = 0.5,
           n_folds = 2, seed = 42, model = "YFB", max_iter = 8)
  direct <- select_K_cv(Y, tv, sv_, K_grid = 3, n_folds = 2, use_1se = FALSE,
             seed = 42, verbose = FALSE, model = "YFB", alpha = 0.5, max_iter = 8)
  assert_near(obj$Score, direct$cv_table$mean_cindex[1], tol = 1e-8,
              "Fractional K=3.4 should round to K=3 and match its select_K_cv() score")
})

# T6: invalid bounds raise an error
run_test("BO-T6: invalid K_bounds raises an error", {
  err <- tryCatch({
    select_k_alpha_bayesopt(Y, tv, sv_, K_bounds = c(5L, 2L), alpha_bounds = c(0.2, 0.8),
      n_folds = 2, init_points = 4, n_iter = 1, max_iter = 8, model = "YFB", verbose = FALSE)
    NULL
  }, error = function(e) conditionMessage(e))
  assert_true(is.character(err), "K_bounds with lower > upper should raise an error")
})

run_test("BO-T7: invalid alpha_bounds (outside [0,1]) raises an error", {
  err <- tryCatch({
    select_k_alpha_bayesopt(Y, tv, sv_, K_bounds = c(2L, 4L), alpha_bounds = c(-0.1, 0.8),
      n_folds = 2, init_points = 4, n_iter = 1, max_iter = 8, model = "YFB", verbose = FALSE)
    NULL
  }, error = function(e) conditionMessage(e))
  assert_true(is.character(err), "alpha_bounds outside [0,1] should raise an error")
})

# ============================================================
# pick_trustworthy_bo_winner() -- guards against a degenerate alpha=0 (or
# alpha=1) "winner" that scores well by incidental unsupervised-reconstruction
# alignment with the outcome, not genuine survival signal (K_eff=0). This is
# a documented real failure mode (DECISIONS.md 2026-05-05: alpha=1.0
# degenerate K-CV selection; the "lucky PCA direction alignment" finding).
#
# These tests inject a stub fit_fn (via pick_trustworthy_bo_winner's fit_fn
# parameter) that deterministically returns EBeta=0 for alpha=0 and a
# non-zero EBeta otherwise -- exercising the selection LOGIC in isolation
# from CAVI's actual numerical behavior (which is highly sensitive to K/DGP
# and not a reliable way to construct guaranteed K_eff=0-or-not fixtures).
# Y/time/status are placeholders here since the stub ignores them.
# ============================================================
stub_fit_fn <- function(Y, time, status, K, alpha, verbose = FALSE, ...) {
  list(EBeta = if (alpha == 0) rep(0, K) else rep(0.5, K))
}

Y_stub  <- matrix(rnorm(10 * 5), 10, 5)
tv_stub <- rexp(10, 0.1)
sv_stub <- rbinom(10, 1, 0.6)

# T8: when the top-scoring candidate already has K_eff > 0, it is returned as-is
run_test("BO-T8: top candidate with K_eff > 0 is returned directly", {
  fake_history <- data.frame(K = c(3, 5, 7), alpha = c(0.5, 0.3, 0.6),
                              Value = c(0.62, 0.58, 0.55))
  res <- pick_trustworthy_bo_winner(fake_history, Y_stub, tv_stub, sv_stub,
           n_candidates = 3, beta_threshold = 0.001, fit_fn = stub_fit_fn)
  assert_equal(res$K, 3, "should pick the top-scoring candidate (K=3) since it has K_eff > 0")
  assert_true(res$k_eff > 0, "winner should have K_eff > 0")
  assert_equal(nrow(res$candidates_checked), 3L, "should have checked all 3 candidates")
})

# T9: when the top-scoring candidate is degenerate (alpha=0, forced K_eff=0),
#     the function skips it and returns the next-best candidate with K_eff > 0
run_test("BO-T9: degenerate top candidate (alpha=0) is skipped for the next trustworthy one", {
  fake_history <- data.frame(K = c(9, 5), alpha = c(0, 0.5),
                              Value = c(0.70, 0.55))  # K=9/alpha=0 scores "best"
  res <- pick_trustworthy_bo_winner(fake_history, Y_stub, tv_stub, sv_stub,
           n_candidates = 2, beta_threshold = 0.001, fit_fn = stub_fit_fn)
  assert_equal(res$K, 5, "should skip the degenerate alpha=0 winner and pick K=5 instead")
  assert_true(res$k_eff > 0, "winner should have K_eff > 0")
  assert_equal(res$candidates_checked$k_eff[1], 0L,
               "the alpha=0 candidate should be recorded with K_eff=0")
})

# T10: if NO candidate among the top n has K_eff > 0, raise an informative error
#      (a genuine finding, not silently swallowed)
run_test("BO-T10: all-degenerate candidates raise an informative error", {
  fake_history <- data.frame(K = c(4, 6), alpha = c(0, 0), Value = c(0.6, 0.55))
  err <- tryCatch({
    pick_trustworthy_bo_winner(fake_history, Y_stub, tv_stub, sv_stub,
      n_candidates = 2, beta_threshold = 0.001, fit_fn = stub_fit_fn)
    NULL
  }, error = function(e) conditionMessage(e))
  assert_true(is.character(err), "all-K_eff=0 candidates should raise an error")
  assert_true(grepl("K_eff", err), "error message should mention K_eff")
})

# T11: candidates_checked has the expected columns and is ordered by cv_value
#      descending (matching the input history's ranking)
run_test("BO-T11: candidates_checked has expected columns, ordered by cv_value desc", {
  fake_history <- data.frame(K = c(3, 5), alpha = c(0.5, 0.5), Value = c(0.5, 0.7))
  res <- pick_trustworthy_bo_winner(fake_history, Y_stub, tv_stub, sv_stub,
           n_candidates = 2, beta_threshold = 0.001, fit_fn = stub_fit_fn)
  required_cols <- c("K", "alpha", "cv_value", "k_eff")
  assert_true(all(required_cols %in% names(res$candidates_checked)),
              "candidates_checked should contain K, alpha, cv_value, k_eff")
  assert_true(all(diff(res$candidates_checked$cv_value) <= 0),
              "candidates_checked should be ordered by cv_value descending")
})

# T12: integration smoke test -- the REAL fit_cox_on_yf() default still works
#      end to end (not just the stub), on a small but genuine signal DGP
#      (same recipe as tests/test_cox_on_yf_smoke.R).
run_test("BO-T12: real fit_cox_on_yf() default recovers a genuine top candidate", {
  set.seed(101); n3 <- 80; p3 <- 150; K3 <- 5
  beta_true <- c(0.8, rep(0, K3 - 1))
  L_true <- matrix(rnorm(n3 * K3), n3, K3)
  F_true <- matrix(0, p3, K3)
  for (k in seq_len(K3)) {
    idx <- sample(seq_len(p3), round(p3 * 0.05))
    F_true[idx, k] <- rnorm(length(idx), 0, 4)
  }
  Y3 <- L_true %*% t(F_true) + matrix(rnorm(n3 * p3), n3, p3)
  Y3 <- sweep(Y3, 2, colMeans(Y3), "-")
  eta3 <- as.vector(L_true %*% beta_true)
  raw_times  <- (-log(runif(n3)) / (0.01 * exp(eta3)))^(1 / 1.5)
  cens_times <- rexp(n3, rate = 1 / 50)
  time3   <- pmin(raw_times, cens_times)
  status3 <- as.integer(raw_times <= cens_times)

  fake_history <- data.frame(K = 9, alpha = 0.5, Value = 0.6)
  res <- pick_trustworthy_bo_winner(fake_history, Y3, time3, status3,
           n_candidates = 1, beta_threshold = 1e-3,
           max_iter = 30, N_burnin = 5, prior_beta = "normal")
  assert_true(res$k_eff > 0, "the real CAVI fit should recover a genuinely non-zero beta")
})
