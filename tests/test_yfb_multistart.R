# ============================================================
# Script:  test_yfb_multistart.R
# Purpose: Tests for fit_cox_on_yf_multistart() (code/fit_modular_multistart.R)
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_yfb_multistart.R  (standalone)
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
if (!exists("fit_cox_on_yf_multistart"))
  source("code/fit_modular_multistart.R")

cat("\n--- test_yfb_multistart.R ---\n")

set.seed(99)
n <- 40; p <- 60; K <- 3
Y   <- matrix(abs(rnorm(n * p, mean = 3)), n, p)
tv  <- rexp(n, 0.1)
sv  <- rbinom(n, 1, 0.6)

# T1: n_init=1 returns a valid best fit
run_test("YFB-MS-T1: n_init=1 returns valid fit", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 1, init_seed_base = 42)
  assert_true(!is.null(ms$best), "best fit should exist")
  assert_equal(ms$best_idx, 1L, "best_idx should be 1 with n_init=1")
  assert_equal(nrow(ms$restarts), 1L, "restarts should have 1 row")
})

# T2: n_init=1 with SVD matches a direct single fit_cox_on_yf() call
run_test("YFB-MS-T2: n_init=1 SVD matches direct single fit", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 1, init_seed_base = 42)
  ref <- fit_cox_on_yf(Y, tv, sv,
           K = K, max_iter = 15, tol = 1e-2,
           prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
           init_method = "svd", verbose = FALSE)
  assert_true(max(abs(ms$best$EBeta - ref$EBeta)) < 1e-8,
              "n_init=1 multistart should match direct single fit EBeta")
})

# T3: n_init=3 returns correctly structured output
run_test("YFB-MS-T3: n_init=3 returns correctly structured output", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 3, init_seed_base = 42)
  assert_equal(nrow(ms$restarts), 3L, "restarts should have 3 rows with n_init=3")
  assert_true(!is.null(ms$best), "best fit should exist")
  assert_true(ms$best_idx %in% 1:3, "best_idx should be in 1:3")
  assert_equal(ms$restarts$init_method[1], "svd", "restart 1 should be svd")
  assert_true(all(ms$restarts$init_method[2:3] == "random"), "restarts 2-3 should be random")
})

# T4: best_idx corresponds to highest final_elbo
run_test("YFB-MS-T4: best_idx corresponds to highest ELBO restart", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 3, init_seed_base = 42)
  expected_best <- which.max(ms$restarts$final_elbo)
  assert_equal(ms$best_idx, expected_best,
               "best_idx should equal which.max(restarts$final_elbo)")
})

# T5: restarts data.frame has required columns
run_test("YFB-MS-T5: restarts data.frame has all required columns", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 2, init_seed_base = 42)
  required_cols <- c("init_id", "init_method", "seed", "final_elbo",
                     "k_eff", "beta_max", "n_iter", "converged", "train_cindex")
  for (col in required_cols) {
    assert_true(col %in% names(ms$restarts),
                sprintf("restarts should contain column '%s'", col))
  }
})

# T6: invalid n_init raises error
run_test("YFB-MS-T6: non-positive n_init raises error", {
  err <- tryCatch(
    fit_cox_on_yf_multistart(Y, tv, sv,
      K = K, max_iter = 10, tol = 1e-2, alpha = 0.5,
      sign_correction = FALSE, n_init = 0),
    error = function(e) conditionMessage(e)
  )
  assert_true(is.character(err), "n_init=0 should raise an error")
})

# T7: train_cindex values are in [0, 1]
run_test("YFB-MS-T7: train_cindex values are in [0, 1]", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 2, init_seed_base = 42)
  c_vals <- ms$restarts$train_cindex
  c_vals <- c_vals[!is.na(c_vals)]
  assert_true(all(c_vals >= 0 & c_vals <= 1),
              "train_cindex values should be in [0, 1]")
})

# T8: train_cindex uses the Cluster B predictor Y·EF·beta (not EL·beta) --
#     regression guard against accidentally reusing the LB formula.
run_test("YFB-MS-T8: train_cindex matches Y-EF-beta predictor, not EL-beta", {
  ms <- fit_cox_on_yf_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
          n_init = 1, init_seed_base = 42)
  f   <- ms$best
  ZF  <- Y %*% sweep(f$EF, 2, f$EF_norms, "/")
  eta_correct <- as.vector(ZF %*% f$EBeta)
  eta_wrong   <- as.vector(f$EL %*% f$EBeta)
  c_correct <- as.numeric(concordance(Surv(tv, sv) ~ eta_correct)$concordance)
  reported  <- ms$restarts$train_cindex[1]
  assert_near(reported, round(c_correct, 4), tol = 1e-4,
              "reported train_cindex should match the ZF-based (not EL-based) predictor")
})
