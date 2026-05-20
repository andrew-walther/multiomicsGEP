# ============================================================
# Script:  test_multistart.R
# Purpose: Tests for fit_supervised_mf_modular_multistart()
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-06
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_multistart.R  (standalone)
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
if (!exists("fit_supervised_mf_modular_multistart"))
  source("code/fit_modular_multistart.R")

cat("\n--- test_multistart.R ---\n")

set.seed(99)
n <- 40; p <- 60; K <- 3
Y   <- matrix(abs(rnorm(n * p, mean = 3)), n, p)
tv  <- rexp(n, 0.1)
sv  <- rbinom(n, 1, 0.6)

# T1: n_init=1 returns a valid best fit (same as single fit)
run_test("MS-T1: n_init=1 returns valid fit", {
  ms <- fit_supervised_mf_modular_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_LF = "point_exponential", prior_beta = "point_normal",
          alpha = 0.5, sign_correction = FALSE,
          n_init = 1, init_seed_base = 42)
  assert_true(!is.null(ms$best), "best fit should exist")
  assert_equal(ms$best_idx, 1L, "best_idx should be 1 with n_init=1")
  assert_equal(nrow(ms$restarts), 1L, "restarts should have 1 row")
})

# T2: n_init=1 with SVD produces same EBeta as direct single fit with same settings
run_test("MS-T2: n_init=1 SVD matches direct single fit", {
  ms <- fit_supervised_mf_modular_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_LF = "point_exponential", prior_beta = "point_normal",
          alpha = 0.5, sign_correction = FALSE,
          n_init = 1, init_seed_base = 42)
  ref <- fit_supervised_mf_modular(Y, tv, sv,
           K = K, max_iter = 15, tol = 1e-2,
           prior_LF = "point_exponential", prior_beta = "point_normal",
           alpha = 0.5, sign_correction = FALSE,
           init_method = "svd", verbose = FALSE)
  # EBeta should be identical (same SVD init, same CAVI path)
  assert_true(max(abs(ms$best$EBeta - ref$EBeta)) < 1e-8,
              "n_init=1 multistart should match direct single fit EBeta")
})

# T3: n_init=3 returns correct structure
run_test("MS-T3: n_init=3 returns correctly structured output", {
  ms <- fit_supervised_mf_modular_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_LF = "point_exponential", prior_beta = "point_normal",
          alpha = 0.5, sign_correction = FALSE,
          n_init = 3, init_seed_base = 42)
  assert_equal(nrow(ms$restarts), 3L, "restarts should have 3 rows with n_init=3")
  assert_true(!is.null(ms$best), "best fit should exist")
  assert_true(ms$best_idx %in% 1:3, "best_idx should be in 1:3")
  # Restart 1 is SVD, rest are random
  assert_equal(ms$restarts$init_method[1], "svd", "restart 1 should be svd")
  assert_true(all(ms$restarts$init_method[2:3] == "random"), "restarts 2-3 should be random")
})

# T4: best_idx selects the restart with highest final_elbo
run_test("MS-T4: best_idx corresponds to highest ELBO restart", {
  ms <- fit_supervised_mf_modular_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_LF = "point_exponential", prior_beta = "point_normal",
          alpha = 0.5, sign_correction = FALSE,
          n_init = 3, init_seed_base = 42)
  expected_best <- which.max(ms$restarts$final_elbo)
  assert_equal(ms$best_idx, expected_best,
               "best_idx should equal which.max(restarts$final_elbo)")
})

# T5: restarts data.frame has required columns
run_test("MS-T5: restarts data.frame has all required columns", {
  ms <- fit_supervised_mf_modular_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_LF = "point_exponential", prior_beta = "point_normal",
          alpha = 0.5, sign_correction = FALSE,
          n_init = 2, init_seed_base = 42)
  required_cols <- c("init_id", "init_method", "seed", "final_elbo",
                     "k_eff", "beta_max", "n_iter", "converged", "train_cindex")
  for (col in required_cols) {
    assert_true(col %in% names(ms$restarts),
                sprintf("restarts should contain column '%s'", col))
  }
})

# T6: invalid n_init raises error
run_test("MS-T6: non-positive n_init raises error", {
  err <- tryCatch(
    fit_supervised_mf_modular_multistart(Y, tv, sv,
      K = K, max_iter = 10, tol = 1e-2, alpha = 0.5,
      sign_correction = FALSE, n_init = 0),
    error = function(e) conditionMessage(e)
  )
  assert_true(is.character(err), "n_init=0 should raise an error")
})

# T7: train_cindex values are in [0, 1]
run_test("MS-T7: train_cindex values are in [0, 1]", {
  ms <- fit_supervised_mf_modular_multistart(Y, tv, sv,
          K = K, max_iter = 15, tol = 1e-2,
          prior_LF = "point_exponential", prior_beta = "point_normal",
          alpha = 0.5, sign_correction = FALSE,
          n_init = 2, init_seed_base = 42)
  c_vals <- ms$restarts$train_cindex
  c_vals <- c_vals[!is.na(c_vals)]
  assert_true(all(c_vals >= 0 & c_vals <= 1),
              "train_cindex values should be in [0, 1]")
})
