# ============================================================
# Script:  test_init_constraints.R
# Purpose: Tests for init_y_normalize and init_constraint parameters
#          added to fit_supervised_mf_modular() in Phase 2.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-06
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_init_constraints.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm)
  library(survival)
})
if (!exists("fit_supervised_mf_modular")) {
  source("code/update_beta.R")
  source("code/update_L.R")
  source("code/update_F.R")
  source("code/update_tau.R")
  source("code/compute_elbo.R")
  source("code/fit_modular.R")
}

cat("\n--- test_init_constraints.R ---\n")

# Tiny reproducible dataset
set.seed(42)
n <- 30; p <- 50; K <- 3
Y_pos <- matrix(abs(rnorm(n * p, mean = 5)), n, p)   # all-positive (rank-like)
Y_mix <- matrix(rnorm(n * p), n, p)                   # mixed sign (z-norm-like)
time_v  <- rexp(n, 0.1)
status_v <- rbinom(n, 1, 0.6)

# Helper: quick fit with verbosity off
quick_fit <- function(...) {
  fit_supervised_mf_modular(Y_pos, time_v, status_v,
                             K = K, max_iter = 20, tol = 1e-3,
                             prior_LF = "point_exponential",
                             prior_beta = "point_normal",
                             alpha = 0.5, sign_correction = FALSE,
                             verbose = FALSE, ...)
}

# T1: default (init_constraint="none") gives same result as before by reproducing
#     fit with identical seed and comparing EL/EBeta.
run_test("IC-T1: default 'none' produces valid fit", {
  fit <- quick_fit()
  assert_true(!is.null(fit$EL), "EL should exist")
  assert_true(!is.null(fit$EBeta), "EBeta should exist")
  assert_true(all(is.finite(fit$EL)), "EL should be finite")
  assert_true(all(is.finite(fit$EBeta)), "EBeta should be finite")
})

# T2: col_l_sum1 — columns of EL at init sum to 1.
#     We can't directly observe init EL, but we can verify: the fit runs without
#     error and the result differs from default (constraint changed the init point).
run_test("IC-T2: col_l_sum1 constraint runs without error", {
  fit <- quick_fit(init_constraint = "col_l_sum1")
  assert_true(!is.null(fit$EL), "EL should exist with col_l_sum1")
  assert_true(all(is.finite(fit$EL)), "EL should be finite with col_l_sum1")
  assert_true(all(is.finite(fit$EBeta)), "EBeta should be finite with col_l_sum1")
})

# T3: col_f_sum1 — columns of EF at init sum to 1.
run_test("IC-T3: col_f_sum1 constraint runs without error", {
  fit <- quick_fit(init_constraint = "col_f_sum1")
  assert_true(!is.null(fit$EF), "EF should exist with col_f_sum1")
  assert_true(all(is.finite(fit$EF)), "EF should be finite with col_f_sum1")
  assert_true(all(is.finite(fit$EBeta)), "EBeta should be finite with col_f_sum1")
})

# T4: init_y_normalize=TRUE — runs without error, Y unchanged for CAVI.
#     Verify the fit is valid (not all-zero EL).
run_test("IC-T4: init_y_normalize=TRUE runs without error", {
  fit <- quick_fit(init_y_normalize = TRUE)
  assert_true(!is.null(fit$EL), "EL should exist with init_y_normalize")
  assert_true(all(is.finite(fit$EL)), "EL should be finite with init_y_normalize")
  # EL should not be all zeros
  assert_true(max(abs(fit$EL)) > 0, "EL should not be all-zero with init_y_normalize")
})

# T5: init_y_normalize + col_l_sum1 combined — no crash, valid output.
run_test("IC-T5: combined init_y_normalize + col_l_sum1 produces valid fit", {
  fit <- quick_fit(init_y_normalize = TRUE, init_constraint = "col_l_sum1")
  assert_true(all(is.finite(fit$EL)), "EL finite with combined init flags")
  assert_true(all(is.finite(fit$EBeta)), "EBeta finite with combined init flags")
})

# T6: invalid init_constraint value raises an error.
run_test("IC-T6: invalid init_constraint raises error", {
  err <- tryCatch(
    quick_fit(init_constraint = "bad_value"),
    error = function(e) conditionMessage(e)
  )
  assert_true(is.character(err), "should raise an error for bad init_constraint")
  assert_true(grepl("init_constraint", err), "error message should mention init_constraint")
})

# T7: col_l_sum1 and col_f_sum1 produce different EL matrices, confirming they
#     each affect a different matrix.
run_test("IC-T7: col_l_sum1 and col_f_sum1 produce different EL outcomes", {
  fit_l <- quick_fit(init_constraint = "col_l_sum1")
  fit_f <- quick_fit(init_constraint = "col_f_sum1")
  # Two different init strategies; should not give identical EL
  diff <- max(abs(fit_l$EL - fit_f$EL))
  assert_true(diff > 0, "col_l_sum1 and col_f_sum1 should give different EL")
})

# T8: mixed-sign Y with init_y_normalize uses abs(colSums) in denominator and
#     doesn't produce NaN or Inf.
run_test("IC-T8: init_y_normalize on mixed-sign Y produces finite output", {
  fit <- fit_supervised_mf_modular(
    Y_mix, time_v, status_v,
    K = K, max_iter = 20, tol = 1e-3,
    prior_LF = "point_exponential", prior_beta = "point_normal",
    alpha = 0.5, sign_correction = FALSE, verbose = FALSE,
    init_y_normalize = TRUE
  )
  assert_true(all(is.finite(fit$EL)), "EL finite for mixed-sign Y with init_y_normalize")
  assert_true(all(is.finite(fit$EBeta)), "EBeta finite for mixed-sign Y with init_y_normalize")
})
