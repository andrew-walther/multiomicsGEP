# ============================================================
# Script:  test_concordance_ci.R
# Purpose: Tests for concordance_ci.R -- bootstrap C-index CI and paired-
#          difference CI (bootstrap_concordance_ci, bootstrap_concordance_diff_ci)
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-16
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_concordance_ci.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages(library(survival))
if (!exists("bootstrap_concordance_ci"))
  source("code/concordance_ci.R")

cat("\n--- test_concordance_ci.R ---\n")

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

# A dataset with a genuinely strong, separable prognostic signal (risk score
# tightly tracks the true hazard) -- CI on its C-index should clear 0.5.
.cci_strong <- local({
  set.seed(11)
  n <- 150
  risk   <- rnorm(n)
  time   <- rexp(n, rate = exp(0.9 * scale(risk)[, 1]))
  status <- rbinom(n, 1, 0.85)
  list(risk = risk, time = time, status = status)
})

# A risk score with NO relationship to survival (pure noise) -- used for the
# paired-difference test against the strong signal above (same time/status).
.cci_noise_risk <- local({
  set.seed(22)
  rnorm(length(.cci_strong$risk))
})

# ============================================================
# bootstrap_concordance_ci() ----
# ============================================================

run_test("CCI-T1: returns estimate/lower/upper/se with lower <= estimate <= upper", {
  d <- .cci_strong
  res <- bootstrap_concordance_ci(d$risk, d$time, d$status, B = 500, seed = 1)
  assert_true(!is.null(res$estimate), "estimate missing")
  assert_true(!is.null(res$lower) && !is.null(res$upper), "CI bounds missing")
  assert_true(res$lower <= res$estimate && res$estimate <= res$upper,
              "CI does not bracket the point estimate")
})

run_test("CCI-T2: a genuinely strong prognostic signal has a CI clearing chance (0.5)", {
  d <- .cci_strong
  res <- bootstrap_concordance_ci(d$risk, d$time, d$status, B = 1000, seed = 1)
  assert_true(res$lower > 0.5, sprintf("expected lower bound > 0.5, got %.4f", res$lower))
})

run_test("CCI-T3: CI width shrinks as sample size grows (more precision, not less)", {
  set.seed(33)
  n_big <- 600
  risk   <- rnorm(n_big)
  time   <- rexp(n_big, rate = exp(0.9 * scale(risk)[, 1]))
  status <- rbinom(n_big, 1, 0.85)
  d <- .cci_strong
  res_small <- bootstrap_concordance_ci(d$risk, d$time, d$status, B = 800, seed = 5)
  res_big   <- bootstrap_concordance_ci(risk, time, status, B = 800, seed = 5)
  width_small <- res_small$upper - res_small$lower
  width_big   <- res_big$upper - res_big$lower
  assert_true(width_big < width_small,
              sprintf("expected narrower CI at n=%d than n=%d (got %.4f vs %.4f)",
                      n_big, length(d$risk), width_big, width_small))
})

run_test("CCI-T4: reproducible -- identical seed gives identical CI", {
  d <- .cci_strong
  r1 <- bootstrap_concordance_ci(d$risk, d$time, d$status, B = 300, seed = 42)
  r2 <- bootstrap_concordance_ci(d$risk, d$time, d$status, B = 300, seed = 42)
  assert_equal(r1$lower, r2$lower, "lower bound not reproducible")
  assert_equal(r1$upper, r2$upper, "upper bound not reproducible")
})

run_test("CCI-T5: rejects mismatched-length inputs (fail loud, not silent recycling)", {
  d <- .cci_strong
  err <- tryCatch({
    bootstrap_concordance_ci(d$risk, d$time[-1], d$status, B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "mismatched lengths should raise an error")
})

run_test("CCI-T6: rejects NA in risk/time/status (fail loud, not silent drop)", {
  d <- .cci_strong
  risk_na <- d$risk; risk_na[1] <- NA
  err <- tryCatch({
    bootstrap_concordance_ci(risk_na, d$time, d$status, B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "NA in risk should raise an error")
})

run_test("CCI-T7: sanity check against survival::concordance()'s analytic SE (same order of magnitude)", {
  d <- .cci_strong
  res <- bootstrap_concordance_ci(d$risk, d$time, d$status, B = 1500, seed = 7)
  analytic <- concordance(Surv(d$time, d$status) ~ d$risk)
  analytic_se <- sqrt(analytic$var)
  ratio <- res$se / analytic_se
  assert_true(ratio > 0.4 && ratio < 2.5,
              sprintf("bootstrap SE (%.4f) not within a sane range of the analytic SE (%.4f), ratio=%.2f",
                      res$se, analytic_se, ratio))
})

run_test("CCI-T7b: rejects too few observations (n < 10)", {
  err <- tryCatch({
    bootstrap_concordance_ci(rnorm(5), rexp(5), rep(1L, 5), B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "n < 10 should raise an error")
})

run_test("CCI-T7c: rejects too few events (< 2) in the full sample", {
  err <- tryCatch({
    bootstrap_concordance_ci(rnorm(20), rexp(20), c(1L, rep(0L, 19)), B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "fewer than 2 events overall should raise an error")
})

run_test("CCI-T7d: a degenerate bootstrap resample (< 2 events) fails loud with an informative message, not a bare quantile() error", {
  # n=10, exactly 2 events; seed=3's first resample draws 0 of the 2 event
  # patients (verified by direct simulation), so this deterministically
  # exercises the per-replicate guard rather than the full-sample guard above.
  risk   <- rnorm(10)
  time   <- rexp(10)
  status <- c(1L, 1L, rep(0L, 8))
  msg <- tryCatch({
    bootstrap_concordance_ci(risk, time, status, B = 50, seed = 3)
    "no error"
  }, error = function(e) conditionMessage(e))
  assert_true(grepl("fewer than 2 events after resampling", msg, fixed = TRUE),
              sprintf("expected an informative degenerate-resample message, got: %s", msg))
})

# ============================================================
# bootstrap_concordance_diff_ci() ----
# ============================================================

run_test("CCI-T8: identical risk vectors give an exact zero-width CI centered at 0", {
  d <- .cci_strong
  res <- bootstrap_concordance_diff_ci(d$risk, d$risk, d$time, d$status, B = 200, seed = 1)
  assert_equal(res$estimate, 0, "identical risk vectors should give an exact zero difference")
  assert_equal(res$lower, 0, "lower bound should be exactly 0")
  assert_equal(res$upper, 0, "upper bound should be exactly 0")
  assert_false(res$significant, "a zero-width CI at 0 should not be reported as excluding 0")
})

run_test("CCI-T9: a strong signal vs. pure noise gives a diff CI excluding 0 (significant=TRUE)", {
  d <- .cci_strong
  res <- bootstrap_concordance_diff_ci(d$risk, .cci_noise_risk, d$time, d$status,
                                        B = 1500, seed = 3)
  assert_true(res$lower > 0,
              sprintf("expected the paired-diff CI to exclude 0, lower=%.4f", res$lower))
  assert_true(res$significant, "expected significant=TRUE when the CI excludes 0")
  assert_true(res$estimate > 0, "expected a positive estimated advantage for the strong signal")
})

run_test("CCI-T10: reproducible -- identical seed gives identical diff CI", {
  d <- .cci_strong
  r1 <- bootstrap_concordance_diff_ci(d$risk, .cci_noise_risk, d$time, d$status, B = 300, seed = 9)
  r2 <- bootstrap_concordance_diff_ci(d$risk, .cci_noise_risk, d$time, d$status, B = 300, seed = 9)
  assert_equal(r1$lower, r2$lower, "diff lower bound not reproducible")
  assert_equal(r1$upper, r2$upper, "diff upper bound not reproducible")
})

run_test("CCI-T11: rejects mismatched-length inputs across risk_a/risk_b/time/status", {
  d <- .cci_strong
  err <- tryCatch({
    bootstrap_concordance_diff_ci(d$risk, .cci_noise_risk[-1], d$time, d$status,
                                    B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "mismatched lengths should raise an error")
})

run_test("CCI-T12: swapping the two risk vectors exactly negates the diff CI (same seed)", {
  d <- .cci_strong
  r_ab <- bootstrap_concordance_diff_ci(d$risk, .cci_noise_risk, d$time, d$status, B = 400, seed = 15)
  r_ba <- bootstrap_concordance_diff_ci(.cci_noise_risk, d$risk, d$time, d$status, B = 400, seed = 15)
  assert_near(r_ab$estimate, -r_ba$estimate, tol = 1e-10, msg = "swap should negate the point estimate")
  assert_near(r_ab$lower, -r_ba$upper, tol = 1e-10, msg = "swap should negate and flip the CI bounds")
  assert_near(r_ab$upper, -r_ba$lower, tol = 1e-10, msg = "swap should negate and flip the CI bounds")
})

run_test("CCI-T12b: rejects too few observations (n < 10)", {
  err <- tryCatch({
    bootstrap_concordance_diff_ci(rnorm(5), rnorm(5), rexp(5), rep(1L, 5), B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "n < 10 should raise an error")
})

run_test("CCI-T12c: rejects too few events (< 2) in the full sample", {
  err <- tryCatch({
    bootstrap_concordance_diff_ci(rnorm(20), rnorm(20), rexp(20), c(1L, rep(0L, 19)),
                                    B = 50, seed = 1)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "fewer than 2 events overall should raise an error")
})

run_test("CCI-T12d: a degenerate bootstrap resample (< 2 events) fails loud with an informative message", {
  risk_a <- rnorm(10); risk_b <- rnorm(10)
  time   <- rexp(10)
  status <- c(1L, 1L, rep(0L, 8))
  msg <- tryCatch({
    bootstrap_concordance_diff_ci(risk_a, risk_b, time, status, B = 50, seed = 3)
    "no error"
  }, error = function(e) conditionMessage(e))
  assert_true(grepl("fewer than 2 events after resampling", msg, fixed = TRUE),
              sprintf("expected an informative degenerate-resample message, got: %s", msg))
})

report_results("test_concordance_ci.R")
