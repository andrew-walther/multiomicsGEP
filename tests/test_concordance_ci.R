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

# ============================================================
# Section: frozen_reverse_cindex() and the flip/flip_a/flip_b override
# (review finding, Step 2, 2026-09-04 -- see DECISIONS.md)
# ============================================================

run_test("CCI-T13: frozen_reverse_cindex matches concordance(reverse=TRUE) exactly", {
  ref <- as.numeric(concordance(Surv(.cci_strong$time, .cci_strong$status) ~ .cci_strong$risk,
                                 reverse = TRUE)$concordance)
  got <- frozen_reverse_cindex(.cci_strong$risk, .cci_strong$time, .cci_strong$status)
  assert_near(got, ref, tol = 1e-12, msg = "must be exactly concordance(reverse=TRUE)")
})

run_test("CCI-T14: frozen_reverse_cindex can report below 0.5 (not masked by max(c,1-c))", {
  # Flip the strong signal's sign: now "larger risk = LONGER survival", the
  # wrong direction for a frozen "higher = worse" convention -- a genuine
  # below-chance finding that max(c, 1-c) would have hidden.
  got <- frozen_reverse_cindex(-.cci_strong$risk, .cci_strong$time, .cci_strong$status)
  assert_true(got < 0.5, msg = sprintf("expected a below-0.5 frozen score, got %.4f", got))
})

run_test("CCI-T15: frozen_reverse_cindex rejects mismatched lengths", {
  err <- tryCatch({
    frozen_reverse_cindex(rnorm(5), rexp(6), rbinom(6, 1, 0.5))
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "mismatched lengths should raise an error")
})

run_test("CCI-T16: bootstrap_concordance_ci(flip=NULL) is bit-for-bit identical to the pre-existing default", {
  r <- .cci_strong
  a <- bootstrap_concordance_ci(r$risk, r$time, r$status, B = 200, seed = 5)
  b <- bootstrap_concordance_ci(r$risk, r$time, r$status, B = 200, seed = 5, flip = NULL)
  assert_near(a$estimate, b$estimate, tol = 1e-12, msg = "flip=NULL must not change estimate")
  assert_near(a$lower, b$lower, tol = 1e-12, msg = "flip=NULL must not change lower")
  assert_near(a$upper, b$upper, tol = 1e-12, msg = "flip=NULL must not change upper")
})

run_test("CCI-T17: bootstrap_concordance_ci(flip=TRUE/FALSE) bypasses the same-sample orientation decision", {
  r <- .cci_strong
  # The strong-signal fixture is already "higher risk = shorter survival"
  # (reverse=TRUE concordance > 0.5); forcing flip=TRUE must give an estimate
  # BELOW 0.5, unlike the flip=NULL default which would re-orient it back up.
  forced <- bootstrap_concordance_ci(r$risk, r$time, r$status, B = 200, seed = 5, flip = TRUE)
  natural <- bootstrap_concordance_ci(r$risk, r$time, r$status, B = 200, seed = 5, flip = NULL)
  assert_true(forced$estimate < 0.5, msg = sprintf("forced flip=TRUE estimate should read <0.5, got %.4f", forced$estimate))
  assert_near(forced$estimate, 1 - natural$estimate, tol = 1e-8,
              msg = "flip=TRUE estimate should be the complement of the unflipped natural estimate")
})

run_test("CCI-T18: bootstrap_concordance_ci rejects a non-logical flip", {
  err <- tryCatch({
    bootstrap_concordance_ci(.cci_strong$risk, .cci_strong$time, .cci_strong$status, B = 20, flip = "yes")
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "non-logical flip should raise an error")
})

run_test("CCI-T19: bootstrap_concordance_diff_ci(flip_a=flip_b=NULL) is bit-for-bit identical to the pre-existing default", {
  r <- .cci_strong
  a <- bootstrap_concordance_diff_ci(r$risk, .cci_noise_risk, r$time, r$status, B = 200, seed = 7)
  b <- bootstrap_concordance_diff_ci(r$risk, .cci_noise_risk, r$time, r$status, B = 200, seed = 7,
                                      flip_a = NULL, flip_b = NULL)
  assert_near(a$estimate, b$estimate, tol = 1e-12, msg = "flip_a/flip_b=NULL must not change estimate")
  assert_near(a$lower, b$lower, tol = 1e-12, msg = "flip_a/flip_b=NULL must not change lower")
  assert_near(a$upper, b$upper, tol = 1e-12, msg = "flip_a/flip_b=NULL must not change upper")
})

run_test("CCI-T20: bootstrap_concordance_diff_ci(flip_a=FALSE, flip_b=FALSE) scores both risk scores as-is (reverse=TRUE, Cox convention)", {
  r <- .cci_strong
  # Force no re-orientation on either score; the estimate must equal the
  # direct reverse=TRUE (Cox risk-score) concordance difference, matching
  # frozen_reverse_cindex()'s convention -- not the same-sample-adjusted,
  # non-reverse arithmetic the flip=NULL path uses.
  ci <- bootstrap_concordance_diff_ci(r$risk, .cci_noise_risk, r$time, r$status, B = 50, seed = 9,
                                       flip_a = FALSE, flip_b = FALSE)
  ref <- frozen_reverse_cindex(r$risk, r$time, r$status) -
         frozen_reverse_cindex(.cci_noise_risk, r$time, r$status)
  assert_near(ci$estimate, ref, tol = 1e-8, msg = "flip_a=flip_b=FALSE should score both risk vectors as-is under reverse=TRUE")
})

# ============================================================
# Section: bootstrap_concordance_diff_ci_stratified() -- paired,
# cohort-equal-weighted bootstrap (review finding, Step 4, 2026-09-04,
# DECISIONS.md -- fixes the pooled-bootstrap estimand mismatch)
# ============================================================

.cci_strata <- local({
  set.seed(31)
  make_stratum <- function(n, seed) {
    set.seed(seed)
    risk_a <- rnorm(n); risk_b <- rnorm(n)
    time   <- rexp(n, rate = exp(0.7 * scale(risk_a)[, 1]))
    status <- rbinom(n, 1, 0.8)
    list(risk_a = risk_a, risk_b = risk_b, time = time, status = status)
  }
  # Deliberately UNEQUAL sizes so a patient-weighted vs. cohort-equal-weighted
  # aggregation would give visibly different results if conflated.
  list(make_stratum(20, 1), make_stratum(20, 2), make_stratum(200, 3))
})

run_test("CCI-T21: bootstrap_concordance_diff_ci_stratified's point estimate is the UNWEIGHTED mean of per-stratum estimates", {
  s <- .cci_strata
  ci <- bootstrap_concordance_diff_ci_stratified(
    lapply(s, `[[`, "risk_a"), lapply(s, `[[`, "risk_b"),
    lapply(s, `[[`, "time"), lapply(s, `[[`, "status"), B = 100, seed = 4
  )
  per_stratum <- vapply(s, function(st) {
    frozen_reverse_cindex(st$risk_a, st$time, st$status) -
      frozen_reverse_cindex(st$risk_b, st$time, st$status)
  }, numeric(1))
  assert_near(ci$estimate, mean(per_stratum), tol = 1e-8,
              "point estimate must be the unweighted mean across strata, not the patient-weighted pooled concordance")
  assert_equal(ci$n_strata, 3L, "n_strata should echo the number of strata supplied")
})

run_test("CCI-T22: bootstrap_concordance_diff_ci_stratified differs from naively concatenating all strata (the bug this fixes)", {
  s <- .cci_strata
  ci_stratified <- bootstrap_concordance_diff_ci_stratified(
    lapply(s, `[[`, "risk_a"), lapply(s, `[[`, "risk_b"),
    lapply(s, `[[`, "time"), lapply(s, `[[`, "status"), B = 100, seed = 4
  )
  ci_pooled <- bootstrap_concordance_diff_ci(
    unlist(lapply(s, `[[`, "risk_a")), unlist(lapply(s, `[[`, "risk_b")),
    unlist(lapply(s, `[[`, "time")), unlist(lapply(s, `[[`, "status")),
    B = 100, seed = 4, flip_a = FALSE, flip_b = FALSE
  )
  assert_true(abs(ci_stratified$estimate - ci_pooled$estimate) > 1e-6,
              "cohort-equal-weighted and patient-weighted pooling should give visibly different estimates under unequal cohort sizes")
})

run_test("CCI-T23: bootstrap_concordance_diff_ci_stratified rejects mismatched list lengths", {
  s <- .cci_strata
  err <- tryCatch({
    bootstrap_concordance_diff_ci_stratified(
      lapply(s, `[[`, "risk_a"), lapply(s, `[[`, "risk_b")[1:2],
      lapply(s, `[[`, "time"), lapply(s, `[[`, "status"), B = 20
    )
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "mismatched number of strata across lists should raise an error")
})

run_test("CCI-T24: bootstrap_concordance_diff_ci_stratified rejects a too-small stratum", {
  s <- .cci_strata
  tiny <- list(risk_a = rnorm(5), risk_b = rnorm(5), time = rexp(5), status = rep(1L, 5))
  err <- tryCatch({
    bootstrap_concordance_diff_ci_stratified(
      c(lapply(s, `[[`, "risk_a"), list(tiny$risk_a)),
      c(lapply(s, `[[`, "risk_b"), list(tiny$risk_b)),
      c(lapply(s, `[[`, "time"), list(tiny$time)),
      c(lapply(s, `[[`, "status"), list(tiny$status)),
      B = 20
    )
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", "a stratum with n < 10 should raise an error")
})

run_test("CCI-T25: bootstrap_concordance_diff_ci_stratified is reproducible given the same seed", {
  s <- .cci_strata
  a <- bootstrap_concordance_diff_ci_stratified(
    lapply(s, `[[`, "risk_a"), lapply(s, `[[`, "risk_b"),
    lapply(s, `[[`, "time"), lapply(s, `[[`, "status"), B = 100, seed = 8
  )
  b <- bootstrap_concordance_diff_ci_stratified(
    lapply(s, `[[`, "risk_a"), lapply(s, `[[`, "risk_b"),
    lapply(s, `[[`, "time"), lapply(s, `[[`, "status"), B = 100, seed = 8
  )
  assert_near(a$lower, b$lower, tol = 1e-12, "identical seed should give identical lower bound")
  assert_near(a$upper, b$upper, tol = 1e-12, "identical seed should give identical upper bound")
})

report_results("test_concordance_ci.R")
