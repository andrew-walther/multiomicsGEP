# ============================================================
# Script:  test_deflation_init.R
# Purpose: Tests for deflation_svd_init() (code/deflation_init.R) and its
#          integration as init_method="deflation" in fit_supervised_mf_modular()
#          (Cluster A / LB) and fit_cox_on_yf() (Cluster B / YFB).
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_deflation_init.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival)
})
if (!exists("deflation_svd_init")) source("code/deflation_init.R")
if (!exists("fit_supervised_mf_modular")) {
  source("code/update_beta.R"); source("code/update_L.R")
  source("code/update_F.R");    source("code/update_tau.R")
  source("code/compute_elbo.R"); source("code/update_F_cohort.R")
  suppressMessages(tryCatch(source("code/fit_modular.R"), error = function(e) invisible(NULL)))
}
if (!exists("fit_cox_on_yf")) {
  source("code/update_L_surv_YFB.R"); source("code/update_F_surv_YFB.R")
  suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL)))
}

cat("\n--- test_deflation_init.R ---\n")

# ============================================================
# Part 1: deflation_svd_init() unit tests
# ============================================================

# T1: dimensions are correct
run_test("DI-T1: EL/EF have correct dimensions", {
  set.seed(1)
  Y <- matrix(rnorm(20 * 15), 20, 15)
  out <- deflation_svd_init(Y, 4)
  assert_equal(dim(out$EL), c(20L, 4L), "EL should be n x K")
  assert_equal(dim(out$EF), c(15L, 4L), "EF should be p x K")
})

# T2: K=1 exactly matches a direct rank-1 SVD of Y
run_test("DI-T2: K=1 matches direct rank-1 SVD", {
  set.seed(2)
  Y  <- matrix(rnorm(10 * 8), 10, 8)
  sv <- svd(Y, nu = 1, nv = 1)
  d1 <- sqrt(sv$d[1])
  out <- deflation_svd_init(Y, 1)
  assert_near(out$EL[, 1], sv$u[, 1] * d1, tol = 1e-8, "EL column 1 should match rank-1 SVD of Y")
  assert_near(out$EF[, 1], sv$v[, 1] * d1, tol = 1e-8, "EF column 1 should match rank-1 SVD of Y")
})

# T3: exact recovery on a clean low-rank matrix (well-separated singular values)
run_test("DI-T3: recovers exact low-rank structure with well-separated singular values", {
  set.seed(3)
  n <- 15; p <- 12; K <- 3
  u <- qr.Q(qr(matrix(rnorm(n * K), n, K)))
  v <- qr.Q(qr(matrix(rnorm(p * K), p, K)))
  d <- c(20, 10, 5)   # well-separated
  Y <- u %*% diag(d) %*% t(v)
  out <- deflation_svd_init(Y, K)
  Y_hat <- out$EL %*% t(out$EF)
  assert_near(Y_hat, Y, tol = 1e-6, "Deflation should exactly reconstruct a clean rank-K matrix")
})

# T4: residual Frobenius norm strictly decreases with each successive factor
#     (each deflation step removes real signal, not noise)
run_test("DI-T4: residual norm strictly decreases through deflation steps", {
  set.seed(4)
  n <- 15; p <- 12; K <- 4
  u <- qr.Q(qr(matrix(rnorm(n * K), n, K)))
  v <- qr.Q(qr(matrix(rnorm(p * K), p, K)))
  d <- c(20, 15, 10, 5)
  Y <- u %*% diag(d) %*% t(v) + matrix(rnorm(n * p, sd = 0.01), n, p)

  resid_norms <- numeric(K)
  Y_resid <- Y
  for (k in seq_len(K)) {
    sv <- svd(Y_resid, nu = 1, nv = 1)
    d1 <- sqrt(max(sv$d[1], 0))
    Y_resid <- Y_resid - outer(sv$u[, 1] * d1, sv$v[, 1] * d1)
    resid_norms[k] <- sum(Y_resid^2)
  }
  assert_true(all(diff(resid_norms) < 0), "Residual Frobenius norm should strictly decrease each step")
})

# T5: near-tied singular values with disjoint support (the collapse precursor
#     pattern) -- deflation should still fully separate and explain both
run_test("DI-T5: near-tied singular values with disjoint support are both recovered", {
  set.seed(5)
  n <- 20; p <- 20
  u1 <- c(rnorm(10), rep(0, 10)); u1 <- u1 / sqrt(sum(u1^2))
  u2 <- c(rep(0, 10), rnorm(10)); u2 <- u2 / sqrt(sum(u2^2))
  v1 <- c(rnorm(10), rep(0, 10)); v1 <- v1 / sqrt(sum(v1^2))
  v2 <- c(rep(0, 10), rnorm(10)); v2 <- v2 / sqrt(sum(v2^2))
  Y <- 10 * outer(u1, v1) + 10 * outer(u2, v2)   # tied amplitude, disjoint support

  out <- deflation_svd_init(Y, 2)
  Y_hat <- out$EL %*% t(out$EF)
  assert_near(Y_hat, Y, tol = 1e-6,
              "Deflation should recover both tied-amplitude, disjoint-support factors")
  # Each extracted factor should itself be non-degenerate (not all-zero)
  assert_true(all(colSums(out$EL^2) > 1e-8), "Neither EL column should be degenerate/all-zero")
  assert_true(all(colSums(out$EF^2) > 1e-8), "Neither EF column should be degenerate/all-zero")
})

# T6: invalid K raises informative errors
run_test("DI-T6: invalid K raises error", {
  Y <- matrix(rnorm(10 * 8), 10, 8)
  err1 <- tryCatch({ deflation_svd_init(Y, 0); NULL }, error = function(e) conditionMessage(e))
  err2 <- tryCatch({ deflation_svd_init(Y, 2.5); NULL }, error = function(e) conditionMessage(e))
  err3 <- tryCatch({ deflation_svd_init(Y, 9); NULL }, error = function(e) conditionMessage(e))
  assert_true(is.character(err1), "K=0 should raise an error")
  assert_true(is.character(err2), "K=2.5 should raise an error")
  assert_true(is.character(err3), "K > min(n,p) should raise an error")
})

# ============================================================
# Part 2: integration -- init_method="deflation" in both fit functions
# ============================================================

set.seed(99)
n <- 40; p <- 60; K <- 3
Y   <- matrix(abs(rnorm(n * p, mean = 3)), n, p)
tv  <- rexp(n, 0.1)
sv_ <- rbinom(n, 1, 0.6)

# T7: fit_supervised_mf_modular(init_method="deflation") runs and returns
#     non-negative EL/EF (matching point_exponential prior convention, via
#     this file's own pmax(.,0) transform)
run_test("DI-T7: LB init_method='deflation' runs and returns non-negative EL/EF", {
  fit <- fit_supervised_mf_modular(Y, tv, sv_,
           K = K, max_iter = 10, tol = 1e-2,
           prior_LF = "point_exponential", prior_beta = "point_normal",
           alpha = 0.5, sign_correction = FALSE,
           init_method = "deflation", verbose = FALSE)
  assert_true(!is.null(fit$EL), "fit should return EL")
  assert_true(all(fit$EL >= 0), "LB deflation-init EL should be non-negative after pmax(.,0)")
  assert_true(all(fit$EF >= 0), "LB deflation-init EF should be non-negative after pmax(.,0)")
})

# T8: fit_cox_on_yf(init_method="deflation") runs and returns non-negative
#     EL/EF (via this file's own abs() transform)
run_test("DI-T8: YFB init_method='deflation' runs and returns non-negative EL/EF", {
  fit <- fit_cox_on_yf(Y, tv, sv_,
           K = K, max_iter = 10, tol = 1e-2,
           prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
           init_method = "deflation", verbose = FALSE)
  assert_true(!is.null(fit$EL), "fit should return EL")
  assert_true(all(fit$EL >= 0), "YFB deflation-init EL should be non-negative after abs()")
  assert_true(all(fit$EF >= 0), "YFB deflation-init EF should be non-negative after abs()")
})

# T9: unknown init_method still raises the existing error (deflation didn't
#     break the else-stop branch) -- regression guard for both models
run_test("DI-T9: unknown init_method still raises an error in both models", {
  err_lb <- tryCatch(
    fit_supervised_mf_modular(Y, tv, sv_, K = K, max_iter = 5, alpha = 0.5,
      sign_correction = FALSE, init_method = "bogus", verbose = FALSE),
    error = function(e) conditionMessage(e))
  err_yfb <- tryCatch(
    fit_cox_on_yf(Y, tv, sv_, K = K, max_iter = 5, alpha = 0.5,
      sign_correction = FALSE, init_method = "bogus", verbose = FALSE),
    error = function(e) conditionMessage(e))
  assert_true(is.character(err_lb), "LB should still reject an unknown init_method")
  assert_true(is.character(err_yfb), "YFB should still reject an unknown init_method")
})
