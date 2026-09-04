# ============================================================
# Script:  test_compute_bic.R
# Purpose: Tests for compute_joint_ll_bic() (code/compute_bic.R) -- joint
#          (genomics + survival) log-likelihood and BIC for a
#          fit_cox_on_yf() (YFB) model.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-27
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_compute_bic.R  (standalone)
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
if (!exists("compute_joint_ll_bic")) source("code/compute_bic.R")

cat("\n--- test_compute_bic.R ---\n")

# --------------------------------------------------------------------------
# Fixture: small synthetic fit (n=60, p=80, K=3, max_iter=10)
# --------------------------------------------------------------------------
set.seed(2027)
n <- 60; p <- 80; K <- 3
Y      <- matrix(abs(rnorm(n * p, mean = 3)), n, p)
time   <- rexp(n, 0.1)
status <- rbinom(n, 1, 0.6)

fit <- fit_cox_on_yf(Y, time, status,
         K = K, max_iter = 10, tol = 1e-2,
         prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
         init_method = "svd", verbose = FALSE)

# ============================================================
# compute_joint_ll_bic() ----
# ============================================================

run_test("BIC-T1: BIC identity bic == -2*loglik_joint + log(n)*df", {
  res <- compute_joint_ll_bic(fit, Y, time, status)
  assert_near(res$bic, -2 * res$loglik_joint + log(res$n) * res$df, tol = 1e-8,
              "BIC identity violated")
})

run_test("BIC-T2: additivity loglik_joint == genomics + survival", {
  res <- compute_joint_ll_bic(fit, Y, time, status)
  assert_near(res$loglik_joint, res$loglik_genomics + res$loglik_survival, tol = 1e-8,
              "loglik_joint is not additive")
})

run_test("BIC-T3: df is exactly K_init*(n+p+1)", {
  res <- compute_joint_ll_bic(fit, Y, time, status)
  assert_equal(res$df, K * (n + p + 1), "df formula mismatch")
})

run_test("BIC-T4: df(K=4) - df(K=3) == n+p+1 (df is linear in K_init)", {
  fit4 <- fit_cox_on_yf(Y, time, status,
            K = K + 1, max_iter = 10, tol = 1e-2,
            prior_beta = "normal", alpha = 0.5, sign_correction = FALSE,
            init_method = "svd", verbose = FALSE)
  res3 <- compute_joint_ll_bic(fit, Y, time, status)
  res4 <- compute_joint_ll_bic(fit4, Y, time, status)
  assert_equal(res4$df - res3$df, n + p + 1, "df step size mismatch")
})

run_test("BIC-T5: df is unchanged when EBeta is mostly pruned (K_init-based, not K_eff-based)", {
  fit_pruned <- fit
  fit_pruned$EBeta[-1] <- 0
  fit_pruned$EBeta2[-1] <- 0
  res_full   <- compute_joint_ll_bic(fit, Y, time, status)
  res_pruned <- compute_joint_ll_bic(fit_pruned, Y, time, status)
  assert_equal(res_pruned$df, res_full$df,
               "df must not shrink when EBeta entries are pruned to zero")
})

#' Reorient eta exactly as compute_joint_ll_bic() now does internally
#' (DECISIONS.md 2026-09-04): a correct-direction (reverse=TRUE) concordance
#' check, flipping eta if below 0.5. Test oracles must apply the same
#' reorientation to compare against compute_joint_ll_bic()'s output.
.oriented_eta <- function(eta, time, status) {
  c_check <- as.numeric(concordance(Surv(time, status) ~ eta, reverse = TRUE)$concordance)
  if (is.finite(c_check) && c_check < 0.5) eta <- -eta
  eta
}

run_test("BIC-T6: survival term matches an independent coxph() oracle", {
  res <- compute_joint_ll_bic(fit, Y, time, status)
  ZF  <- Y %*% sweep(fit$EF, 2, fit$EF_norms, "/")
  eta <- .oriented_eta(as.vector(ZF %*% fit$EBeta), time, status)
  oracle <- coxph(Surv(time, status) ~ offset(eta), ties = "breslow")
  assert_near(res$loglik_survival, oracle$loglik[1], tol = 1e-6,
              "survival term does not match coxph() oracle")
})

run_test("BIC-T7: survival term is NOT the stale in-loop logPL", {
  res <- compute_joint_ll_bic(fit, Y, time, status)
  # A fit with sign_correction=TRUE (default) exercises the post-loop EBeta
  # flip path; the recomputed term must reflect the RETURNED EBeta, not
  # whatever logPL happened to be tracked mid-loop.
  fit_sc <- fit_cox_on_yf(Y, time, status,
              K = K, max_iter = 10, tol = 1e-2,
              prior_beta = "normal", alpha = 0.5, sign_correction = TRUE,
              init_method = "svd", verbose = FALSE)
  ZF  <- Y %*% sweep(fit_sc$EF, 2, fit_sc$EF_norms, "/")
  eta <- .oriented_eta(as.vector(ZF %*% fit_sc$EBeta), time, status)
  oracle <- coxph(Surv(time, status) ~ offset(eta), ties = "breslow")
  res_sc <- compute_joint_ll_bic(fit_sc, Y, time, status)
  assert_near(res_sc$loglik_survival, oracle$loglik[1], tol = 1e-6,
              "survival term must match the RETURNED (post sign-correction) EBeta")
})

run_test("BIC-T16: survival term is invariant to fit_cox_on_yf()'s own (unreliable) EBeta sign", {
  # Regression test for the orientation bug documented in DECISIONS.md
  # 2026-09-04: fit_cox_on_yf()'s Phase C sign-correction omits reverse=TRUE
  # and can return an anti-hazard-oriented EBeta while believing it is
  # correctly oriented. compute_joint_ll_bic() must not simply trust
  # fit$EBeta's sign -- negating it should not change the reported
  # loglik_survival, because the function re-derives the correct
  # orientation itself regardless of the input sign.
  res_asis    <- compute_joint_ll_bic(fit, Y, time, status)
  fit_flipped <- fit
  fit_flipped$EBeta  <- -fit$EBeta
  fit_flipped$EBeta2 <- fit$EBeta2  # unchanged: EBeta^2 is sign-invariant
  res_flipped <- compute_joint_ll_bic(fit_flipped, Y, time, status)
  assert_near(res_asis$loglik_survival, res_flipped$loglik_survival, tol = 1e-8,
              "loglik_survival must be invariant to the input EBeta's sign")
})

run_test("BIC-T8: genomics term equals elbo_proxy[n_iter] minus the Gaussian constant", {
  res <- compute_joint_ll_bic(fit, Y, time, status, include_gaussian_const = TRUE)
  raw <- fit$history$elbo_proxy[fit$history$n_iter]
  const <- -(n * p / 2) * log(2 * pi)
  assert_near(res$loglik_genomics, raw + const, tol = 1e-8,
              "Gaussian constant not applied correctly")
})

run_test("BIC-T9: include_gaussian_const=FALSE equals the raw elbo_proxy element exactly", {
  res <- compute_joint_ll_bic(fit, Y, time, status, include_gaussian_const = FALSE)
  raw <- fit$history$elbo_proxy[fit$history$n_iter]
  assert_equal(res$loglik_genomics, raw,
               "raw genomics term should equal elbo_proxy[n_iter] with tol=0")
})

run_test("BIC-T10: row-permutation invariance (Y, time, status permuted together)", {
  res_orig <- compute_joint_ll_bic(fit, Y, time, status)
  perm <- sample(n)
  res_perm <- compute_joint_ll_bic(fit, Y[perm, , drop = FALSE], time[perm], status[perm])
  assert_near(res_perm$loglik_survival, res_orig$loglik_survival, tol = 1e-6,
              "survival term not invariant to a consistent row permutation")
  assert_near(res_perm$loglik_genomics, res_orig$loglik_genomics, tol = 1e-8,
              "genomics term (read from fit history) should be unaffected by permutation")
})

run_test("BIC-T11: degenerate fit (zero EF column, zero EBeta) returns finite BIC", {
  fit_degen <- fit
  fit_degen$EF[, 1] <- 0
  fit_degen$EF_norms[1] <- 1e-10  # matches the fit's own floor for an all-zero column
  fit_degen$EBeta[1] <- 0
  res <- compute_joint_ll_bic(fit_degen, Y, time, status)
  assert_finite(res$bic, "BIC should stay finite for a degenerate (ARD-pruned) column")
  assert_finite(res$loglik_survival, "survival term should stay finite (no NaN from 0/0)")
})

run_test("BIC-T12: stop()s (not a wrong number) on a dimension mismatch", {
  err <- tryCatch({
    compute_joint_ll_bic(fit, Y[1:10, , drop = FALSE], time, status)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error on Y/time/status dimension mismatch")
})

run_test("BIC-T13: stop()s on a fit missing required fields", {
  err <- tryCatch({
    compute_joint_ll_bic(list(EF = fit$EF), Y, time, status)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error on a fit missing EF_norms/EBeta/history")
})

run_test("BIC-T14: stop()s on a short elbo_proxy history", {
  fit_short <- fit
  fit_short$history$elbo_proxy <- fit_short$history$elbo_proxy[1:2]
  fit_short$history$n_iter <- 10L
  err <- tryCatch({
    compute_joint_ll_bic(fit_short, Y, time, status)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error on a short elbo_proxy history")
})

run_test("BIC-T15: stop()s on a non-finite elbo_proxy entry", {
  fit_nonfinite <- fit
  fit_nonfinite$history$elbo_proxy[fit_nonfinite$history$n_iter] <- NA_real_
  err <- tryCatch({
    compute_joint_ll_bic(fit_nonfinite, Y, time, status)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error on a non-finite genomics term")
})

report_results("test_compute_bic.R")
