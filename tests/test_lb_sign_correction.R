# ============================================================
# Script:  test_lb_sign_correction.R
# Purpose: Check fit_ssbmf_lb()'s post-convergence sign check (Phase C):
#          it must orient eta = L beta as a Cox risk score and must keep
#          the posterior second moment of beta when it flips the sign.
# Author: Andrew Walther
# Created: 2026-11-06
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_lb_sign_correction.R  (standalone)
# Dependencies: survival, ebnm
# ============================================================
#
# Background. survival::concordance(Surv(time, status) ~ x) treats a larger x
# as a longer survival time. A Cox risk score has the opposite direction, so
# its concordance needs reverse = TRUE. The LB check omitted it, so it flipped
# correctly oriented fits (the same bug was fixed in fit_ssbmf_yfb() on
# 2026-09-04). A flip of beta -> -beta also leaves E[beta^2] unchanged; the
# old code replaced it with EBeta^2 and so discarded Var(beta).
#
# CAVI fits beta to the Cox likelihood given L, so the unflipped fits below
# are already oriented (training C > 0.5). The correct check therefore leaves
# them unchanged, while the old check flipped them and erased Var(beta).

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival)
})
if (!exists("fit_ssbmf_lb")) {
  suppressMessages(tryCatch(source("code/fit_ssbmf_lb.R"),
                            error = function(e) invisible(NULL)))
}

cat("\n--- test_lb_sign_correction.R ---\n")

## Simulated LB data with a planted survival signal ----
# eta_true = L beta with beta = (1, -0.8, 0); exponential event times.
.lbsc_sim <- function(seed, n = 120, p = 100, K = 3, beta = c(1, -0.8, 0)) {
  set.seed(seed)
  L  <- matrix(rexp(n * K), n, K)
  Fm <- matrix(0, p, K)
  for (k in seq_len(K)) Fm[sample.int(p, 15), k] <- 1
  Y   <- L %*% t(Fm) + matrix(rnorm(n * p, sd = 0.5), n, p)
  eta <- as.vector(L %*% beta)
  ev  <- -log(runif(n)) / (0.05 * exp(eta))
  cen <- rexp(n, 1) * median(ev) * 2
  list(Y = Y, time = pmin(ev, cen), status = as.integer(ev <= cen))
}

# Cox-direction concordance: larger eta = higher hazard = shorter time.
.lbsc_cindex <- function(eta, time, status) {
  as.numeric(concordance(Surv(time, status) ~ eta, reverse = TRUE)$concordance)
}

.lbsc_fit <- function(d, sign_correction) {
  suppressWarnings(fit_ssbmf_lb(d$Y, d$time, d$status, K = 3, max_iter = 60,
                                tol = 1e-4, sign_correction = sign_correction,
                                verbose = FALSE))
}

for (.seed in 1:2) {
  .d   <- .lbsc_sim(.seed)
  .off <- .lbsc_fit(.d, FALSE)
  .on  <- .lbsc_fit(.d, TRUE)
  .c_off <- .lbsc_cindex(as.vector(.off$EL %*% .off$EBeta), .d$time, .d$status)
  .c_on  <- .lbsc_cindex(as.vector(.on$EL %*% .on$EBeta), .d$time, .d$status)

  run_test(sprintf("LBSC-T1 (seed %d): corrected fit has Cox-direction training C >= 0.5", .seed), {
    assert_true(.c_on >= 0.5,
                sprintf("training C (reverse = TRUE) = %.3f after sign correction", .c_on))
  })

  run_test(sprintf("LBSC-T2 (seed %d): sign check flips exactly when the uncorrected C < 0.5", .seed), {
    # Same SVD init and CAVI path, so the two fits differ only by Phase C.
    s <- if (.c_off < 0.5) -1 else 1
    assert_near(.on$EBeta, s * .off$EBeta, tol = 1e-10,
                msg = "EBeta equals the uncorrected EBeta times the expected sign")
    assert_near(.on$EL, .off$EL, tol = 1e-10, msg = "loadings are untouched by Phase C")
  })

  run_test(sprintf("LBSC-T3 (seed %d): E[beta^2] (and so Var(beta)) is preserved", .seed), {
    assert_near(.on$EBeta2, .off$EBeta2, tol = 1e-10,
                msg = "a sign flip leaves the second moment unchanged")
    assert_true(any(.off$EBeta2 - .off$EBeta^2 > 1e-4),
                "fixture has nonzero posterior variance, so T3 can detect its loss")
  })
}
