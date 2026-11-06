# ============================================================
# Script:  test_lb_orientation.R
# Purpose: Check that converged fit_ssbmf_lb() fits are oriented as Cox risk
#          scores (training C >= 0.5 with reverse = TRUE) without any
#          post-fit sign flip.
# Author: Andrew Walther
# Created: 2026-11-06
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_lb_orientation.R  (standalone)
# Dependencies: survival, ebnm
# ============================================================
#
# The LB sign check was first fixed (reverse = TRUE) and then removed with
# sign_correction on 2026-11-06 (DECISIONS.md): CAVI fits beta to the Cox
# likelihood given L, so a converged fit is already oriented. This test
# checks that property directly.

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival)
})
if (!exists("fit_ssbmf_lb")) {
  suppressMessages(tryCatch(source("code/fit_ssbmf_lb.R"),
                            error = function(e) invisible(NULL)))
}

cat("\n--- test_lb_orientation.R ---\n")

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
  list(Y = Y, time = pmin(ev, cen), status = as.integer(ev <= cen), eta = eta)
}

# Cox-direction concordance: larger eta = higher hazard = shorter time.
.lbsc_cindex <- function(eta, time, status) {
  as.numeric(concordance(Surv(time, status) ~ eta, reverse = TRUE)$concordance)
}

.lbsc_fit <- function(d) {
  suppressWarnings(fit_ssbmf_lb(d$Y, d$time, d$status, K = 3, max_iter = 300,
                                tol = 1e-4, verbose = FALSE))
}

for (.seed in 1:3) {
  run_test(sprintf("LBO-T1 (seed %d): converged LB fit has Cox-direction training C near the oracle", .seed), {
    d   <- .lbsc_sim(.seed)
    fit <- .lbsc_fit(d)
    assert_true(isTRUE(fit$history$converged), "fit converged")
    c_fit    <- .lbsc_cindex(as.vector(fit$EL %*% fit$EBeta), d$time, d$status)
    c_oracle <- .lbsc_cindex(d$eta, d$time, d$status)
    assert_true(c_fit >= 0.5, sprintf("training C (reverse = TRUE) = %.3f", c_fit))
    assert_true(c_fit >= c_oracle - 0.1,
                sprintf("training C %.3f is more than 0.1 below the oracle %.3f", c_fit, c_oracle))
  })
}
