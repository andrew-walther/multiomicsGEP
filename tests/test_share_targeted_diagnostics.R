# ============================================================
# Script: test_share_targeted_diagnostics.R
# Purpose: Verify opt-in S-16 instrumentation preserves the fit and reports exact shares.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: survival, ebnm
# ============================================================
source("code/fit_ssbmf_yfb.R")
set.seed(42)
Y <- matrix(rnorm(240), 30, 8)
time <- rexp(30); status <- rep(c(0, 1), 15)
args <- list(Y = Y, time = time, status = status, K = 2, max_iter = 3,
             canonicalize = TRUE, w_S = 2, prior_LF = "point_laplace",
             prior_beta = "normal", survival_block_weight = TRUE, verbose = FALSE)
set.seed(55); plain <- do.call(fit_ssbmf_yfb, args)
set.seed(55); traced <- do.call(fit_ssbmf_yfb, c(args, list(trace_s16 = TRUE)))
stopifnot(is.null(plain$history$s16), is.null(plain$CoxW))
traced$history$s16 <- NULL; traced$CoxW <- NULL
stopifnot(identical(plain, traced))
set.seed(55); traced <- do.call(fit_ssbmf_yfb, c(args, list(trace_s16 = TRUE)))
k <- 1; f <- traced$EF[, k] / sqrt(sum(traced$EF[, k]^2))
a <- sum(traced$EL2[, k]) * sum(traced$Tau * f^2)
b <- 2 * traced$EBeta2[k] * sum(traced$CoxW * as.vector(Y %*% f)^2)
stopifnot(isTRUE(all.equal(tail(traced$history$s16$directional_share[traced$history$s16$program == k], 1), b / (a + b))))
cat("S-16 instrumentation checks passed\n")
# Ramp follows the registered schedule and cannot stop while its objective changes.
set.seed(55); zero <- do.call(fit_ssbmf_yfb, c(args, list(survival_weight_ramp = 0L)))
stopifnot(identical(plain, zero))
ramp_args <- args; ramp_args$max_iter <- 52L; ramp_args$tol <- 1e9
set.seed(55); ramped <- do.call(fit_ssbmf_yfb, c(ramp_args, list(trace_s16 = TRUE,
                                                               survival_weight_ramp = 50L)))
r <- ramped$history$s16[ramped$history$s16$program == 1, ]
stopifnot(ramped$history$n_iter == 51L, r$weight[1] == 1,
          isTRUE(all.equal(r$weight[25], exp(log(2) * 24 / 49))),
          r$weight[50] == 2, r$weight[51] == 2,
          all(r$ramp[1:50]), !r$ramp[51])
cat("S-16 ramp checks passed\n")
# Compare against a fixture produced by the pre-change fitter at w_S = 1.
source("tests/fixtures/s16_w1_baseline.R")
set.seed(4242); baseline_Y <- matrix(rnorm(160), 20, 8)
baseline_time <- rexp(20); baseline_status <- rep(c(0, 1), 10)
set.seed(42)
current_w1 <- fit_ssbmf_yfb(baseline_Y, baseline_time, baseline_status, K = 2,
  max_iter = 60, canonicalize = TRUE, w_S = 1, prior_LF = "point_laplace",
  prior_beta = "normal", survival_block_weight = TRUE, verbose = FALSE)
stopifnot(identical(current_w1, s16_w1_baseline))
cat("S-16 pre-change w=1 baseline check passed\n")
