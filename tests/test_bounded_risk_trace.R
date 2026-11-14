# ============================================================
# Script: test_bounded_risk_trace.R
# Purpose: Check the per-program canonicalization trace for S-17.
# Author: Andrew Walther
# Created: 2026-11-14
# Dependencies: survival, ebnm
# ============================================================

source("code/fit_ssbmf_yfb.R")
set.seed(42)
Y <- matrix(rnorm(240), 30, 8)
time <- rexp(30)
status <- rep(c(0, 1), 15)
args <- list(Y = Y, time = time, status = status, K = 2, max_iter = 3,
             canonicalize = TRUE, w_S = 2, prior_LF = "point_laplace",
             prior_beta = "normal", survival_block_weight = TRUE, verbose = FALSE)
set.seed(55)
fit <- do.call(fit_ssbmf_yfb, c(args, list(trace_s16 = TRUE)))
trace <- fit$history$s16
stopifnot(all(c("score_sd_before", "canonical_scale", "beta_before", "beta_after") %in% names(trace)))
stopifnot(all(is.finite(trace$score_sd_before)), all(is.finite(trace$beta_before)),
          any(abs(trace$canonical_scale - 1) > 1e-4),
          any(abs(trace$beta_before) > 1e-4),
          isTRUE(all.equal(trace$beta_after, trace$canonical_scale * trace$beta_before)))
last <- trace[trace$iteration == fit$history$n_iter, ]
stopifnot(isTRUE(all.equal(last$beta_after, as.vector(fit$EBeta))),
          isTRUE(all.equal(last$score_sd_before,
                           apply(Y %*% fit$EF, 2, sd) * last$canonical_scale)))
set.seed(55)
plain <- do.call(fit_ssbmf_yfb, args)
fit$history$s16 <- NULL
fit$CoxW <- NULL
stopifnot(identical(plain, fit))
cat("S-17 canonicalization trace checks passed\n")

# A loose ELBO tolerance stops the native fit early. The opt-in gate must
# collect the entire registered window before reporting stable convergence.
window_args <- args
set.seed(91)
signal <- rnorm(80)
window_args$Y <- scale(outer(signal, rnorm(12)) + matrix(rnorm(960, sd = .3), 80, 12))
window_args$time <- rexp(80, exp(signal))
window_args$status <- rep(1, 80)
window_args$max_iter <- 55L
window_args$tol <- 1e9
window_args$prior_LF <- "normal"
window_args$beta_prior_mode <- "shared"
window_args$beta_prior_fix <- TRUE
window_args$beta_prior_init <- ashr::normalmix(1, 0, .5)
set.seed(55)
native <- do.call(fit_ssbmf_yfb, c(window_args, list(trace_s16 = TRUE)))
set.seed(55)
stable <- do.call(fit_ssbmf_yfb, c(window_args, list(trace_s16 = TRUE,
                                                      stable_convergence_window = 50L)))
stopifnot(native$history$converged, native$history$n_iter < 50L,
          stable$history$n_iter >= 50L, isTRUE(stable$history$native_elbo_converged),
          any(abs(stable$EBeta) > 1e-4))
one <- stable$history$s16[stable$history$s16$program == 1, ]
last50 <- tail(one, 50)
window_pass <- all(is.finite(last50$eta_rel_change)) &&
  all(last50$eta_rel_change < 1e-3) &&
  is.finite(last50$eta_start_max[1]) && is.finite(tail(last50$eta_end_max, 1)) &&
  abs(tail(last50$eta_end_max, 1) - last50$eta_start_max[1]) /
    max(last50$eta_start_max[1], .Machine$double.eps) < .01
stopifnot(identical(stable$history$converged, window_pass))
window_args$max_iter <- 90L
set.seed(55)
settled <- do.call(fit_ssbmf_yfb, c(window_args, list(trace_s16 = TRUE,
                                                       stable_convergence_window = 50L)))
settled_trace <- settled$history$s16[settled$history$s16$program == 1, ]
settled_last <- tail(settled_trace, 50)
stopifnot(settled$history$converged, settled$history$n_iter >= 50L,
          settled$history$n_iter < 90L, settled$history$native_elbo_converged,
          all(settled_last$eta_rel_change < 1e-3),
          abs(tail(settled_last$eta_end_max, 1) - settled_last$eta_start_max[1]) /
            max(settled_last$eta_start_max[1], .Machine$double.eps) < .01,
          any(abs(settled$EBeta) > 1e-4))
bad <- try(do.call(fit_ssbmf_yfb, c(args, list(stable_convergence_window = 50L))),
           silent = TRUE)
stopifnot(inherits(bad, "try-error"))
cat("S-17 stability window checks passed\n")
