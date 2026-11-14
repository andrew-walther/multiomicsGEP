# ============================================================
# Script: test_bounded_risk_controls.R
# Purpose: Check opt-in S-17 controls and calibration pass retention.
# Author: Andrew Walther
# Created: 2026-11-14
# Dependencies: base R
# ============================================================
local({
source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
cfg <- list(benchmark = list(k_merged_yfb_desurv = 1, alpha = .5),
            cavi = list(max_iter = 3, tol = 1e-6))
calls <- list()
fit_ssbmf_yfb <- function(...) {
  args <- list(...); calls[[length(calls) + 1L]] <<- args
  list(posterior_state = list(g_beta = "fixed", id = length(calls)),
       history = list(converged = TRUE, trace = length(calls)))
}
Y <- matrix(1:4, 2)
fit <- s16_fit(Y, 1:2, c(1, 0), 1, cfg,
               fit_control = list(prior_beta = "point_laplace", trace_s16 = TRUE))
stopifnot(calls[[1]]$prior_beta == "point_laplace", calls[[1]]$trace_s16)
share_targeted_ws <- function(fit, Y, time, status, target, w_current, pilot_loading = NULL) {
  list(no_signal_guard = FALSE, pilot_loading = 1, w_S = w_current * 2,
       capped = FALSE, ratio = if (w_current < 4) 0.01 else target / (1 - target),
       lost_activity = FALSE, matching_correlation = 1,
       diagnostics = list(per_gene_median_share = .1), selected_program = 1)
}
cal <- s16_calibrate(fit$fit, Y, 1:2, c(1, 0), .1, cfg,
                     fit_control = list(prior_beta = "point_laplace"),
                     retain_pass_fits = TRUE)
stopifnot(length(cal$pass_fits) == 2L,
          identical(cal$pass_fits[[1]]$history$trace, 2L),
          identical(calls[[3]]$posterior_init, cal$pass_fits[[1]]$posterior_state),
          identical(calls[[3]]$posterior_init$g_beta, "fixed"))
cal_default <- s16_calibrate(fit$fit, Y, 1:2, c(1, 0), .1, cfg)
stopifnot(is.null(cal_default$pass_fits), calls[[4]]$prior_beta == "normal")
cat("Bounded risk controls: passed\n")
})

# The real fitter must retain the same fixed shared Laplace prior on warm start.
if (!exists("fit_ssbmf_yfb", mode = "function")) source("code/fit_ssbmf_yfb.R")
local({
  source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
  set.seed(432)
  Y <- scale(matrix(rnorm(30 * 12), 30, 12))
  time <- rexp(30); status <- rep(c(1L, 1L, 0L), 10)
  cfg <- list(benchmark = list(k_merged_yfb_desurv = 2L, alpha = .5),
              cavi = list(max_iter = 2L, tol = 1e-6))
  prior <- ebnm::laplacemix(1, 0, .25 / sqrt(2))
  controls <- list(prior_beta = "point_laplace", beta_prior_mode = "shared",
                   beta_prior_fix = TRUE, beta_prior_init = prior)
  first <- s16_fit(Y, time, status, 1, cfg, fit_control = controls)$fit
  second <- s16_fit(Y, time, status, 2, cfg, initial = first$posterior_state,
                    fit_control = controls)$fit
  stopifnot(identical(first$posterior_state$g_beta, prior),
            identical(second$posterior_state$g_beta, prior))
})
cat("Fixed shared Laplace warm start: passed\n")

# Verify the dataset dispatcher routes the two fits to their requested controls.
local({
  source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
  calls <- list()
  s16_fit <- function(Y, time, status, weight, cfg, initial = NULL, fit_control = list()) {
    calls$pilot <<- fit_control
    list(fit = list(history = list(converged = TRUE), tag = "pilot"))
  }
  s16_twostep <- function(Y, time, status, cfg, fit_control = list()) {
    calls$twostep <<- fit_control
    list(fit = list(history = list(converged = TRUE), tag = "twostep"))
  }
  s16_checkpoint <- function(fit, Y, time, status) list(tag = fit$tag)
  s16_sim_score <- function(fit, d, null) list(tag = fit$tag)
  d <- list(Ytr = matrix(1:4, 2), time = 1:2, event = c(1, 0))
  fixed <- list(beta_prior_mode = "shared", beta_prior_fix = TRUE)
  reference <- list(beta_prior_mode = "per_program")
  result <- s16_dataset(d, list(), targets = numeric(), fit_control = fixed,
                        reference_control = reference)
  stopifnot(identical(calls$pilot, fixed), identical(calls$twostep, reference),
            identical(result$pilot$metrics$tag, "pilot"),
            identical(result$twostep$metrics$tag, "twostep"))
  calls <- list()
  s16_dataset(d, list(), targets = numeric(), fit_control = fixed)
  stopifnot(identical(calls$pilot, fixed), identical(calls$twostep, fixed))
})
cat("Reference-control routing: passed\n")
