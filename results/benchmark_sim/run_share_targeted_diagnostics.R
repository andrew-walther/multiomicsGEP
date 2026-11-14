# ============================================================
# Script: run_share_targeted_diagnostics.R
# Purpose: Trace S-16 raw-weight mechanism diagnostics on registered training data.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: survival, ebnm, yaml, S-9 runner
# ============================================================
# One worker and one BLAS thread. On Longleaf allocate one CPU and enough memory
# for the training matrix and posterior moments; checkpoints contain patient-level fits.
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
s9_load()
cfg <- yaml::read_yaml("config/globals.yml")
cfg$cavi$max_iter <- 500L
out <- "results/benchmark_sim/outputs/share_targeted_ws_2026-11-13"
private <- "data/share_targeted_ws_2026-11-13/diagnostics"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(private, recursive = TRUE, showWarnings = FALSE)

#' Run one instrumented A1 fit and save the complete state and warnings.
#' @param Y Centered training matrix. @param time,status Survival outcomes.
#' @param weight Survival weight. @param K Number of programs. @param name Checkpoint label.
#' @return Full fit invisibly. An error is saved and rethrown.
#' @examples
#' # s16_diagnostic_fit(Y, time, status, 10, 5, "example")
s16_diagnostic_fit <- function(Y, time, status, weight, K, name, ramp = 0L) {
  path <- file.path(private, paste0(name, ".rds"))
  if (file.exists(path)) stop("Checkpoint already exists: ", path)
  warnings <- character(); start <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers(
    fit_ssbmf_yfb(Y, time, status, K = K, max_iter = cfg$cavi$max_iter,
                  tol = cfg$cavi$tol, prior_LF = "point_laplace", prior_beta = "normal",
                  canonicalize = TRUE, w_S = weight, survival_block_weight = TRUE,
                  beta_posterior = "meanfield", alpha = cfg$benchmark$alpha,
                  trace_s16 = TRUE, survival_weight_ramp = ramp, verbose = TRUE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(result, "error")) {
    saveRDS(list(error = conditionMessage(result), warnings = warnings, elapsed = elapsed), path)
    stop("Diagnostic ", name, " failed: ", conditionMessage(result))
  }
  saveRDS(list(fit = result, center = attr(Y, "scaled:center"), warnings = warnings,
               elapsed = elapsed, weight = weight), path)
  trace <- result$history$s16
  trace$dataset <- name; trace$weight <- weight
  write.csv(trace, file.path(out, paste0(name, "_trace.csv")), row.names = FALSE)
  cat(sprintf("%s: %.1f seconds, %d sweeps, converged=%s, warnings=%d\n",
              name, elapsed, result$history$n_iter, result$history$converged, length(warnings)))
  invisible(result)
}

if (sys.nframe() == 0L && identical(commandArgs(TRUE)[1], "fallback_real")) {
  source("code/share_targeted_ws.R")
  z <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  Y <- z$Y[, match(z$sets$D4, z$genes), drop = FALSE]
  if (any(!is.finite(Y))) stop("D4 features missing or non-finite")
  Y <- scale(Y, center = TRUE, scale = FALSE)
  pilot <- readRDS(file.path(private, "real_pilot_w1.rds"))$fit
  for (target in c(.01, .05, .1)) {
    calibration <- share_targeted_ws(pilot, Y, z$time, z$status, target)
    if (calibration$no_signal_guard) stop("Real pilot has no survival-active program")
    saveRDS(calibration, file.path(private, sprintf("real_target_%s_calibration.rds", target)))
    s16_diagnostic_fit(Y, z$time, z$status, calibration$w_S, cfg$benchmark$k_working_yfb,
                       sprintf("real_target_%s", target))
  }
} else if (sys.nframe() == 0L && identical(commandArgs(TRUE)[1], "ramp_real")) {
  z <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  Y <- z$Y[, match(z$sets$D4, z$genes), drop = FALSE]
  if (any(!is.finite(Y))) stop("D4 features missing or non-finite")
  Y <- scale(Y, center = TRUE, scale = FALSE)
  for (target in c(.25, .5)) {
    calibration <- readRDS(file.path(private, sprintf("real_target_%s_calibration.rds", target)))
    s16_diagnostic_fit(Y, z$time, z$status, calibration$w_S, cfg$benchmark$k_working_yfb,
                       sprintf("real_target_%s_ramp50", target), ramp = 50L)
  }
} else if (sys.nframe() == 0L && identical(commandArgs(TRUE)[1], "ramp_sim")) {
  d <- yfb_decision_data(1, 2000L, .01, truth_mode = "shared")
  Y <- s9_zscore(d$Y)
  for (w in c(10, 40, 100, 1000))
    s16_diagnostic_fit(Y, d$time, d$event, w, cfg$benchmark$k_merged_yfb_desurv,
                       sprintf("simulation_w%g_ramp50", w), ramp = 50L)
} else if (sys.nframe() == 0L && identical(commandArgs(TRUE)[1], "real")) {
  source("code/share_targeted_ws.R")
  z <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  Y <- z$Y[, match(z$sets$D4, z$genes), drop = FALSE]
  if (anyNA(Y)) stop("D4 features missing or non-finite")
  Y <- scale(Y, center = TRUE, scale = FALSE)
  pilot <- s16_diagnostic_fit(Y, z$time, z$status, 1, cfg$benchmark$k_working_yfb,
                              "real_pilot_w1")
  for (target in c(.25, .5)) {
    calibration <- share_targeted_ws(pilot, Y, z$time, z$status, target)
    if (calibration$no_signal_guard) stop("Real pilot has no survival-active program")
    saveRDS(calibration, file.path(private, sprintf("real_target_%s_calibration.rds", target)))
    s16_diagnostic_fit(Y, z$time, z$status, calibration$w_S, cfg$benchmark$k_working_yfb,
                       sprintf("real_target_%s", target))
  }
} else if (sys.nframe() == 0L) {
  d <- yfb_decision_data(1, 2000L, .01, truth_mode = "shared")
  Y <- s9_zscore(d$Y)
  for (w in c(10, 40, 100, 1000))
    s16_diagnostic_fit(Y, d$time, d$event, w, cfg$benchmark$k_merged_yfb_desurv,
                       sprintf("simulation_w%g", w))
}
