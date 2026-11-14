# ============================================================
# Script: run_bounded_risk_screen.R
# Purpose: Run the registered S-17 fixed-prior convergence screen.
# Author: Andrew Walther
# Created: 2026-11-14
# Dependencies: survival, ebnm, ashr, yaml; S-16 runners
# ============================================================
# Serial execution uses one worker. Set BLAS threads before starting R.
# Checkpoints contain patient-level fits and must remain under data/.

#' Assess the registered post-fit stability rule
#' @param fit Fitted model with a per-program S-16 trace.
#' @param sigma Fixed prior standard deviation.
#' @return One-row data frame with convergence, stability and coefficient checks.
#' @examples
#' # s17_stability(fit, 0.5)
s17_stability <- function(fit, sigma) {
  tr <- fit$history$s16
  required <- c("iteration", "eta_rel_change", "eta_start_max", "eta_end_max", "beta_max")
  if (!is.data.frame(tr) || !all(required %in% names(tr)) || anyNA(tr$iteration))
    stop("Missing or malformed S-16 trace")
  # All program rows within a sweep must agree on the sweep-level quantities.
  groups <- split(tr, tr$iteration)
  fields <- required[-1L]
  if (any(vapply(groups, function(g) any(vapply(fields, function(f)
    length(unique(g[[f]])) != 1L, logical(1))), logical(1))))
    stop("Inconsistent sweep-level trace values")
  sweep <- tr[!duplicated(tr$iteration), , drop = FALSE]
  sweep <- sweep[order(sweep$iteration), , drop = FALSE]
  n <- nrow(sweep); assessable <- n >= 50L
  window <- if (assessable) tail(sweep, 50L) else sweep
  valid <- assessable && all(is.finite(as.matrix(window[, fields])))
  drift <- if (valid) abs(tail(window$eta_end_max, 1) - window$eta_start_max[1]) /
    max(window$eta_start_max[1], .Machine$double.eps) else NA_real_
  max_deviation <- if (valid) max(abs(window$eta_end_max - window$eta_start_max[1])) /
    max(window$eta_start_max[1], .Machine$double.eps) else NA_real_
  max_beta <- if (length(fit$EBeta) && all(is.finite(fit$EBeta))) max(abs(fit$EBeta)) else NA_real_
  converged <- isTRUE(fit$history$converged)
  native <- isTRUE(fit$history$native_elbo_converged)
  stable <- valid && converged && native && all(window$eta_rel_change < 1e-3) && drift < .01
  beta_pass <- is.finite(max_beta) && max_beta <= 4 * sigma
  data.frame(sweeps = n, assessable = assessable, fitter_converged = converged, native_elbo_converged = native,
    final_eta_max = tail(sweep$eta_end_max, 1), peak_eta_max = max(sweep$eta_end_max),
    peak_beta_max = max(sweep$beta_max),
    max_eta_rel_change = if (valid) max(window$eta_rel_change) else NA_real_,
    eta_endpoint_drift = drift, eta_max_deviation = max_deviation,
    stable = stable, max_beta = max_beta, beta_limit = 4 * sigma,
    beta_pass = beta_pass, pass = stable && beta_pass)
}

#' Construct a registered fixed-variance coefficient prior
#' @param family Normal or Laplace.
#' @param sigma Positive prior standard deviation.
#' @return Named opt-in fitter controls.
#' @examples
#' # s17_controls("normal", 0.5)
s17_controls <- function(family, sigma) {
  stopifnot(family %in% c("normal", "laplace"), length(sigma) == 1L,
            is.finite(sigma), sigma > 0)
  prior <- if (family == "normal") ashr::normalmix(1, 0, sigma) else
    ebnm::laplacemix(1, 0, sigma / sqrt(2))
  list(prior_beta = if (family == "normal") "normal" else "point_laplace",
       beta_prior_mode = "shared", beta_prior_fix = TRUE,
       beta_prior_init = prior, trace_s16 = TRUE, stable_convergence_window = 50L)
}

#' Run and checkpoint one dataset and prior combination
#' @param dataset Dataset label: real, s9, or weak.
#' @param family,sigma Registered prior family and standard deviation.
#' @param d List with Y, time and status for this training dataset.
#' @param cfg Fitter configuration, including the appropriate rank.
#' @param signature Source, input and configuration signature.
#' @param private_dir Directory for patient-level checkpoint files.
#' @return Completed checkpoint, including all fits and calibration histories.
#' @examples
#' # s17_case("real", "normal", 0.5, d, cfg, signature, "data/s17")
s17_case <- function(dataset, family, sigma, d, cfg, signature, private_dir) {
  key <- sprintf("%s_%s_%s", dataset, family, format(sigma, trim = TRUE, scientific = FALSE))
  path <- file.path(private_dir, paste0(key, ".rds"))
  arms <- list(); prior_elapsed <- 0
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old$schema, 1L) || !identical(old$signature, signature) ||
        !identical(old$dataset, dataset) || !identical(old$family, family) ||
        !identical(old$sigma, sigma)) stop("Stale checkpoint: ", path)
    if (!is.null(old$error)) stop("Checkpoint records an unresolved error: ", old$error)
    if (isTRUE(old$complete)) return(old)
    arms <- old$arms; prior_elapsed <- old$elapsed
  }
  controls <- s17_controls(family, sigma)
  start <- proc.time()[["elapsed"]]
  #' Save all completed arms atomically, including any explicit failure
  #' @param complete Whether all arms have been processed.
  #' @param error Unexpected error message, or NULL.
  #' @return Checkpoint object.
  #' @examples
  #' # save_progress(FALSE)
  save_progress <- function(complete, error = NULL) {
    result <- list(schema = 1L, complete = complete, signature = signature, dataset = dataset,
      family = family, sigma = sigma, elapsed = prior_elapsed + proc.time()[["elapsed"]] - start,
      arms = arms, error = error)
    s16_atomic_save(result, path)
    result
  }
  tryCatch({
    if (dataset == "s9") {
      for (weight in c(40, 100, 1000)) {
        name <- paste0("w", weight)
        if (!is.null(arms[[name]])) next
        run <- s16_fit(d$Y, d$time, d$status, weight, cfg, fit_control = controls)
        arms[[name]] <- list(fits = list(run$fit), weights = weight, capped = FALSE,
          warnings = list(run$warnings), elapsed = run$elapsed, guard = FALSE,
          checkpoints = list(s16_checkpoint(run$fit, d$Y, d$time, d$status)))
        save_progress(FALSE)
      }
    } else {
      if (is.null(arms$pilot)) {
        run <- s16_fit(d$Y, d$time, d$status, 1, cfg, fit_control = controls)
        arms$pilot <- list(fits = list(run$fit), weights = 1, capped = FALSE,
          warnings = list(run$warnings), elapsed = run$elapsed, guard = FALSE,
          checkpoints = list(s16_checkpoint(run$fit, d$Y, d$time, d$status)))
        save_progress(FALSE)
      }
      pilot <- arms$pilot$fits[[1]]
      for (target in c(.05, .25, .5)) {
        name <- paste0("s", target)
        if (!is.null(arms[[name]])) next
        if (!isTRUE(pilot$history$converged)) {
          arms[[name]] <- list(fits = list(), guard = FALSE, error = "pilot did not converge")
        } else {
          cal <- s16_calibrate(pilot, d$Y, d$time, d$status, target, cfg,
                               fit_control = controls, retain_pass_fits = TRUE)
          guarded <- isTRUE(cal$no_signal_guard)
          fits <- if (guarded) list(pilot) else cal$pass_fits
          arms[[name]] <- list(fits = fits,
            weights = if (guarded) 1 else vapply(cal$history, `[[`, numeric(1), "weight"),
            capped = if (guarded) FALSE else vapply(cal$history, `[[`, logical(1), "capped"),
            warnings = if (guarded) list(character()) else lapply(cal$history, `[[`, "warnings"),
            elapsed = if (guarded) 0 else vapply(cal$history, `[[`, numeric(1), "elapsed"),
            guard = guarded, calibration = cal[setdiff(names(cal), c("fit", "pass_fits"))],
            checkpoints = lapply(fits, s16_checkpoint, Y = d$Y, time = d$time, status = d$status))
        }
        save_progress(FALSE)
      }
    }
    save_progress(TRUE)
  }, error = function(e) {
    save_progress(FALSE, conditionMessage(e))
    stop("S-17 case ", key, " failed: ", conditionMessage(e))
  })
}

#' Export trace, per-program shares and convergence rows for one checkpoint
#' @param object Completed checkpoint from s17_case.
#' @param out_dir Aggregate output directory.
#' @return One row per fit or guarded target.
#' @examples
#' # s17_export(checkpoint, "results/benchmark_sim/outputs/example")
s17_export <- function(object, out_dir) {
  rows <- list(); shares <- list()
  for (arm_name in names(object$arms)) {
    arm <- object$arms[[arm_name]]
    if (!length(arm$fits)) {
      rows[[length(rows) + 1L]] <- data.frame(dataset = object$dataset, family = object$family,
        sigma = object$sigma, arm = arm_name, pass_number = 0L, final = TRUE,
        weight = NA_real_, capped = FALSE, guard = arm$guard, elapsed = 0, warnings = "",
        error = if (is.null(arm$error)) "no fit" else arm$error,
        sweeps = 0L, assessable = FALSE, fitter_converged = FALSE, native_elbo_converged = FALSE,
        final_eta_max = NA_real_, peak_eta_max = NA_real_, peak_beta_max = NA_real_,
        max_eta_rel_change = NA_real_, eta_endpoint_drift = NA_real_, eta_max_deviation = NA_real_,
        stable = FALSE, max_beta = NA_real_, beta_limit = 4 * object$sigma,
        beta_pass = FALSE, pass = FALSE)
      next
    }
    for (i in seq_along(arm$fits)) {
      fit <- arm$fits[[i]]
      prefix <- sprintf("%s_%s_%s_%s_pass%d", object$dataset, object$family,
                        format(object$sigma, trim = TRUE, scientific = FALSE), arm_name, i)
      utils::write.csv(fit$history$s16, file.path(out_dir, paste0(prefix, "_trace.csv")), row.names = FALSE)
      p <- arm$checkpoints[[i]]$shares
      p$dataset <- object$dataset; p$family <- object$family; p$sigma <- object$sigma
      p$arm <- arm_name; p$pass_number <- i; p$weight <- arm$weights[i]
      shares[[length(shares) + 1L]] <- p
      rows[[length(rows) + 1L]] <- cbind(data.frame(dataset = object$dataset,
        family = object$family, sigma = object$sigma, arm = arm_name, pass_number = i,
        final = i == length(arm$fits), weight = arm$weights[i], capped = arm$capped[i], guard = arm$guard,
        elapsed = arm$elapsed[i], warnings = paste(arm$warnings[[i]], collapse = " | "), error = ""),
        s17_stability(fit, object$sigma))
    }
  }
  list(rows = do.call(rbind, rows), shares = if (length(shares)) do.call(rbind, shares) else NULL)
}

#' Execute the complete serial S-17 convergence screen
#' @param family One of normal or laplace.
#' @param stage Real data first, then simulations.
#' @return Invisibly, the case results for that family.
#' @examples
#' # s17_screen("normal", "real")
s17_screen <- function(family, stage) {
  stopifnot(family %in% c("normal", "laplace"), stage %in% c("real", "sim"))
  source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
  source("code/share_targeted_ws.R")
  source("results/benchmark_sim/run_yfb_supervision_strength.R")
  source("results/benchmark_sim/run_yfb_supervised_real_refit.R"); s11_load()
  source("results/benchmark_sim/run_share_targeted_weight.R")
  cfg <- yaml::read_yaml("config/globals.yml"); cfg$cavi$max_iter <- 500L
  stopifnot(cfg$benchmark$k_working_yfb == 5L, cfg$benchmark$k_merged_yfb_desurv == 7L)
  files <- unique(c(list.files("code", pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE),
    list.files("results/benchmark_sim", pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE),
    list.files("config", recursive = TRUE, full.names = TRUE), LOCAL))
  signature <- list(files = tools::md5sum(files), targets = c(.05, .25, .5),
                    weights = c(40, 100, 1000), max_iter = 500L, seed = 1L,
                    families = c("normal", "laplace"), sigma = c(.1, .25, .5, 1),
                    packages = vapply(c("ashr", "ebnm", "survival", "yaml"),
                      function(p) as.character(utils::packageVersion(p)), character(1)))
  if (anyNA(signature$files)) stop("Cannot hash required source or input files")
  private_dir <- "data/bounded_risk_2026-11-14/screen"
  out_dir <- "results/benchmark_sim/outputs/bounded_risk_2026-11-14"
  dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (stage == "sim") {
    real_grid <- expand.grid(family = c("normal", "laplace"), sigma = c(.1, .25, .5, 1),
                             stringsAsFactors = FALSE)
    real_paths <- file.path(private_dir, sprintf("real_%s_%s.rds", real_grid$family,
      vapply(real_grid$sigma, function(x) format(x, trim = TRUE, scientific = FALSE), character(1))))
    if (!all(file.exists(real_paths))) stop("All real-data cases must finish before simulation")
    real_cases <- lapply(real_paths, readRDS)
    if (!all(vapply(real_cases, function(x) isTRUE(x$complete) &&
                    identical(x$schema, 1L) && identical(x$signature, signature), logical(1))))
      stop("Real-data cases do not match the current source and settings")
  }
  results <- list()
  for (dataset in if (stage == "real") "real" else c("s9", "weak")) {
    if (dataset == "real") {
      raw <- s15_data()
      if (nrow(raw$Y) != 273L || ncol(raw$Y) != 2218L || sum(raw$status) != 139L ||
          any(!is.finite(raw$Y))) stop("Real training data do not match registration")
      d <- list(Y = scale(raw$Y, center = TRUE, scale = FALSE), time = raw$time, status = raw$status)
      case_cfg <- cfg; case_cfg$benchmark$k_merged_yfb_desurv <- cfg$benchmark$k_working_yfb
    } else if (dataset == "s9") {
      raw <- yfb_decision_data(1, 2000L, .01, truth_mode = "shared")
      d <- list(Y = s9_zscore(raw$Y), time = raw$time, status = raw$event); case_cfg <- cfg
    } else {
      raw <- s13_data(1, "weak", FALSE)
      d <- list(Y = raw$Ytr, time = raw$time, status = raw$event); case_cfg <- cfg
    }
    for (sigma in c(.1, .25, .5, 1)) {
      message(sprintf("S-17 %s %s sigma=%g", dataset, family, sigma))
      x <- s17_case(dataset, family, sigma, d, case_cfg, signature, private_dir)
      export <- s17_export(x, out_dir)
      s16_atomic_save(export, file.path(private_dir, sprintf("%s_%s_%s_export.rds",
        dataset, family, format(sigma, trim = TRUE, scientific = FALSE))))
      utils::write.csv(export$rows, file.path(out_dir, sprintf("%s_%s_%s_fits.csv",
        dataset, family, format(sigma, trim = TRUE, scientific = FALSE))), row.names = FALSE)
      if (!is.null(export$shares)) utils::write.csv(export$shares,
        file.path(out_dir, sprintf("%s_%s_%s_shares.csv", dataset, family,
          format(sigma, trim = TRUE, scientific = FALSE))), row.names = FALSE)
      results[[length(results) + 1L]] <- export
    }
  }
  invisible(results)
}

#' Assemble all eight family and prior cases after both partitions finish
#' @return Invisibly, fit-level rows and the real-data pair gate.
#' @examples
#' # s17_summarize()
s17_summarize <- function() {
  private_dir <- "data/bounded_risk_2026-11-14/screen"
  out_dir <- "results/benchmark_sim/outputs/bounded_risk_2026-11-14"
  keys <- expand.grid(dataset = c("real", "s9", "weak"), family = c("normal", "laplace"),
                      sigma = c(.1, .25, .5, 1), stringsAsFactors = FALSE)
  paths <- file.path(private_dir, sprintf("%s_%s_%s.rds", keys$dataset, keys$family,
    vapply(keys$sigma, function(x) format(x, trim = TRUE, scientific = FALSE), character(1))))
  if (!all(file.exists(paths))) stop("Missing case checkpoints: ", paste(paths[!file.exists(paths)], collapse = ", "))
  cases <- lapply(paths, readRDS)
  if (!all(vapply(cases, function(x) isTRUE(x$complete) && identical(x$schema, 1L), logical(1))) ||
      !all(vapply(cases, function(x) identical(x$signature, cases[[1]]$signature), logical(1))))
    stop("Checkpoints have incomplete or inconsistent signatures")
  hashes <- cases[[1]]$signature$files
  if (!identical(hashes, tools::md5sum(names(hashes)))) stop("Sources or input changed after fitting")
  results <- lapply(cases, s17_export, out_dir = out_dir)
  rows <- do.call(rbind, lapply(results, `[[`, "rows"))
  shares <- do.call(rbind, lapply(results, `[[`, "shares"))
  real <- rows[rows$dataset == "real" & rows$final, , drop = FALSE]
  gate <- do.call(rbind, lapply(split(real, interaction(real$family, real$sigma, drop = TRUE)), function(z)
    data.frame(family = z$family[1], sigma = z$sigma[1],
      pass = nrow(z) == 4L && setequal(z$arm, c("pilot", "s0.05", "s0.25", "s0.5")) && all(z$pass))))
  utils::write.csv(rows, file.path(out_dir, "screen_fits.csv"), row.names = FALSE)
  utils::write.csv(shares, file.path(out_dir, "screen_program_shares.csv"), row.names = FALSE)
  utils::write.csv(gate, file.path(out_dir, "screen_real_pair_gate.csv"), row.names = FALSE)
  invisible(list(rows = rows, shares = shares, gate = gate))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "summarize")) {
    s17_summarize()
  } else if (length(args) == 2L && args[1] %in% c("normal", "laplace") &&
             args[2] %in% c("real", "sim")) {
    s17_screen(args[1], args[2])
  } else stop("Usage: Rscript results/benchmark_sim/run_bounded_risk_screen.R normal|laplace real|sim OR summarize")
}
