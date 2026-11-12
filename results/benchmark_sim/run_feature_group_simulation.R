# ============================================================
# Script: run_feature_group_simulation.R
# Purpose: Evaluate grouped YFB recovery and modality precision in the fixed Task D simulation.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, ebnm, survival, existing YFB fitter and multimodal simulator
# ============================================================
# Run from the repository root. Sequential fits use one CPU worker. Set BLAS
# thread counts to 1. Small synthetic matrices and private fit checkpoints are
# kept separately from the aggregate CSVs. No patient data are required.
if (!exists("fit_ssbmf_yfb")) source("code/fit_ssbmf_yfb.R")
if (!exists("predict_ssbmf_yfb")) source("code/predict_ssbmf_yfb.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")
if (!exists("frozen_reverse_cindex")) source("code/concordance_ci.R")

#' Standardize each block using its training rows only
#' @param training,validation Matching named lists of numeric matrices.
#' @return Training/validation standardized blocks and training centers/SDs.
#' @examples
#' fgsim_scale(list(a = matrix(1:6, 3)), list(a = matrix(c(4, 7), 1)))
fgsim_scale <- function(training, validation) {
  if (!identical(names(training), names(validation))) stop("Modality names differ.")
  center <- lapply(training, colMeans)
  sd <- lapply(training, function(Y) apply(Y, 2, stats::sd))
  if (any(vapply(sd, function(s) any(!is.finite(s) | s <= 0), logical(1))))
    stop("A simulated training feature has zero or nonfinite SD; no feature was dropped.")
  #' Apply stored feature centers and SDs to one block
  #' @param Y Numeric patient-by-feature matrix.
  #' @param center,sd Training feature vectors.
  #' @return Standardized matrix.
  #' @examples
  #' # scale(training[[1]], center[[1]], sd[[1]])
  scale <- function(Y, center, sd) sweep(sweep(Y, 2, center), 2, sd, "/")
  list(training = Map(scale, training, center, sd),
       validation = Map(scale, validation, center, sd), center = center, sd = sd)
}

#' Compute each modality's share of the shared-score precision
#' @param fit Returned YFB moments and feature precisions.
#' @param groups One modality label per feature.
#' @return Long data frame of Σ_j τ_mj E[F²_mjk] and its share for each program.
#' @examples
#' fgsim_precision(list(Tau = c(1, 2), EF2 = matrix(c(2, 3), 2)), c("a", "b"))
fgsim_precision <- function(fit, groups) {
  rows <- lapply(seq_len(ncol(fit$EF2)), function(k) {
    precision <- vapply(unique(groups), function(g)
      sum(fit$Tau[groups == g] * fit$EF2[groups == g, k]), numeric(1))
    data.frame(program = k, modality = names(precision), precision = unname(precision),
      share = if (sum(precision) > 0) unname(precision) / sum(precision) else NA_real_)
  })
  do.call(rbind, rows)
}

#' Score recovery of the third, shared prognostic program
#'
#' Uses the existing multimodal comparator's |loading correlation| ≥ 0.7 and
#' survival-activity rule. Reconstruction truth is F_m / training_sd_m after
#' z-scoring. Raw-truth correlations are also retained to expose scaling effects.
#' @param fit Returned posterior moments.
#' @param truth Named list of true loading matrices in generator units.
#' @param sd Training feature SD vectors, with matching modality names.
#' @param z_threshold Configured coefficient activity threshold.
#' @return Recovery flag, match/correlations, false-active count and program table.
#' @examples
#' # fgsim_recovery(fit, data$truth$F, scaled$sd, cfg$multimodal_yfb$prognostic_z_threshold)
fgsim_recovery <- function(fit, truth, sd, z_threshold) {
  standardized <- do.call(rbind, Map(function(F, s) sweep(F, 1, s, "/"), truth, sd))
  raw <- do.call(rbind, truth)
  correlation <- raw_correlation <- rep(NA_real_, ncol(fit$EF))
  for (k in seq_len(ncol(fit$EF))) {
    if (stats::sd(fit$EF[, k]) > 0) {
      correlation[k] <- abs(stats::cor(standardized[, 3], fit$EF[, k]))
      raw_correlation[k] <- abs(stats::cor(raw[, 3], fit$EF[, k]))
    }
  }
  beta_sd <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  active <- abs(fit$EBeta) > z_threshold * beta_sd
  matched <- if (any(is.finite(correlation))) which.max(replace(correlation, is.na(correlation), -Inf)) else NA_integer_
  loading_hit <- !is.na(matched) && correlation[matched] >= 0.7
  recovered <- loading_hit && active[matched]
  true_matches <- if (loading_hit) matched else integer()
  list(recovered = recovered, matched_program = matched,
    loading_correlation = if (is.na(matched)) NA_real_ else correlation[matched],
    raw_loading_correlation = if (any(is.finite(raw_correlation))) max(raw_correlation, na.rm = TRUE) else NA_real_,
    false_active = as.integer(sum(active[setdiff(seq_along(active), true_matches)])),
    programs = data.frame(program = seq_along(active), beta = fit$EBeta,
      beta_sd = beta_sd, active = active, loading_correlation = correlation,
      raw_loading_correlation = raw_correlation))
}

#' Run one prespecified seed/feature-ratio cell and save its private checkpoint
#' @param seed One of 1 through 10.
#' @param ratio Methylation feature count divided by expression count, 1 or 5.
#' @param cfg Parsed config/globals.yml.
#' @param private Directory under ignored data/ for synthetic fits and warnings.
#' @return Aggregate metrics, program and precision tables, and warning records.
#' @examples
#' # fgsim_one(1, 1, yaml::read_yaml("config/globals.yml"), "data/feature_group_simulation_2026-11-12")
fgsim_one <- function(seed, ratio, cfg, private) {
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  p <- eval(formals(simulate_multimodal_yfb_data)$p_expression)
  base <- data.frame(seed = seed, feature_ratio = ratio, p_expression = p,
                     p_methylation = p * ratio, K = cfg$benchmark$k_working_yfb)
  result <- tryCatch(withCallingHandlers({
    data <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
      p_expression = p, p_methylation = p * ratio, seed = seed)
    scaled <- fgsim_scale(data$training$Y, data$validation$Y)
    Y <- do.call(cbind, scaled$training)
    groups <- rep(names(scaled$training), vapply(scaled$training, ncol, integer(1)))
    fit <- fit_ssbmf_yfb(Y, data$training$time, data$training$event,
      K = cfg$benchmark$k_working_yfb, max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
      prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
      w_S = 1, alpha = cfg$benchmark$alpha, feature_group = groups, verbose = FALSE)
    recovery <- fgsim_recovery(fit, data$truth$F, scaled$sd,
                               cfg$multimodal_yfb$prognostic_z_threshold)
    precision <- fgsim_precision(fit, groups)
    new_Y <- do.call(cbind, scaled$validation)
    score <- predict_ssbmf_yfb(new_Y, fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
    c_index <- frozen_reverse_cindex(score, data$validation$time, data$validation$event)
    saveRDS(list(fit = fit, data = data, scaled = scaled),
            file.path(private, sprintf("ratio%d_seed%d.rds", ratio, seed)))
    list(metrics = cbind(base, data.frame(error = "", recovered = recovery$recovered,
      matched_program = recovery$matched_program, loading_correlation = recovery$loading_correlation,
      raw_loading_correlation = recovery$raw_loading_correlation,
      false_active = recovery$false_active, converged = isTRUE(fit$history$converged),
      n_iter = fit$history$n_iter, c_validation = c_index,
      elbo = tail(fit$history$elbo_full, 1))),
      programs = cbind(base[rep(1, nrow(recovery$programs)), ], recovery$programs),
      precision = cbind(base[rep(1, nrow(precision)), ], precision))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    # Warnings stay visible in the log and are retained individually below.
  }), error = function(e) list(metrics = cbind(base, data.frame(error = conditionMessage(e),
    recovered = FALSE, matched_program = NA_integer_, loading_correlation = NA_real_,
    raw_loading_correlation = NA_real_, false_active = NA_integer_, converged = FALSE,
    n_iter = NA_integer_, c_validation = NA_real_, elbo = NA_real_)),
    programs = NULL, precision = NULL))
  result$metrics$elapsed_seconds <- proc.time()[["elapsed"]] - start
  result$metrics$n_warnings <- length(warnings)
  result$warnings <- if (length(warnings)) cbind(base[rep(1, length(warnings)), ],
    data.frame(warning = warnings)) else NULL
  saveRDS(result, file.path(private, sprintf("ratio%d_seed%d_aggregate.rds", ratio, seed)))
  result
}

#' Run the fixed 20-cell study sequentially and write aggregate results
#' @param out Aggregate output directory. Existing directories are refused.
#' @param private Ignored checkpoint directory. Existing directories are refused.
#' @return Invisibly the per-cell aggregate results. Any failed cell yields a failing gate.
#' @examples
#' # fgsim_run("results/benchmark_sim/outputs/feature_group_simulation_2026-11-12",
#' #           "data/feature_group_simulation_2026-11-12")
fgsim_run <- function(out, private) {
  if (dir.exists(out) || dir.exists(private)) stop("Refusing to overwrite an existing study.")
  cfg <- yaml::read_yaml("config/globals.yml")
  if (cfg$benchmark$k_working_yfb != 5L) stop("Task D requires K = 5.")
  dir.create(out, recursive = TRUE); dir.create(private, recursive = TRUE)
  sources <- c("config/globals.yml", "code/fit_ssbmf_yfb.R", "code/feature_groups_yfb.R",
    "code/multimodal_yfb_updates.R", "code/loading_prior_fixed_q.R",
    "code/update_L_surv_YFB.R", "code/update_F_surv_YFB.R", "code/update_beta.R",
    "code/update_tau.R", "code/compute_elbo.R", "code/simulate_multimodal_yfb.R",
    "results/benchmark_sim/run_feature_group_simulation.R")
  write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
            file.path(out, "source_hashes.csv"), row.names = FALSE)
  writeLines(c("scenario=low_variance_prognostic_factor", "seeds=1:10", "feature_ratios=1:1,1:5",
    "generator=defaults (120 training, 120 validation, 50 expression, noise_scale=1, truncate=TRUE)",
    "scaling=training-only z-score within modality", "prior_LF=point_laplace",
    "prior_beta=normal per-program", "canonicalize=TRUE", "w_S=1", "prune=FALSE",
    "strata_id=NULL", "center_projection=FALSE", "workers=1",
    paste0("K=", cfg$benchmark$k_working_yfb), paste0("max_iter=", cfg$cavi$max_iter),
    paste0("tol=", cfg$cavi$tol), "loading_recovery=absolute correlation >=0.7 on standardized reconstruction truth",
    paste0("activity_z=", cfg$multimodal_yfb$prognostic_z_threshold)), file.path(out, "settings.txt"))
  results <- list()
  for (ratio in c(1L, 5L)) for (seed in 1:10) {
    cell <- fgsim_one(seed, ratio, cfg, private)
    results[[length(results) + 1L]] <- cell
    for (type in c("metrics", "programs", "precision", "warnings")) {
      rows <- do.call(rbind, lapply(results, `[[`, type))
      if (!is.null(rows)) write.csv(rows, file.path(out, paste0(type, ".csv")), row.names = FALSE)
    }
    cat(sprintf("Completed ratio 1:%d seed %d/10: recovered=%s converged=%s error=%s\n",
      ratio, seed, cell$metrics$recovered, cell$metrics$converged, cell$metrics$error))
    flush.console()
  }
  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  summary <- do.call(rbind, lapply(c(1L, 5L), function(ratio) {
    x <- metrics[metrics$feature_ratio == ratio, ]
    data.frame(feature_ratio = ratio, seeds = nrow(x), failed = sum(nzchar(x$error)),
      recovered = sum(x$recovered), mean_false_active = mean(x$false_active),
      converged = sum(x$converged),
      passed = !any(nzchar(x$error)) && sum(x$recovered) >= 8 &&
        mean(x$false_active) <= 1 && sum(x$converged) >= 9)
  }))
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  print(summary)
  invisible(results)
}

if (sys.nframe() == 0L) {
  stop("Task D simulation is paused pending Andrew's choice of recovery target: prognostic D F or reconstruction D^-1 F.")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  run_date <- format(Sys.Date(), "%Y-%m-%d")
  fgsim_run(paste0("results/benchmark_sim/outputs/feature_group_simulation_", run_date),
             paste0("data/feature_group_simulation_", run_date))
}
