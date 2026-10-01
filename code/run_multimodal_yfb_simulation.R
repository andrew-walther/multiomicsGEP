# ============================================================
# Script: run_multimodal_yfb_simulation.R
# Purpose: Run the Stage 4 multimodal YFB smoke-study benchmark.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R; survival; flashier; ebnm; multimodal YFB modules
# ============================================================

source("code/simulate_multimodal_yfb.R")
source("code/fit_multimodal_yfb.R")
source("code/predict_multimodal_yfb.R")
suppressMessages(tryCatch(
  source("code/fit_cox_on_yf.R"),
  error = function(e) invisible(NULL)
))
source("code/predict_cox_on_yf.R")

#' Return a correlation only when both quantities have nonzero variation
#'
#' @param x,y Numeric vectors of equal length.
#' @return Pearson correlation, or `NA_real_` when the correlation is undefined.
#' @examples
#' multimodal_yfb_safe_correlation(c(1, 1), c(0, 0))
#' @family multimodal_yfb_simulation
multimodal_yfb_safe_correlation <- function(x, y) {
  if (length(x) != length(y) || any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must be finite vectors of equal length.")
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}

#' Fit a nonnegative two-step EBMF-to-Cox benchmark
#'
#' EBMF is deliberately survival-blind. Cox regression is fitted only after
#' the training projection is fixed, so validation outcomes cannot influence
#' either stage.
#'
#' @param Y Named training molecular blocks.
#' @param time,event Training survival outcomes.
#' @param K Maximum EBMF factor count.
#' @return Fitted stacked loadings, Cox coefficients, and training risk scores.
#' @examples
#' # fit_multimodal_yfb_ebmf_cox(Y, time, event, K = 3)
#' @family multimodal_yfb_simulation
fit_multimodal_yfb_ebmf_cox <- function(Y, time, event, K) {
  if (!requireNamespace("flashier", quietly = TRUE) ||
      !requireNamespace("ebnm", quietly = TRUE) ||
      !requireNamespace("survival", quietly = TRUE)) {
    stop("flashier, ebnm, and survival are required for the EBMF-to-Cox baseline.")
  }
  stacked <- do.call(cbind, Y)
  flash_fit <- flashier::flash(
    stacked, var_type = 2, greedy_Kmax = K, backfit = TRUE, verbose = 0,
    ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential)
  )
  loadings <- as.matrix(flashier::ldf(flash_fit, type = "2")$F)
  scores <- stacked %*% loadings %*%
    solve(crossprod(loadings) + diag(1e-8, ncol(loadings)))
  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ scores),
    error = function(e) NULL
  )
  coefficient <- if (is.null(cox_fit)) rep(0, ncol(scores)) else as.numeric(stats::coef(cox_fit))
  coefficient[!is.finite(coefficient)] <- 0
  list(loadings = loadings, coefficient = coefficient,
       training_risk = as.vector(scores %*% coefficient),
       n_factors = ncol(loadings))
}

#' Score molecular blocks with a frozen EBMF-to-Cox benchmark
#'
#' @param fit Result from [fit_multimodal_yfb_ebmf_cox()].
#' @param Y Named validation molecular blocks in training feature order.
#' @return Numeric validation risk score.
#' @examples
#' # predict_multimodal_yfb_ebmf_cox(fit, Y_validation)
#' @family multimodal_yfb_simulation
predict_multimodal_yfb_ebmf_cox <- function(fit, Y) {
  stacked <- do.call(cbind, Y)
  scores <- stacked %*% fit$loadings %*%
    solve(crossprod(fit$loadings) + diag(1e-8, ncol(fit$loadings)))
  as.vector(scores %*% fit$coefficient)
}

#' Fit a training-only stacked modality-risk Cox baseline
#'
#' Each modality first receives its own survival-blind EBMF-to-Cox score. A
#' final Cox regression combines those two *training-derived* scores; no
#' validation labels or recalibration enter this baseline.
#'
#' @param Y,time,event,K As in [fit_multimodal_yfb_ebmf_cox()].
#' @return Component EBMF fits and two Cox stacking coefficients.
#' @examples
#' # fit_multimodal_yfb_stacked_risk_cox(Y, time, event, K = 3)
#' @family multimodal_yfb_simulation
fit_multimodal_yfb_stacked_risk_cox <- function(Y, time, event, K) {
  components <- lapply(Y, function(block) {
    fit_multimodal_yfb_ebmf_cox(list(modality = block), time, event, K)
  })
  training_scores <- do.call(cbind, lapply(components, `[[`, "training_risk"))
  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ training_scores),
    error = function(e) NULL
  )
  coefficient <- if (is.null(cox_fit)) rep(0, ncol(training_scores)) else as.numeric(stats::coef(cox_fit))
  coefficient[!is.finite(coefficient)] <- 0
  list(components = components, coefficient = coefficient)
}

#' Score a frozen training-only stacked modality-risk Cox baseline
#'
#' @param fit Result from [fit_multimodal_yfb_stacked_risk_cox()].
#' @param Y Named validation molecular blocks.
#' @return Numeric validation risk score.
#' @examples
#' # predict_multimodal_yfb_stacked_risk_cox(fit, Y_validation)
#' @family multimodal_yfb_simulation
predict_multimodal_yfb_stacked_risk_cox <- function(fit, Y) {
  scores <- mapply(function(component, block) {
    predict_multimodal_yfb_ebmf_cox(component, list(modality = block))
  }, fit$components, Y, SIMPLIFY = FALSE)
  as.vector(do.call(cbind, scores) %*% fit$coefficient)
}

#' Fit the established single-modality YFB benchmark
#'
#' This benchmark uses the existing single-modality implementation directly.
#' It does not construct a degenerate placeholder block for the unused modality.
#'
#' @param Y,time,event As in [fit_multimodal_yfb()].
#' @param modality Either `"expression"` or `"methylation"`.
#' @param K Initial factor count.
#' @param control Fitter control list.
#' @return An established single-modality YFB fit carrying the selected modality.
#' @examples
#' # fit_multimodal_yfb_single_modality(Y, time, event, "expression", 3)
#' @family multimodal_yfb_simulation
fit_multimodal_yfb_single_modality <- function(Y, time, event, modality, K,
                                               control = list()) {
  if (!modality %in% names(Y)) stop("modality must be present in Y.")
  fit <- fit_cox_on_yf(
    Y[[modality]], time, event, K = K,
    max_iter = control$max_outer %||% 100L, tol = control$tol %||% 1e-5,
    prior_LF = "point_exponential", prior_beta = "normal", verbose = FALSE
  )
  attr(fit, "selected_modality") <- modality
  fit
}

#' Predict from a single-modality raw-YFB benchmark
#'
#' @param fit Result from [fit_multimodal_yfb_single_modality()].
#' @param Y Named validation molecular blocks.
#' @return Named risk scores and ignored validation features.
#' @examples
#' # predict_multimodal_yfb_single_modality(fit, Y_validation)
#' @family multimodal_yfb_simulation
predict_multimodal_yfb_single_modality <- function(fit, Y) {
  selected <- attr(fit, "selected_modality")
  predict_cox_on_yf(Y[[selected]], fit$EF, fit$EBeta, fit$EF_norms)
}

#' Calculate recovery, sparsity, and factor-selection metrics for one fit
#'
#' @param fit Fitted multimodal YFB object.
#' @param data Simulation object from [simulate_multimodal_yfb_data()].
#' @return Named scalar metrics.
#' @examples
#' # multimodal_yfb_fit_metrics(fit, simulated_data)
#' @family multimodal_yfb_simulation
multimodal_yfb_fit_metrics <- function(fit, data) {
  true_loadings <- do.call(rbind, data$truth$F)
  estimated_loadings <- do.call(rbind, fit$EF)
  matches <- multimodal_yfb_match_factors(true_loadings, estimated_loadings)
  tau_metrics <- unlist(Map(function(estimate, truth) {
    c(cor = stats::cor(estimate, truth),
      rmse = sqrt(mean((estimate - truth)^2)))
  }, fit$Tau, data$truth$Tau), use.names = TRUE)
  true_prognostic <- data$truth$beta != 0
  matched_estimates <- integer(length(true_prognostic))
  matched_estimates[matches$true_factor] <- matches$estimated_factor
  detected <- rep(FALSE, length(true_prognostic))
  matched_beta <- rep(NA_real_, length(true_prognostic))
  matched <- matched_estimates > 0L
  detected[matched] <- fit$diagnostics$prognostic[matched_estimates[matched]]
  matched_beta[matched] <- fit$EBeta[matched_estimates[matched]]
  list(
    loading_recovery = mean(matches$abs_correlation),
    factor_recovery = mean(matches$abs_correlation >= 0.70),
    tau_correlation = mean(tau_metrics[grep("cor", names(tau_metrics))]),
    tau_rmse = mean(tau_metrics[grep("rmse", names(tau_metrics))]),
    sparsity = mean(unlist(lapply(fit$EF, function(x) x == 0))),
    K_eff = fit$diagnostics$K_eff,
    K_eff_reconstruction = fit$diagnostics$K_eff_reconstruction,
    K_eff_survival = fit$diagnostics$K_eff_survival,
    K_eff_retained = fit$diagnostics$K_eff_retained,
    prognostic_recall = if (any(true_prognostic)) {
      mean(detected[true_prognostic])
    } else {
      mean(!detected)
    },
    prognostic_direction_recall = if (any(true_prognostic)) {
      mean(detected[true_prognostic] &
        sign(matched_beta[true_prognostic]) == sign(data$truth$beta[true_prognostic]))
    } else {
      NA_real_
    },
    converged = fit$diagnostics$converged
  )
}

#' Run one scenario/replicate across all Stage 4 benchmark arms
#'
#' @param scenario Scenario name.
#' @param replicate_id Replicate identifier.
#' @param seed RNG seed.
#' @param n_train,n_validation,p_expression,p_methylation Simulation dimensions.
#' @param K_init Over-specified fitting factor count.
#' @param control Control passed to raw-YFB fitters.
#' @return Data frame with one row per benchmark arm and required metrics.
#' @examples
#' # run_multimodal_yfb_replicate("both_informative", 1, 20260918, 40, 40, 20, 20, 4)
#' @family multimodal_yfb_simulation
run_multimodal_yfb_replicate <- function(scenario, replicate_id, seed,
                                         n_train = 120L, n_validation = 120L,
                                         p_expression = 50L, p_methylation = 50L,
                                         K_init = 5L,
                                         control = list()) {
  data <- simulate_multimodal_yfb_data(
    scenario, n_train, n_validation, p_expression, p_methylation, seed = seed
  )
  evaluate <- function(name, fit, prediction, metrics = list()) {
    survival <- multimodal_yfb_survival_metrics(prediction, data$validation$time,
                                                 data$validation$event)
    runtime <- attr(fit, "runtime_seconds")
    out <- data.frame(
      scenario = scenario, replicate_id = replicate_id, seed = seed, method = name,
      censoring_train = data$training$censoring_rate,
      censoring_validation = data$validation$censoring_rate,
      reconstruction_expression = NA_real_, reconstruction_methylation = NA_real_,
      factor_loading_recovery = NA_real_, factor_recovery = NA_real_,
      tau_correlation = NA_real_, tau_rmse = NA_real_, sparsity = NA_real_,
      K_eff = NA_real_, K_eff_reconstruction = NA_real_,
      K_eff_survival = NA_real_, K_eff_retained = NA_real_,
      prognostic_recall = NA_real_, prognostic_direction_recall = NA_real_, converged = NA,
      runtime_seconds = runtime, memory_mb = sum(gc()[, 2]), c_index = survival$c_index,
      heldout_partial_log_likelihood = survival$partial_log_likelihood,
      true_risk_correlation = multimodal_yfb_safe_correlation(
        prediction, data$validation$eta
      )
    )
    for (metric in names(metrics)) out[[metric]] <- metrics[[metric]]
    out
  }

  joint_start <- proc.time()[["elapsed"]]
  joint <- fit_multimodal_yfb(data$training$Y, data$training$time, data$training$event,
                              K_init, control)
  attr(joint, "runtime_seconds") <- proc.time()[["elapsed"]] - joint_start
  joint_prediction <- predict_multimodal_yfb(joint, data$validation$Y)$risk_scores
  joint_metrics <- multimodal_yfb_fit_metrics(joint, data)
  joint_metrics$factor_loading_recovery <- joint_metrics$loading_recovery
  joint_metrics$loading_recovery <- NULL
  reconstructed <- lapply(names(data$training$Y), function(modality) {
    fitted <- joint$EL %*% t(joint$EF[[modality]])
    1 - sum((data$training$Y[[modality]] - fitted)^2) / sum(data$training$Y[[modality]]^2)
  })
  names(reconstructed) <- names(data$training$Y)
  joint_metrics$reconstruction_expression <- reconstructed$expression
  joint_metrics$reconstruction_methylation <- reconstructed$methylation
  rows <- list(evaluate("joint_multimodal_yfb", joint, joint_prediction, joint_metrics))

  ebmf_start <- proc.time()[["elapsed"]]
  ebmf <- fit_multimodal_yfb_ebmf_cox(data$training$Y, data$training$time,
                                       data$training$event, K_init)
  attr(ebmf, "runtime_seconds") <- proc.time()[["elapsed"]] - ebmf_start
  ebmf_prediction <- predict_multimodal_yfb_ebmf_cox(ebmf, data$validation$Y)
  ebmf_match <- multimodal_yfb_match_factors(do.call(rbind, data$truth$F), ebmf$loadings)
  rows[[length(rows) + 1L]] <- evaluate("two_step_ebmf_cox", ebmf, ebmf_prediction,
    list(factor_loading_recovery = mean(ebmf_match$abs_correlation),
         factor_recovery = mean(ebmf_match$abs_correlation >= 0.70), K_eff = ebmf$n_factors))

  for (modality in names(data$training$Y)) {
    single_start <- proc.time()[["elapsed"]]
    single <- fit_multimodal_yfb_single_modality(
      data$training$Y, data$training$time, data$training$event, modality, K_init, control
    )
    attr(single, "runtime_seconds") <- proc.time()[["elapsed"]] - single_start
    prediction <- predict_multimodal_yfb_single_modality(single, data$validation$Y)$risk_scores
    single_pve <- single$history$factor_pve[single$history$n_iter, ]
    rows[[length(rows) + 1L]] <- evaluate(paste0(modality, "_only_yfb"), single, prediction,
      list(K_eff = sum(single_pve >= 0.01), converged = single$history$converged))
  }
  stacked_start <- proc.time()[["elapsed"]]
  stacked <- fit_multimodal_yfb_stacked_risk_cox(
    data$training$Y, data$training$time, data$training$event, K_init
  )
  attr(stacked, "runtime_seconds") <- proc.time()[["elapsed"]] - stacked_start
  stacked_prediction <- predict_multimodal_yfb_stacked_risk_cox(stacked, data$validation$Y)
  rows[[length(rows) + 1L]] <- evaluate("training_only_stacked_risk_cox", stacked,
                                        stacked_prediction)
  do.call(rbind, rows)
}

#' Run the approved ten-replicate Stage 4 smoke study
#'
#' @param output_path CSV path to receive per-arm metrics.
#' @param replicates Number of smoke-study replicates, default 10.
#' @param seed_base Seed offset; scenario/replicate seeds are deterministic.
#' @param ... Passed to [run_multimodal_yfb_replicate()].
#' @return Data frame written to `output_path`.
#' @examples
#' # run_multimodal_yfb_smoke("results/multimodal_yfb_smoke.csv", replicates = 10)
#' @family multimodal_yfb_simulation
run_multimodal_yfb_smoke <- function(output_path, replicates = 10L,
                                     seed_base = 20260917L, ...) {
  rows <- list()
  scenarios <- names(multimodal_yfb_scenarios())
  index <- 1L
  for (scenario in scenarios) for (replicate_id in seq_len(replicates)) {
    rows[[index]] <- run_multimodal_yfb_replicate(
      scenario, replicate_id, seed_base + replicate_id, ...
    )
    index <- index + 1L
  }
  results <- do.call(rbind, rows)
  directory <- dirname(output_path)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  utils::write.csv(results, output_path, row.names = FALSE)
  results
}

#' Run a converged K-init sensitivity analysis for the multimodal smoke design
#'
#' Every K value sees the same scenario/replicate seed, so its comparisons are
#' paired at the simulated-data level. The default controls are read from the
#' production multimodal-YFB configuration rather than the short smoke cap.
#'
#' @param output_path CSV path to receive per-arm metrics.
#' @param K_values Positive integer K-init values to compare.
#' @param replicates Number of independent replicates per scenario and K.
#' @param scenarios Scenario names to include; defaults to all approved cases.
#' @param seed_base Base seed; replicate `r` uses `seed_base + r` for all K.
#' @param n_train,n_validation,p_expression,p_methylation Simulation dimensions.
#' @param control Named overrides to configured production fitting controls.
#' @return Per-arm benchmark data frame with a `K_init` column.
#' @examples
#' # run_multimodal_yfb_k_sensitivity("results/k_sensitivity.csv", c(3, 5, 7))
#' @family multimodal_yfb_simulation
run_multimodal_yfb_k_sensitivity <- function(
    output_path, K_values = 7:15, replicates = 10L,
    scenarios = names(multimodal_yfb_scenarios()), seed_base = 20260917L,
    n_train = 120L, n_validation = 120L, p_expression = 50L,
    p_methylation = 50L, control = list()) {
  if (!is.numeric(K_values) || length(K_values) == 0L ||
      any(!is.finite(K_values)) || any(K_values < 1L) ||
      any(K_values != as.integer(K_values))) {
    stop("K_values must be a nonempty vector of positive integers.")
  }
  if (!is.numeric(replicates) || length(replicates) != 1L || replicates < 1L ||
      replicates != as.integer(replicates)) {
    stop("replicates must be a positive integer.")
  }
  allowed_scenarios <- names(multimodal_yfb_scenarios())
  if (!is.character(scenarios) || length(scenarios) == 0L ||
      any(!scenarios %in% allowed_scenarios)) {
    stop("scenarios must be nonempty and drawn from the approved scenario names.")
  }
  configured <- multimodal_yfb_default_control()
  production_control <- modifyList(
    list(max_outer = configured$max_outer, max_inner = configured$max_inner,
         tol = configured$tol), control
  )
  rows <- list()
  index <- 1L
  for (scenario in scenarios) for (replicate_id in seq_len(as.integer(replicates))) {
    for (K_init in sort(unique(as.integer(K_values)))) {
      result <- run_multimodal_yfb_replicate(
        scenario = scenario, replicate_id = replicate_id,
        seed = as.integer(seed_base + replicate_id), n_train = n_train,
        n_validation = n_validation, p_expression = p_expression,
        p_methylation = p_methylation, K_init = K_init,
        control = production_control
      )
      result$K_init <- K_init
      rows[[index]] <- result
      index <- index + 1L
    }
  }
  results <- do.call(rbind, rows)
  directory <- dirname(output_path)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  utils::write.csv(results, output_path, row.names = FALSE)
  results
}

#' Select a parsimonious K-init from paired multimodal simulation results
#'
#' Selects the smallest K whose mean held-out C-index lies within one standard
#' error of the best mean.  The result is restricted to converged joint-model
#' fits; reconstruction-active counts remain reported diagnostics, not the
#' selection criterion.
#'
#' @param results Per-arm output from [run_multimodal_yfb_k_sensitivity()].
#' @return List containing the selected K and per-K summary data frame.
#' @examples
#' summarize_multimodal_yfb_k_selection(data.frame(
#'   method = "joint_multimodal_yfb", K_init = 7, c_index = 0.6,
#'   heldout_partial_log_likelihood = -1, converged = TRUE))
#' @family multimodal_yfb_simulation
summarize_multimodal_yfb_k_selection <- function(results) {
  required <- c("method", "K_init", "c_index", "heldout_partial_log_likelihood", "converged")
  if (!is.data.frame(results) || !all(required %in% names(results))) {
    stop("results must contain method, K_init, C-index, held-out likelihood, and convergence columns.")
  }
  joint <- results[results$method == "joint_multimodal_yfb" & results$converged, , drop = FALSE]
  if (nrow(joint) == 0L || any(!is.finite(joint$K_init)) || any(!is.finite(joint$c_index))) {
    stop("At least one converged joint multimodal result with finite K and C-index is required.")
  }
  K_values <- sort(unique(as.integer(joint$K_init)))
  summary <- do.call(rbind, lapply(K_values, function(K) {
    rows <- joint[joint$K_init == K, , drop = FALSE]
    data.frame(method = "joint_multimodal_yfb", K_init = K,
               n_replicates = nrow(rows), mean_c_index = mean(rows$c_index),
               se_c_index = if (nrow(rows) > 1L) stats::sd(rows$c_index) / sqrt(nrow(rows)) else 0,
               mean_heldout_partial_log_likelihood = mean(rows$heldout_partial_log_likelihood))
  }))
  best <- max(summary$mean_c_index)
  best_se <- summary$se_c_index[which.max(summary$mean_c_index)][1]
  eligible <- summary$K_init[summary$mean_c_index >= best - best_se]
  list(selected_K = min(eligible), summary = summary)
}

#' Run paired joint-only fits for multimodal K-init selection
#'
#' This diagnostic deliberately avoids the benchmark comparators: all compute is
#' spent on the joint model whose K-init is being selected. Each K receives the
#' identical simulated training/validation draw within a replicate.
#'
#' @param output_path CSV path to receive per-fit results.
#' @param K_values Positive integer K-init values to compare.
#' @param replicates Number of independent paired simulation replicates.
#' @param scenario One approved scenario name.
#' @param seed_base Base seed; replicate `r` uses `seed_base + r` for all K.
#' @param n_train,n_validation,p_expression,p_methylation Simulation dimensions.
#' @param control Named fitter-control overrides.
#' @param verbose Whether to report and checkpoint each completed fit.
#' @return Joint-model results, written to `output_path`.
#' @examples
#' # run_multimodal_yfb_k_selection("results/k_select.csv", 3:7, replicates = 3)
#' @family multimodal_yfb_simulation
run_multimodal_yfb_k_selection <- function(
    output_path, K_values = 3:15, replicates = 10L, scenario = "both_informative",
    seed_base = 20260917L, n_train = 120L, n_validation = 120L,
    p_expression = 50L, p_methylation = 50L, control = list(), verbose = TRUE) {
  if (!is.numeric(K_values) || length(K_values) == 0L || any(!is.finite(K_values)) ||
      any(K_values < 1L) || any(K_values != as.integer(K_values))) {
    stop("K_values must be a nonempty vector of positive integers.")
  }
  if (!is.numeric(replicates) || length(replicates) != 1L || replicates < 1L ||
      replicates != as.integer(replicates) || !is.character(scenario) ||
      length(scenario) != 1L || !scenario %in% names(multimodal_yfb_scenarios())) {
    stop("replicates must be positive and scenario must be an approved name.")
  }
  configured <- multimodal_yfb_default_control()
  production_control <- modifyList(
    list(max_outer = configured$max_outer, max_inner = configured$max_inner,
         tol = configured$tol), control
  )
  directory <- dirname(output_path)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  rows <- list()
  index <- 1L
  for (replicate_id in seq_len(as.integer(replicates))) {
    seed <- as.integer(seed_base + replicate_id)
    data <- simulate_multimodal_yfb_data(
      scenario = scenario, n_train = n_train, n_validation = n_validation,
      p_expression = p_expression, p_methylation = p_methylation, seed = seed
    )
    for (K_init in sort(unique(as.integer(K_values)))) {
      start <- proc.time()[["elapsed"]]
      fit <- fit_multimodal_yfb(data$training$Y, data$training$time,
                                 data$training$event, K = K_init,
                                 control = production_control)
      prediction <- predict_multimodal_yfb(fit, data$validation$Y)$risk_scores
      survival <- multimodal_yfb_survival_metrics(
        prediction, data$validation$time, data$validation$event
      )
      recovery <- multimodal_yfb_fit_metrics(fit, data)
      rows[[index]] <- data.frame(
        scenario = scenario, replicate_id = replicate_id, seed = seed,
        method = "joint_multimodal_yfb", K_init = K_init,
        converged = fit$diagnostics$converged,
        K_eff_reconstruction = fit$diagnostics$K_eff_reconstruction,
        K_eff_survival = fit$diagnostics$K_eff_survival,
        K_eff_retained = fit$diagnostics$K_eff_retained,
        c_index = survival$c_index,
        heldout_partial_log_likelihood = survival$partial_log_likelihood,
        training_partial_log_likelihood = tail(fit$diagnostics$objective, 1),
        true_risk_correlation = multimodal_yfb_safe_correlation(
          prediction, data$validation$eta
        ),
        prognostic_recall = recovery$prognostic_recall,
        prognostic_direction_recall = recovery$prognostic_direction_recall,
        runtime_seconds = proc.time()[["elapsed"]] - start
      )
      utils::write.csv(do.call(rbind, rows), output_path, row.names = FALSE)
      if (isTRUE(verbose)) {
        message(sprintf(
          "Completed replicate %d/%d, K_init=%d: converged=%s, K_rec=%d, K_surv=%d, C-index=%.3f",
          replicate_id, replicates, K_init, fit$diagnostics$converged,
          fit$diagnostics$K_eff_reconstruction, fit$diagnostics$K_eff_survival,
          survival$c_index
        ))
      }
      index <- index + 1L
    }
  }
  results <- do.call(rbind, rows)
  results
}
