# ============================================================
# Script: simulate_multimodal_yfb.R
# Purpose: Generate matched multimodal raw-YFB training and validation data.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R; survival; code/multimodal_yfb_helpers.R
# ============================================================

if (!exists("multimodal_yfb_cox_working")) {
  source("code/multimodal_yfb_helpers.R")
}

#' Calibrate independent censoring to a target empirical rate
#'
#' @param event_time Simulated event times.
#' @param censor_base Positive base censoring draws.
#' @param target_censoring Target probability of censoring.
#' @param iterations Number of bisection iterations.
#' @return Positive multiplier for `censor_base`.
#' @examples
#' multimodal_yfb_calibrate_censoring(c(1, 2), c(1, 1), 0.5)
#' @family multimodal_yfb_simulation
multimodal_yfb_calibrate_censoring <- function(event_time, censor_base,
                                               target_censoring = 0.30,
                                               iterations = 50L) {
  if (length(event_time) == 0L || length(censor_base) != length(event_time) ||
      any(!is.finite(event_time)) || any(event_time <= 0) ||
      any(!is.finite(censor_base)) || any(censor_base <= 0) ||
      length(target_censoring) != 1L || target_censoring <= 0 ||
      target_censoring >= 1) {
    stop("Event/censoring times must be positive and target_censoring must be in (0, 1).")
  }
  lower <- 1e-8
  upper <- max(event_time / censor_base) * 100
  for (iteration in seq_len(iterations)) {
    middle <- sqrt(lower * upper)
    if (mean(censor_base * middle < event_time) > target_censoring) {
      lower <- middle
    } else {
      upper <- middle
    }
  }
  upper
}

#' Return the seven approved multimodal simulation scenarios
#'
#' @return Named list of scenario-specific parameter overrides.
#' @examples
#' names(multimodal_yfb_scenarios())
#' @family multimodal_yfb_simulation
multimodal_yfb_scenarios <- function() {
  list(
    both_informative = list(active_modalities = c(TRUE, TRUE), beta = c(0.70, 0.35, 0)),
    adverse_protective = list(active_modalities = c(TRUE, TRUE), beta = c(0.70, -0.70, 0)),
    expression_only = list(active_modalities = c(TRUE, FALSE), beta = c(0.70, 0.35, 0)),
    methylation_only = list(active_modalities = c(FALSE, TRUE), beta = c(0.70, 0.35, 0)),
    weak_methylation = list(active_modalities = c(TRUE, TRUE), methylation_scale = 0.25,
                            beta = c(0.70, 0.35, 0)),
    unequal_methylation_features = list(active_modalities = c(TRUE, TRUE), p_methylation = 75L,
                                        beta = c(0.70, 0.35, 0)),
    low_variance_prognostic_factor = list(active_modalities = c(TRUE, TRUE),
                                           score_scale = c(1, 1, 0.15),
                                           beta = c(0, 0, 1.10)),
    null_survival = list(active_modalities = c(TRUE, TRUE), beta = c(0, 0, 0))
  )
}

#' Generate sparse nonnegative modality-specific loading matrices
#'
#' @param p Number of features.
#' @param K Number of factors.
#' @param active_fraction Fraction of nonzero loadings per factor.
#' @param scale Overall loading scale.
#' @return `p` by `K` nonnegative matrix with unit-norm columns before scaling.
#' @examples
#' simulate_multimodal_yfb_loadings(10, 2)
#' @family multimodal_yfb_simulation
simulate_multimodal_yfb_loadings <- function(p, K, active_fraction = 0.15,
                                             scale = 1) {
  F <- matrix(0, p, K)
  for (k in seq_len(K)) {
    active <- sample.int(p, max(1L, ceiling(p * active_fraction)))
    F[active, k] <- rexp(length(active))
    F[, k] <- F[, k] / sqrt(sum(F[, k]^2)) * scale
  }
  F
}

#' Simulate one matched cohort from fixed multimodal YFB parameters
#'
#' @param parameters True loadings, precisions, survival coefficients, and
#'   Weibull/censoring settings created by [simulate_multimodal_yfb_data()].
#' @param n Number of subjects.
#' @param id_prefix Prefix used for unique subject identifiers.
#' @return Named molecular blocks, survival outcomes, true scores, and risk.
#' @examples
#' # Used internally by simulate_multimodal_yfb_data().
#' @family multimodal_yfb_simulation
simulate_multimodal_yfb_cohort <- function(parameters, n, id_prefix) {
  K <- length(parameters$beta)
  L <- sweep(matrix(rexp(n * K), n, K), 2, parameters$score_scale, "*")
  Y <- lapply(names(parameters$F), function(modality) {
    signal <- L %*% t(parameters$F[[modality]])
    noise <- sweep(matrix(rnorm(n * nrow(parameters$F[[modality]])), n), 2,
                   sqrt(parameters$Tau[[modality]]), "/") * (parameters$noise_scale %||% 1)
    # The production fitter requires raw inputs >= 0. Truncation is explicit
    # here; Tau remains the pre-truncation Gaussian precision for recovery.
    if (isFALSE(parameters$truncate)) signal + noise else pmax(signal + noise, 0)
  })
  names(Y) <- names(parameters$F)
  ids <- paste0(id_prefix, seq_len(n))
  Y <- Map(function(block, feature_names) {
    dimnames(block) <- list(ids, feature_names)
    block
  }, Y, parameters$feature_names)
  eta <- Reduce(`+`, Map(`%*%`, Y, parameters$F)) %*% parameters$beta
  event_time <- (-log(stats::runif(n)) /
    (parameters$weibull_scale * exp(as.vector(eta))))^(1 / parameters$weibull_shape)
  censor_base <- stats::rexp(n)
  censor_scale <- multimodal_yfb_calibrate_censoring(
    event_time, censor_base, parameters$target_censoring
  )
  censor_time <- censor_base * censor_scale
  list(Y = Y, time = stats::setNames(pmin(event_time, censor_time), ids),
       event = stats::setNames(as.integer(event_time <= censor_time), ids),
       L = L, eta = as.vector(eta), censoring_rate = mean(event_time > censor_time))
}

#' Simulate independent matched training and validation multimodal YFB data
#'
#' Survival follows Weibull proportional hazards using the same raw joint
#' predictor fitted by the proposed model, `η = (Y_expression F_expression +
#' Y_methylation F_methylation)β`.
#'
#' @param scenario One name from [multimodal_yfb_scenarios()].
#' @param n_train,n_validation Independent training and validation sample sizes.
#' @param p_expression,p_methylation Feature counts before scenario overrides.
#' @param K_true Number of true factors.
#' @param target_censoring Target censoring fraction in each cohort.
#' @param seed RNG seed.
#' @param noise_scale Multiplier on every feature's noise SD (default 1, the
#'   original design, where each true factor explains < 1% of the variance).
#'   0.25 gives a moderate signal-to-noise ratio (~25% variance explained).
#' @param truncate If TRUE (default), truncate the molecular data at zero so it
#'   meets the point-exponential (nonnegative) input requirement. FALSE keeps
#'   the Gaussian data as generated, for signed loading priors and centering.
#' @return Training/validation data plus true model parameters and scenario name.
#' @examples
#' simulate_multimodal_yfb_data("both_informative", n_train = 20, n_validation = 10)
#' @family multimodal_yfb_simulation
simulate_multimodal_yfb_data <- function(scenario = "both_informative",
                                         n_train = 120L, n_validation = 120L,
                                         p_expression = 50L, p_methylation = 50L,
                                         K_true = 3L, target_censoring = 0.30,
                                         seed = 1L, truncate = TRUE, noise_scale = 1) {
  scenarios <- multimodal_yfb_scenarios()
  if (!scenario %in% names(scenarios)) {
    stop("scenario must be one of: ", paste(names(scenarios), collapse = ", "), ".")
  }
  if (K_true != 3L) stop("The current approved scenario definitions require K_true = 3.")
  if (n_train < 2L || n_validation < 2L || p_expression < 2L || p_methylation < 2L) {
    stop("Training/validation sample sizes and feature counts must be at least two.")
  }
  setting <- scenarios[[scenario]]
  p_methylation <- setting$p_methylation %||% p_methylation
  score_scale <- setting$score_scale %||% rep(1, K_true)
  methylation_scale <- setting$methylation_scale %||% 1
  set.seed(seed)
  F <- list(
    expression = simulate_multimodal_yfb_loadings(p_expression, K_true),
    methylation = simulate_multimodal_yfb_loadings(p_methylation, K_true,
                                                    scale = methylation_scale)
  )
  for (m in seq_along(F)) {
    if (!setting$active_modalities[m]) F[[m]][, which(setting$beta != 0)] <- 0
  }
  feature_names <- list(expression = paste0("g", seq_len(p_expression)),
                        methylation = paste0("cg", seq_len(p_methylation)))
  parameters <- list(
    F = F, Tau = lapply(F, function(x) stats::rgamma(nrow(x), shape = 2, rate = 2)),
    beta = setting$beta, score_scale = score_scale, feature_names = feature_names,
    weibull_shape = 1.5, weibull_scale = 0.01, target_censoring = target_censoring,
    truncate = truncate, noise_scale = noise_scale
  )
  training <- simulate_multimodal_yfb_cohort(parameters, n_train, "train_")
  validation <- simulate_multimodal_yfb_cohort(parameters, n_validation, "validation_")
  list(training = training, validation = validation, truth = parameters,
       scenario = scenario, seed = seed)
}

#' Greedily match estimated factors to true factors by absolute correlation
#'
#' @param truth True feature-by-factor matrix.
#' @param estimate Estimated feature-by-factor matrix.
#' @return Data frame with one non-reused estimated match per true factor.
#' @examples
#' multimodal_yfb_match_factors(diag(2), diag(2))
#' @family multimodal_yfb_simulation
multimodal_yfb_match_factors <- function(truth, estimate) {
  correlation <- stats::cor(truth, estimate)
  score <- abs(correlation)
  matches <- vector("list", min(ncol(truth), ncol(estimate)))
  for (i in seq_along(matches)) {
    selected <- which(score == max(score), arr.ind = TRUE)[1, ]
    matches[[i]] <- data.frame(true_factor = selected[1], estimated_factor = selected[2],
                               correlation = correlation[selected[1], selected[2]],
                               abs_correlation = score[selected[1], selected[2]])
    score[selected[1], ] <- -Inf
    score[, selected[2]] <- -Inf
  }
  do.call(rbind, matches)
}

#' Compute held-out survival metrics without validation reorientation
#'
#' @param risk Frozen validation risk score.
#' @param time,event Validation outcomes.
#' @return C-index and Breslow partial log likelihood.
#' @examples
#' multimodal_yfb_survival_metrics(c(0, 1), c(1, 2), c(1, 0))
#' @family multimodal_yfb_simulation
multimodal_yfb_survival_metrics <- function(risk, time, event) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival is required.")
  list(
    c_index = as.numeric(survival::concordance(
      survival::Surv(time, event) ~ risk, reverse = TRUE
    )$concordance),
    partial_log_likelihood = multimodal_yfb_cox_working(time, event, risk)$log_partial_likelihood,
    true_risk_correlation = NA_real_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
