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

#' Regenerate survival from the registered standardized reconstruction directions
#' @param data Unchanged generator output, including Weibull settings and raw η.
#' @param seed Simulation seed. Survival uses the fixed offset seed + 12000.
#' @param programs 3L for the historical low-variance target, or 1:2 for the
#'   both_informative directions weighted by their true survival coefficients.
#' @return Scaled matrices, unit direction(s) g, multiplier a, outcomes and training statistics.
#' @examples
#' # x <- fgsim_prepare(simulate_multimodal_yfb_data(), 1)
fgsim_prepare <- function(data, seed, programs = 3L) {
  if (!identical(programs, 3L) && !identical(programs, 1:2))
    stop("Use the registered program3 target or programs1:2.")
  scaled <- fgsim_scale(data$training$Y, data$validation$Y)
  Y <- do.call(cbind, scaled$training)
  Y_validation <- do.call(cbind, scaled$validation)
  if (length(programs) == 1L) {
    # Preserve the earlier low-variance path and its exact arithmetic.
    g <- unlist(Map(function(F, s) F[, 3] / s, data$truth$F, scaled$sd), use.names = FALSE)
    g <- g / sqrt(sum(g^2))
    direction <- g
  } else {
    g <- do.call(rbind, Map(function(F, s) sweep(F[, programs, drop = FALSE], 1, s, "/"),
                           data$truth$F, scaled$sd))
    g <- sweep(g, 2, sqrt(colSums(g^2)), "/")
    colnames(g) <- paste0("program", programs)
    direction <- drop(g %*% data$truth$beta[programs])
  }
  training_center <- colMeans(Y)
  projected <- drop(sweep(Y, 2, training_center) %*% direction)
  a <- stats::sd(data$training$eta) / stats::sd(projected)
  if (!is.finite(a) || a <= 0) stop("Cannot preserve the registered survival signal SD.")
  eta <- projected * a
  eta_validation <- drop(sweep(Y_validation, 2, training_center) %*% direction) * a
  set.seed(seed + 12000L)
  #' Draw Weibull outcomes and calibrate censoring as in the generator
  #' @param risk Numeric Cox predictor, with larger values meaning greater hazard.
  #' @return Time and event vectors.
  #' @examples
  #' # draw(eta)
  draw <- function(risk) {
    et <- (-log(runif(length(risk))) /
      (data$truth$weibull_scale * exp(risk)))^(1 / data$truth$weibull_shape)
    cb <- rexp(length(risk))
    cs <- multimodal_yfb_calibrate_censoring(et, cb, data$truth$target_censoring)
    list(time = pmin(et, cb * cs), event = as.integer(et <= cb * cs))
  }
  tr <- draw(eta); va <- draw(eta_validation)
  list(Y = Y, Y_validation = Y_validation, g = g, a = a, eta = eta,
    eta_validation = eta_validation, time = tr$time, event = tr$event,
    time_validation = va$time, event_validation = va$event, scaled = scaled,
    training_center = training_center,
    groups = rep(names(scaled$training), vapply(scaled$training, ncol, integer(1))))
}

#' Require identical observed data, outcomes and initialization in paired arms
#' @param a,b Lists of prepared matrices, survival outcomes and initial moments.
#' @return Invisible TRUE, or an error naming the mismatched field.
#' @examples
#' # fgsim_assert_inputs(grouped_input, pooled_input)
fgsim_assert_inputs <- function(a, b) {
  for (field in c("Y", "Y_validation", "time", "event", "time_validation",
                  "event_validation", "EL_init", "EF_init")) {
    if (is.null(a[[field]]) || !identical(a[[field]], b[[field]]))
      stop("Paired arm input differs: ", field)
  }
  invisible(TRUE)
}

#' Build the paired joint arms and their expression-only diagnostic
#' @param x Prepared joint matrices, truth and survival outcomes.
#' @param K Number of common signed-SVD starting programs.
#' @return Three input lists. Expression uses matching rows of the joint F start.
#' @examples
#' # inputs <- fgsim_arm_inputs(fgsim_prepare(data, 1), 5)
fgsim_arm_inputs <- function(x, K) {
  s <- svd(x$Y, nu = K, nv = K)
  x$EL_init <- s$u %*% diag(sqrt(s$d[seq_len(K)]), K)
  x$EF_init <- s$v %*% diag(sqrt(s$d[seq_len(K)]), K)
  j <- which(x$groups == "expression")
  expression <- x
  expression$Y <- x$Y[, j, drop = FALSE]
  expression$Y_validation <- x$Y_validation[, j, drop = FALSE]
  expression$EF_init <- x$EF_init[j, , drop = FALSE]
  expression$training_center <- x$training_center[j]
  expression$g <- if (is.matrix(x$g))
    sweep(x$g[j, , drop = FALSE], 2, sqrt(colSums(x$g[j, , drop = FALSE]^2)), "/") else
      x$g[j] / sqrt(sum(x$g[j]^2))
  expression$groups <- x$groups[j]
  expression$scaled <- lapply(x$scaled, function(blocks) blocks["expression"])
  list(grouped = x, pooled = x, expression_only = expression)
}

#' Require the expression arm to use joint survival and matching molecular inputs
#' @param joint Prepared concatenated input with common starting matrices.
#' @param expression Expression-only input.
#' @return Invisible TRUE, or an error identifying a mismatched field.
#' @examples
#' # fgsim_assert_expression_inputs(inputs$grouped, inputs$expression_only)
fgsim_assert_expression_inputs <- function(joint, expression) {
  j <- which(joint$groups == "expression")
  expected <- list(Y = joint$Y[, j, drop = FALSE],
    Y_validation = joint$Y_validation[, j, drop = FALSE],
    EF_init = joint$EF_init[j, , drop = FALSE], EL_init = joint$EL_init,
    training_center = joint$training_center[j], groups = joint$groups[j],
    g = if (is.matrix(joint$g))
      sweep(joint$g[j, , drop = FALSE], 2, sqrt(colSums(joint$g[j, , drop = FALSE]^2)), "/") else
        joint$g[j] / sqrt(sum(joint$g[j]^2)))
  for (field in c("time", "event", "time_validation", "event_validation", "eta", "eta_validation"))
    expected[[field]] <- joint[[field]]
  for (field in names(expected)) {
    if (is.null(expression[[field]]) || !identical(expression[[field]], expected[[field]]))
      stop("Expression arm input differs: ", field)
  }
  invisible(TRUE)
}

#' Evaluate the known risk used to draw validation outcomes without reorientation
#' @param x Prepared input containing eta_validation and validation outcomes.
#' @return Validation concordance, or NA if the true predictor is constant.
#' @examples
#' # fgsim_true_risk(inputs$grouped)
fgsim_true_risk <- function(x) {
  frozen_reverse_cindex(x$eta_validation, x$time_validation, x$event_validation)
}

#' Require a separate F prior per modality and program
#' @param fit Grouped fit object.
#' @param K Registered program count.
#' @return Invisible TRUE or an explicit metadata error.
#' @examples
#' # fgsim_assert_priors(fit, 5)
fgsim_assert_priors <- function(fit, K) {
  if (!identical(names(fit$prior_F), c("expression", "methylation")) ||
      any(lengths(fit$prior_F) != K)) stop("Grouped fit lacks K priors per modality.")
  for (priors in fit$prior_F) for (prior in priors) {
    if (!identical(prior$family, "point_laplace") || is.null(prior$pi) ||
        length(prior$pi) != 1L || !is.finite(prior$pi) || prior$pi < 0 || prior$pi > 1 ||
        (!isTRUE(prior$point_mass) &&
         (is.null(prior$rate) || length(prior$rate) != 1L ||
          !is.finite(prior$rate) || prior$rate <= 0)))
      stop("Grouped fit contains an invalid or missing fitted point-Laplace prior.")
  }
  invisible(TRUE)
}

#' Fit a fixed arm and observe final pooled priors without changing the fitter
#' @param inputs Prepared data and common EL_init/EF_init matrices.
#' @param groups NULL for pooled, otherwise modality labels per feature.
#' @param cfg Parsed configuration with rank and stopping settings.
#' @return A list containing the unchanged fit and its final canonical F priors.
#' @examples
#' # answer <- fgsim_fit(inputs, NULL, cfg)
fgsim_fit <- function(inputs, groups, cfg) {
  # A child environment confines the observer to this call. It delegates all
  # arithmetic and captures the priors after their canonical scale adjustment.
  observer <- new.env(parent = environment(fit_ssbmf_yfb))
  original <- canonicalize_yf_factors
  last_priors <- NULL
  #' Observe canonicalization while preserving its returned value
  #' @param ... Unchanged canonicalization arguments.
  #' @return The original canonicalization result.
  #' @examples
  #' # Used only by the local fitter binding.
  observer$canonicalize_yf_factors <- function(...) {
    answer <- original(...)
    last_priors <<- answer$prior_F
    answer
  }
  fitter <- fit_ssbmf_yfb
  environment(fitter) <- observer
  fit <- fitter(inputs$Y, inputs$time, inputs$event,
    K = cfg$benchmark$k_working_yfb, max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
    prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
    w_S = 1, alpha = cfg$benchmark$alpha, feature_group = groups,
    EL_init = inputs$EL_init, EF_init = inputs$EF_init, verbose = FALSE)
  if (!is.null(groups)) fgsim_assert_priors(fit, cfg$benchmark$k_working_yfb)
  list(fit = fit, priors = if (is.null(groups)) last_priors else fit$prior_F)
}

#' Score a shared prognostic direction with the registered sign and activity rule
#' @param fit Returned posterior moments, using VBeta when available.
#' @param g Unit standardized reconstruction truth.
#' @param groups Modality label for each feature.
#' @param z_threshold Configured coefficient activity threshold.
#' @return Recovery, matched program, false activity, program and block tables.
#' @examples
#' # fgsim_recovery(fit, arm_x$g, arm_x$groups, 1.96)
fgsim_recovery <- function(fit, g, groups, z_threshold) {
  correlation <- vapply(seq_len(ncol(fit$EF)), function(k) {
    if (stats::sd(fit$EF[, k]) > 0) stats::cor(g, fit$EF[, k]) else NA_real_
  }, numeric(1))
  variance <- if (!is.null(fit$VBeta)) fit$VBeta else pmax(fit$EBeta2 - fit$EBeta^2, 0)
  beta_sd <- sqrt(variance)
  active <- abs(fit$EBeta) >= z_threshold * beta_sd & beta_sd > 0
  matched <- if (any(is.finite(correlation)))
    which.max(replace(abs(correlation), !is.finite(correlation), -Inf)) else NA_integer_
  recovered <- !is.na(matched) && abs(correlation[matched]) >= 0.7 &&
    active[matched] && fit$EBeta[matched] * correlation[matched] > 0
  blocks <- do.call(rbind, lapply(unique(groups), function(m) {
    j <- groups == m
    r <- if (!is.na(matched) && stats::sd(fit$EF[j, matched]) > 0 && stats::sd(g[j]) > 0)
      stats::cor(g[j], fit$EF[j, matched]) else NA_real_
    data.frame(modality = m, correlation = r)
  }))
  list(recovered = recovered, matched_program = matched,
    loading_correlation = if (is.na(matched)) NA_real_ else abs(correlation[matched]),
    false_active = as.integer(sum(active) - as.integer(recovered)), blocks = blocks,
    programs = data.frame(program = seq_along(active), beta = fit$EBeta,
      beta_sd = beta_sd, active = active, loading_correlation = correlation))
}

#' Score both informative truth directions with independent loading matches
#' @param fit Fitter posterior moments.
#' @param g Feature-by-two unit truth directions.
#' @param groups Modality label per fitted feature.
#' @param z_threshold Posterior coefficient activity threshold.
#' @return Primary and secondary recovery, threshold-qualified false activity,
#'   matching diagnostics, program moments and block correlations.
#' @examples
#' # fgsim_recovery_two(fit, inputs$g, inputs$groups, 1.96)
fgsim_recovery_two <- function(fit, g, groups, z_threshold) {
  if (!is.matrix(g) || ncol(g) != 2L) stop("Informative scoring needs two truth directions.")
  first <- fgsim_recovery(fit, g[, 1], groups, z_threshold)
  second <- fgsim_recovery(fit, g[, 2], groups, z_threshold)
  matches <- integer()
  if (isTRUE(first$loading_correlation >= 0.7)) matches <- c(matches, first$matched_program)
  if (isTRUE(second$loading_correlation >= 0.7)) matches <- c(matches, second$matched_program)
  # Andrew's 2026-11-12 clarification: S-9 loading matches are exempt even
  # with a wrong coefficient sign. Wrong signs still fail recovery. Reused
  # fitted columns count once, so independent matching cannot exempt twice.
  first$false_active <- as.integer(sum(first$programs$active[
    setdiff(seq_len(nrow(first$programs)), unique(matches))]))
  first$recovered2 <- second$recovered
  first$matched_program2 <- second$matched_program
  first$loading_correlation2 <- second$loading_correlation
  first$match_collision <- !is.na(first$matched_program) && !is.na(second$matched_program) &&
    first$matched_program == second$matched_program
  first$programs$loading_correlation2 <- second$programs$loading_correlation
  first$blocks$truth_program <- 1L
  second$blocks$truth_program <- 2L
  first$blocks <- rbind(first$blocks, second$blocks)
  first
}

#' Apply the registered informative-scenario separation rule cell by cell
#' @param summary Gate table containing the three arms in every cell.
#' @return Expression-gate eligibility, joint-arm gates and descriptive conclusion.
#' @examples
#' # fgsim_separation(fgsim_gate(metrics))
fgsim_separation <- function(summary) {
  cells <- unique(summary[c("feature_ratio", "methylation_active_fraction")])
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    x <- summary[summary$feature_ratio == cells$feature_ratio[i] &
      summary$methylation_active_fraction == cells$methylation_active_fraction[i], ]
    if (nrow(x) != 3L || !setequal(x$arm, c("grouped", "pooled", "expression_only")))
      stop("Separation rule requires all three arms exactly once per cell.")
    expression <- x$passed[x$arm == "expression_only"]
    grouped <- x$passed[x$arm == "grouped"]
    pooled <- x$passed[x$arm == "pooled"]
    conclusion <- if (!expression) "not testable" else if (grouped && !pooled)
      "separate priors helped" else if (pooled && !grouped) "separate priors hurt" else
        if (grouped) "both priors pass" else "both priors fail"
    cbind(cells[i, ], data.frame(testable = expression, grouped_passed = grouped,
      pooled_passed = pooled, conclusion = conclusion))
  }))
}

#' Evaluate the registered gate in each cell and arm, retaining failed seeds
#' @param metrics Per-fit aggregate rows.
#' @return Cell/arm recovery, false activity, convergence, defined C count and gate.
#' @examples
#' # fgsim_gate(metrics)
fgsim_gate <- function(metrics) {
  cells <- split(metrics, interaction(metrics$feature_ratio,
    metrics$methylation_active_fraction, metrics$arm, drop = TRUE))
  do.call(rbind, lapply(cells, function(x) {
    false <- mean(x$false_active)
    data.frame(x[1, c("feature_ratio", "methylation_active_fraction", "arm")],
      seeds = nrow(x), failed = sum(nzchar(x$error)), recovered = sum(x$recovered),
      mean_false_active = false, converged = sum(x$converged),
      total_empty_blocks = sum(x$empty_blocks, na.rm = TRUE),
      empty_blocks_defined = sum(!is.na(x$empty_blocks)),
      true_c_defined = sum(is.finite(x$c_true_risk)),
      mean_true_c = if (any(is.finite(x$c_true_risk))) mean(x$c_true_risk[is.finite(x$c_true_risk)]) else NA_real_,
      c_defined = sum(is.finite(x$c_validation)),
      mean_c = if (any(is.finite(x$c_validation))) mean(x$c_validation[is.finite(x$c_validation)]) else NA_real_,
      passed = nrow(x) == 10L && !any(nzchar(x$error)) && sum(x$recovered) >= 8L &&
        is.finite(false) && false <= 1 && sum(x$converged) >= 9L)
  }))
}

#' Apply the rerun prior preference and selected-arm gate rule
#' @param summary Cell/arm gate table.
#' @return Selected prior, step4_allowed and a reason for the decision.
#' @examples
#' # fgsim_select(fgsim_gate(metrics))
fgsim_select <- function(summary) {
  grouped <- summary[summary$arm == "grouped", ]
  pooled <- summary[summary$arm == "pooled", ]
  #' Form a cell key from the registered design levels
  #' @param x Cell summary rows.
  #' @return Character keys.
  #' @examples
  #' # key(grouped)
  key <- function(x) paste(x$feature_ratio, x$methylation_active_fraction)
  pooled <- pooled[match(key(grouped), key(pooled)), ]
  if (anyNA(pooled$passed)) stop("Missing paired arm gate.")
  selected <- if (any(pooled$passed & !grouped$passed)) "pooled" else "grouped"
  allowed <- any(summary$passed[summary$arm == selected])
  list(selected_prior = selected, step4_allowed = allowed,
    reason = if (!allowed) "Selected prior passes no cell: stop before step 4."
      else if (selected == "pooled") "Pooled passes a cell where grouped fails."
      else "Grouped passes a cell; no cell passes pooled while failing grouped.")
}

#' Fit three registered arms on one shared seed/cell and save private synthetic fits
#' @param seed Simulation seed, 1 through 10.
#' @param ratio Methylation/expression feature ratio, 1 or 5.
#' @param fraction Methylation nonzero loading fraction, 0.15 or 0.03.
#' @param cfg Parsed global configuration.
#' @param scenario Registered low-variance or both_informative generator scenario.
#' @param private Ignored checkpoint directory.
#' @return Aggregate metrics, programs, blocks, precision, priors and warnings.
#' @examples
#' # fgsim_one(1, 1, 0.15, cfg, "data/feature_group_simulation_smoke")
fgsim_one <- function(seed, ratio, fraction, cfg, private,
                      scenario = "low_variance_prognostic_factor") {
  informative <- identical(scenario, "both_informative")
  p <- eval(formals(simulate_multimodal_yfb_data)$p_expression)
  data <- simulate_multimodal_yfb_data(scenario,
    p_expression = p, p_methylation = p * ratio, seed = seed,
    methylation_active_fraction = fraction)
  x <- fgsim_prepare(data, seed, programs = if (informative) 1:2 else 3L)
  K <- cfg$benchmark$k_working_yfb
  inputs <- fgsim_arm_inputs(x, K)
  fgsim_assert_inputs(inputs$grouped, inputs$pooled)
  fgsim_assert_expression_inputs(inputs$grouped, inputs$expression_only)
  c_true <- fgsim_true_risk(x)
  results <- list()
  for (arm in names(inputs)) {
    arm_x <- inputs[[arm]]
    warnings <- character()
    start <- proc.time()[["elapsed"]]
    base <- data.frame(seed = seed, scenario = scenario, feature_ratio = ratio,
      methylation_active_fraction = fraction, arm = arm,
      p_expression = p, p_methylation = p * ratio, K = K)
    result <- tryCatch(withCallingHandlers({
      answer <- fgsim_fit(arm_x, if (arm == "grouped") arm_x$groups else NULL, cfg)
      fit <- answer$fit
      recovery <- if (informative) fgsim_recovery_two(fit, arm_x$g, arm_x$groups,
        cfg$multimodal_yfb$prognostic_z_threshold) else
          fgsim_recovery(fit, arm_x$g, arm_x$groups, cfg$multimodal_yfb$prognostic_z_threshold)
      precision <- fgsim_precision(fit, arm_x$groups)
      precision$feature_share <- vapply(precision$modality, function(m) mean(arm_x$groups == m), numeric(1))
      precision$live <- colSums(fit$EF^2)[precision$program] > 0
      precision$matched <- precision$program == recovery$matched_program
      precision$dominant <- precision$share > precision$feature_share
      if (informative) precision$matched2 <- precision$program == recovery$matched_program2
      prior_rows <- do.call(rbind, lapply(seq_len(K), function(k) {
        do.call(rbind, lapply(unique(arm_x$groups), function(m) {
          prior <- if (arm == "grouped") answer$priors[[m]][[k]] else answer$priors[[k]]
          if (is.null(prior$pi)) stop("Missing final fitted F prior.")
          data.frame(program = k, modality = m, prior_scope = if (arm == "grouped") m else "pooled",
            point_mass_weight = 1 - prior$pi,
            laplace_scale = if (isTRUE(prior$point_mass)) NA_real_ else 1 / prior$rate,
            true_nonzero_fraction = mean(data$truth$F[[m]][, 3] != 0))
        }))
      }))
      score <- predict_ssbmf_yfb(arm_x$Y_validation, fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
      c_index <- frozen_reverse_cindex(score, x$time_validation, x$event_validation)
      # Private checkpoints retain the common input and the returned moments.
      saveRDS(list(fit = fit, priors = answer$priors, data = data, inputs = arm_x),
        file.path(private, sprintf("ratio%d_fraction%s_seed%d_%s.rds", ratio, fraction, seed, arm)))
      list(metrics = cbind(base, data.frame(error = "", recovered = recovery$recovered,
        matched_program = recovery$matched_program, loading_correlation = recovery$loading_correlation,
        expression_loading_correlation = recovery$blocks$correlation[recovery$blocks$modality == "expression" &
          if (informative) recovery$blocks$truth_program == 1L else TRUE],
        empty_blocks = if (is.null(fit$empty_blocks)) 0L else sum(fit$empty_blocks),
        false_active = recovery$false_active, converged = isTRUE(fit$history$converged),
        n_iter = fit$history$n_iter, c_validation = c_index,
        elbo = tail(fit$history$elbo_full, 1))),
        programs = recovery$programs, blocks = recovery$blocks,
        empty_blocks = if (!is.null(fit$empty_blocks))
          data.frame(modality = rep(rownames(fit$empty_blocks), times = K),
            program = rep(seq_len(K), each = nrow(fit$empty_blocks)),
            events = as.vector(fit$empty_blocks)) else NULL,
        precision = precision, priors = prior_rows)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      # Keep each warning, even if repeated. Muffle only R's duplicate display.
      invokeRestart("muffleWarning")
    }), error = function(e) list(metrics = cbind(base, data.frame(error = conditionMessage(e),
      recovered = FALSE, matched_program = NA_integer_, loading_correlation = NA_real_,
      expression_loading_correlation = NA_real_, empty_blocks = NA_integer_,
      false_active = NA_integer_, converged = FALSE, n_iter = NA_integer_,
      c_validation = NA_real_, elbo = NA_real_))))
    for (type in c("programs", "blocks", "precision", "priors", "empty_blocks")) {
      if (!is.null(result[[type]])) result[[type]] <- cbind(base[rep(1, nrow(result[[type]])), ], result[[type]])
    }
    if (informative) {
      failed <- nzchar(result$metrics$error)
      result$metrics$recovered2 <- if (failed) FALSE else recovery$recovered2
      result$metrics$matched_program2 <- if (failed) NA_integer_ else recovery$matched_program2
      result$metrics$loading_correlation2 <- if (failed) NA_real_ else recovery$loading_correlation2
      result$metrics$match_collision <- if (failed) NA else recovery$match_collision
    }
    result$metrics$c_true_risk <- c_true
    result$metrics$elapsed_seconds <- proc.time()[["elapsed"]] - start
    result$metrics$n_warnings <- length(warnings)
    result$warnings <- if (length(warnings)) cbind(base[rep(1, length(warnings)), ],
      data.frame(warning_number = seq_along(warnings), warning = warnings)) else NULL
    if (length(warnings)) cat(sprintf("%s: %d warnings retained individually.\n", arm, length(warnings)))
    cat(sprintf("ratio 1:%d methylation %.2f seed %d %s: recovered=%s converged=%s error=%s\n",
      ratio, fraction, seed, arm, result$metrics$recovered, result$metrics$converged, result$metrics$error))
    flush.console()
    results[[arm]] <- result
  }
  # Recheck after all three calls so accidental input mutation also fails loudly.
  fgsim_assert_inputs(inputs$grouped, inputs$pooled)
  fgsim_assert_expression_inputs(inputs$grouped, inputs$expression_only)
  lapply(c("metrics", "programs", "blocks", "precision", "priors", "empty_blocks", "warnings"), function(type)
    do.call(rbind, lapply(results, `[[`, type))) |>
    stats::setNames(c("metrics", "programs", "blocks", "precision", "priors", "empty_blocks", "warnings"))
}

#' Run the registered 120 fits, or a twelve-fit smoke check, sequentially
#' @param out Aggregate output directory, which must not exist.
#' @param private Ignored checkpoint directory, which must not exist.
#' @param seeds 1:10 for the study, or 1 for the smoke check.
#' @param scenario Registered generator scenario.
#' @return Invisible gate summary and selection. Smoke results cannot pass the gate.
#' @examples
#' # fgsim_run("results/benchmark_sim/outputs/feature_group_simulation_smoke",
#' #           "data/feature_group_simulation_smoke", seeds = 1L)
fgsim_run <- function(out, private, seeds = 1:10,
                      scenario = "low_variance_prognostic_factor") {
  if (!scenario %in% c("low_variance_prognostic_factor", "both_informative"))
    stop("Scenario is outside the registered studies.")
  informative <- identical(scenario, "both_informative")
  if (dir.exists(out) || dir.exists(private)) stop("Refusing to overwrite an existing study.")
  if (!identical(seeds, 1:10) && !identical(seeds, 1L)) stop("Use all registered seeds or smoke seed 1.")
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
  writeLines(c(paste0("registration=docs/prompts/", if (informative)
    "2026-11-12-codex-d-next-prompt.md" else "2026-11-12-codex-d-rerun-prompt.md"),
    paste0("seeds=", paste(seeds, collapse = ",")), paste0("scenario=", scenario),
    "feature_ratios=50:50,50:250", "methylation_active_fraction=0.15,0.03",
    "expression_active_fraction=generator default 0.15", "arms=grouped,pooled,expression_only",
    "generator=defaults (120 training, 120 validation, noise_scale=1, truncate=TRUE)",
    "scaling=training-only z-score within modality",
    if (informative) "truth=g_k=unit stacked D^-1 F[,k], k=1,2; beta=(0.70,0.35) from generator" else
      "truth=g=unit stacked D^-1 F[,3]",
    if (informative) "eta=a * centered standardized X * (g1 beta1 + g2 beta2); a preserves original training eta SD" else
      "eta=centered standardized X times g times a; a preserves original training eta SD",
    "survival_seed_offset=12000", "initialization=joint signed SVD; expression-only uses matching F rows and same L",
    "true_risk_C=frozen validation C of the joint eta used to draw all outcomes",
    "expression_truth=each expression part of joint g renormalized to unit norm",
    "empty_blocks=all refitted point-mass block resets, including repeated empty resets",
    "prior_LF=point_laplace", "prior_beta=normal per-program", "canonicalize=TRUE",
    "w_S=1", "prune=FALSE", "strata_id=NULL", "center_projection=FALSE", "workers=1",
    paste0("K=", cfg$benchmark$k_working_yfb), paste0("max_iter=", cfg$cavi$max_iter),
    paste0("tol=", cfg$cavi$tol), paste0("alpha=", cfg$benchmark$alpha),
    paste0("activity_z=", cfg$multimodal_yfb$prognostic_z_threshold),
    "recovery=abs(cor(g,F_k))>=0.7 and active and beta_k*cor(g,F_k)>0",
    if (informative) "gate_target=program1; program2 recovery descriptive; independent best loading matches" else
      "gate_target=program3",
    if (informative) "false_active=active programs outside unique abs-loading-correlation>=0.7 matches to truths1,2 regardless coefficient sign (Andrew confirmed)" else
      "false_active=active programs other than a recovered match",
    "gate=recovery>=8/10; mean false_active<=1; converged>=9/10; no failed fits",
    if (informative) "separation=compare joint arms only in cells where expression-only passes the full gate; PartA stays descriptive" else
      "selection=grouped unless pooled passes a cell grouped fails; selected arm must pass a cell",
    "prior_pi=slab probability; reported point_mass_weight=1-pi",
    "objectives=historical pooled incoming-prior KL vs grouped returned-prior KL; not comparable"),
    file.path(out, "settings.txt"))
  results <- list()
  for (ratio in c(1L, 5L)) for (fraction in c(0.15, 0.03)) for (seed in seeds) {
    results[[length(results) + 1L]] <- fgsim_one(seed, ratio, fraction, cfg, private, scenario)
    for (type in c("metrics", "programs", "blocks", "precision", "priors", "empty_blocks", "warnings")) {
      rows <- do.call(rbind, lapply(results, `[[`, type))
      if (!is.null(rows)) write.csv(rows, file.path(out, paste0(type, ".csv")), row.names = FALSE)
    }
  }
  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  summary <- fgsim_gate(metrics)
  if (informative) {
    summary$recovered2 <- vapply(seq_len(nrow(summary)), function(i) {
      x <- metrics[metrics$feature_ratio == summary$feature_ratio[i] &
        metrics$methylation_active_fraction == summary$methylation_active_fraction[i] &
        metrics$arm == summary$arm[i], ]
      as.integer(sum(x$recovered2))
    }, integer(1))
  }
  selection <- if (informative) fgsim_separation(summary) else fgsim_select(summary)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  write.csv(as.data.frame(selection), file.path(out, if (informative) "separation.csv" else "selection.csv"), row.names = FALSE)
  print(summary); print(selection)
  invisible(list(summary = summary, selection = selection))
}

if (sys.nframe() == 0L) {
  smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)
  informative <- "--informative" %in% commandArgs(trailingOnly = TRUE)
  suffix <- if (informative) {
    if (smoke) "2026-11-12_informative_smoke" else "2026-11-12_informative"
  } else if (smoke) "2026-11-12_rerun_smoke" else "2026-11-12_rerun"
  fgsim_run(paste0("results/benchmark_sim/outputs/feature_group_simulation_", suffix),
    paste0("data/feature_group_simulation_", suffix), seeds = if (smoke) 1L else 1:10,
    scenario = if (informative) "both_informative" else "low_variance_prognostic_factor")
}
