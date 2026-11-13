# ============================================================
# Script: run_modality_weight_simulation.R
# Purpose: Run registered balance and joint/exclusive modality-weight simulations.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, ebnm, survival and existing feature-group simulation helpers
# ============================================================
# Sequential worker, one numerical thread. Synthetic fits stay under ignored
# data/. Aggregate checkpoints retain failed fits and every warning. No tuning.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_feature_group_simulation.R")

#' Compute weighted modality precision, final emptiness and repeated reset counts
#' @param fit Fitter result.
#' @param groups Modality labels in feature order.
#' @return Long program-by-modality data frame. Undefined zero-total shares stay NA.
#' @examples
#' # mwsim_precision(fit, labels)
mwsim_precision <- function(fit, groups) {
  weights <- if (is.null(fit$modality_weight)) rep(1, length(groups)) else
    unname(fit$modality_weight[groups])
  do.call(rbind, lapply(seq_len(ncol(fit$EF)), function(k) {
    precision <- vapply(unique(groups), function(m) {
      j <- groups == m
      sum(weights[j] * fit$Tau[j] * fit$EF2[j, k])
    }, numeric(1))
    data.frame(program = k, modality = names(precision), precision = unname(precision),
      share = if (sum(precision) > 0) unname(precision) / sum(precision) else NA_real_,
      empty = vapply(unique(groups), function(m) all(fit$EF2[groups == m, k] == 0), logical(1)),
      point_mass_resets = if (is.null(fit$empty_blocks)) NA_integer_ else fit$empty_blocks[, k])
  }))
}

#' Score independently matched loading directions and all three modality patterns
#' @param fit Returned posterior moments.
#' @param g Feature-by-three standardized, unit-norm reconstruction truth.
#' @param groups Modality labels in feature order.
#' @param z_threshold Posterior coefficient activity threshold.
#' @return One row per truth, with loading/pattern and separate survival recovery.
#' @examples
#' # mwsim_patterns(fit, truth, groups, 1.96)
mwsim_patterns <- function(fit, g, groups, z_threshold) {
  precision <- mwsim_precision(fit, groups)
  do.call(rbind, lapply(1:3, function(t) {
    r <- vapply(seq_len(ncol(fit$EF)), function(k)
      if (sd(fit$EF[, k]) > 0) cor(g[, t], fit$EF[, k]) else NA_real_, numeric(1))
    k <- if (any(is.finite(r))) which.max(replace(abs(r), !is.finite(r), -Inf)) else NA_integer_
    matched <- !is.na(k)
    p <- precision[precision$program == k & !is.na(k), ]
    shares <- stats::setNames(p$share, p$modality)
    empty <- stats::setNames(p$empty, p$modality)
    pattern <- if (!matched) FALSE else if (t < 3L) {
      absent <- if (t == 1L) "methylation" else "expression"
      isTRUE(empty[[absent]]) || isTRUE(shares[[absent]] < 0.05)
    } else all(is.finite(shares)) && all(shares >= 0.05)
    recovered <- matched && abs(r[k]) >= 0.7
    variance <- if (is.null(fit$VBeta)) pmax(fit$EBeta2 - fit$EBeta^2, 0) else fit$VBeta
    active <- matched && variance[k] > 0 && abs(fit$EBeta[k]) >= z_threshold * sqrt(variance[k])
    data.frame(truth_program = t, matched_program = k,
      loading_correlation = if (matched) r[k] else NA_real_,
      loading_recovered = recovered, pattern_correct = pattern,
      expression_share = if (matched) shares[["expression"]] else NA_real_,
      methylation_share = if (matched) shares[["methylation"]] else NA_real_,
      active = active, correct_sign = matched && fit$EBeta[k] * r[k] > 0,
      survival_recovered = t < 3L && recovered && active && fit$EBeta[k] * r[k] > 0)
  }))
}

#' Apply the preregistered gate and tie breakers without external outcomes
#' @param A Block A cell/arm summary including gate and primary recovery counts.
#' @param B Block B summary including full pattern success counts.
#' @return Single-row data frame with effectiveness and selected rule.
#' @examples
#' # mwsim_select(summary_A, summary_B)
mwsim_select <- function(A, B) {
  effective <- stats::setNames(logical(2), c("feature_count", "eigenvalue"))
  arms <- c(feature_count = "grouped_count", eigenvalue = "grouped_eigenvalue")
  for (rule in names(arms)) {
    x <- A[A$arm == arms[[rule]], ]
    baseline <- A[A$arm == "grouped_unweighted", ]
    key <- paste(x$feature_ratio, x$methylation_active_fraction)
    baseline <- baseline[match(key, paste(baseline$feature_ratio, baseline$methylation_active_fraction)), ]
    control <- x$feature_ratio == 1 & x$methylation_active_fraction == .03
    effective[rule] <- sum(control) == 1 && isTRUE(x$passed[control]) &&
      any(x$feature_ratio == 5 & x$passed & !baseline$passed)
  }
  candidates <- names(effective)[effective]
  selected <- "none"
  if (length(candidates) == 1L) selected <- candidates
  if (length(candidates) == 2L) {
    successes <- vapply(arms, function(arm) sum(B$full_pattern_successes[B$arm == arm]), numeric(1))
    recovery <- vapply(arms, function(arm) sum(A$recovered[A$arm == arm & A$feature_ratio == 5]), numeric(1))
    selected <- names(arms)[order(-successes, -recovery, seq_along(arms))[1L]]
  }
  data.frame(count_effective = effective[[1]], eigenvalue_effective = effective[[2]],
    selected_rule = selected, part4_allowed = selected != "none")
}

#' Fit one registered arm from fixed data and common initial moments
#' @param x Prepared input with Y, outcomes and common starts.
#' @param weights NULL or named reconstruction powers.
#' @param pooled Whether to pool F priors across modalities.
#' @param cfg Global configuration.
#' @return Fitter result.
#' @examples
#' # mwsim_fit(x, rules$feature_count, FALSE, cfg)
mwsim_fit <- function(x, weights, pooled, cfg) {
  fit_ssbmf_yfb(x$Y, x$time, x$event, K = cfg$benchmark$k_working_yfb,
    max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol, alpha = cfg$benchmark$alpha,
    prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
    w_S = 1, feature_group = x$groups, modality_weight = weights,
    F_prior_by_group = !pooled, EL_init = x$EL_init, EF_init = x$EF_init, verbose = FALSE)
}

#' Run and checkpoint the complete registered study or a one-seed smoke study
#' @param smoke Use seed 1 only. Its aggregates cannot select a weight rule.
#' @return Invisible A/B summaries and selection, also written to new directories.
#' @examples
#' # mwsim_run(TRUE)
mwsim_run <- function(smoke = FALSE) {
  suffix <- if (smoke) "2026-11-12_weights_smoke" else "2026-11-12_weights"
  out <- paste0("results/benchmark_sim/outputs/feature_group_simulation_", suffix)
  private <- paste0("data/feature_group_simulation_", suffix)
  if (dir.exists(out) || dir.exists(private)) stop("Refusing to overwrite a registered run.")
  dir.create(out, recursive = TRUE); dir.create(private, recursive = TRUE)
  cfg <- yaml::read_yaml("config/globals.yml")
  if (cfg$benchmark$k_working_yfb != 5L) stop("The registered working rank is five.")
  seeds <- if (smoke) 1L else 1:10
  # Freeze every rule before the first scientific fit.
  writeLines(c("registration=docs/prompts/2026-11-12-codex-d-modality-weights-prompt.md",
    paste0("seeds=", paste(seeds, collapse = ",")), "n_train=120; n_validation=120; K=5",
    "A_cells=50:250/15%,50:250/3%,50:50/3%",
    "A_arms=grouped_unweighted,grouped_count,grouped_eigenvalue,expression_only",
    "B_cells=50:50/15%,50:250/15%",
    "B_arms=grouped_count,grouped_eigenvalue,pooled_count,grouped_unweighted",
    "generator=both_informative; existing defaults; exclusive option only for B",
    "B_truth=expression_only1,methylation_only2,joint3; beta=0.70,0.35,0",
    "scaling=per-feature training-only z-score; weights from training rows only",
    "survival=fgsim_prepare programs1:2; seed+12000; preserve original training eta SD",
    "initialization=common joint signed SVD; expression uses corresponding F rows",
    "fit=point_laplace L/F; per-program Normal beta; canonicalize TRUE; wS1; no pruning",
    paste0("max_iter=", cfg$cavi$max_iter, "; tol=", cfg$cavi$tol, "; alpha=", cfg$benchmark$alpha),
    "A_recovery=program1 loading abs-cor>=0.7, posterior activity, correct beta sign",
    "A_false_active=active outside unique abs-cor>=0.7 matches to truths1,2, regardless beta sign",
    "gate=recovery>=8/10; mean false-active<=1; convergence>=9/10; no errors",
    "effective=passes a 50:250 cell unweighted fails AND passes 50:50/3% control",
    "selection=only effective rule; else B grouped full successes, then A 50:250 primary recoveries, then count",
    "B_matching=independent largest abs loading correlation per truth",
    "B_pattern=absent modality empty or precision share<0.05; joint both>=0.05",
    "B_full_success=all three loading recoveries and correct modality patterns; survival separately reported",
    "empty=final all-zero second moments per block; reset events separately, NA when pooled not instrumented",
    "warnings=retain each event; failed fits retained; zero-total precision shares NA",
    "workers=1; numerical_threads=1"), file.path(out, "settings.txt"))
  sources <- c("code/fit_ssbmf_yfb.R", "code/modality_weights_yfb.R", "code/feature_groups_yfb.R",
    "code/simulate_multimodal_yfb.R", "code/multimodal_yfb_updates.R", "config/globals.yml",
    "results/benchmark_sim/run_feature_group_simulation.R",
    "results/benchmark_sim/run_modality_weight_simulation.R")
  write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
    file.path(out, "source_hashes.csv"), row.names = FALSE)
  tables <- list(metrics = NULL, precision = NULL, patterns = NULL, warnings = NULL)
  cells <- rbind(data.frame(block = "A", ratio = c(5L, 5L, 1L), fraction = c(.15, .03, .03)),
                 data.frame(block = "B", ratio = c(1L, 5L), fraction = .15))
  for (cell in seq_len(nrow(cells))) for (seed in seeds) {
    row <- cells[cell, ]; exclusive <- row$block == "B"
    data <- simulate_multimodal_yfb_data("both_informative", p_methylation = 50L * row$ratio,
      seed = seed, methylation_active_fraction = row$fraction,
      modality_exclusive_programs = exclusive)
    x <- fgsim_prepare(data, seed, programs = 1:2)
    arms <- if (exclusive) c("grouped_count", "grouped_eigenvalue", "pooled_count", "grouped_unweighted") else
      c("grouped_unweighted", "grouped_count", "grouped_eigenvalue", "expression_only")
    inputs <- fgsim_arm_inputs(x, cfg$benchmark$k_working_yfb)
    rules <- yf_modality_weight_rules(x$Y, x$groups)
    g3 <- unlist(Map(function(F, s) F[, 3] / s, data$truth$F, x$scaled$sd), use.names = FALSE)
    g3 <- g3 / sqrt(sum(g3^2))
    for (arm in arms) {
      arm_x <- if (arm == "expression_only") inputs$expression_only else inputs$grouped
      weights <- if (grepl("_count$", arm)) rules$feature_count else
        if (arm == "grouped_eigenvalue") rules$eigenvalue else NULL
      base <- data.frame(block = row$block, feature_ratio = row$ratio,
        methylation_active_fraction = row$fraction, seed = seed, arm = arm)
      start <- proc.time()[["elapsed"]]; warnings <- character()
      answer <- tryCatch(withCallingHandlers({
        fit <- mwsim_fit(arm_x, weights, arm == "pooled_count", cfg)
        recovery <- fgsim_recovery_two(fit, arm_x$g, arm_x$groups, cfg$multimodal_yfb$prognostic_z_threshold)
        precision <- mwsim_precision(fit, arm_x$groups)
        precision$weight <- if (is.null(weights)) 1 else unname(weights[precision$modality])
        pattern <- if (exclusive) mwsim_patterns(fit, cbind(x$g, g3), x$groups,
          cfg$multimodal_yfb$prognostic_z_threshold) else NULL
        score <- predict_ssbmf_yfb(arm_x$Y_validation, fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
        saveRDS(list(fit = fit, inputs = arm_x, truth = data$truth), file.path(private,
          sprintf("%s_ratio%d_fraction%s_seed%d_%s.rds", row$block, row$ratio, row$fraction, seed, arm)))
        list(metrics = data.frame(error = "", recovered = recovery$recovered,
          recovered2 = recovery$recovered2, false_active = recovery$false_active,
          converged = isTRUE(fit$history$converged), n_iter = fit$history$n_iter,
          empty_blocks = sum(precision$empty),
          point_mass_resets = if (is.null(fit$empty_blocks)) NA_integer_ else sum(fit$empty_blocks),
          full_pattern_success = if (exclusive) all(pattern$loading_recovered & pattern$pattern_correct) else NA,
          c_validation = frozen_reverse_cindex(score, arm_x$time_validation, arm_x$event_validation)),
          precision = precision, patterns = pattern)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) list(metrics = data.frame(error = conditionMessage(e),
        recovered = FALSE, recovered2 = FALSE, false_active = NA_real_, converged = FALSE,
        n_iter = NA_integer_, empty_blocks = NA_integer_, point_mass_resets = NA_integer_,
        full_pattern_success = if (exclusive) FALSE else NA, c_validation = NA_real_)))
      answer$metrics$elapsed_seconds <- proc.time()[["elapsed"]] - start
      answer$metrics$n_warnings <- length(warnings)
      answer$warnings <- if (length(warnings)) data.frame(warning_number = seq_along(warnings), warning = warnings) else NULL
      for (type in names(tables)) if (!is.null(answer[[type]])) {
        rows <- cbind(base[rep(1, nrow(answer[[type]])), ], answer[[type]])
        tables[[type]] <- rbind(tables[[type]], rows)
        write.csv(tables[[type]], file.path(out, paste0(type, ".csv")), row.names = FALSE)
      }
      cat(sprintf("%s ratio%d fraction%.2f seed%d %s recovered=%s converged=%s error=%s\n",
        row$block, row$ratio, row$fraction, seed, arm, answer$metrics$recovered,
        answer$metrics$converged, answer$metrics$error)); flush.console()
    }
  }
  A <- fgsim_gate(tables$metrics[tables$metrics$block == "A", ])
  bmetrics <- tables$metrics[tables$metrics$block == "B", ]
  B <- do.call(rbind, lapply(split(bmetrics, interaction(bmetrics$feature_ratio, bmetrics$arm)), function(x)
    data.frame(feature_ratio = x$feature_ratio[1], arm = x$arm[1], seeds = nrow(x),
      failed = sum(nzchar(x$error)), full_pattern_successes = sum(x$full_pattern_success),
      converged = sum(x$converged), empty_blocks = sum(x$empty_blocks), mean_false_active = mean(x$false_active))))
  selection <- if (smoke) data.frame(count_effective = FALSE, eigenvalue_effective = FALSE,
    selected_rule = "smoke_only", part4_allowed = FALSE) else mwsim_select(A, B)
  write.csv(A, file.path(out, "summary_A.csv"), row.names = FALSE)
  write.csv(B, file.path(out, "summary_B.csv"), row.names = FALSE)
  write.csv(selection, file.path(out, "selection.csv"), row.names = FALSE)
  print(A); print(B); print(selection)
  invisible(list(A = A, B = B, selection = selection))
}
if (sys.nframe() == 0L) mwsim_run("--smoke" %in% commandArgs(trailingOnly = TRUE))
