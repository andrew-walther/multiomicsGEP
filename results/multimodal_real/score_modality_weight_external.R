# ============================================================
# Script: score_modality_weight_external.R
# Purpose: Score fixed Task C/D TCGA fits once on primary matched ICGC donors.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, survival, impute and existing YFB/data/concordance helpers
# ============================================================
# One CPU worker, one numerical thread. No fitting, tuning or orientation change.
# Only aggregates are written outside ignored data/. The primary subset is set
# by the established cohort loader, before within-series feature z-scoring.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
if (!exists("predict_ssbmf_yfb")) source("code/predict_ssbmf_yfb.R")
if (!exists("frozen_reverse_cindex")) source("code/concordance_ci.R")

#' Z-score an external series, align fixed training features and fill missing zeros
#' @param Y Finite numeric donor-by-feature matrix in original external units.
#' @param features Unique training feature names in fitted F row order.
#' @return List with aligned Y, absent-feature count and constant-feature count.
#' @examples
#' mwext_align(cbind(a = 1:3), c("a", "missing"))
mwext_align <- function(Y, features) {
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y)) ||
      is.null(colnames(Y)) || anyDuplicated(colnames(Y)) || anyDuplicated(features))
    stop("External alignment requires finite named features without duplicates.")
  center <- colMeans(Y); sd <- apply(Y, 2, stats::sd)
  constant <- sd == 0
  if (any(!is.finite(sd))) stop("External SDs are undefined.")
  sd[constant] <- 1
  Z <- sweep(sweep(Y, 2, center), 2, sd, "/")
  aligned <- matrix(0, nrow(Y), length(features), dimnames = list(rownames(Y), features))
  present <- features %in% colnames(Z)
  aligned[, present] <- Z[, features[present], drop = FALSE]
  list(Y = aligned, missing = sum(!present),
       constant = sum(features[present] %in% names(sd)[constant]))
}

#' Freeze program matches on TCGA and describe their score directions on ICGC
#' @param tcga_joint,tcga_expression Training patient-by-program projections.
#' @param icgc_joint,icgc_expression External donor-by-program projections.
#' @return Program matches, training/external correlations and direction agreement.
#' @examples
#' # mwext_direction(Z_joint_train, Z_expr_train, Z_joint_icgc, Z_expr_icgc)
mwext_direction <- function(tcga_joint, tcga_expression, icgc_joint, icgc_expression) {
  if (nrow(tcga_joint) != nrow(tcga_expression) || nrow(icgc_joint) != nrow(icgc_expression))
    stop("Program score rows are not aligned.")
  do.call(rbind, lapply(seq_len(ncol(tcga_joint)), function(k) {
    correlations <- vapply(seq_len(ncol(tcga_expression)), function(j)
      if (sd(tcga_joint[, k]) > 0 && sd(tcga_expression[, j]) > 0)
        cor(tcga_joint[, k], tcga_expression[, j]) else NA_real_, numeric(1))
    matched <- any(is.finite(correlations))
    j <- if (matched) which.max(replace(abs(correlations), !is.finite(correlations), -Inf)) else NA_integer_
    tr <- if (matched) correlations[j] else NA_real_
    ex <- if (matched && sd(icgc_joint[, k]) > 0 && sd(icgc_expression[, j]) > 0)
      cor(icgc_joint[, k], icgc_expression[, j]) else NA_real_
    assessable <- matched && abs(tr) >= 0.5 && is.finite(ex) && ex != 0
    data.frame(program = k, expression_program = j, tcga_correlation = tr,
      matched_on_tcga = matched && abs(tr) >= 0.5, icgc_correlation = ex,
      kept_direction = if (assessable) sign(tr) == sign(ex) else NA)
  }))
}

#' Score every existing full-data model after applying the registered primary filter
#' @return Invisible aggregate CIs, paired differences and program directions.
#' @examples
#' # mwext_run()
mwext_run <- function() {
  source("code/load_multiomics_data.R")
  source("code/fit_ssbmf_yfb.R")
  cfg <- yaml::read_yaml("config/globals.yml")$multiomics_data
  selection <- read.csv("results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weights/selection.csv")
  if (selection$part4_allowed || selection$selected_rule != "none")
    stop("This fixed scoring manifest requires the completed no-weight-selected study.")
  out <- "results/multimodal_real/outputs/modality_weights_icgc_2026-11-12"
  private <- "data/modality_weights_icgc_2026-11-12"
  if (dir.exists(out) || dir.exists(private)) stop("Refusing to rescore existing fixed models.")
  dir.create(out, recursive = TRUE); dir.create(private, recursive = TRUE)
  C <- "data/methylation_single_modality_2026-11-12/K05"
  paths <- c(expression_only = file.path(C, "expression_2218_PL_f0.rds"),
    methylation_10k = file.path(C, "methylation_10000_PL_f0.rds"),
    methylation_100k = file.path(C, "methylation_100000_PL_f0.rds"),
    joint_unweighted_10k = "data/multimodal_tcga_comparison_2026-11-12/joint_grouped_fold0.rds")
  models <- lapply(paths, readRDS)
  design_path <- file.path(C, "design.rds")
  design <- readRDS(design_path)
  fixed_inputs <- c(paths, design = design_path,
    external = file.path(cfg$root, "icgc/icgc_matched_data.rds"))
  before <- unname(tools::md5sum(fixed_inputs))
  write.csv(data.frame(input = names(fixed_inputs), md5 = before),
    file.path(out, "input_hashes.csv"), row.names = FALSE)
  sources <- c("results/multimodal_real/score_modality_weight_external.R",
    "code/load_multiomics_data.R", "code/predict_ssbmf_yfb.R", "code/concordance_ci.R", "config/globals.yml")
  write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
    file.path(out, "source_hashes.csv"), row.names = FALSE)
  writeLines(c("registration=docs/prompts/2026-11-12-codex-d-modality-weights-prompt.md Part5",
    "models=fixed K5 point-Laplace full fits; no refitting, tuning, sign or rank changes",
    "population=build_multiomics_cohorts validation_primary; one record per donor",
    "scaling=primary external series mean and sample SD per feature, before training-order alignment",
    "missing=zero after scaling; constant external features also zero, counted separately",
    "projection=sum_m Y_m F_m / training EF_norms; frozen training beta and Cox orientation",
    "interval=donor percentile bootstrap; B=2000; seed=1; flip=FALSE",
    "paired=joint minus expression_only, identical donor draws; flip_a=flip_b=FALSE",
    "direction=largest absolute TCGA patient-score correlation to expression; match abs-r>=0.5",
    "direction_readout=same correlation sign on ICGC; no orientation change",
    "status=descriptive; joint fit not validated by simulation; external chooses nothing",
    "imputation=existing ICGC reader on all available Task C design CpGs before primary restriction",
    "workers=1; numerical_threads=1"), file.path(out, "settings.txt"))
  # Use the established loader to select the donors, not an independently defined
  # clinical filter. Read the full frozen feature universe separately so its
  # outcome-free feature screen cannot discard a fixed model's input features.
  cohort <- build_multiomics_cohorts()
  ids <- names(cohort$validation_primary$time)
  reference_time <- cohort$validation_primary$time
  reference_event <- cohort$validation_primary$event
  rm(cohort); gc(FALSE)
  raw <- readRDS(file.path(cfg$root, "icgc/icgc_matched_data.rds"))
  before_filter <- ncol(raw$expr)
  cpgs <- intersect(colnames(design$methylation), rownames(raw$meth))
  rm(raw); gc(FALSE)
  external <- read_icgc_multiomics(cfg$root, cpgs,
    cfg$icgc_primary_histology, cfg$icgc_primary_sample_type)
  if (anyDuplicated(ids) || !all(ids %in% names(external$time)) ||
      !identical(unname(external$time[ids]), unname(reference_time)) ||
      !identical(unname(external$event[ids]), unname(reference_event)))
    stop("Primary donor survival does not reproduce the established loader.")
  time <- unname(reference_time); event <- unname(reference_event)
  expression <- t(external$expression[, ids, drop = FALSE])
  methylation <- t(external$methylation[, ids, drop = FALSE])
  n_usable <- length(external$time)
  rm(external); gc(FALSE)
  write.csv(data.frame(before_filter = before_filter, survival_eligible = n_usable,
    primary_donors = length(ids), events = sum(event), unique_donors = length(unique(ids))),
    file.path(out, "donors.csv"), row.names = FALSE)
  scores <- train_scores <- risks <- list(); summaries <- list()
  for (name in names(models)) {
    model <- models[[name]]; fit <- model$fit
    if (ncol(fit$EF) != 5L || !isTRUE(fit$history$converged))
      stop("Fixed model is not the registered converged K5 fit: ", name)
    groups <- if (name == "joint_unweighted_10k") model$groups else
      rep(if (name == "expression_only") "expression" else "methylation", length(model$features))
    if (length(groups) != nrow(fit$EF) || length(model$features) != nrow(fit$EF))
      stop("Frozen feature manifest does not match F.")
    blocks <- lapply(unique(groups), function(m)
      mwext_align(if (m == "expression") expression else methylation, model$features[groups == m]))
    names(blocks) <- unique(groups)
    Y <- matrix(0, length(ids), nrow(fit$EF))
    for (m in names(blocks)) Y[, groups == m] <- blocks[[m]]$Y
    # One projection of each fixed fit. Bootstrap operates on these frozen scores.
    predicted <- predict_ssbmf_yfb(Y, fit$EF, fit$EBeta, fit$EF_norms,
      projection_center = fit$projection_center)
    scores[[name]] <- predicted$ZF_test; risks[[name]] <- predicted$risk_scores
    if (name == "joint_unweighted_10k") {
      tr <- lapply(unique(groups), function(m) {
        X <- design[[m]][, model$features[groups == m], drop = FALSE]
        scale(X)
      })
      Ytr <- matrix(0, nrow(design$expression), nrow(fit$EF))
      for (i in seq_along(tr)) Ytr[, groups == unique(groups)[i]] <- tr[[i]]
      train_scores[[name]] <- predict_ssbmf_yfb(Ytr, fit$EF, fit$EBeta,
        fit$EF_norms, projection_center = fit$projection_center)$ZF_test
    } else {
      if (!identical(model$tr, seq_len(nrow(design$expression))))
        stop("Task C full-data training rows do not align with the saved design.")
      train_scores[[name]] <- model$scores
    }
    ci <- bootstrap_concordance_ci(risks[[name]], time, event, B = 2000, seed = 1, flip = FALSE)
    observed <- frozen_reverse_cindex(risks[[name]], time, event)
    if (!isTRUE(all.equal(ci$estimate, observed))) stop("External CI orientation mismatch.")
    summaries[[name]] <- data.frame(model = name, donors = length(ids), events = sum(event),
      C = observed, lower = ci$lower, upper = ci$upper, bootstrap_se = ci$se,
      missing_genes = if ("expression" %in% names(blocks)) blocks$expression$missing else 0L,
      missing_cpgs = if ("methylation" %in% names(blocks)) blocks$methylation$missing else 0L,
      constant_genes = if ("expression" %in% names(blocks)) blocks$expression$constant else 0L,
      constant_cpgs = if ("methylation" %in% names(blocks)) blocks$methylation$constant else 0L,
      simulation_validated = FALSE, converged = TRUE, rank = 5L)
    cat(name, "C =", observed, "\n"); flush.console()
  }
  joint <- "joint_unweighted_10k"
  ci <- bootstrap_concordance_diff_ci(risks[[joint]], risks$expression_only, time, event,
    B = 2000, seed = 1, flip_a = FALSE, flip_b = FALSE)
  difference <- data.frame(model = joint, reference = "expression_only", donors = length(ids),
    events = sum(event), difference = ci$estimate, lower = ci$lower, upper = ci$upper, bootstrap_se = ci$se)
  direction <- mwext_direction(train_scores[[joint]], train_scores$expression_only,
    scores[[joint]], scores$expression_only)
  direction$model <- joint
  if (!identical(before, unname(tools::md5sum(fixed_inputs))))
    stop("An external input or fixed fit changed during scoring.")
  saveRDS(list(donors = ids, time = time, event = event, scores = scores,
    risks = risks, train_scores = train_scores), file.path(private, "fixed_scores.rds"))
  summary <- do.call(rbind, summaries)
  write.csv(summary, file.path(out, "concordance.csv"), row.names = FALSE)
  write.csv(difference, file.path(out, "paired_differences.csv"), row.names = FALSE)
  write.csv(direction, file.path(out, "program_directions.csv"), row.names = FALSE)
  print(summary); print(difference); print(direction)
  invisible(list(summary = summary, difference = difference, direction = direction))
}
if (sys.nframe() == 0L) mwext_run()
