# ============================================================
# Script: run_multimodal_tcga_comparison.R
# Purpose: Describe matched-TCGA single and joint modality fits.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, survival, ebnm, Rcpp and existing YFB helpers
# ============================================================
# Run from the repository root: Rscript <script> [--smoke]
# One sequential worker, one BLAS thread. Fits and scores remain under data/.
# This descriptive comparison is not validated by simulation.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_feature_group_simulation.R")

#' Prepare a saved fold with training-only modality scaling
#' @param design Task C design containing matrices, ranked CpGs, outcomes and folds.
#' @param fold Zero for all patients or a saved held-out fold.
#' @param n_cpg Number of ranked CpGs, fixed at 10000 for the study.
#' @return Training and test indices, frozen fold labels, and standardized matrices.
#' @examples
#' # x <- mtcga_prepare_fold(design, 1L)
mtcga_prepare_fold <- function(design, fold, n_cpg = 10000L) {
  n <- length(design$time)
  if (length(design$fold) != n || !all(design$fold %in% 1:5) ||
      length(design$ranked) < n_cpg) stop("Invalid frozen design or CpG count.")
  tr <- if (fold == 0L) seq_len(n) else which(design$fold != fold)
  te <- if (fold == 0L) integer() else which(design$fold == fold)
  raw <- list(expression = design$expression,
              methylation = design$methylation[, design$ranked[seq_len(n_cpg)], drop = FALSE])
  scaled <- fgsim_scale(lapply(raw, function(x) x[tr, , drop = FALSE]),
                         lapply(raw, function(x) x[te, , drop = FALSE]))
  list(tr = tr, te = te, fold = design$fold, training = scaled$training,
       test_data = scaled$validation, center = scaled$center, sd = scaled$sd)
}

#' Check fitted grouped priors using the registered assertion
#' @param fit Grouped YFB fit.
#' @param K Number of programs.
#' @return Invisible TRUE when valid; otherwise an error.
#' @examples
#' # mtcga_assert_priors(fit, 5L)
mtcga_assert_priors <- function(fit, K) fgsim_assert_priors(fit, K)

#' Calculate paired five-fold mean and standard error
#' @param joint,expression Fold-ordered concordance vectors.
#' @return Mean difference and SE, or both NA when any fold is undefined.
#' @examples
#' mtcga_paired(rep(.6, 5), rep(.5, 5))
mtcga_paired <- function(joint, expression) {
  stopifnot(length(joint) == 5L, length(expression) == 5L)
  d <- joint - expression
  if (any(!is.finite(d))) return(c(mean_difference = NA_real_, se = NA_real_))
  c(mean_difference = mean(d), se = stats::sd(d) / sqrt(5))
}

#' Identify joint programs unmatched to either single-modality fit
#' @param joint,expression,methylation Feature-by-program loading matrices.
#' @param groups Modality label for each joint feature, in joint feature order.
#' @param threshold Absolute loading-correlation match threshold.
#' @return Data frame with maximum correlations and joint-only classification.
#' @examples
#' # mtcga_joint_only(joint, expression, methylation, groups)
mtcga_joint_only <- function(joint, expression, methylation, groups, threshold = 0.5) {
  singles <- list(expression = expression, methylation = methylation)
  correlations <- lapply(names(singles), function(m) {
    block <- joint[groups == m, , drop = FALSE]
    if (nrow(block) != nrow(singles[[m]])) stop("Loading blocks do not align.")
    suppressWarnings(abs(stats::cor(block, singles[[m]])))
  })
  maxima <- sapply(correlations, function(r) apply(r, 1, function(x)
    if (all(!is.finite(x))) NA_real_ else max(x[is.finite(x)])))
  counts <- sapply(correlations, function(r) rowSums(is.finite(r)))
  live <- colSums(joint^2) > 0
  defined <- counts[, 1] > 0 & counts[, 2] > 0
  data.frame(program = seq_len(ncol(joint)), expression_correlation = maxima[, 1],
    methylation_correlation = maxima[, 2], expression_defined = counts[, 1],
    methylation_defined = counts[, 2], live = live, correlations_defined = defined,
    joint_only = live & defined & ifelse(defined, maxima[, 1] < threshold & maxima[, 2] < threshold, FALSE))
}

#' Initialize one matrix with the same rank-K SVD for both joint arms
#' @param Y Standardized patient-by-feature matrix.
#' @param K Number of programs.
#' @return List of loading and feature initial moments.
#' @examples
#' # mtcga_init(Y, 5L)
mtcga_init <- function(Y, K) {
  s <- svd(Y, nu = K, nv = K)
  list(EL_init = s$u[, seq_len(K), drop = FALSE] %*% diag(sqrt(s$d[seq_len(K)])),
       EF_init = s$v[, seq_len(K), drop = FALSE] %*% diag(sqrt(s$d[seq_len(K)])))
}

#' Fit one arm and retain private results and individual warnings
#' @param design Frozen input; x prepared fold; arm fit name; cfg configuration.
#' @param fold Fold number, zero for full fit; private private output directory.
#' @return Checkpoint list with aggregate row and private fitted state.
#' @examples
#' # mtcga_fit(design, x, "expression", cfg, 1L, "data/private")
mtcga_fit <- function(design, x, arm, cfg, fold, private) {
  key <- sprintf("%s_fold%d", arm, fold)
  path <- file.path(private, paste0(key, ".rds"))
  if (file.exists(path)) stop("Checkpoint already exists: ", path)
  modalities <- switch(arm, expression = "expression", methylation = "methylation",
                        joint_grouped = c("expression", "methylation"),
                        joint_pooled = c("expression", "methylation"))
  Y <- do.call(cbind, x$training[modalities])
  test <- do.call(cbind, x$test_data[modalities])
  groups <- rep(modalities, vapply(x$training[modalities], ncol, integer(1)))
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  answer <- tryCatch(withCallingHandlers({
    init <- mtcga_init(Y, 5L)
    inputs <- c(list(Y = Y, time = design$time[x$tr], event = design$event[x$tr]), init)
    fitted <- fgsim_fit(inputs, if (arm == "joint_grouped") groups else NULL, cfg)
    fit <- fitted$fit
    if (arm == "joint_grouped") mtcga_assert_priors(fit, 5L)
    scoring <- if (fold == 0L) Y else test
    pred <- predict_ssbmf_yfb(scoring, fit$EF, fit$EBeta, fit$EF_norms,
                              projection_center = fit$projection_center)$risk_scores
    ii <- if (fold == 0L) x$tr else x$te
    concordance <- if (stats::sd(pred) > 1e-8)
      frozen_reverse_cindex(pred, design$time[ii], design$event[ii]) else NA_real_
    list(row = data.frame(arm, fold, status = if (isTRUE(fit$history$converged)) "converged" else "capped",
      converged = isTRUE(fit$history$converged), iterations = fit$history$n_iter,
      seconds = proc.time()[["elapsed"]] - started, C = concordance, error = NA_character_),
      fit = fit, priors = fitted$priors, groups = groups, features = colnames(Y), scores = pred)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(row = data.frame(arm, fold, status = "failed", converged = FALSE,
    iterations = NA_integer_, seconds = proc.time()[["elapsed"]] - started,
    C = NA_real_, error = conditionMessage(e))))
  answer$warnings <- warnings
  saveRDS(answer, path)
  answer
}

#' Write descriptive aggregate tables for all completed checkpoints
#' @param records List of private checkpoint results.
#' @param out Aggregate output directory.
#' @return NULL invisibly; writes CSV tables and settings.
#' @examples
#' # mtcga_summarize(records, "results/output")
mtcga_summarize <- function(records, out) {
  label <- "not validated by simulation"
  #' Write one aggregate table with its validation label
  #' @param x Aggregate data frame.
  #' @param name Output file stem.
  #' @return NULL invisibly after writing the CSV.
  #' @examples
  #' # write_table(fits, "fits")
  write_table <- function(x, name) {
    x$validation <- rep(label, nrow(x))
    utils::write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
  }
  rows <- do.call(rbind, lapply(records, `[[`, "row"))
  write_table(rows, "fits")
  warnings <- do.call(rbind, lapply(records, function(r) data.frame(arm = rep(r$row$arm, length(r$warnings)),
    fold = rep(r$row$fold, length(r$warnings)), warning = r$warnings)))
  write_table(warnings, "warnings")
  cv <- rows[rows$fold > 0, , drop = FALSE]
  summary <- do.call(rbind, lapply(unique(rows$arm), function(a) {
    v <- cv[cv$arm == a, ]; e <- cv[cv$arm == "expression", ]
    all_folds <- identical(sort(v$fold), 1:5) && identical(sort(e$fold), 1:5)
    paired <- if (all_folds && a %in% c("joint_grouped", "joint_pooled"))
      mtcga_paired(v$C[match(1:5, v$fold)], e$C[match(1:5, e$fold)]) else c(NA_real_, NA_real_)
    data.frame(arm = a, folds = nrow(v), mean_C = if (all_folds && all(is.finite(v$C))) mean(v$C) else NA_real_,
      mean_difference = paired[1], paired_SE = paired[2])
  }))
  write_table(summary, "cv_summary")
  joint <- Filter(function(r) r$row$arm %in% c("joint_grouped", "joint_pooled") && !is.null(r$fit), records)
  precision <- do.call(rbind, lapply(joint, function(r) {
    p <- fgsim_precision(r$fit, r$groups)
    p <- p[colSums(r$fit$EF^2)[p$program] > 0, , drop = FALSE]
    cbind(arm = rep(r$row$arm, nrow(p)), fold = rep(r$row$fold, nrow(p)), p)
  }))
  if (is.null(precision)) precision <- data.frame(arm = character(), fold = integer(), program = integer(), modality = character(), precision = numeric(), share = numeric())
  write_table(precision, "precision")
  blocks <- do.call(rbind, lapply(Filter(function(r) r$row$arm == "joint_grouped" && !is.null(r$fit), records), function(r) {
    b <- r$fit$empty_blocks
    if (is.null(b)) stop("Grouped fit lacks empty-block counts.")
    data.frame(arm = r$row$arm, fold = r$row$fold, modality = rep(rownames(b), times = ncol(b)),
      program = rep(seq_len(ncol(b)), each = nrow(b)), events = as.vector(b))
  }))
  if (is.null(blocks)) blocks <- data.frame(arm = character(), fold = integer(), modality = character(), program = integer(), events = integer())
  write_table(blocks, "empty_blocks")
  full <- Filter(function(r) r$row$fold == 0L && !is.null(r$fit), records)
  names(full) <- vapply(full, function(r) r$row$arm, character(1))
  priors <- data.frame(arm = character(), program = integer(), modality = character(), point_mass_weight = numeric(), laplace_scale = numeric())
  for (arm in intersect(c("joint_grouped", "joint_pooled"), names(full))) for (k in 1:5) for (m in c("expression", "methylation")) {
    p <- if (arm == "joint_grouped") full[[arm]]$priors[[m]][[k]] else full[[arm]]$priors[[k]]
    if (is.null(p$pi)) stop("Missing fitted prior.")
    priors <- rbind(priors, data.frame(arm = arm, program = k, modality = m,
      point_mass_weight = 1 - p$pi, laplace_scale = if (isTRUE(p$point_mass)) NA_real_ else 1 / p$rate))
  }
  write_table(priors, "full_priors")
  matches <- data.frame(arm = character(), program = integer(), expression_correlation = numeric(), methylation_correlation = numeric(), expression_defined = integer(), methylation_defined = integer(), live = logical(), correlations_defined = logical(), joint_only = logical())
  top <- data.frame(arm = character(), program = integer(), modality = character(), rank = integer(), feature = character(), loading = numeric())
  if (all(c("expression", "methylation") %in% names(full))) for (arm in intersect(c("joint_grouped", "joint_pooled"), names(full))) {
    r <- full[[arm]]
    m <- mtcga_joint_only(r$fit$EF, full$expression$fit$EF, full$methylation$fit$EF, r$groups)
    matches <- rbind(matches, cbind(arm = arm, m))
    for (k in m$program[!is.na(m$joint_only) & m$joint_only]) for (modality in c("expression", "methylation")) {
      idx <- which(r$groups == modality); idx <- idx[order(abs(r$fit$EF[idx, k]), decreasing = TRUE)][seq_len(min(20L, length(idx)))]
      top <- rbind(top, data.frame(arm = arm, program = k, modality = modality,
        rank = seq_along(idx), feature = r$features[idx], loading = r$fit$EF[idx, k]))
    }
  }
  write_table(matches, "full_program_matches")
  write_table(top, "joint_only_top_features")
  invisible(NULL)
}

#' Run the registered descriptive comparison
#' @param smoke TRUE runs one fold in separate directories.
#' @return NULL invisibly; failures are recorded and result in an error.
#' @examples
#' # mtcga_run(TRUE)
mtcga_run <- function(smoke = FALSE) {
  suffix <- if (smoke) "2026-11-12_smoke" else "2026-11-12"
  private <- paste0("data/multimodal_tcga_comparison_", suffix)
  out <- paste0("results/multimodal_real/outputs/multimodal_tcga_comparison_", suffix)
  dir.create(private, recursive = TRUE, showWarnings = FALSE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  design <- readRDS("data/methylation_single_modality_2026-11-12/design.rds")
  if (length(design$time) != 144L || !identical(sort(unique(design$fold)), 1:5)) stop("Frozen design does not match registration.")
  cfg <- yaml::read_yaml("config/globals.yml")
  if (cfg$benchmark$k_working_yfb != 5L) stop("Configuration rank differs from registered settings.")
  if (!multimodal_yfb_use_cpp()) stop("Compiled sweep is required.")
  sources <- c("results/multimodal_real/run_multimodal_tcga_comparison.R",
    "results/benchmark_sim/run_feature_group_simulation.R",
    list.files("code", pattern = "\\.(R|cpp)$", full.names = TRUE), "config/globals.yml")
  input <- "data/methylation_single_modality_2026-11-12/design.rds"
  before <- unname(tools::md5sum(input))
  manifest <- data.frame(path = sources, md5 = unname(tools::md5sum(sources)),
                         validation = "not validated by simulation")
  utils::write.csv(manifest, file.path(out, "source_hashes.csv"), row.names = FALSE)
  utils::write.csv(data.frame(path = input, md5_before = before, md5_after = NA_character_,
    validation = "not validated by simulation"), file.path(out, "input_hashes.csv"), row.names = FALSE)
  packages <- c("yaml", "survival", "ebnm", "Rcpp")
  utils::write.csv(data.frame(package = packages,
    version = vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)),
    validation = "not validated by simulation"), file.path(out, "package_versions.csv"), row.names = FALSE)
  writeLines(c("not validated by simulation", paste("K=5; max_iter=", cfg$cavi$max_iter,
    "; tol=", cfg$cavi$tol, "; alpha=", cfg$benchmark$alpha),
    "point-Laplace L/F; Normal beta; w_S=1; canonicalize=TRUE; no strata or pruning.",
    "Task C saved folds and D4 genes; 10000 ranked CpGs; training-only z-scores; one worker."),
    file.path(out, "settings.txt"))
  # The smoke fold may be reused only when its exact code and frozen input agree.
  smoke_private <- "data/multimodal_tcga_comparison_2026-11-12_smoke"
  smoke_out <- "results/multimodal_real/outputs/multimodal_tcga_comparison_2026-11-12_smoke"
  reuse <- FALSE
  if (!smoke && file.exists(file.path(smoke_out, "source_hashes.csv")) &&
      file.exists(file.path(smoke_out, "input_hashes.csv"))) {
    old_sources <- utils::read.csv(file.path(smoke_out, "source_hashes.csv"))
    old_input <- utils::read.csv(file.path(smoke_out, "input_hashes.csv"))
    reuse <- identical(old_sources, manifest) && identical(old_input$md5_before, before) &&
      identical(old_input$md5_after, before)
    if (!reuse) stop("Smoke sources or input differ; cannot reuse its checkpoints.")
  }
  records <- list()
  for (fold in if (smoke) 1L else 0:5) {
    x <- mtcga_prepare_fold(design, fold)
    for (arm in c("expression", "methylation", "joint_grouped", "joint_pooled")) {
      checkpoint <- file.path(private, sprintf("%s_fold%d.rds", arm, fold))
      old <- file.path(smoke_private, sprintf("%s_fold%d.rds", arm, fold))
      if (reuse && fold == 1L && file.exists(old)) {
        if (file.exists(checkpoint) || !file.copy(old, checkpoint)) stop("Cannot reuse smoke checkpoint safely.")
        r <- readRDS(checkpoint)
        if (r$row$status == "failed") stop("Smoke checkpoint failed; full run cannot proceed.")
      } else r <- mtcga_fit(design, x, arm, cfg, fold, private)
      records[[length(records) + 1L]] <- r
      mtcga_summarize(records, out)
      if (!identical(unname(tools::md5sum(sources)), manifest$md5))
        stop("A source file changed during execution.")
      after <- unname(tools::md5sum(input))
      utils::write.csv(data.frame(path = input, md5_before = before, md5_after = after,
        validation = "not validated by simulation"), file.path(out, "input_hashes.csv"), row.names = FALSE)
      if (!identical(before, after)) stop("Frozen task C input changed during execution.")
      message(arm, " fold ", fold, ": ", r$row$status)
    }
  }
  mtcga_summarize(records, out)
  if (any(vapply(records, function(r) r$row$status == "failed", logical(1)))) stop("One or more fits failed; see aggregate diagnostics.")
  invisible(NULL)
}

if (sys.nframe() == 0L) mtcga_run("--smoke" %in% commandArgs(trailingOnly = TRUE))
