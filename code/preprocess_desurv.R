# ============================================================
# Script: preprocess_desurv.R
# Purpose: DeSurv-aligned preprocessing helpers for cross-cohort
#          benchmarking on PDAC data.
# ============================================================

validate_expression_inputs <- function(Y, gene_names) {
  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix.")
  if (anyNA(Y))
    stop("Y must not contain missing values.")
  if (length(gene_names) != ncol(Y))
    stop("gene_names length must equal ncol(Y).")
  invisible(NULL)
}

log2_plus1_transform <- function(Y) {
  if (any(Y < 0, na.rm = TRUE))
    stop("log2_plus1_transform requires non-negative values.")
  log2(Y + 1)
}

select_top_variable_genes <- function(Y, gene_names, top_n = 2000) {
  validate_expression_inputs(Y, gene_names)

  if (is.null(top_n) || top_n >= ncol(Y)) {
    return(list(Y = Y, gene_names = gene_names, gene_var = apply(Y, 2, stats::var)))
  }
  if (length(top_n) != 1 || !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n))
    stop("top_n must be NULL or an integer >= 1.")

  gene_var <- apply(Y, 2, stats::var)
  ord <- order(gene_var, decreasing = TRUE, na.last = NA)
  keep <- ord[seq_len(min(as.integer(top_n), length(ord)))]

  list(
    Y = Y[, keep, drop = FALSE],
    gene_names = gene_names[keep],
    gene_var = gene_var[keep]
  )
}

rank_transform_subjects <- function(Y, ties_method = "average") {
  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix.")
  if (ncol(Y) == 0)
    stop("Y must have at least one gene/column.")

  ranked <- t(apply(Y, 1, function(x) rank(x, ties.method = ties_method)))
  storage.mode(ranked) <- "double"
  rownames(ranked) <- rownames(Y)
  colnames(ranked) <- colnames(Y)
  ranked
}

preprocess_desurv_cohort <- function(Y, gene_names,
                                     top_n = 2000,
                                     log_transform = TRUE,
                                     ties_method = "average",
                                     cohort_name = NULL) {
  validate_expression_inputs(Y, gene_names)

  Y_proc <- Y
  if (log_transform)
    Y_proc <- log2_plus1_transform(Y_proc)

  filtered <- select_top_variable_genes(Y_proc, gene_names, top_n = top_n)
  Y_ranked <- rank_transform_subjects(filtered$Y, ties_method = ties_method)
  colnames(Y_ranked) <- filtered$gene_names

  list(
    Y = Y_ranked,
    gene_names = filtered$gene_names,
    n = nrow(Y_ranked),
    p = ncol(Y_ranked),
    cohort_name = cohort_name,
    top_n = top_n,
    log_transform = log_transform,
    ties_method = ties_method
  )
}

intersect_preprocessed_cohorts <- function(cohort_list, reference = 1) {
  if (!is.list(cohort_list) || length(cohort_list) < 1)
    stop("cohort_list must be a non-empty list.")
  if (length(reference) != 1 || !is.finite(reference) || reference < 1 ||
      reference > length(cohort_list) || reference != as.integer(reference)) {
    stop("reference must be a valid cohort index.")
  }

  for (cohort in cohort_list) {
    if (is.null(cohort$Y) || is.null(cohort$gene_names))
      stop("Each cohort must contain Y and gene_names.")
  }

  common_genes <- Reduce(intersect, lapply(cohort_list, function(x) x$gene_names))
  if (length(common_genes) == 0)
    stop("No common genes across the supplied cohorts.")

  ref_genes <- cohort_list[[reference]]$gene_names
  common_genes <- ref_genes[ref_genes %in% common_genes]

  lapply(cohort_list, function(cohort) {
    idx <- match(common_genes, cohort$gene_names)
    cohort$Y <- cohort$Y[, idx, drop = FALSE]
    cohort$gene_names <- common_genes
    cohort$p <- length(common_genes)
    cohort
  })
}

# ============================================================
# v2 merged-cohort preprocessing (reordered pipeline)
# ============================================================

#' Z-standardize each cohort's gene columns independently.
#'
#' Applied before the merged row-bind so that each platform's per-gene
#' mean and SD are removed prior to joint quantile normalisation. This
#' addresses the A_surv/A_gen ≈ 10⁻³ imbalance observed in YFB merged
#' training: cross-platform QN equalises marginal distributions but not
#' per-platform variance contributions to A_k = Σᵢ wᵢ · ZF_ik².
#'
#' Operates independently per cohort — cohorts may have different nrow.
#' Only ncol must be identical (common gene set already intersected).
#'
#' @param cohort_matrices named list of numeric matrices (each n_c × p),
#'   already log-transformed and subsetted to the common gene universe.
#' @return named list of the same structure; each matrix has per-gene
#'   colMean ≈ 0 and colSD ≈ 1 (SD floored at 1e-8 to handle constant
#'   gene columns without crashing).
#' @family v2 preprocessing
per_platform_standardize_cohorts <- function(cohort_matrices) {
  lapply(cohort_matrices, function(Y) {
    gene_means <- colMeans(Y)
    gene_sds   <- apply(Y, 2, stats::sd)
    gene_sds   <- pmax(gene_sds, 1e-8)   # floor: constant genes stay finite
    sweep(sweep(Y, 2, gene_means, "-"), 2, gene_sds, "/")
  })
}

#' Quantile-normalize a merged expression matrix (genes × samples).
#'
#' Replaces each sample's gene distribution with the average quantile
#' distribution computed across all samples. Operates on the full merged
#' matrix blindly — no cohort labels required — so the correction does not
#' bake in batch-group information that would prevent generalisation to new
#' cohorts at prediction time.
#'
#' Thin wrapper around \code{preprocessCore::normalize.quantiles()}, which
#' expects a genes × samples matrix (transposed relative to the n × p
#' convention used elsewhere in the codebase).
#'
#' @param Y numeric matrix (n × p): samples in rows, genes in columns.
#' @return numeric matrix (n × p) with quantile-normalised values; row/col
#'   names preserved.
#' @family v2 preprocessing
quantile_normalize_merged <- function(Y) {
  if (!requireNamespace("preprocessCore", quietly = TRUE))
    stop("Package 'preprocessCore' is required. Install with:\n",
         "  BiocManager::install('preprocessCore')")

  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix (n x p).")

  # preprocessCore expects genes × samples (p × n); transpose in/out
  Y_qn <- preprocessCore::normalize.quantiles(t(Y))
  Y_out <- t(Y_qn)
  rownames(Y_out) <- rownames(Y)
  colnames(Y_out) <- colnames(Y)
  Y_out
}

#' v2 merged-cohort preprocessing: intersect first, then normalise jointly.
#'
#' Corrects the gene-selection order bug in the v1 pipeline. v1 runs
#' \code{preprocess_desurv_cohort()} per cohort (log2 → top-N by per-cohort
#' variance → rank), then intersects. Because TCGA top-2000 and CPTAC
#' top-2000 are computed on different assay types, only ~838 genes overlap —
#' and those genes are dominated by platform-level variation rather than
#' biology.
#'
#' v2 pipeline order (7 steps):
#' \enumerate{
#'   \item Intersect raw gene universes across all training cohorts.
#'   \item Log2(x + 1) transform per cohort (platform-aware: RNA-seq only).
#'   \item Row-bind into a single merged matrix.
#'   \item Quantile-normalize across all merged samples jointly
#'         (\code{preprocessCore::normalize.quantiles}).
#'   \item Compute per-gene variance across the full merged matrix.
#'   \item Select the top \code{top_n} most-variable genes.
#'   \item Rank-transform each subject within the selected gene set.
#' }
#'
#' Steps 4–7 operate on the merged matrix, so variance reflects combined
#' biological + residual batch variation rather than per-platform variance.
#' Quantile normalisation (step 4) reduces platform-scale differences without
#' requiring explicit cohort labels, preserving generalisability.
#'
#' @param cohort_raw_list   named list of raw cohort objects as returned by
#'   \code{load_pdac_raw()}: each must contain \code{$Y} (n × p numeric
#'   matrix) and \code{$gene_names} (character vector, length p).
#' @param log_transform_flags named logical vector; TRUE entries trigger
#'   log2(x + 1) for that cohort. Names must match \code{cohort_raw_list}.
#' @param top_n             integer; number of most-variable genes to retain
#'   after quantile normalisation. Default 2000.
#' @param ties_method       ties method passed to \code{rank()}. Default
#'   "average".
#' @param rank_transform    logical; if TRUE (default) apply per-subject rank
#'   transform (Step 7). Set FALSE for the no-rank sensitivity run: subjects
#'   retain the raw quantile-normalised values on the selected gene set.
#'   Rank transform forces all genes onto a uniform ordinal scale within each
#'   subject; disabling it preserves absolute QN expression differences, which
#'   may carry additional signal that rank-normalisation discards.
#' @return list with components:
#'   \describe{
#'     \item{Y}{numeric matrix (n_total × top_n) — QN-transformed, and
#'       rank-transformed if \code{rank_transform = TRUE}.}
#'     \item{gene_names}{character vector of retained gene names.}
#'     \item{n}{total number of training samples.}
#'     \item{p}{number of retained genes (= top_n or fewer if universe is
#'       smaller).}
#'     \item{dataset_labels}{factor of length n indicating cohort membership.}
#'     \item{n_raw_intersect}{number of genes in the raw intersection (before
#'       top-N selection) — use to verify Step 1 recovers substantially more
#'       than 838.}
#'     \item{rank_transform}{logical; echoes the input flag so downstream
#'       callers can record which preprocessing variant was used.}
#'   }
#' @family v2 preprocessing
#' @seealso \code{\link{preprocess_desurv_cohort}} (v1 single-cohort path),
#'   \code{\link{quantile_normalize_merged}}
preprocess_merged_cohorts <- function(cohort_raw_list,
                                      log_transform_flags,
                                      top_n                    = 2000,
                                      ties_method              = "average",
                                      rank_transform           = TRUE,
                                      per_platform_standardize = FALSE,
                                      normalize_method         = c("quantile", "z_score", "none")) {
  normalize_method <- match.arg(normalize_method)
  cohort_names <- names(cohort_raw_list)
  stopifnot(!is.null(cohort_names), all(cohort_names %in% names(log_transform_flags)))

  # Step 1: intersect raw gene universes (no preprocessing yet)
  gene_lists   <- lapply(cohort_raw_list, function(x) x$gene_names)
  common_genes <- Reduce(intersect, gene_lists)
  if (length(common_genes) == 0)
    stop("No common genes found across cohorts — check gene_names fields.")
  cat(sprintf("  [v2] Raw gene intersection: %d genes across %s\n",
              length(common_genes), paste(cohort_names, collapse = " + ")))

  # Steps 2–3: log2 transform per cohort (platform-aware), then subset to
  #            common genes and row-bind into a single n_total × p matrix.
  cohort_matrices <- lapply(cohort_names, function(ds) {
    raw <- cohort_raw_list[[ds]]
    idx <- match(common_genes, raw$gene_names)
    Y   <- raw$Y[, idx, drop = FALSE]
    if (log_transform_flags[[ds]])
      Y <- log2_plus1_transform(Y)
    Y
  })
  names(cohort_matrices) <- cohort_names

  # Optional: per-platform z-standardisation before row-bind.
  # Removes cohort-specific gene mean and SD so that no single platform's
  # variance scale dominates the merged survival precision A_k.
  if (per_platform_standardize) {
    cat("  [v2] Per-platform z-standardization (colMean=0, colSD=1 per cohort) ...\n")
    cohort_matrices <- per_platform_standardize_cohorts(cohort_matrices)
  }

  Y_merged <- do.call(rbind, cohort_matrices)
  colnames(Y_merged) <- common_genes

  dataset_labels <- factor(
    rep(cohort_names, vapply(cohort_raw_list, function(x) nrow(x$Y), integer(1))),
    levels = cohort_names
  )

  # Step 4: distribution normalization across merged samples (method-dependent)
  if (normalize_method == "quantile") {
    cat(sprintf("  [v2] Quantile normalising merged matrix (%d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- quantile_normalize_merged(Y_merged)
  } else if (normalize_method == "z_score") {
    cat(sprintf("  [v2] Joint z-standardizing merged matrix (%d x %d, colMean=0, colSD=1) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- scale(Y_merged, center = TRUE, scale = TRUE)
    colnames(Y_norm) <- common_genes
  } else {
    # normalize_method == "none": skip normalization; pass log-transformed matrix directly.
    # Gene-level mean differences across platforms are NOT removed here.
    # Downstream column-centering inside fit_modular/fit_cox_on_yf provides
    # partial correction, but platform-scale differences remain.
    cat(sprintf("  [v2] Skipping normalization (log-transform only; %d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- Y_merged
  }

  # Steps 5–6: per-gene variance on the normalized matrix → top-N selection
  selected <- select_top_variable_genes(Y_norm, common_genes, top_n = top_n)
  cat(sprintf("  [v2] Genes retained after top-%d variance filter: %d\n",
              top_n, length(selected$gene_names)))

  # Step 7 (optional): rank-transform each subject within the selected gene set.
  # Disabled when rank_transform = FALSE for the no-rank sensitivity run.
  if (rank_transform) {
    Y_final <- rank_transform_subjects(selected$Y, ties_method = ties_method)
    cat("  [v2] Rank-transforming subjects (Step 7).\n")
  } else {
    Y_final <- selected$Y
    cat("  [v2] Skipping rank transform (rank_transform = FALSE).\n")
  }
  colnames(Y_final) <- selected$gene_names

  list(
    Y                        = Y_final,
    gene_names               = selected$gene_names,
    n                        = nrow(Y_final),
    p                        = ncol(Y_final),
    dataset_labels           = dataset_labels,
    n_raw_intersect          = length(common_genes),
    rank_transform           = rank_transform,
    per_platform_standardize = per_platform_standardize,
    normalize_method         = normalize_method
  )
}

merge_preprocessed_cohorts <- function(cohort_list, dataset_labels = names(cohort_list)) {
  if (!is.list(cohort_list) || length(cohort_list) < 1)
    stop("cohort_list must be a non-empty list.")

  if (is.null(dataset_labels) || length(dataset_labels) != length(cohort_list)) {
    dataset_labels <- paste0("cohort_", seq_along(cohort_list))
  }

  first_genes <- cohort_list[[1]]$gene_names
  if (any(vapply(cohort_list, function(x) !identical(x$gene_names, first_genes), logical(1)))) {
    stop("All cohorts must share an identical, ordered gene list before merging.")
  }

  merged_Y <- do.call(rbind, lapply(cohort_list, function(x) x$Y))
  merged_labels <- factor(rep(dataset_labels, vapply(cohort_list, function(x) nrow(x$Y), integer(1))),
                          levels = dataset_labels)

  list(
    Y = merged_Y,
    gene_names = first_genes,
    n = nrow(merged_Y),
    p = ncol(merged_Y),
    dataset_labels = merged_labels
  )
}
