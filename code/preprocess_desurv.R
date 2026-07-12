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

#' Select the top-N most informative genes from an expression matrix.
#'
#' @param Y          numeric matrix (n × p), rows = subjects, columns = genes.
#' @param gene_names character vector of length p; gene identifiers.
#' @param top_n      integer or NULL. If NULL or >= ncol(Y), all genes are returned.
#' @param method     "variance" (default) ranks by variance only.
#'                   "combined_rank" ranks each gene by mean expression and by
#'                   variance independently (highest = rank 1), sums the two ranks,
#'                   and retains the top_n genes with the smallest rank sum.
#'                   This matches the DeSurv (Young et al. 2025) gene selection
#'                   criterion: genes that are both highly expressed AND highly
#'                   variable receive the lowest combined rank.
#' @return list with components Y (n × top_n), gene_names (length top_n), gene_var.
#' @family v2 preprocessing
select_top_variable_genes <- function(Y, gene_names, top_n = 2000,
                                       method = c("variance", "combined_rank")) {
  validate_expression_inputs(Y, gene_names)
  method <- match.arg(method)

  if (is.null(top_n) || top_n >= ncol(Y)) {
    return(list(Y = Y, gene_names = gene_names, gene_var = apply(Y, 2, stats::var)))
  }
  if (length(top_n) != 1 || !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n))
    stop("top_n must be NULL or an integer >= 1.")

  if (method == "variance") {
    # Original behaviour: order by variance descending, keep top top_n.
    gene_var <- apply(Y, 2, stats::var)
    ord      <- order(gene_var, decreasing = TRUE, na.last = NA)
  } else {
    # combined_rank (DeSurv criterion):
    # rank_mean = rank of negative mean (rank 1 = highest mean expression).
    # rank_var  = rank of negative variance (rank 1 = highest variance).
    # Select genes with the smallest combined rank_mean + rank_var.
    gene_mean <- colMeans(Y)
    gene_var  <- apply(Y, 2, stats::var)
    rank_mean <- rank(-gene_mean, ties.method = "average", na.last = "keep")
    rank_var  <- rank(-gene_var,  ties.method = "average", na.last = "keep")
    rank_sum  <- rank_mean + rank_var
    # order ascending: rank_sum=2 (rank 1 in both) is the best gene.
    ord <- order(rank_sum, na.last = NA)
  }

  keep         <- ord[seq_len(min(as.integer(top_n), length(ord)))]
  Y_keep       <- Y[, keep, drop = FALSE]
  gene_var_out <- apply(Y_keep, 2, stats::var)

  list(
    Y          = Y_keep,
    gene_names = gene_names[keep],
    gene_var   = gene_var_out
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

#' @param rank_transform logical; apply the per-subject rank transform (default
#'   TRUE, backward compatible). Phase 1c (see DECISIONS.md): external-cohort
#'   preprocessing must match whatever training used -- pass FALSE together
#'   with per_platform_standardize=TRUE to reproduce the DeSurv-aligned
#'   training pipeline's per-platform z-standardization instead.
#' @param per_platform_standardize logical; z-standardize each gene column
#'   (colMean=0, colSD=1) instead of/alongside the rank transform (default
#'   FALSE, backward compatible). Mirrors
#'   \code{\link{per_platform_standardize_cohorts}}'s per-cohort treatment in
#'   the training (\code{preprocess_merged_cohorts}) pipeline.
preprocess_desurv_cohort <- function(Y, gene_names,
                                     top_n = 2000,
                                     log_transform = TRUE,
                                     ties_method = "average",
                                     cohort_name = NULL,
                                     rank_transform = TRUE,
                                     per_platform_standardize = FALSE) {
  validate_expression_inputs(Y, gene_names)

  Y_proc <- Y
  if (log_transform)
    Y_proc <- log2_plus1_transform(Y_proc)

  filtered <- select_top_variable_genes(Y_proc, gene_names, top_n = top_n)
  Y_out <- filtered$Y
  colnames(Y_out) <- filtered$gene_names

  # Per-platform z-standardization and the per-subject rank transform are
  # alternative, mutually-exclusive-in-practice normalization strategies (see
  # code/preprocess_desurv.R's preprocess_merged_cohorts, which the training
  # pipeline uses): apply whichever the caller selects, matching training.
  if (per_platform_standardize) {
    Y_out <- per_platform_standardize_cohorts(list(Y_out))[[1]]
  }
  if (rank_transform) {
    Y_out <- rank_transform_subjects(Y_out, ties_method = ties_method)
    colnames(Y_out) <- filtered$gene_names
  }

  list(
    Y = Y_out,
    gene_names = filtered$gene_names,
    n = nrow(Y_out),
    p = ncol(Y_out),
    cohort_name = cohort_name,
    top_n = top_n,
    log_transform = log_transform,
    ties_method = ties_method,
    rank_transform = rank_transform,
    per_platform_standardize = per_platform_standardize
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
#' Now supports two gene selection strategies:
#'
#' \strong{selection_per_cohort = FALSE (default, original behavior):}
#' Gene selection runs AFTER normalization on the merged matrix (Steps 5-6).
#' top_n genes are selected by \code{selection_method} on the normalized data.
#'
#' \strong{selection_per_cohort = TRUE (DeSurv-aligned):}
#' Gene selection runs PER COHORT before normalization (Step 2b).
#' Each cohort independently selects its top-\code{top_n} genes, then the
#' intersection of those gene sets is carried forward. This avoids selecting
#' genes on data where variance has been equalized by z-standardization
#' (which defeats variance-based filtering).
#'
#' DeSurv-aligned pipeline (selection_per_cohort=TRUE):
#' \enumerate{
#'   \item Intersect raw gene universes across all training cohorts.
#'   \item Log2(x+1) transform per cohort (RNA-seq only).
#'   \item Per-cohort: select top-N genes by combined_rank, then intersect.
#'   \item Per-platform z-standardization on the intersected gene set.
#'   \item Row-bind into merged matrix.
#'   \item (Optional) rank-transform each subject.
#' }
#'
#' @param cohort_raw_list        named list of raw cohort objects with \code{$Y}
#'   and \code{$gene_names}.
#' @param log_transform_flags    named logical vector; TRUE = apply log2(x+1).
#' @param top_n                  integer; number of genes to select. Default 2000.
#' @param ties_method            ties method for rank(). Default "average".
#' @param rank_transform         logical; apply per-subject rank transform. Default TRUE.
#' @param per_platform_standardize logical; z-standardize each cohort before merging. Default FALSE.
#' @param normalize_method       "quantile", "z_score", or "none". Default "quantile".
#' @param selection_per_cohort   logical; if TRUE, gene selection runs per cohort
#'   before normalization (DeSurv-aligned). Default FALSE.
#' @param selection_method       "variance" or "combined_rank"; passed to
#'   \code{select_top_variable_genes()}. Default "variance".
#' @return list with Y, gene_names, n, p, dataset_labels, n_raw_intersect,
#'   rank_transform, per_platform_standardize, normalize_method,
#'   selection_per_cohort, selection_method.
#' @family v2 preprocessing
#' @seealso \code{\link{select_top_variable_genes}}
preprocess_merged_cohorts <- function(cohort_raw_list,
                                      log_transform_flags,
                                      top_n                    = 2000,
                                      ties_method              = "average",
                                      rank_transform           = TRUE,
                                      per_platform_standardize = FALSE,
                                      normalize_method         = c("quantile", "z_score", "none"),
                                      selection_per_cohort     = FALSE,
                                      selection_method         = c("variance", "combined_rank")) {
  normalize_method <- match.arg(normalize_method)
  selection_method <- match.arg(selection_method)
  cohort_names <- names(cohort_raw_list)
  stopifnot(!is.null(cohort_names), all(cohort_names %in% names(log_transform_flags)))

  # Step 1: intersect raw gene universes (no preprocessing yet)
  gene_lists   <- lapply(cohort_raw_list, function(x) x$gene_names)
  common_genes    <- Reduce(intersect, gene_lists)
  n_raw_intersect <- length(common_genes)
  if (n_raw_intersect == 0)
    stop("No common genes found across cohorts — check gene_names fields.")
  cat(sprintf("  [v2] Raw gene intersection: %d genes across %s\n",
              n_raw_intersect, paste(cohort_names, collapse = " + ")))

  # Steps 2–3: log2 transform per cohort (platform-aware), subset to common genes.
  cohort_matrices <- lapply(cohort_names, function(ds) {
    raw <- cohort_raw_list[[ds]]
    idx <- match(common_genes, raw$gene_names)
    Y   <- raw$Y[, idx, drop = FALSE]
    if (log_transform_flags[[ds]])
      Y <- log2_plus1_transform(Y)
    Y
  })
  names(cohort_matrices) <- cohort_names

  # Step 2b (DeSurv-aligned, optional): per-cohort gene selection BEFORE normalization.
  # Selects top-N genes within each cohort independently on the log-transformed data,
  # then intersects the selected sets. This prevents variance-equalization from
  # neutralizing the variance signal used for gene selection.
  if (selection_per_cohort) {
    cat(sprintf("  [v2] Per-cohort gene selection (top-%d, method=%s) before normalization ...\n",
                top_n, selection_method))
    selected_per_cohort <- lapply(cohort_names, function(ds) {
      sel <- select_top_variable_genes(cohort_matrices[[ds]], common_genes,
                                       top_n = top_n, method = selection_method)
      sel$gene_names
    })
    names(selected_per_cohort) <- cohort_names
    selected_genes <- Reduce(intersect, selected_per_cohort)
    if (length(selected_genes) == 0)
      stop("Per-cohort gene selection produced no common genes. Reduce top_n or check inputs.")
    cat(sprintf("  [v2] Genes after per-cohort selection and intersection: %d\n",
                length(selected_genes)))
    # Subset cohort matrices to the intersected gene set
    cohort_matrices <- lapply(cohort_names, function(ds) {
      idx <- match(selected_genes, common_genes)
      cohort_matrices[[ds]][, idx, drop = FALSE]
    })
    names(cohort_matrices) <- cohort_names
    common_genes <- selected_genes
  }

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

  # Step 4: distribution normalization across merged samples (method-dependent).
  if (normalize_method == "quantile") {
    cat(sprintf("  [v2] Quantile normalising merged matrix (%d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- quantile_normalize_merged(Y_merged)
  } else if (normalize_method == "z_score") {
    cat(sprintf("  [v2] Joint z-standardizing merged matrix (%d x %d, colMean=0, colSD=1) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- scale(Y_merged, center = TRUE, scale = TRUE)
    n_nan  <- sum(is.nan(Y_norm))
    if (n_nan > 0)
      cat(sprintf("  [v2] Note: %d NaN entries from zero-variance genes replaced with 0.\n", n_nan))
    Y_norm[is.nan(Y_norm)] <- 0
    colnames(Y_norm) <- common_genes
  } else {
    cat(sprintf("  [v2] Skipping normalization (log-transform only; %d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- Y_merged
  }

  # Steps 5–6: post-normalization gene selection (only when NOT using per-cohort selection).
  # When selection_per_cohort=TRUE the gene set was fixed in Step 2b; skip here.
  if (!selection_per_cohort) {
    selected <- select_top_variable_genes(Y_norm, common_genes, top_n = top_n,
                                          method = selection_method)
    cat(sprintf("  [v2] Genes retained after top-%d %s filter: %d\n",
                top_n, selection_method, length(selected$gene_names)))
  } else {
    selected <- list(Y = Y_norm, gene_names = common_genes)
  }

  # Step 7 (optional): rank-transform each subject within the selected gene set.
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
    n_raw_intersect          = n_raw_intersect,
    rank_transform           = rank_transform,
    per_platform_standardize = per_platform_standardize,
    normalize_method         = normalize_method,
    selection_per_cohort     = selection_per_cohort,
    selection_method         = selection_method
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
