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
