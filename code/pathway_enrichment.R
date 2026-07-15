# ============================================================
# Script: pathway_enrichment.R
# Purpose: Pathway / gene-set enrichment on the recommended D4 model's
#          survival-active gene expression programs (Program 3, Program 7).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-14
# Dependencies: fgsea, clusterProfiler, msigdbr, org.Hs.eg.db
# ============================================================

# Section: Loading & program labeling ----

#' Load the D4 fit's gene weights and program labels.
#'
#' Reads the recommended (D4) DeSurv-comparison fit and its saved gene-symbol
#' vector, and attaches the direction-corrected program labels. Program 7 is
#' labeled "Adverse" and Program 3 "Protective" by the *marginal* (YF)-projection
#' survival direction (DECISIONS.md 2026-06-16) -- the joint posterior mean
#' beta-hat signs for these two programs are opposite (suppression among
#' correlated programs) and must NOT be used to derive the label. All other
#' programs are labeled "Inactive" (EBeta approx 0 in the recommended config).
#'
#' @return Named list:
#'   $EF             numeric matrix, 2064 genes x 7 programs, rows named by gene symbol
#'   $EBeta          numeric vector, length 7, posterior mean survival coefficients
#'   $EL             numeric matrix, n patients x 7 programs (posterior mean loadings)
#'   $gene_names     character vector, length 2064, gene symbols (same order as EF rows)
#'   $program_labels named list keyed by program index as a string ("1".."7"),
#'                    values in {"Adverse", "Protective", "Inactive"}
load_d4_weights <- function() {
  fits <- readRDS("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds")
  d4 <- fits[["D4"]]
  gene_names <- readRDS("results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")

  if (nrow(d4$EF) != length(gene_names)) {
    stop(sprintf("EF has %d rows but gene_names has length %d", nrow(d4$EF), length(gene_names)))
  }
  n_dup <- sum(duplicated(gene_names))
  if (n_dup > 0) {
    stop(sprintf("gene_names contains %d duplicate symbol(s); de-dup rule not implemented", n_dup))
  }

  EF <- d4$EF
  rownames(EF) <- gene_names

  program_labels <- as.list(rep("Inactive", ncol(EF)))
  names(program_labels) <- as.character(seq_len(ncol(EF)))
  program_labels[["7"]] <- "Adverse"
  program_labels[["3"]] <- "Protective"

  list(
    EF = EF,
    EBeta = d4$EBeta,
    EL = d4$EL,
    gene_names = gene_names,
    program_labels = program_labels
  )
}
