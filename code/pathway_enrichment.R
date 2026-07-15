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

# Section: Enrichment engines ----

#' Rank-based gene-set enrichment (fgsea) for one gene-weight vector.
#'
#' Because the point-exponential F prior makes all gene weights >= 0, the
#' weight vector is a natural one-sided continuous ranking statistic --
#' enrichment here is one-sided by construction (only positive NES is
#' meaningful; "depletion" has no interpretation for a non-negative weight).
#'
#' @param weights_k  named numeric vector of gene weights (names = gene symbols)
#' @param genesets   named list of character vectors (gene symbols per set)
#' @param seed       integer seed set before fgsea's permutation step (reproducibility)
#' @param minSize    minimum gene-set size after intersecting with weights_k's names
#' @param maxSize    maximum gene-set size after intersecting with weights_k's names
#' @param collection optional label identifying which gene-set collection this came from
#' @param scoreType  fgsea scoreType; default "pos" matches the non-negative
#'                    point-exponential weight vector (one-sided ranking, no "depletion")
#'
#' @return data.frame with columns: collection, set, size, NES, pval, padj, leading_edge
#'   (leading_edge genes joined by ";"), sorted by padj ascending.
run_fgsea_program <- function(weights_k, genesets, seed = 1, minSize = 10, maxSize = 500,
                               collection = NA_character_, scoreType = "pos") {
  if (is.null(names(weights_k)) || any(names(weights_k) == "")) {
    stop("run_fgsea_program: weights_k must be a named numeric vector (names = gene symbols)")
  }

  set.seed(seed)
  res <- fgsea::fgsea(pathways = genesets, stats = weights_k, minSize = minSize, maxSize = maxSize,
                       scoreType = scoreType)

  if (nrow(res) == 0 || all(is.na(res$pval))) {
    stop("run_fgsea_program: fgsea returned zero results or all-NA p-values -- ",
         "check gene-symbol overlap between weights_k and genesets")
  }

  df <- data.frame(
    collection = collection,
    set = res$pathway,
    size = res$size,
    NES = res$NES,
    pval = res$pval,
    padj = res$padj,
    leading_edge = vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  df[order(df$padj), ]
}

#' Over-representation analysis (hypergeometric test) on top-N weighted genes.
#'
#' Confirmatory cross-check for run_fgsea_program(): tests whether the top-N
#' genes by weight are enriched in each gene set relative to a specified
#' background (the 2064 selected genes, NOT the whole genome -- using the
#' genome as background inflates significance because the 2064 genes are
#' already survival/variance-selected).
#'
#' @param top_genes  character vector of gene symbols (top-N weighted genes for one program)
#' @param background character vector of gene symbols (the full selected-gene universe)
#' @param genesets   named list of character vectors (gene symbols per set)
#' @param collection optional label identifying which gene-set collection this came from
#'
#' @return data.frame with columns: collection, set, size, GeneRatio, pval, padj, leading_edge
#'   (leading_edge = overlapping gene symbols, "/"-joined by clusterProfiler)
run_ora_program <- function(top_genes, background, genesets, collection = NA_character_) {
  overlap_total <- sum(vapply(genesets, function(s) length(intersect(s, top_genes)), integer(1)))
  if (overlap_total == 0) {
    stop("run_ora_program: zero overlap between top_genes and every gene set -- ",
         "check gene-symbol case/mapping")
  }

  term2gene <- do.call(rbind, lapply(names(genesets), function(nm) {
    data.frame(term = nm, gene = genesets[[nm]], stringsAsFactors = FALSE)
  }))

  res <- clusterProfiler::enricher(
    gene = top_genes, universe = background, TERM2GENE = term2gene,
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = 100000
  )

  if (is.null(res) || nrow(res@result) == 0) {
    stop("run_ora_program: enricher() returned zero results")
  }

  rdf <- res@result
  data.frame(
    collection = collection,
    set = rdf$ID,
    size = rdf$Count,
    GeneRatio = rdf$GeneRatio,
    pval = rdf$pvalue,
    padj = rdf$p.adjust,
    leading_edge = rdf$geneID,
    stringsAsFactors = FALSE
  )
}
