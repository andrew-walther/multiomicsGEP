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
#' correlated programs) and must NOT be used to derive the label.
#'
#' The remaining 5 programs are split three ways rather than lumped as
#' "Inactive" (DECISIONS.md 2026-08-19, `classify_factors()` in select_K.R):
#' programs 5 and 6 are **genomics-only** (real gene-expression programs with
#' PVE > 1% but no survival coefficient), while 1, 2 and 4 are **dead** (fully
#' pruned by ARD, nothing retained). The 4 "kept" factors are therefore
#' {3, 5, 6, 7}. Superseding the earlier blanket "Inactive" label matters
#' because it changes which programs are worth characterizing: the pathway
#' enrichment for programs 5 and 6 was computed under the old convention but
#' filtered out of every reported table.
#'
#' Survival-active membership is *verified* here against |EBeta| > beta_thresh
#' rather than trusted, so a changed fit fails loudly instead of silently
#' mislabeling. The genomics-only/dead split is taken from the verified
#' `classify_factors()` result, since recomputing PVE requires the training
#' matrix Y, which this loader does not read.
#'
#' @param beta_thresh numeric: |EBeta_k| threshold for survival-active, matching
#'                    `config/globals.yml` `k_selection$beta_threshold` (default 0.001)
#'
#' @return Named list:
#'   $EF              numeric matrix, 2064 genes x 7 programs, rows named by gene symbol
#'   $EBeta           numeric vector, length 7, posterior mean survival coefficients
#'   $EL              numeric matrix, n patients x 7 programs (posterior mean loadings)
#'   $gene_names      character vector, length 2064, gene symbols (same order as EF rows)
#'   $program_labels  named list keyed by program index as a string ("1".."7"),
#'                     values in {"Adverse", "Protective", "Genomics-only", "Dead"}
#'   $program_class   named list, same keys, values in
#'                     {"survival_active", "genomics_only", "dead"}
#'   $survival_active integer vector: program indices with |EBeta| > beta_thresh
#'   $genomics_only   integer vector: genomics-only program indices
#'   $kept_factors    integer vector: survival_active + genomics_only, sorted
#'
#' @seealso \code{\link{classify_factors}} in `code/select_K.R`, the source of
#'   the three-way split.
load_d4_weights <- function(beta_thresh = 0.001) {
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

  # Verified from the fit: which programs carry survival weight.
  surv_active <- sort(which(abs(d4$EBeta) > beta_thresh))
  if (!identical(as.integer(surv_active), c(3L, 7L))) {
    stop(sprintf(paste0("survival-active programs are {%s} but the labeling convention in this ",
                        "function (Program 3 Protective / Program 7 Adverse, DECISIONS.md ",
                        "2026-06-16) assumes {3, 7}. The fit appears to have changed; re-derive ",
                        "the marginal survival directions before relabeling."),
                 paste(surv_active, collapse = ", ")))
  }

  # From the verified classify_factors() result on this fit (DECISIONS.md 2026-08-19).
  # PVE is not recomputed here because that needs the training matrix Y.
  genomics_only <- c(5L, 6L)

  program_class <- as.list(rep("dead", ncol(EF)))
  names(program_class) <- as.character(seq_len(ncol(EF)))
  for (k in surv_active)   program_class[[as.character(k)]] <- "survival_active"
  for (k in genomics_only) program_class[[as.character(k)]] <- "genomics_only"

  program_labels <- lapply(program_class, function(cl) {
    switch(cl, genomics_only = "Genomics-only", dead = "Dead", cl)
  })
  program_labels[["7"]] <- "Adverse"
  program_labels[["3"]] <- "Protective"

  list(
    EF = EF,
    EBeta = d4$EBeta,
    EL = d4$EL,
    gene_names = gene_names,
    program_labels = program_labels,
    program_class = program_class,
    survival_active = as.integer(surv_active),
    genomics_only = genomics_only,
    kept_factors = sort(c(as.integer(surv_active), genomics_only))
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

# Section: Custom PDAC gene sets (Moffitt / Bailey / DeSurv) ----

#' Parse gene-symbol tables out of pdftotext -layout output (DeSurv SI appendix style).
#'
#' Pure text-parsing helper, factored out of extract_desurv_genesets() so it can
#' be unit-tested without the actual PDF. Each "table" is a block of rows with a
#' fixed number of whitespace-separated gene-symbol columns; this function
#' strips the surrounding column-range header row (e.g. "1-45  46-90 ..."),
#' blank lines, and bare page-number footer lines, and flattens the remaining
#' rows into one gene vector per factor (row-major within each row, tables
#' read left-to-right top-to-bottom as in the source PDF).
#'
#' @param lines           character vector, one element per PDF text line
#'                        (as returned by `pdftotext -layout <pdf> -`)
#' @param factor_patterns named character vector of regex patterns marking the
#'                        start of each factor's table (e.g. c(D1 = "Table S7:"))
#' @param end_pattern     regex marking the end of the last table (first line
#'                        of the section that follows, e.g. the references list)
#' @param expected_n      if not NULL, stop() if any factor's gene count differs
#'                        from this value (guards against silent PDF-layout drift)
#'
#' @return named list (same names as factor_patterns), each a character vector of gene symbols
parse_desurv_si_text <- function(lines, factor_patterns, end_pattern, expected_n = NULL) {
  lines <- gsub("\f", "", lines, fixed = TRUE)

  starts <- vapply(factor_patterns, function(p) grep(p, lines)[1], integer(1))
  if (any(is.na(starts))) {
    stop("parse_desurv_si_text: could not locate start pattern(s): ",
         paste(names(factor_patterns)[is.na(starts)], collapse = ", "))
  }
  ends <- c(starts[-1], grep(end_pattern, lines)[1])
  if (is.na(ends[length(ends)])) {
    stop("parse_desurv_si_text: could not locate end_pattern: ", end_pattern)
  }

  parse_block <- function(block) {
    keep <- vapply(block, function(r) {
      tr <- trimws(r)
      if (nchar(tr) == 0 || grepl("^[0-9]+$", tr)) return(FALSE)
      toks <- strsplit(tr, "[ ]{2,}")[[1]]
      length(toks) > 0 && grepl("^[A-Za-z0-9.-]+$", toks[1])
    }, logical(1))
    data_rows <- block[keep]
    unlist(lapply(data_rows, function(r) strsplit(trimws(r), "[ ]{2,}")[[1]]))
  }

  result <- Map(function(s, e) parse_block(lines[(s + 1):(e - 1)]), starts, ends)
  names(result) <- names(factor_patterns)

  if (!is.null(expected_n)) {
    bad <- names(result)[vapply(result, length, integer(1)) != expected_n]
    if (length(bad) > 0) {
      stop(sprintf("parse_desurv_si_text: factor(s) %s did not yield exactly %d genes -- ",
                    paste(bad, collapse = ", "), expected_n),
           "check for SI appendix layout changes or a parsing regression")
    }
  }

  result
}

#' Extract the DeSurv D1/D2/D3 factor-specific gene lists from the SI appendix PDF.
#'
#' Per DECISIONS.md 2026-06-16 (D3): DeSurv program gene lists are extracted
#' from the paper's SI appendix (Tables S7-S9, top 270 genes per factor by
#' factor-specificity score). Falls back with an informative error (not a
#' silent empty result) if the `pdftotext` system binary is unavailable or the
#' expected tables aren't found -- the documented fallback in that case is the
#' DeSurv GitHub repo's own gene lists (not implemented here; only needed if
#' this path fails).
#'
#' @param si_pdf_path path to si_appendix.pdf (Young et al., DeSurv)
#' @return named list: D1 (Classical tumor), D2 (stromal/immune), D3 (Basal-like tumor);
#'   each a character vector of 270 gene symbols
extract_desurv_genesets <- function(si_pdf_path) {
  if (!file.exists(si_pdf_path)) {
    stop("extract_desurv_genesets: si_appendix.pdf not found at ", si_pdf_path,
         " -- fall back to the DeSurv GitHub repo's gene lists (D3 fallback path)")
  }
  if (nchar(Sys.which("pdftotext")) == 0) {
    stop("extract_desurv_genesets: 'pdftotext' (poppler-utils) not found on PATH -- ",
         "fall back to the DeSurv GitHub repo's gene lists (D3 fallback path)")
  }

  lines <- system2("pdftotext", c("-layout", shQuote(si_pdf_path), "-"), stdout = TRUE)

  parse_desurv_si_text(
    lines,
    factor_patterns = c(D1 = "Table S7:", D2 = "Table S8:", D3 = "Table S9:"),
    end_pattern = "^Bailey, Peter",
    expected_n = 270
  )
}

#' Load Moffitt (basal/classical) and Bailey (4-subtype) gene signatures from local reference data.
#'
#' Both signatures live in `cmbSubtypes.RData` (schemaList / subtypeGeneList),
#' an ancillary reference file already present in PDAC_DATA_ROOT for subtype
#' classification elsewhere in this project. "MT" (Moffitt Tumor) is the
#' 25-basal + 25-classical classifier from Moffitt et al. 2015; "Bailey" is
#' the 612-gene, 4-subtype signature from Bailey et al. 2016 (Squamous,
#' Immunogenic, Pancreatic Progenitor, ADEX; genes flagged NotUnique in the
#' source table are not assigned a single subtype and are excluded here).
#'
#' @param pdac_data_root path to PDAC_DATA_ROOT (must contain original/cmbSubtypes.RData)
#' @return named list of character vectors: Moffitt_BasalLike, Moffitt_Classical,
#'   Bailey_Squamous, Bailey_Immunogenic, Bailey_PancreaticProgenitor, Bailey_ADEX
load_moffitt_bailey_genesets <- function(pdac_data_root) {
  cmb_path <- file.path(pdac_data_root, "original", "cmbSubtypes.RData")
  if (!file.exists(cmb_path)) {
    stop("load_moffitt_bailey_genesets: cmbSubtypes.RData not found at ", cmb_path)
  }

  e <- new.env()
  load(cmb_path, envir = e)
  if (!all(c("schemaList", "subtypeGeneList") %in% ls(e))) {
    stop("load_moffitt_bailey_genesets: cmbSubtypes.RData does not contain ",
         "schemaList/subtypeGeneList as expected")
  }

  idx_mt <- which(e$schemaList == "MT")
  idx_bailey <- which(e$schemaList == "Bailey")
  if (length(idx_mt) != 1 || length(idx_bailey) != 1) {
    stop("load_moffitt_bailey_genesets: expected exactly one 'MT' and one 'Bailey' ",
         "schema in cmbSubtypes.RData$schemaList")
  }

  mt <- e$subtypeGeneList[[idx_mt]]
  bailey <- e$subtypeGeneList[[idx_bailey]]

  list(
    Moffitt_BasalLike = mt$geneSymbol[mt[["Basal-like"]]],
    Moffitt_Classical = mt$geneSymbol[mt[["Classical"]]],
    Bailey_Squamous = bailey$geneSymbol[bailey$Squamous],
    Bailey_Immunogenic = bailey$geneSymbol[bailey$Immunogenic],
    Bailey_PancreaticProgenitor = bailey$geneSymbol[bailey$PancreaticProgenitor],
    Bailey_ADEX = bailey$geneSymbol[bailey$ADEX]
  )
}

#' Assemble the full custom PDAC gene-set collection (Moffitt + Bailey + DeSurv).
#'
#' Orchestrates load_moffitt_bailey_genesets() and extract_desurv_genesets(),
#' checks every set is non-empty and that its gene symbols substantially
#' overlap the MSigDB human gene-symbol universe (a parsing-sanity check, not
#' a filter -- genes are kept even if unmapped, but unmapped counts are logged,
#' never silently dropped), writes the combined list plus a manifest recording
#' the source citation per set.
#'
#' @param pdac_data_root path to PDAC_DATA_ROOT
#' @param si_pdf_path    path to the DeSurv si_appendix.pdf
#' @param output_dir     directory to write pdac_genesets.rds + genesets_manifest.txt
#' @return named list of all 9 gene sets (invisibly; also written to disk)
build_pdac_genesets <- function(pdac_data_root, si_pdf_path, output_dir) {
  moffitt_bailey <- load_moffitt_bailey_genesets(pdac_data_root)
  desurv <- extract_desurv_genesets(si_pdf_path)
  desurv_labeled <- list(
    DeSurv_D1_ClassicalTumor = desurv$D1,
    DeSurv_D2_StromalImmune = desurv$D2,
    DeSurv_D3_BasalLikeTumor = desurv$D3
  )
  genesets <- c(moffitt_bailey, desurv_labeled)

  empty <- names(genesets)[vapply(genesets, length, integer(1)) == 0]
  if (length(empty) > 0) {
    stop("build_pdac_genesets: empty gene set(s): ", paste(empty, collapse = ", "))
  }

  msigdb_universe <- unique(msigdbr::msigdbr(species = "Homo sapiens")$gene_symbol)
  mapping_lines <- vapply(names(genesets), function(nm) {
    g <- genesets[[nm]]
    n_mapped <- sum(g %in% msigdb_universe)
    if (n_mapped == 0) {
      stop("build_pdac_genesets: gene set '", nm, "' has ZERO overlap with the ",
           "MSigDB gene-symbol universe -- likely a parsing failure, not real biology")
    }
    sprintf("%s (n=%d): %d/%d (%.0f%%) mapped to MSigDB universe",
            nm, length(g), n_mapped, length(g), 100 * n_mapped / length(g))
  }, character(1))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(genesets, file.path(output_dir, "pdac_genesets.rds"))

  citations <- c(
    Moffitt_BasalLike = "Moffitt RA et al. 2015. Nat Med 21(10):1168-78. (cmbSubtypes.RData, schema 'MT', column 'Basal-like')",
    Moffitt_Classical = "Moffitt RA et al. 2015. Nat Med 21(10):1168-78. (cmbSubtypes.RData, schema 'MT', column 'Classical')",
    Bailey_Squamous = "Bailey P et al. 2016. Nature 531(7592):47-52. (cmbSubtypes.RData, schema 'Bailey', column 'Squamous')",
    Bailey_Immunogenic = "Bailey P et al. 2016. Nature 531(7592):47-52. (cmbSubtypes.RData, schema 'Bailey', column 'Immunogenic')",
    Bailey_PancreaticProgenitor = "Bailey P et al. 2016. Nature 531(7592):47-52. (cmbSubtypes.RData, schema 'Bailey', column 'PancreaticProgenitor')",
    Bailey_ADEX = "Bailey P et al. 2016. Nature 531(7592):47-52. (cmbSubtypes.RData, schema 'Bailey', column 'ADEX')",
    DeSurv_D1_ClassicalTumor = "Young AM et al. DeSurv (Young et al., PNAS 2026). si_appendix.pdf, Table S7.",
    DeSurv_D2_StromalImmune = "Young AM et al. DeSurv (Young et al., PNAS 2026). si_appendix.pdf, Table S8.",
    DeSurv_D3_BasalLikeTumor = "Young AM et al. DeSurv (Young et al., PNAS 2026). si_appendix.pdf, Table S9."
  )

  manifest <- c(
    "PDAC custom gene-set collection -- pdac_genesets.rds",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("R version:", R.version.string),
    paste("msigdbr version:", as.character(utils::packageVersion("msigdbr"))),
    "",
    "Per-set source citation and MSigDB-universe mapping check:",
    paste0(" - ", citations[names(genesets)]),
    "",
    mapping_lines
  )
  writeLines(manifest, file.path(output_dir, "genesets_manifest.txt"))

  invisible(genesets)
}

# Section: MSigDB collections ----

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Fetch MSigDB gene-set collections as gene-symbol lists suitable for fgsea/ORA.
#'
#' Default collections match the plan's methodology: Hallmark (compact,
#' interpretable), Reactome and KEGG (pathway-level detail), GO Biological
#' Process (broad coverage). Uses msigdbr's `collection`/`subcollection`
#' arguments (msigdbr >= 25: `category`/`subcategory` are deprecated aliases).
#'
#' @param species     passed to msigdbr() (default "human")
#' @param collections named list of msigdbr collection specs, each a list with
#'                     `collection` and optionally `subcollection`
#' @return named list (same names as `collections`), each element a named list
#'   mapping gene-set name -> character vector of gene symbols
get_msigdb_collections <- function(species = "human", collections = list(
  Hallmark = list(collection = "H"),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  GO_BP = list(collection = "C5", subcollection = "GO:BP")
)) {
  lapply(collections, function(spec) {
    df <- if (is.null(spec$subcollection)) {
      msigdbr::msigdbr(species = species, collection = spec$collection)
    } else {
      msigdbr::msigdbr(species = species, collection = spec$collection,
                        subcollection = spec$subcollection)
    }
    if (nrow(df) == 0) {
      stop("get_msigdb_collections: zero gene sets returned for collection=",
           spec$collection, " subcollection=", spec$subcollection %||% "NULL")
    }
    split(df$gene_symbol, df$gs_name)
  })
}

#' Top-N weighted genes per program, formalized as a tidy table (T2).
#'
#' @param EF          gene x program weight matrix (rows named by gene symbol)
#' @param programs    integer vector of program indices to include
#' @param program_labels named list mapping program index (as string) to label
#' @param n           number of top genes to report per program
#' @return data.frame: program, label, rank, gene, weight
top_n_genes_table <- function(EF, programs, program_labels, n = 50) {
  rows <- lapply(programs, function(k) {
    w <- EF[, k]
    ord <- order(w, decreasing = TRUE)[seq_len(min(n, length(w)))]
    data.frame(
      program = k,
      label = program_labels[[as.character(k)]],
      rank = seq_along(ord),
      gene = rownames(EF)[ord],
      weight = w[ord],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Section: Figures (F1-F3) ----

#' Prepare top-set-per-program data for the F1 enrichment dot plot.
#'
#' Pure data-prep step, factored out of plot_enrichment_dotplot() so the
#' filtering/ranking logic is unit-testable without rendering a plot.
#'
#' @param fgsea_results data.frame as returned by run_fgsea_program() (with a
#'   `program` column added, e.g. by run_pathway_enrichment.R)
#' @param programs      integer vector of program indices to include
#' @param top_n         number of top (lowest-padj) sets to keep per program
#' @return data.frame with an added neglog10padj column, one row per (program, set)
prepare_dotplot_data <- function(fgsea_results, programs, top_n = 10) {
  sub <- fgsea_results[fgsea_results$program %in% programs, ]
  sub <- sub[order(sub$program, sub$padj), ]
  sub <- do.call(rbind, lapply(programs, function(k) {
    rows_k <- sub[sub$program == k, ]
    rows_k[seq_len(min(top_n, nrow(rows_k))), ]
  }))
  sub$neglog10padj <- -log10(pmax(sub$padj, .Machine$double.xmin))
  sub
}

#' F1: enrichment dot plot, top sets per program (x=NES, size=set size, color=-log10(padj)).
#'
#' @param dotplot_data output of prepare_dotplot_data()
#' @return a ggplot object
plot_enrichment_dotplot <- function(dotplot_data) {
  dotplot_data$set_label <- paste0(dotplot_data$set, " (", dotplot_data$collection, ")")
  dotplot_data$set_label <- factor(dotplot_data$set_label,
                                    levels = rev(unique(dotplot_data$set_label)))
  ggplot2::ggplot(dotplot_data, ggplot2::aes(x = NES, y = set_label, size = size, color = neglog10padj)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~label, scales = "free_y") +
    ggplot2::scale_color_viridis_c(name = "-log10(padj)") +
    ggplot2::labs(x = "NES", y = NULL, size = "Set size",
                  title = "Top enriched gene sets by program") +
    ggplot2::theme_bw()
}

#' F2: fgsea running-enrichment-score plot for one headline gene set.
#'
#' @param weights_k named numeric vector of gene weights for one program
#' @param geneset   character vector of gene symbols (the gene set to plot)
#' @param title     plot title (typically "<set name> (Program <k>, <label>)")
#' @return a ggplot object (fgsea::plotEnrichment output)
plot_running_es <- function(weights_k, geneset, title = "") {
  fgsea::plotEnrichment(geneset, weights_k) + ggplot2::labs(title = title)
}

#' F3: gene-weight heatmap for a set of genes across all 7 programs.
#'
#' @param EF    gene x program weight matrix (rows named by gene symbol)
#' @param genes character vector of gene symbols to include (e.g. leading-edge
#'              genes from headline sets); rows not found in EF are dropped
#'              (logged via message(), not silently ignored)
#' @param program_labels named list mapping program index (as string) to label
#' @param filename if provided, written directly to this path (pheatmap's own
#'                 file-output mechanism); otherwise the pheatmap object is returned
#' @param cellheight_in row height in inches, scales the saved figure's height
#'                       and row-label font size so gene names stay readable
#'                       and non-overlapping regardless of how many genes are
#'                       plotted (rather than a fixed figure size for all N)
#' @param title    plot title; default matches the original all-7-programs caller.
#'                 Callers plotting a subset of programs should override this so
#'                 the figure doesn't claim to show more columns than it does.
#' @return the pheatmap grob (invisibly if filename is given)
plot_geneweight_heatmap <- function(EF, genes, program_labels, filename = NA, cellheight_in = 0.14,
                                     title = "Gene weights (EF) across all 7 programs") {
  present <- intersect(genes, rownames(EF))
  missing <- setdiff(genes, rownames(EF))
  if (length(missing) > 0) {
    message(sprintf("plot_geneweight_heatmap: %d gene(s) not found in EF and dropped: %s",
                     length(missing), paste(missing, collapse = ", ")))
  }
  if (length(present) == 0) {
    stop("plot_geneweight_heatmap: none of the requested genes were found in EF")
  }

  mat <- EF[present, , drop = FALSE]
  colnames(mat) <- paste0("P", seq_len(ncol(mat)), "_", vapply(seq_len(ncol(mat)),
                                                                function(k) program_labels[[as.character(k)]],
                                                                character(1)))

  dims <- heatmap_dimensions(nrow(mat), ncol(mat), cellheight_in)

  pheatmap::pheatmap(mat, cluster_cols = FALSE, filename = filename,
                      fontsize_row = dims$fontsize_row,
                      width = dims$width_in, height = dims$height_in,
                      main = title)
}

#' Compute row-label font size and saved-figure dimensions for the F3 heatmap.
#'
#' Pure sizing logic, factored out of plot_geneweight_heatmap() so it's
#' unit-testable without rendering: more genes -> smaller font, taller figure,
#' so labels stay legible and non-overlapping regardless of gene count
#' (a fixed figure size was found to squeeze >100 gene labels illegibly).
#'
#' @param n_genes integer, number of rows (genes) to plot
#' @param n_cols  integer, number of columns (programs)
#' @param cellheight_in row height in inches
#' @return list(fontsize_row, width_in, height_in)
heatmap_dimensions <- function(n_genes, n_cols, cellheight_in = 0.14) {
  list(
    fontsize_row = max(4, min(10, 500 / n_genes)),
    width_in = max(6, n_cols * 1.1 + 3),
    height_in = max(4, n_genes * cellheight_in + 2)
  )
}

# Section: PDAC subtype concordance (Step 7) ----

#' Merge D4 patient loadings (EL) with per-sample Moffitt/PurIST subtype calls.
#'
#' Uses PurIST (categorical Basal-like/Classical) and PurIST.prob (continuous
#' basal-likelihood score), NOT "MS"/"MS_K2" -- those columns in the local
#' reference data represent the Moffitt STROMA activation axis (Activated vs
#' Normal), a different biological question, not the tumor basal/classical
#' axis this concordance check targets. (This corrects an assumption in the
#' original 2026-06-16 plan draft, which conflated the two.)
#'
#' @param EL          patient x program loading matrix (rows in the same order
#'                    as `sample_ids`, e.g. d4$EL[1:n_tcga, ] for the TCGA-first
#'                    block of the pooled D4 training set)
#' @param sample_ids  character vector, TCGA barcodes, same length/order as EL's rows
#' @param subtype_df  data.frame with columns sampID, PurIST, PurIST.prob
#'                    (e.g. readRDS(".../TCGA_PAAD.caf_subtype.rds")$Subtype)
#' @param min_match_frac minimum fraction of sample_ids that must match subtype_df$sampID;
#'                        stop() if not met (fail loud, per project convention)
#' @return data.frame: sampID, one column per EL program (named EL_<k>), PurIST, PurIST.prob
merge_loadings_with_subtype <- function(EL, sample_ids, subtype_df, min_match_frac = 0.80) {
  if (nrow(EL) != length(sample_ids)) {
    stop("merge_loadings_with_subtype: nrow(EL) must equal length(sample_ids)")
  }
  el_df <- as.data.frame(EL)
  colnames(el_df) <- paste0("EL_", seq_len(ncol(EL)))
  el_df$sampID <- sample_ids

  merged <- merge(el_df, subtype_df[, c("sampID", "PurIST", "PurIST.prob")], by = "sampID")

  match_frac <- nrow(merged) / length(sample_ids)
  cat(sprintf("merge_loadings_with_subtype: matched %d/%d samples (%.1f%%)\n",
              nrow(merged), length(sample_ids), 100 * match_frac))
  if (match_frac < min_match_frac) {
    stop(sprintf("merge_loadings_with_subtype: only %.1f%% of samples matched a subtype ",
                  100 * match_frac), "call -- below the ", 100 * min_match_frac, "% threshold")
  }
  merged
}

#' Spearman correlation (continuous) + Kruskal-Wallis test (categorical) of
#' program loadings against the PurIST basal/classical axis.
#'
#' @param matched_df output of merge_loadings_with_subtype()
#' @param programs   integer vector of program indices (matching EL_<k> column names)
#' @return data.frame (T3): program, spearman_rho, spearman_p, kruskal_stat, kruskal_p, n
compute_subtype_concordance <- function(matched_df, programs) {
  rows <- lapply(programs, function(k) {
    col <- paste0("EL_", k)
    sp <- suppressWarnings(cor.test(matched_df[[col]], matched_df$PurIST.prob, method = "spearman"))
    kw <- kruskal.test(matched_df[[col]], as.factor(matched_df$PurIST))
    data.frame(
      program = k,
      spearman_rho = unname(sp$estimate),
      spearman_p = sp$p.value,
      kruskal_stat = unname(kw$statistic),
      kruskal_p = kw$p.value,
      n = nrow(matched_df)
    )
  })
  do.call(rbind, rows)
}

#' F4: violin plot of a program's patient loading by PurIST subtype, with test p-value.
#'
#' @param matched_df output of merge_loadings_with_subtype()
#' @param program    integer program index
#' @param program_label character label (e.g. "Adverse") for the plot title
#' @return a ggplot object
plot_loading_by_subtype <- function(matched_df, program, program_label = "") {
  col <- paste0("EL_", program)
  kw <- kruskal.test(matched_df[[col]], as.factor(matched_df$PurIST))
  ggplot2::ggplot(matched_df, ggplot2::aes(x = PurIST, y = .data[[col]], fill = PurIST)) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.4) +
    ggplot2::labs(x = "PurIST subtype", y = sprintf("Program %d loading", program),
                  title = sprintf("Program %d (%s) vs. PurIST (Kruskal-Wallis p=%.3g)",
                                   program, program_label, kw$p.value)) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")
}

# Section: External cohort robustness (Step 8) ----

#' Score a program's leading-edge gene signature in one cohort's expression data.
#'
#' Signature score = mean of per-gene z-scores across the leading-edge genes,
#' evaluated within this cohort (per-cohort standardization, matching the
#' project's per-platform z-std convention -- each cohort is its own platform).
#' Missing genes are logged (never silently zero-filled or dropped without
#' a trace), and the function fails loud only if NONE of the requested genes
#' are found (a genuine data problem, not a bad prediction).
#'
#' @param Y     cohort expression matrix, n patients x p genes
#' @param gene_names character vector, length p, matching Y's columns
#' @param leading_edge_genes character vector of gene symbols to score
#' @return list(score = numeric vector length n, n_genes_used, missing_genes)
score_leading_edge_signature <- function(Y, gene_names, leading_edge_genes) {
  colnames(Y) <- gene_names
  present <- intersect(leading_edge_genes, gene_names)
  missing <- setdiff(leading_edge_genes, gene_names)
  if (length(missing) > 0) {
    message(sprintf("score_leading_edge_signature: %d/%d leading-edge gene(s) not found in this cohort: %s",
                     length(missing), length(leading_edge_genes), paste(missing, collapse = ", ")))
  }
  if (length(present) == 0) {
    stop("score_leading_edge_signature: none of the leading-edge genes were found in this cohort")
  }
  Z <- scale(Y[, present, drop = FALSE])
  list(score = rowMeans(Z), n_genes_used = length(present), missing_genes = missing)
}

#' Cox model of survival on a single signature score, plus an oriented C-index.
#'
#' survival::concordance() treats a HIGHER predictor as predicting LONGER
#' survival by default (verified empirically: an unambiguously high-risk
#' synthetic score gave C-index ~0.01 unnegated vs ~0.99 negated). A fixed
#' negation is wrong here because this function is called for both Adverse
#' (higher score = higher risk) and Protective (higher score = lower risk)
#' program signatures -- negating unconditionally silently reports the
#' *wrong-direction* C-index for the Protective case (caught in code review:
#' a constructed strong, significant protective signature, HR=0.25,
#' p=1e-18, reported cindex=0.22, i.e. "worse than chance", when its true
#' discriminative accuracy in the correct direction is 1-0.22=0.78).
#'
#' Matches this project's existing oriented_cindex() convention (used in
#' run_desurv_comparison.R, run_k_parsimony_curve.R, and 6 other benchmark
#' scripts): report max(c_raw, 1-c_raw), the C-index in whichever direction
#' is more discriminative, which sidesteps needing to know a priori which
#' sign convention a given score/program follows.
#'
#' @param score  numeric vector, length n (signature score per patient)
#' @param time   numeric vector, length n
#' @param status integer vector, length n, in {0,1}
#' @return list(HR, p, cindex, n)
cohort_signature_cox <- function(score, time, status) {
  fit <- survival::coxph(survival::Surv(time, status) ~ score)
  s <- summary(fit)
  c_raw <- as.numeric(survival::concordance(survival::Surv(time, status) ~ score)$concordance)
  list(
    HR = unname(exp(coef(fit))),
    p = unname(s$coefficients[1, "Pr(>|z|)"]),
    cindex = max(c_raw, 1 - c_raw),
    n = length(score)
  )
}

# Section: SBMF vs DeSurv overlap (Step 9) ----

#' Jaccard overlap + hypergeometric enrichment p-value between two gene sets.
#'
#' Both sets are first restricted to `background` -- the hypergeometric model
#' assumes both are drawn from that universe, which does not hold in general
#' (e.g. DeSurv's own factor gene lists are not a subset of SBMF's 2064-gene
#' selected universe: ~245-259 of each 270-gene DeSurv list actually fall
#' inside it). Using the un-restricted set sizes as phyper()'s `m`/`n`/`k`
#' understates significance (verified: recomputing with m=249 instead of the
#' unrestricted 270 for Program 3 vs. DeSurv D1 changes p from ~2.4e-6 to
#' ~1.1e-7, a ~22x difference) -- caught in code review, fixed here.
#'
#' @param set_a, set_b character vectors of gene symbols
#' @param background   character vector of all gene symbols in the universe
#'                      both sets are conceptually drawn from (e.g. the 2064
#'                      D4 selected genes)
#' @return list(overlap_n, jaccard, hyper_p) -- computed on each set restricted to `background`
compute_geneset_overlap <- function(set_a, set_b, background) {
  set_a <- intersect(set_a, background)
  set_b <- intersect(set_b, background)
  overlap_n <- length(intersect(set_a, set_b))
  union_n <- length(union(set_a, set_b))
  jaccard <- if (union_n == 0) 0 else overlap_n / union_n
  # P(X >= overlap_n) under the hypergeometric null: drawing length(set_a)
  # genes without replacement from a background universe containing
  # length(set_b) "successes".
  hyper_p <- stats::phyper(overlap_n - 1, length(set_b), length(background) - length(set_b),
                            length(set_a), lower.tail = FALSE)
  list(overlap_n = overlap_n, jaccard = jaccard, hyper_p = hyper_p)
}

#' SBMF Programs 3 & 7 vs. DeSurv D1/D2/D3 gene-list overlap table (T4).
#'
#' Compares each active program's own top-N weighted genes (independent of
#' any particular enrichment test result, to avoid circularity with Step 6's
#' fgsea/ORA output) against each of the 3 DeSurv factor gene lists.
#'
#' @param EF          gene x program weight matrix (rows named by gene symbol)
#' @param program_labels named list mapping program index (as string) to label
#' @param desurv_genesets named list with DeSurv_D1_ClassicalTumor,
#'                         DeSurv_D2_StromalImmune, DeSurv_D3_BasalLikeTumor
#'                         (e.g. the relevant 3 elements of build_pdac_genesets()'s output)
#' @param programs    integer vector of program indices to compare
#' @param top_n       number of top-weighted genes per program to compare (270,
#'                     matching DeSurv's own top-N-per-factor convention, for
#'                     an apples-to-apples comparison)
#' @param background  character vector of all gene symbols in the universe
#'                     (default: rownames(EF), the full 2064-gene D4 universe)
#' @return data.frame (T4): program, label, desurv_factor, overlap_n, jaccard, hyper_p
sbmf_desurv_overlap_table <- function(EF, program_labels, desurv_genesets, programs,
                                       top_n = 270, background = rownames(EF)) {
  rows <- lapply(programs, function(k) {
    top_genes <- top_n_genes_table(EF, k, program_labels, n = top_n)$gene
    do.call(rbind, lapply(names(desurv_genesets), function(dname) {
      ov <- compute_geneset_overlap(top_genes, desurv_genesets[[dname]], background)
      data.frame(
        program = k, label = program_labels[[as.character(k)]], desurv_factor = dname,
        overlap_n = ov$overlap_n, jaccard = ov$jaccard, hyper_p = ov$hyper_p,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

#' F5: SBMF x DeSurv leading-edge/top-gene Jaccard overlap heatmap.
#'
#' @param overlap_table output of sbmf_desurv_overlap_table()
#' @return a ggplot object
plot_sbmf_desurv_overlap <- function(overlap_table) {
  overlap_table$program_label <- sprintf("P%d (%s)", overlap_table$program, overlap_table$label)
  ggplot2::ggplot(overlap_table, ggplot2::aes(x = desurv_factor, y = program_label, fill = jaccard)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f\np=%.1e", jaccard, hyper_p)), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, NA)) +
    ggplot2::labs(x = "DeSurv factor", y = NULL, fill = "Jaccard",
                  title = "SBMF vs. DeSurv gene-list overlap") +
    ggplot2::theme_bw()
}
