# ============================================================
# Script:  results/benchmark_sim/generate_pathway_enrichment_figures.R
# Purpose: Generate F1-F3 from a completed run_pathway_enrichment.R run's
#          cached fgsea_results_all.rds + pdac_genesets.rds, without
#          re-running the expensive fgsea/ORA computation. Re-fetches the
#          (fast) MSigDB collections needed for F2's raw gene-set lookups.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-14
# Usage:   Rscript results/benchmark_sim/generate_pathway_enrichment_figures.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")
suppressPackageStartupMessages(library(ggplot2))

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
ACTIVE_PROGRAMS <- c(3, 7)

d4 <- load_d4_weights()
fgsea_results <- readRDS(file.path(OUT_DIR, "fgsea_results_all.rds"))
pdac_genesets <- readRDS(file.path(OUT_DIR, "pdac_genesets.rds"))
msigdb_collections <- get_msigdb_collections()
collections <- c(msigdb_collections, list(PDAC_custom = pdac_genesets))

t1 <- fgsea_results[fgsea_results$program %in% ACTIVE_PROGRAMS & fgsea_results$padj < 0.10, ]
t1 <- t1[order(t1$program, t1$padj),
         c("program", "label", "collection", "set", "size", "NES", "pval", "padj", "leading_edge")]

# F1
f1_data <- prepare_dotplot_data(fgsea_results, programs = ACTIVE_PROGRAMS, top_n = 10)
if (nrow(f1_data) > 0) {
  f1 <- plot_enrichment_dotplot(f1_data)
  ggsave(file.path(OUT_DIR, "F1_enrichment_dotplot.png"), f1, width = 16, height = 8, dpi = 150)
  ggsave(file.path(OUT_DIR, "F1_enrichment_dotplot.pdf"), f1, width = 16, height = 8)
  cat("F1 written.\n")
}

# F2: top 2 headline sets per active program
for (k in ACTIVE_PROGRAMS) {
  top_rows <- head(t1[t1$program == k, ], 2)
  if (nrow(top_rows) == 0) {
    message(sprintf("no headline sets for program %d -- skipping F2", k))
    next
  }
  for (i in seq_len(nrow(top_rows))) {
    row <- top_rows[i, ]
    geneset <- collections[[row$collection]][[row$set]]
    if (is.null(geneset)) {
      message(sprintf("F2: geneset '%s' not found in collection '%s' -- skipping", row$set, row$collection))
      next
    }
    f2 <- plot_running_es(d4$EF[, k], geneset,
                           title = sprintf("%s (Program %d, %s)", row$set, k, row$label))
    fname <- gsub("[^A-Za-z0-9_-]", "_", sprintf("F2_running_es_program%d_%s", k, row$set))
    ggsave(file.path(OUT_DIR, paste0(fname, ".png")), f2, width = 7, height = 4.5, dpi = 120)
    cat(sprintf("F2 written for program %d (%s, %s).\n", k, row$set, row$collection))
  }
}

# F3: union of top-20-per-set leading-edge genes from Programs 3 & 7's headline sets
f3_genes <- unique(unlist(lapply(strsplit(t1$leading_edge[t1$program %in% ACTIVE_PROGRAMS], ";"),
                                  function(g) g[seq_len(min(20, length(g)))])))
if (length(f3_genes) > 0) {
  plot_geneweight_heatmap(d4$EF, f3_genes, d4$program_labels,
                           filename = file.path(OUT_DIR, "F3_geneweight_heatmap.png"))
  cat(sprintf("F3 written (%d genes).\n", length(f3_genes)))
}

cat("Done.\n")
