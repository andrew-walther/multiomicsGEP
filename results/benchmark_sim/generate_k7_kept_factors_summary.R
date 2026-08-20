# ============================================================
# Script:  results/benchmark_sim/generate_k7_kept_factors_summary.R
# Purpose: Gene-level summary of the recommended D4 model's 4 "kept" factors
#          at K=7 -- the 2 survival-active (Program 3 Protective, Program 7
#          Adverse) plus the 2 genomics-only factors (real gene-expression
#          programs, no survival effect) identified by classify_factors()
#          (DECISIONS.md 2026-08-19). Requested for the 8/21 progress-book
#          chapter: a gene table plus a heatmap in the style of Figure 3 of
#          docs/reports/pathway_enrichment_report_07_15_26.pdf, but scoped
#          to these 4 factors specifically (not all 7).
#
#          Kept-factor indices {3,5,6,7} are hardcoded from the verified
#          classify_factors() result on the D4 fit (same convention as
#          generate_pathway_enrichment_figures.R's ACTIVE_PROGRAMS <- c(3,7)):
#            factor 3: survival_active, EBeta=+0.0115  (Program 3, Protective)
#            factor 5: genomics_only,   EBeta~0
#            factor 6: genomics_only,   EBeta~0
#            factor 7: survival_active, EBeta=-0.0405  (Program 7, Adverse)
#
#          Top-N genes per factor are ranked by raw EF weight (point_exponential
#          prior => EF >= 0, so |EF|=EF) -- no fresh enrichment analysis for
#          the 2 genomics-only factors (kept as future work, see DECISIONS.md).
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#             K7_kept_factors_geneweight_heatmap.png
#             K7_kept_factors_top_genes.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   Rscript results/benchmark_sim/generate_k7_kept_factors_summary.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")   # load_d4_weights(), plot_geneweight_heatmap()

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N <- 20L
KEPT_FACTORS <- c(3L, 5L, 6L, 7L)
KEPT_LABELS  <- c("3" = "Program 3 (Protective)",
                   "5" = "Genomics-only (A)",
                   "6" = "Genomics-only (B)",
                   "7" = "Program 7 (Adverse)")

d4 <- load_d4_weights()

# --------------------------------------------------------------------------
# 1. Top-N genes per kept factor, by raw EF weight
# --------------------------------------------------------------------------
top_gene_rows <- list()
for (k in KEPT_FACTORS) {
  w <- d4$EF[, k]
  ord <- order(w, decreasing = TRUE)[seq_len(TOP_N)]
  top_gene_rows[[length(top_gene_rows) + 1]] <- data.frame(
    factor      = k,
    label       = KEPT_LABELS[[as.character(k)]],
    rank        = seq_len(TOP_N),
    gene        = d4$gene_names[ord],
    weight      = round(w[ord], 5),
    stringsAsFactors = FALSE
  )
}
top_genes <- do.call(rbind, top_gene_rows)
write.csv(top_genes, file.path(OUT_DIR, "K7_kept_factors_top_genes.csv"), row.names = FALSE)
cat(sprintf("Top-%d genes per factor written: %s\n", TOP_N,
            file.path(OUT_DIR, "K7_kept_factors_top_genes.csv")))

# --------------------------------------------------------------------------
# 2. Heatmap: union of top-N genes across the 4 kept factors, all 4 columns
#    (subset EF/program_labels to just the kept factors, remap to 1..4 so
#    plot_geneweight_heatmap()'s column-labeling logic applies cleanly)
# --------------------------------------------------------------------------
heatmap_genes <- unique(top_genes$gene)

EF_sub <- d4$EF[, KEPT_FACTORS, drop = FALSE]
labels_sub <- as.list(KEPT_LABELS[as.character(KEPT_FACTORS)])
names(labels_sub) <- as.character(seq_along(KEPT_FACTORS))

plot_geneweight_heatmap(
  EF_sub, heatmap_genes, labels_sub,
  filename = file.path(OUT_DIR, "K7_kept_factors_geneweight_heatmap.png"),
  title = sprintf("K=7 fit's 4 kept factors: top-%d genes each (3 fully-pruned factors omitted)", TOP_N)
)
cat(sprintf("Heatmap written (%d genes x %d factors): %s\n",
            length(heatmap_genes), length(KEPT_FACTORS),
            file.path(OUT_DIR, "K7_kept_factors_geneweight_heatmap.png")))
