# ============================================================
# Script:  results/benchmark_sim/run_sbmf_desurv_overlap.R
# Purpose: Step 9 (pathway enrichment plan, F5/T4) -- Jaccard overlap and
#          hypergeometric enrichment between SBMF's Programs 3 & 7 top-270
#          weighted genes and each of DeSurv's D1/D2/D3 factor gene lists
#          (270 genes each, from the SI appendix). Uses each program's own
#          top-N weighted genes (independent of Step 6's fgsea/ORA results)
#          to avoid circularity with that analysis.
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#     T4_sbmf_desurv_overlap.csv
#     F5_sbmf_vs_desurv.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-15
# Usage:   Rscript results/benchmark_sim/run_sbmf_desurv_overlap.R
#          (requires results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds
#           from a prior run_pathway_enrichment.R run)
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")
suppressPackageStartupMessages(library(ggplot2))

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"

genesets_path <- file.path(OUT_DIR, "pdac_genesets.rds")
if (!file.exists(genesets_path)) {
  stop("run_sbmf_desurv_overlap: ", genesets_path, " not found -- run ",
       "run_pathway_enrichment.R first (it builds and saves the DeSurv gene lists)")
}

d4 <- load_d4_weights()

# All 4 kept factors: the genomics-only pair is included because both showed
# coherent enrichment against DeSurv's StromalImmune factor (DECISIONS.md 2026-08-19).
ACTIVE_PROGRAMS <- d4$kept_factors
pdac_genesets <- readRDS(genesets_path)
desurv <- pdac_genesets[c("DeSurv_D1_ClassicalTumor", "DeSurv_D2_StromalImmune", "DeSurv_D3_BasalLikeTumor")]

# background = d4$gene_names (the 2064-gene D4 selected universe); DeSurv's own
# 270-gene-per-factor lists are NOT a subset of it (~245-259 of each fall inside),
# so compute_geneset_overlap() restricts both sets to this background internally.
t4 <- sbmf_desurv_overlap_table(d4$EF, d4$program_labels, desurv, programs = ACTIVE_PROGRAMS,
                                 top_n = 270, background = d4$gene_names)
write.csv(t4, file.path(OUT_DIR, "T4_sbmf_desurv_overlap.csv"), row.names = FALSE)
cat("T4 (SBMF vs DeSurv gene-list overlap):\n")
print(t4)

f5 <- plot_sbmf_desurv_overlap(t4)
ggsave(file.path(OUT_DIR, "F5_sbmf_vs_desurv.png"), f5, width = 8, height = 5, dpi = 120)
cat("F5 written.\n")

cat("\nrun_sbmf_desurv_overlap.R complete. Outputs written to", OUT_DIR, "\n")
