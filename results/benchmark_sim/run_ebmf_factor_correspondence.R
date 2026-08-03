# ============================================================
# Script:  results/benchmark_sim/run_ebmf_factor_correspondence.R
# Purpose: ROADMAP.md "A/B comparison: SSBMF vs unsupervised EBMF" (line 359),
#          Part 2 -- do the same 2 survival-active programs (Program 3
#          Protective, Program 7 Adverse) emerge from EBMF run completely
#          unsupervised, with no survival term at all? Part 1 (external
#          C-index, YFB vs EBMF->Cox) was already done in
#          run_ebmf_cox_external.R; this reuses that fit's gene-weight
#          matrix rather than re-fitting.
#
#          Both fits select genes via the identical DeSurv-aligned
#          combined_rank/top-3000-per-cohort procedure on the same merged
#          TCGA_PAAD + CPTAC training data, so their gene universes are
#          identical and in the same order -- this is checked, not assumed.
#          Correlates each of EBMF's 20 factors' gene loadings against D4's
#          Program 3 and Program 7 loadings; reports the best-matching EBMF
#          factor for each.
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#     T5_ebmf_factor_correspondence.csv       (best EBMF match per program)
#     T5_ebmf_factor_correlation_full.csv     (full 20 x 2 correlation matrix)
#     F6_ebmf_factor_correspondence_heatmap.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-03
# Usage:   Rscript results/benchmark_sim/run_ebmf_factor_correspondence.R
#          (requires results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit.rds
#           from a prior run_ebmf_cox_external.R run, and
#           results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds
#           from a prior run_desurv_comparison.R run)
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")
suppressPackageStartupMessages(library(ggplot2))

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
ACTIVE_PROGRAMS <- c(3, 7)

ebmf_path <- "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit.rds"
if (!file.exists(ebmf_path)) {
  stop("run_ebmf_factor_correspondence: ", ebmf_path, " not found -- run ",
       "run_ebmf_cox_external.R first (it fits and saves the unsupervised EBMF model)")
}

ebmf <- readRDS(ebmf_path)
d4 <- load_d4_weights()

# Section: verify identical gene universe/order (not assumed) ----

if (!identical(ebmf$train_genes, d4$gene_names)) {
  stop("run_ebmf_factor_correspondence: EBMF's train_genes and D4's gene_names ",
       "differ (in content or order) -- the two fits' gene-selection runs have ",
       "diverged since this script was written; correlating F_ebmf against D4$EF ",
       "directly would silently misalign genes. Re-derive the alignment before proceeding.")
}

# Section: correlate every EBMF factor against Programs 3 & 7 ----

full_cor <- stats::cor(ebmf$F_ebmf, d4$EF[, ACTIVE_PROGRAMS])
colnames(full_cor) <- sprintf("P%d_%s", ACTIVE_PROGRAMS,
                               vapply(ACTIVE_PROGRAMS, function(k) d4$program_labels[[as.character(k)]], character(1)))
rownames(full_cor) <- sprintf("EBMF_F%d", seq_len(nrow(full_cor)))

full_df <- data.frame(ebmf_factor = rownames(full_cor), full_cor, check.names = FALSE,
                       row.names = NULL)
write.csv(full_df, file.path(OUT_DIR, "T5_ebmf_factor_correlation_full.csv"), row.names = FALSE)

best <- do.call(rbind, lapply(seq_along(ACTIVE_PROGRAMS), function(i) {
  k <- ACTIVE_PROGRAMS[i]
  r <- full_cor[, i]
  j <- which.max(abs(r))
  data.frame(program = k, label = d4$program_labels[[as.character(k)]],
             best_ebmf_factor = j, r = r[j],
             second_best_abs_r = sort(abs(r), decreasing = TRUE)[2],
             stringsAsFactors = FALSE)
}))
write.csv(best, file.path(OUT_DIR, "T5_ebmf_factor_correspondence.csv"), row.names = FALSE)

cat("T5 (best-matching unsupervised EBMF factor per survival-active program):\n")
print(best)

# Section: F6 -- correlation heatmap, all 20 EBMF factors x Programs 3 & 7 ----

plot_df <- data.frame(
  ebmf_factor = factor(rownames(full_cor), levels = rownames(full_cor)),
  program = rep(colnames(full_cor), each = nrow(full_cor)),
  r = as.vector(full_cor)
)
f6 <- ggplot(plot_df, aes(x = program, y = ebmf_factor, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = "Unsupervised EBMF factor", fill = "r",
       title = "Unsupervised EBMF factors vs. SBMF's survival-active programs") +
  theme_bw()
ggsave(file.path(OUT_DIR, "F6_ebmf_factor_correspondence_heatmap.png"), f6, width = 6, height = 8, dpi = 120)
cat("F6 written.\n")

cat("\nrun_ebmf_factor_correspondence.R complete. Outputs written to", OUT_DIR, "\n")
