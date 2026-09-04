# ============================================================
# Script:  results/benchmark_sim/run_factor_comparison.R
# Purpose: Stage 3b of the 9/4 plan: supervised-vs-unsupervised factor
#          comparison, no performance metric involved. Compares the
#          recommended model's (D4, K_init=7) 4 ARD-kept programs against
#          the unsupervised EBMF K=40 factors on gene-list overlap --
#          motivated by DeSurv (Young et al., PNAS) running the same
#          comparison, and covering manuscript Key Result #1 ("SSBMF vs
#          EBMF for recovering gene programs").
#
#          Unsupervised EBMF/NMF is NOT scored on C-index anywhere in this
#          project's method-comparison work (run_cohort_beta_comparison.R):
#          it has no survival coefficient. This script is the SEPARATE
#          analysis where it belongs -- a factor-identity comparison.
#
#          No re-fitting: reuses the cached D4 K=7 fit
#          (desurv_comparison_fits.rds) and the cached EBMF K=40 fit
#          (ebmf_cox_external_fit_k40.rds, the two-step baseline's own
#          stage-1 factorization).
#
#   Output: results/benchmark_sim/outputs/factor_comparison/
#             factor_comparison_overlap.csv
#             factor_comparison_summary.csv
#   Figure: docs/progress_book/figs/2026-09-04_factor_comparison_heatmap.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/benchmark_sim/run_factor_comparison.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(ggplot2) })
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")
source("code/select_K.R")
source("code/pathway_enrichment.R")   # compute_geneset_overlap(), top_n_genes_table()

OUT_DIR <- "results/benchmark_sim/outputs/factor_comparison"
FIG_DIR <- "docs/progress_book/figs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load the recommended D4 K=7 fit; identify the 4 ARD-kept programs.
# --------------------------------------------------------------------------
d4_fits <- readRDS("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds")
d4 <- d4_fits[["D4"]]
d4_gene_names <- readRDS("results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")
if (nrow(d4$EF) != length(d4_gene_names)) stop("D4 EF/gene_names length mismatch.")
EF_d4 <- d4$EF
rownames(EF_d4) <- d4_gene_names

BETA_THRESH <- cfg$k_selection$beta_threshold
PVE_THRESH  <- cfg$k_selection$pve_threshold

# classify_factors()/compute_pve() need the actual D4 training matrix (not
# derivable from EL/EF alone, since total_var = sum(Y^2) includes residual
# noise) -- rebuild it with the same D4 procedure used everywhere else in
# this session, no re-fitting.
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
pp <- preprocess_merged_cohorts(
  cohort_raw_list = train_raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank"
)
Y_train <- pp$Y
pve_d4 <- compute_pve(d4, Y_train)
ab_beta_d4 <- abs(d4$EBeta)
surv_active <- ab_beta_d4 > BETA_THRESH
geno_active <- pve_d4 > PVE_THRESH
category <- ifelse(surv_active, "survival_active", ifelse(geno_active, "genomics_only", "dead"))
kept_programs <- which(category != "dead")
program_labels <- setNames(as.list(category[kept_programs]), as.character(kept_programs))

cat(sprintf("D4 K_init=7 fit: %d ARD-kept programs (of %d): %s\n",
            length(kept_programs), ncol(d4$EF),
            paste(sprintf("P%d(%s)", kept_programs, category[kept_programs]), collapse = ", ")))

# --------------------------------------------------------------------------
# 2. Load the unsupervised EBMF K=40 fit (two-step baseline's own stage-1
#    factorization) -- same D4 gene universe (2064 genes, D4 preprocessing).
# --------------------------------------------------------------------------
ebmf <- readRDS("results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit_k40.rds")
EF_ebmf <- ebmf$F_ebmf
rownames(EF_ebmf) <- ebmf$train_genes
pve_ebmf <- ebmf$flash_fit$pve
variance_rank_ebmf <- rank(-pve_ebmf, ties.method = "min")   # 1 = highest variance

# --------------------------------------------------------------------------
# 3. Overlap: each kept SSBMF program's top-270 genes vs. each of the 40
#    EBMF factors' top-270 genes (same top_n as the existing DeSurv-overlap
#    precedent, run_sbmf_desurv_overlap.R).
# --------------------------------------------------------------------------
TOP_N <- 270
background <- intersect(rownames(EF_d4), rownames(EF_ebmf))
cat(sprintf("Common gene background: %d genes\n\n", length(background)))

overlap_rows <- list()
for (k in kept_programs) {
  top_sbmf <- top_n_genes_table(EF_d4, k, program_labels, n = TOP_N)$gene
  for (j in seq_len(ncol(EF_ebmf))) {
    top_e <- top_n_genes_table(EF_ebmf, j, setNames(as.list(rep("EBMF", ncol(EF_ebmf))), as.character(seq_len(ncol(EF_ebmf)))), n = TOP_N)$gene
    ov <- compute_geneset_overlap(top_sbmf, top_e, background)
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      sbmf_program = k, sbmf_category = category[k], ebmf_factor = j,
      ebmf_variance_rank = variance_rank_ebmf[j], ebmf_pve = round(pve_ebmf[j], 5),
      overlap_n = ov$overlap_n, jaccard = round(ov$jaccard, 4), hyper_p = ov$hyper_p,
      stringsAsFactors = FALSE
    )
  }
}
overlap <- do.call(rbind, overlap_rows)
write.csv(overlap, file.path(OUT_DIR, "factor_comparison_overlap.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# 4. Summary: each SSBMF program's BEST-matching EBMF factor, and its
#    variance rank. Also: which EBMF factors have no SSBMF counterpart at
#    all (i.e. never anyone's best match).
# --------------------------------------------------------------------------
summary_rows <- list()
best_matches <- integer(0)
for (k in kept_programs) {
  sub <- overlap[overlap$sbmf_program == k, ]
  best <- sub[which.max(sub$jaccard), ]
  best_matches <- c(best_matches, best$ebmf_factor)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    sbmf_program = k, sbmf_category = category[k],
    best_ebmf_factor = best$ebmf_factor, best_jaccard = best$jaccard,
    best_hyper_p = best$hyper_p, ebmf_variance_rank = best$ebmf_variance_rank,
    ebmf_variance_rank_of_40 = sprintf("%d/%d", best$ebmf_variance_rank, ncol(EF_ebmf)),
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "factor_comparison_summary.csv"), row.names = FALSE)

cat("=== Each kept SSBMF program's best-matching EBMF factor ===\n")
print(summary_df)

unmatched_ebmf <- setdiff(seq_len(ncol(EF_ebmf)), best_matches)
cat(sprintf("\nEBMF factors with NO SSBMF program's best match (%d/%d): %s\n",
            length(unmatched_ebmf), ncol(EF_ebmf),
            paste(head(sort(unmatched_ebmf), 10), collapse = ", ")))

# --------------------------------------------------------------------------
# 5. Heatmap figure: kept SSBMF programs x all 40 EBMF factors, Jaccard,
#    annotated with each EBMF factor's variance rank.
# --------------------------------------------------------------------------
overlap$program_label <- sprintf("P%d (%s)", overlap$sbmf_program,
                                  ifelse(overlap$sbmf_category == "survival_active", "surv", "geno"))
p <- ggplot(overlap, aes(x = reorder(factor(ebmf_factor), ebmf_variance_rank), y = program_label, fill = jaccard)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, NA)) +
  labs(x = "EBMF factor (ordered by variance rank, left = highest PVE)", y = NULL,
       fill = "Jaccard",
       title = "SSBMF kept programs vs. unsupervised EBMF (K=40) factors",
       caption = sprintf("Top-%d genes per factor; background = %d common genes.", TOP_N, length(background))) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
ggsave(file.path(FIG_DIR, "2026-09-04_factor_comparison_heatmap.png"), p, width = 12, height = 4.5, dpi = 170, bg = "white")

cat(sprintf("\nOutputs: %s, %s\n", file.path(OUT_DIR, "factor_comparison_overlap.csv"),
            file.path(OUT_DIR, "factor_comparison_summary.csv")))
