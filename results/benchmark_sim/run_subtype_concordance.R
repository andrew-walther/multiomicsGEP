# ============================================================
# Script:  results/benchmark_sim/run_subtype_concordance.R
# Purpose: Step 7 (pathway enrichment plan) -- concordance between the D4
#          model's Program 3/7 patient loadings and the Moffitt/PurIST tumor
#          basal-classical axis, in the TCGA_PAAD training cohort.
#
#          Uses PurIST (categorical Basal-like/Classical) and PurIST.prob
#          (continuous), NOT "MS"/"MS_K2" -- confirmed those columns in the
#          local reference data are the Moffitt STROMA activation axis
#          (Activated/Normal), not the tumor subtype axis this check targets.
#          This corrects an assumption in the original 2026-06-16 plan draft.
#
#          Bailey per-sample labels are NOT available (only the gene-set
#          schema is registered locally, no classifier output) -- deferred,
#          per the plan's own decision (D2).
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#     T3_concordance_stats.csv
#     F4_loading_vs_subtype_program3.png
#     F4_loading_vs_subtype_program7.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-14
# Usage:   export PDAC_DATA_ROOT=...
#          Rscript results/benchmark_sim/run_subtype_concordance.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")
suppressPackageStartupMessages(library(ggplot2))

pdac_root <- Sys.getenv("PDAC_DATA_ROOT")
if (nchar(pdac_root) == 0 || !dir.exists(pdac_root)) {
  stop("run_subtype_concordance: PDAC_DATA_ROOT not set or not found -- this ",
       "step requires local real data (TCGA_PAAD raw + subtype calls)")
}

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

d4 <- load_d4_weights()

# Section: recover TCGA_PAAD's row block of the pooled D4 training set ----
# D4's EL rows are cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))
# (results/benchmark_sim/run_desurv_comparison.R) -- TCGA_PAAD is listed first
# in cfg$pdac$training_cohorts, so EL[1:n_tcga, ] is the TCGA_PAAD block, in
# load_pdac_raw()'s own sample order.

train_cohorts <- cfg$pdac$training_cohorts
if (train_cohorts[1] != "TCGA_PAAD") {
  stop("run_subtype_concordance: expected TCGA_PAAD first in cfg$pdac$training_cohorts, got: ",
       paste(train_cohorts, collapse = ", "))
}

tcga_raw <- load_pdac_raw("TCGA_PAAD", pdac_root)
n_tcga <- tcga_raw$n
cat(sprintf("TCGA_PAAD: n=%d (training-set block, rows 1:%d of D4's pooled EL)\n", n_tcga, n_tcga))

if (n_tcga > nrow(d4$EL)) {
  stop(sprintf("run_subtype_concordance: n_tcga (%d) exceeds nrow(EL) (%d) -- ",
               n_tcga, nrow(d4$EL)), "training-set assembly assumption is violated")
}
EL_tcga <- d4$EL[seq_len(n_tcga), , drop = FALSE]

# Section: load PurIST subtype calls and merge ----

subtype_obj <- readRDS(file.path(pdac_root, "original", "TCGA_PAAD.caf_subtype.rds"))
subtype_df <- subtype_obj$Subtype

matched <- merge_loadings_with_subtype(EL_tcga, tcga_raw$sampID, subtype_df, min_match_frac = 0.80)

# Section: T3 -- concordance stats, Programs 3 & 7 ----

# Extended to all 4 kept factors (survival-active + genomics-only, DECISIONS.md
# 2026-08-19). For the genomics-only pair this asks whether they track the
# classical/basal subtype axis at all, which the 07/15 report never tested.
ACTIVE_PROGRAMS <- d4$kept_factors
t3 <- compute_subtype_concordance(matched, programs = ACTIVE_PROGRAMS)
t3$label <- vapply(t3$program, function(k) d4$program_labels[[as.character(k)]], character(1))
write.csv(t3, file.path(OUT_DIR, "T3_concordance_stats.csv"), row.names = FALSE)
cat("T3 (concordance with PurIST basal/classical axis):\n")
print(t3)

# Directional prediction (stated before computing, per plan Step 7): Program 7
# (Adverse) should track the basal-like axis; Program 3 (Protective) the
# classical axis. PurIST.prob is a basal-likelihood score, so a positive
# spearman_rho for Program 7 and/or negative for Program 3 would confirm this.

# Section: F4 -- violin plots ----

for (k in ACTIVE_PROGRAMS) {
  label_k <- d4$program_labels[[as.character(k)]]
  f4 <- plot_loading_by_subtype(matched, program = k, program_label = label_k)
  ggsave(file.path(OUT_DIR, sprintf("F4_loading_vs_subtype_program%d.png", k)),
         f4, width = 6, height = 5, dpi = 120)
  cat(sprintf("F4 written for program %d (%s).\n", k, label_k))
}

cat("\nrun_subtype_concordance.R complete. Outputs written to", OUT_DIR, "\n")
