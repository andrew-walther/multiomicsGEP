# ============================================================
# Script:  results/benchmark_sim/run_external_cohort_robustness.R
# Purpose: Step 8 (pathway enrichment plan, comparison C2) -- score each
#          active program's headline leading-edge gene signature (mean
#          within-cohort z-score) in each of the 5 held-out external PDAC
#          cohorts, and relate the signature score to survival (Cox HR,
#          C-index) per cohort.
#
#          Verify (per plan): adverse-program signature -> HR>1 (worse
#          survival) in a majority of cohorts; protective -> HR<1. Missing
#          genes per cohort are logged, never silently zero-filled.
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#     C2_external_cohort_robustness.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-14
# Usage:   export PDAC_DATA_ROOT=...
#          Rscript results/benchmark_sim/run_external_cohort_robustness.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")

pdac_root <- Sys.getenv("PDAC_DATA_ROOT")
if (nchar(pdac_root) == 0 || !dir.exists(pdac_root)) {
  stop("run_external_cohort_robustness: PDAC_DATA_ROOT not set or not found -- ",
       "this step requires local real data (5 external cohorts)")
}

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

d4 <- load_d4_weights()
ACTIVE_PROGRAMS <- c(3, 7)
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
cat(sprintf("External cohorts: %s\n", paste(EXTERNAL_COHORTS, collapse = ", ")))

# Section: headline leading-edge signature per active program ----
# Reuses Step 5/6's T1 (already-computed fgsea results at padj<0.10);
# takes each program's single lowest-padj set as its headline signature.

fgsea_results <- readRDS(file.path(OUT_DIR, "fgsea_results_all.rds"))
t1 <- fgsea_results[fgsea_results$program %in% ACTIVE_PROGRAMS & fgsea_results$padj < 0.10, ]
t1 <- t1[order(t1$program, t1$padj), ]

headline <- lapply(ACTIVE_PROGRAMS, function(k) {
  row <- t1[t1$program == k, ][1, ]
  list(program = k, label = d4$program_labels[[as.character(k)]],
       set = row$set, collection = row$collection,
       leading_edge = strsplit(row$leading_edge, ";")[[1]])
})
names(headline) <- paste0("program", ACTIVE_PROGRAMS)
for (h in headline) {
  cat(sprintf("Program %d (%s) headline signature: %s (%s), %d leading-edge genes\n",
              h$program, h$label, h$set, h$collection, length(h$leading_edge)))
}

# Section: score + Cox fit per cohort per program ----

results <- do.call(rbind, lapply(headline, function(h) {
  do.call(rbind, lapply(EXTERNAL_COHORTS, function(cohort) {
    tryCatch({
      raw <- load_pdac_raw(cohort, pdac_root)
      sig <- score_leading_edge_signature(raw$Y, raw$gene_names, h$leading_edge)
      cox <- cohort_signature_cox(sig$score, raw$time, raw$status)
      data.frame(
        program = h$program, label = h$label, set = h$set, cohort = cohort,
        n = cox$n, n_genes_used = sig$n_genes_used,
        n_genes_missing = length(sig$missing_genes),
        HR = cox$HR, p = cox$p, cindex = cox$cindex,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message(sprintf("run_external_cohort_robustness: cohort %s skipped for program %d (%s)",
                       cohort, h$program, conditionMessage(e)))
      NULL
    })
  }))
}))

write.csv(results, file.path(OUT_DIR, "C2_external_cohort_robustness.csv"), row.names = FALSE)
cat("\nC2 (external cohort robustness):\n")
print(results)

# Section: directional summary ----

for (h in headline) {
  sub <- results[results$program == h$program, ]
  if (h$label == "Adverse") {
    n_correct <- sum(sub$HR > 1)
  } else if (h$label == "Protective") {
    n_correct <- sum(sub$HR < 1)
  } else {
    n_correct <- NA
  }
  cat(sprintf("Program %d (%s): %d/%d cohorts have HR in the predicted direction\n",
              h$program, h$label, n_correct, nrow(sub)))
}

cat("\nrun_external_cohort_robustness.R complete. Outputs written to", OUT_DIR, "\n")
