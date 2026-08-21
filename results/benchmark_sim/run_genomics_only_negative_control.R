# ============================================================
# Script:  results/benchmark_sim/run_genomics_only_negative_control.R
# Purpose: Disambiguate the external-cohort negative control for the 2
#          genomics-only programs (5, 6) of the recommended D4 K=7 fit.
#
#          run_external_cohort_robustness.R scores each program's *headline*
#          signature, defined as its single lowest-padj enriched set. For
#          Program 6 that set is DeSurv_D1_ClassicalTumor -- the same
#          tumor-subtype gene list that serves as Program 3's (protective)
#          signature. Scoring it therefore measures the classical-tumor axis's
#          prognostic value, NOT whether Program 6 itself carries survival
#          signal, and Program 6 duly came back protective-leaning (HR<1 in 4/5
#          cohorts, 2 significant). That is a confounded negative control.
#
#          This script re-runs the control for the genomics-only programs using
#          each one's top enriched set *excluding* the PDAC_custom collection,
#          i.e. excluding DeSurv's three published factor gene lists, which are
#          the known prognostic/subtype axes. What remains is the program's own
#          pathway character (Hallmark/Reactome/KEGG/GO-BP), which is what the
#          negative control is supposed to test.
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#             C4_genomics_only_negative_control.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-21
# Usage:   PDAC_DATA_ROOT=... Rscript results/benchmark_sim/run_genomics_only_negative_control.R
# Requires: local real data (5 external cohorts) via PDAC_DATA_ROOT
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R")) setwd("../../")

pdac_root <- Sys.getenv("PDAC_DATA_ROOT")
if (nchar(pdac_root) == 0 || !dir.exists(pdac_root)) {
  stop("run_genomics_only_negative_control: PDAC_DATA_ROOT unset or missing -- ",
       "this step requires local real data (5 external cohorts)")
}

# cfg must exist before benchmark_helpers.R is sourced (it reads cfg at load time).
cfg     <- yaml::read_yaml("config/globals.yml")
OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"

source("code/pathway_enrichment.R")
source("results/benchmark_sim/benchmark_helpers.R")

d4               <- load_d4_weights()
GENOMICS_ONLY    <- d4$genomics_only            # {5, 6}
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

# Collections whose sets are the known prognostic/subtype axes, excluded here so
# a genomics-only program is not scored using a borrowed tumor-subtype signature.
EXCLUDE_COLLECTIONS <- "PDAC_custom"

fgsea_results <- readRDS(file.path(OUT_DIR, "fgsea_results_all.rds"))
sig <- fgsea_results[fgsea_results$padj < 0.10 &
                     !(fgsea_results$collection %in% EXCLUDE_COLLECTIONS), ]

headline <- lapply(GENOMICS_ONLY, function(k) {
  rows <- sig[sig$program == k, ]
  rows <- rows[order(rows$padj), ]
  if (!nrow(rows)) {
    stop(sprintf(paste0("run_genomics_only_negative_control: program %d has no enriched set ",
                        "outside the excluded collection(s) {%s} -- cannot build an unconfounded ",
                        "signature for it"), k, paste(EXCLUDE_COLLECTIONS, collapse = ", ")))
  }
  row <- rows[1, ]
  list(program = k, label = d4$program_labels[[as.character(k)]],
       set = row$set, collection = row$collection, padj = row$padj,
       leading_edge = strsplit(row$leading_edge, ";")[[1]])
})

for (h in headline) {
  cat(sprintf("Program %d (%s): %s [%s], padj=%.4g, %d leading-edge genes\n",
              h$program, h$label, h$set, h$collection, h$padj, length(h$leading_edge)))
}

results <- do.call(rbind, lapply(headline, function(h) {
  do.call(rbind, lapply(EXTERNAL_COHORTS, function(cohort) {
    tryCatch({
      raw   <- load_pdac_raw(cohort, pdac_root)
      score <- score_leading_edge_signature(raw$Y, raw$gene_names, h$leading_edge)
      cox   <- cohort_signature_cox(score$score, raw$time, raw$status)
      data.frame(program = h$program, label = h$label, set = h$set,
                 collection = h$collection, cohort = cohort, n = cox$n,
                 n_genes_used = score$n_genes_used,
                 HR = cox$HR, p = cox$p, cindex = cox$cindex)
    }, error = function(e) {
      message(sprintf("run_genomics_only_negative_control: skipped program %d / %s (%s)",
                      h$program, cohort, conditionMessage(e)))
      NULL
    })
  }))
}))

write.csv(results, file.path(OUT_DIR, "C4_genomics_only_negative_control.csv"), row.names = FALSE)

cat("\nC4 (genomics-only negative control, DeSurv subtype sets excluded):\n")
print(results, row.names = FALSE, digits = 3)

# Fail loud if a cohort silently dropped out.
expected <- length(GENOMICS_ONLY) * length(EXTERNAL_COHORTS)
if (nrow(results) < expected) {
  warning(sprintf("run_genomics_only_negative_control: %d of %d program x cohort fits produced no row",
                  expected - nrow(results), expected))
}

for (k in GENOMICS_ONLY) {
  sub <- results[results$program == k, ]
  cat(sprintf("\nProgram %d: %d/%d cohorts significant at p<0.05; HR range %.2f-%.2f; %d/%d with HR>1\n",
              k, sum(sub$p < 0.05), nrow(sub), min(sub$HR), max(sub$HR),
              sum(sub$HR > 1), nrow(sub)))
}

cat("\nrun_genomics_only_negative_control.R complete.\n")
