# ============================================================
# Script:  results/multi_cohort_sim/run_top2_match_diagnostic.R
# Purpose: Stage 5 leftover. match_factors() (sim_scoring.R) reports each
#          TRUE factor's single BEST-matching estimated factor -- exactly
#          what the recovery-correlation summary needs, but it hides
#          merging: if two true factors both correlate highly with the
#          SAME one estimated factor, match_factors() only shows that from
#          the true-factor side (each true factor still gets its own
#          "best" score), not from the estimated-factor side. The
#          advisors' question ("do multiple true factors merge into one
#          estimated factor?") is answered by looking at each ESTIMATED
#          factor's top-2 correlations with the true factors -- if an
#          estimated factor correlates highly with two different true
#          factors at once, that is a merge signature.
#
#          No re-fitting: reuses the first-seed example fits already saved
#          by run_multicohort_sim.R (results/multi_cohort_sim/outputs/
#          multicohort_sim_example.rds), for the YFB_base arm at every
#          K_init setting and scenario -- merging is a property of the
#          fit's own factor structure, not something the two-step/
#          unsupervised EBMF arm needs a separate check for here.
#
#   Output: results/multi_cohort_sim/outputs/top2_match_diagnostic.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/multi_cohort_sim/run_top2_match_diagnostic.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

EX_RDS <- "results/multi_cohort_sim/outputs/multicohort_sim_example.rds"
if (!file.exists(EX_RDS)) stop("Missing ", EX_RDS, " -- run run_multicohort_sim.R first.")
ex <- readRDS(EX_RDS)

OUT_DIR <- "results/multi_cohort_sim/outputs"

rows <- list()
for (k_setting in names(ex)) {
  for (scenario in names(ex[[k_setting]])) {
    entry <- ex[[k_setting]][[scenario]]
    if (is.null(entry$YFB_base) || is.null(entry$data)) next

    EF     <- entry$YFB_base$EF          # p x K_est
    F_true <- entry$data$F_true          # p x K_true
    labels_true <- entry$data$factor_labels

    Cmat <- abs(suppressWarnings(cor(EF, F_true)))  # K_est x K_true
    Cmat[is.na(Cmat)] <- 0

    for (est_k in seq_len(nrow(Cmat))) {
      ord <- order(Cmat[est_k, ], decreasing = TRUE)
      top1_true <- ord[1]; top2_true <- ord[2]
      top1_cor  <- Cmat[est_k, top1_true]
      top2_cor  <- Cmat[est_k, top2_true]
      # A merge signature: the estimated factor's SECOND-best true-factor
      # correlation is still substantial, not just noise trailing the best.
      merge_flag <- (top2_cor > 0.5) && (top1_true != top2_true)

      rows[[length(rows) + 1]] <- data.frame(
        k_setting = k_setting, scenario = scenario, estimated_factor = est_k,
        best_true_factor = labels_true[top1_true], best_cor = round(top1_cor, 3),
        second_true_factor = labels_true[top2_true], second_cor = round(top2_cor, 3),
        possible_merge = merge_flag,
        stringsAsFactors = FALSE
      )
    }
  }
}

result <- do.call(rbind, rows)
write.csv(result, file.path(OUT_DIR, "top2_match_diagnostic.csv"), row.names = FALSE)

cat("=== Under-specified K_init settings: possible merges (second_cor > 0.5) ===\n")
under_k <- result[result$k_setting %in% c("under_k2", "under_k3", "under_k4"), ]
print(under_k[under_k$possible_merge, ])
cat(sprintf("\n%d / %d estimated factors (under-specified K settings) show a possible merge signature.\n",
            sum(under_k$possible_merge), nrow(under_k)))

cat("\n=== Fully-specified/over-specified settings (oracle_k6, ard_k12, ard_k20): any merges? ===\n")
full_k <- result[result$k_setting %in% c("oracle_k6", "ard_k12", "ard_k20"), ]
print(full_k[full_k$possible_merge, ])
cat(sprintf("%d / %d show a possible merge signature.\n", sum(full_k$possible_merge), nrow(full_k)))

cat(sprintf("\nOutput: %s\n", file.path(OUT_DIR, "top2_match_diagnostic.csv")))
