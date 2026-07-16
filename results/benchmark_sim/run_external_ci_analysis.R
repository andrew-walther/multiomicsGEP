# ============================================================
# Script:  results/benchmark_sim/run_external_ci_analysis.R
# Purpose: Bootstrap confidence intervals for the recommended model's (YFB,
#          DeSurv-aligned gene selection, K=7, no cohort indicator -- "D4")
#          external C-index per held-out cohort, and a paired-bootstrap test
#          of whether it is significantly more concordant than the EBMF->Cox
#          two-step baseline on the SAME patients.
#
#          Uses ALREADY-CACHED per-patient risk scores (no re-fitting):
#            results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_riskscores.rds
#            results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_riskscores.rds
#          Both were regenerated 2026-07-16 under the corrected (post-Phase-1c)
#          preprocessing -- see DECISIONS.md.
#
#          Per cohort: (1) a percentile-bootstrap CI on each model's own
#          C-index, (2) a PAIRED percentile-bootstrap CI on the difference
#          (YFB - EBMF->Cox), since both models are scored on the identical
#          patients in each cohort -- this correctly accounts for the
#          correlation between the two models' errors on the same patients,
#          unlike comparing two independent CIs. Also reports a fixed-effect
#          (inverse-variance-weighted) pooled difference across the 5
#          cohorts, clearly labeled as a simplifying assumption (no formal
#          cross-cohort heterogeneity test is performed here).
#
#   Output: results/benchmark_sim/outputs/desurv_comparison/
#     external_cindex_ci.csv        (per-model, per-cohort C-index + CI)
#     external_paired_diff_ci.csv   (per-cohort + pooled paired-diff CI)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-16
# Dependencies: survival
# Usage:   Rscript results/benchmark_sim/run_external_ci_analysis.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages(library(survival))
source("code/concordance_ci.R")

B    <- 2000L
SEED <- 1L

YFB_MODEL <- "D4"  # recommended configuration
OUT_DIR   <- "results/benchmark_sim/outputs/desurv_comparison"

yfb_scores  <- readRDS(file.path(OUT_DIR, "desurv_comparison_riskscores.rds"))[[YFB_MODEL]]
ebmf_scores <- readRDS("results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_riskscores.rds")

cohorts <- intersect(names(yfb_scores), names(ebmf_scores))
if (length(cohorts) == 0)
  stop("No cohorts in common between the YFB and EBMF->Cox risk-score caches.")

cat(sprintf("--- Bootstrap C-index CIs (B=%d, seed=%d) ---\n", B, SEED))

ci_rows   <- list()
diff_rows <- list()

for (co in cohorts) {
  y <- yfb_scores[[co]]
  e <- ebmf_scores[[co]]

  if (!isTRUE(all.equal(y$time, e$time)) || !isTRUE(all.equal(y$status, e$status)))
    stop(sprintf("%s: YFB and EBMF->Cox time/status do not match -- ", co),
         "these must be the identical patients for a paired comparison.")

  cat(sprintf("  %s (n=%d) ...\n", co, y$n))

  ci_y <- bootstrap_concordance_ci(y$risk, y$time, y$status, B = B, seed = SEED)
  ci_e <- bootstrap_concordance_ci(e$risk, e$time, e$status, B = B, seed = SEED)

  ci_rows[[length(ci_rows) + 1]] <- data.frame(
    cohort = co, model = "YFB (recommended)",
    c_index = round(ci_y$estimate, 4),
    ci_lo = round(ci_y$lower, 4), ci_hi = round(ci_y$upper, 4),
    se = round(ci_y$se, 4), n = y$n, stringsAsFactors = FALSE)
  ci_rows[[length(ci_rows) + 1]] <- data.frame(
    cohort = co, model = "EBMF->Cox (2-step)",
    c_index = round(ci_e$estimate, 4),
    ci_lo = round(ci_e$lower, 4), ci_hi = round(ci_e$upper, 4),
    se = round(ci_e$se, 4), n = e$n, stringsAsFactors = FALSE)

  diff <- bootstrap_concordance_diff_ci(y$risk, e$risk, y$time, y$status,
                                         B = B, seed = SEED)
  diff_rows[[length(diff_rows) + 1]] <- data.frame(
    cohort = co, diff_estimate = round(diff$estimate, 4),
    diff_ci_lo = round(diff$lower, 4), diff_ci_hi = round(diff$upper, 4),
    se = round(diff$se, 4), significant = diff$significant, n = y$n,
    stringsAsFactors = FALSE)
}

ci_table   <- do.call(rbind, ci_rows)
diff_table <- do.call(rbind, diff_rows)

# Fixed-effect (inverse-variance-weighted) pooled difference across cohorts.
# NOTE: this assumes a single common effect size across cohorts -- a
# simplifying assumption, not a formal test of cross-cohort heterogeneity
# (e.g. Cochran's Q). Reported as a clearly-labeled summary line alongside
# the per-cohort CIs, not in place of them.
w        <- 1 / diff_table$se^2
pooled_d <- sum(w * diff_table$diff_estimate) / sum(w)
pooled_se <- sqrt(1 / sum(w))
pooled_lo <- pooled_d - 1.96 * pooled_se
pooled_hi <- pooled_d + 1.96 * pooled_se
diff_table <- rbind(diff_table, data.frame(
  cohort = "POOLED (fixed-effect, weighted by 1/SE^2)",
  diff_estimate = round(pooled_d, 4),
  diff_ci_lo = round(pooled_lo, 4), diff_ci_hi = round(pooled_hi, 4),
  se = round(pooled_se, 4), significant = (pooled_lo > 0) || (pooled_hi < 0),
  n = sum(diff_table$n), stringsAsFactors = FALSE))

write.csv(ci_table, file.path(OUT_DIR, "external_cindex_ci.csv"), row.names = FALSE)
write.csv(diff_table, file.path(OUT_DIR, "external_paired_diff_ci.csv"), row.names = FALSE)

cat("\n=== Per-model C-index with 95% bootstrap CI ===\n")
print(ci_table, row.names = FALSE)

cat("\n=== Paired difference (YFB - EBMF->Cox), 95% bootstrap CI ===\n")
print(diff_table, row.names = FALSE)

cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "external_cindex_ci.csv")))
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "external_paired_diff_ci.csv")))
