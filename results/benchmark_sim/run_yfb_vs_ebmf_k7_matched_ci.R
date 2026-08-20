# ============================================================
# Script:  results/benchmark_sim/run_yfb_vs_ebmf_k7_matched_ci.R
# Purpose: Paired bootstrap CI on YFB (K=7) vs. two-step EBMF+Cox, matched at
#          the SAME K=7 (not the K=20 used in the original 2026-07-16
#          comparison -- an unmatched-K comparison that inflated the apparent
#          advantage to a significant +0.042; DECISIONS.md 2026-08-20).
#          No re-fitting: reuses cached per-cohort risk scores from both
#          run_desurv_comparison.R (YFB, D4) and
#          run_ebmf_cox_external.R --k 7 (the new K=7-matched EBMF+Cox run).
#
#   Output: results/benchmark_sim/outputs/ebmf_cox_external/yfb_vs_ebmf_k7_matched_ci.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   Rscript results/benchmark_sim/run_yfb_vs_ebmf_k7_matched_ci.R
# Requires: results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_riskscores_k7.rds
#           (run `Rscript results/benchmark_sim/run_ebmf_cox_external.R --k 7` first if missing)
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/concordance_ci.R")

d4   <- readRDS("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_riskscores.rds")[["D4"]]
ebmf <- readRDS("results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_riskscores_k7.rds")

cohorts <- names(d4)
rows <- list()
pooled_risk_yfb <- pooled_risk_ebmf <- pooled_time <- pooled_status <- c()

for (co in cohorts) {
  stopifnot(identical(d4[[co]]$time, ebmf[[co]]$time))
  ci <- bootstrap_concordance_diff_ci(d4[[co]]$risk, ebmf[[co]]$risk, d4[[co]]$time, d4[[co]]$status,
                                       B = 2000, seed = 1)
  cat(sprintf("%-20s n=%d diff(YFB-EBMF_K7)=%.4f [%.4f, %.4f] sig=%s\n",
              co, d4[[co]]$n, ci$estimate, ci$lower, ci$upper, ci$significant))

  rows[[length(rows) + 1]] <- data.frame(
    cohort = co, n = d4[[co]]$n, diff_estimate = round(ci$estimate, 4),
    diff_lower = round(ci$lower, 4), diff_upper = round(ci$upper, 4),
    significant = ci$significant, stringsAsFactors = FALSE
  )

  pooled_risk_yfb  <- c(pooled_risk_yfb, d4[[co]]$risk)
  pooled_risk_ebmf <- c(pooled_risk_ebmf, ebmf[[co]]$risk)
  pooled_time      <- c(pooled_time, d4[[co]]$time)
  pooled_status    <- c(pooled_status, d4[[co]]$status)
}

ci_pooled <- bootstrap_concordance_diff_ci(pooled_risk_yfb, pooled_risk_ebmf, pooled_time, pooled_status,
                                            B = 2000, seed = 1)
cat(sprintf("\nPOOLED (n=%d): diff(YFB-EBMF_K7)=%.4f [%.4f, %.4f] sig=%s\n",
            length(pooled_time), ci_pooled$estimate, ci_pooled$lower, ci_pooled$upper, ci_pooled$significant))

rows[[length(rows) + 1]] <- data.frame(
  cohort = "POOLED", n = length(pooled_time), diff_estimate = round(ci_pooled$estimate, 4),
  diff_lower = round(ci_pooled$lower, 4), diff_upper = round(ci_pooled$upper, 4),
  significant = ci_pooled$significant, stringsAsFactors = FALSE
)

results <- do.call(rbind, rows)
out_csv <- "results/benchmark_sim/outputs/ebmf_cox_external/yfb_vs_ebmf_k7_matched_ci.csv"
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nResults saved: %s\n", out_csv))
