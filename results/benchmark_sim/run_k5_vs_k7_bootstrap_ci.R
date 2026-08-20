# ============================================================
# Script:  results/benchmark_sim/run_k5_vs_k7_bootstrap_ci.R
# Purpose: Bootstrap CI on the K=5-vs-K=7 external C-index gap flagged in
#          Analysis A (DECISIONS.md 2026-08-19): K=5's mean external C-index
#          (0.596) is far below K=7's (0.627) and consistent in direction
#          across all 5 held-out cohorts, but no formal significance test
#          had been run on this specific comparison. Reuses the
#          already-fitted best-of-multistart K=5/K=7 models (no re-fitting)
#          and code/concordance_ci.R's bootstrap_concordance_diff_ci()
#          (already used for a similar paired comparison, DECISIONS.md
#          2026-07-16).
#
#   Output: results/benchmark_sim/outputs/k_init_sweep/k5_vs_k7_bootstrap_ci.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   Rscript results/benchmark_sim/run_k5_vs_k7_bootstrap_ci.R
# Requires: PDAC_DATA_ROOT set (real data not in git).
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })
cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/predict_cox_on_yf.R")
source("code/preprocess_desurv.R")
source("code/concordance_ci.R")   # bootstrap_concordance_diff_ci()

TOP_N_DESURV <- cfg$preprocessing$top_n_genes_desurv

# --------------------------------------------------------------------------
# 1. Load the already-fitted K=5 and K=7 best-of-multistart models, and the
#    same training gene set used to fit them (recomputed once, deterministic
#    -- preprocess_merged_cohorts() re-derives the identical gene list from
#    raw data, no re-fitting involved).
# --------------------------------------------------------------------------
fits <- readRDS("results/benchmark_sim/outputs/k_init_sweep/k_init_multistart_best_fits.rds")
fit5 <- fits[["5"]]
fit7 <- fits[["7"]]

TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
pp <- preprocess_merged_cohorts(
  cohort_raw_list          = train_raw,
  log_transform_flags      = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n                    = TOP_N_DESURV,
  rank_transform           = FALSE,
  per_platform_standardize = TRUE,
  normalize_method         = "none",
  selection_per_cohort     = TRUE,
  selection_method         = "combined_rank"
)
train_genes <- pp$gene_names

# --------------------------------------------------------------------------
# 2. External cohorts: compute risk scores for both fits, per cohort
# --------------------------------------------------------------------------
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
rows <- list()
pooled_risk5 <- pooled_risk7 <- pooled_time <- pooled_status <- c()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y = raw_ext$Y, gene_names = raw_ext$gene_names, top_n = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]], cohort_name = ext_cohort,
    rank_transform = FALSE, per_platform_standardize = TRUE
  )
  common    <- intersect(train_genes, pre_ext$gene_names)
  Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
  train_idx <- match(common, train_genes)

  EF5_sub <- fit5$EF[train_idx, , drop = FALSE]
  EF7_sub <- fit7$EF[train_idx, , drop = FALSE]
  risk5 <- predict_cox_on_yf(Y_ext, EF5_sub, fit5$EBeta, EF_norms = fit5$EF_norms)$risk_scores
  risk7 <- predict_cox_on_yf(Y_ext, EF7_sub, fit7$EBeta, EF_norms = fit7$EF_norms)$risk_scores

  ci <- bootstrap_concordance_diff_ci(risk5, risk7, raw_ext$time, raw_ext$status, B = 2000, seed = 1)
  cat(sprintf("  %s: diff (K5-K7) = %.4f [%.4f, %.4f]\n",
              ext_cohort, ci$estimate, ci$lower, ci$upper))

  rows[[length(rows) + 1]] <- data.frame(
    cohort = ext_cohort, n = length(raw_ext$time), n_events = sum(raw_ext$status),
    diff_estimate = round(ci$estimate, 4), diff_lower = round(ci$lower, 4),
    diff_upper = round(ci$upper, 4), significant = ci$significant,
    stringsAsFactors = FALSE
  )

  pooled_risk5  <- c(pooled_risk5, risk5)
  pooled_risk7  <- c(pooled_risk7, risk7)
  pooled_time   <- c(pooled_time, raw_ext$time)
  pooled_status <- c(pooled_status, raw_ext$status)
}

# --------------------------------------------------------------------------
# 3. Pooled (all 5 cohorts stacked) comparison
# --------------------------------------------------------------------------
ci_pooled <- bootstrap_concordance_diff_ci(pooled_risk5, pooled_risk7, pooled_time, pooled_status,
                                            B = 2000, seed = 1)
cat(sprintf("\nPooled (n=%d): diff (K5-K7) = %.4f [%.4f, %.4f]\n",
            length(pooled_time), ci_pooled$estimate, ci_pooled$lower, ci_pooled$upper))

rows[[length(rows) + 1]] <- data.frame(
  cohort = "POOLED", n = length(pooled_time), n_events = sum(pooled_status),
  diff_estimate = round(ci_pooled$estimate, 4), diff_lower = round(ci_pooled$lower, 4),
  diff_upper = round(ci_pooled$upper, 4),
  significant = ci_pooled$significant,
  stringsAsFactors = FALSE
)

results <- do.call(rbind, rows)
out_csv <- "results/benchmark_sim/outputs/k_init_sweep/k5_vs_k7_bootstrap_ci.csv"
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nResults saved: %s\n", out_csv))
print(results)
