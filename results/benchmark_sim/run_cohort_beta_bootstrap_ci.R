# ============================================================
# Script:  results/benchmark_sim/run_cohort_beta_bootstrap_ci.R
# Purpose: Stage 3 (full): paired bootstrap CIs on the pairwise external-C
#          differences between the cohort-beta comparison arms
#          (run_cohort_beta_comparison.R) and joint_yfb (the current
#          recommended model), and between each cohort-aware arm and the
#          two-step EBMF->Cox baseline. No re-fitting: reuses the cached
#          fits from run_cohort_beta_comparison.R and the cached two-step
#          risk scores from run_ebmf_cox_regularized.R.
#
#          Per-cohort AND pooled (patient-aligned risk scores concatenated
#          across the 5 external cohorts) bootstrap CIs, matching the
#          established pattern in run_yfb_vs_ebmf_k7_matched_ci.R
#          (bootstrap_concordance_diff_ci(), B=2000, seed=1).
#
#   Output: results/benchmark_sim/outputs/cohort_beta_comparison/
#             cohort_beta_bootstrap_ci.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/benchmark_sim/run_cohort_beta_bootstrap_ci.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_beta_cohort.R")
source("code/update_L.R"); source("code/update_F.R"); source("code/update_tau.R")
source("code/compute_elbo.R"); source("code/update_F_cohort.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")
source("code/concordance_ci.R")

OUT_DIR <- "results/benchmark_sim/outputs/cohort_beta_comparison"
FITS_RDS <- file.path(OUT_DIR, "cohort_beta_comparison_fits.rds")
if (!file.exists(FITS_RDS)) stop("Missing ", FITS_RDS, " -- run run_cohort_beta_comparison.R first.")
fits <- readRDS(FITS_RDS)

EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
TRAIN_COHORTS    <- cfg$pdac$training_cohorts

# --------------------------------------------------------------------------
# Rebuild train_genes (needed to align external cohort genes to fit$EF rows)
# and reload external cohorts once -- same D4 procedure as the comparison
# script, no re-fitting.
# --------------------------------------------------------------------------
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
pp <- preprocess_merged_cohorts(
  cohort_raw_list = train_raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank"
)
train_genes <- pp$gene_names

ext_data <- list()
for (ec in EXTERNAL_COHORTS) {
  raw_ext <- load_pdac_raw(ec, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(raw_ext$Y, raw_ext$gene_names, top_n = NULL,
               log_transform = PLATFORM_LOG_TRANSFORM[[ec]], cohort_name = ec,
               rank_transform = FALSE, per_platform_standardize = TRUE)
  common <- intersect(train_genes, pre_ext$gene_names)
  ext_data[[ec]] <- list(
    Y_ext = pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE],
    train_idx = match(common, train_genes), time = raw_ext$time, status = raw_ext$status
  )
}

# --------------------------------------------------------------------------
# Per-cohort, per-arm risk scores -- oriented ONCE per arm per cohort (the
# same fixed-sign convention bootstrap_concordance_diff_ci() itself applies
# internally), from the fits already saved by run_cohort_beta_comparison.R.
# --------------------------------------------------------------------------
score_arm <- function(fit, use_pooled_beta) {
  out <- list()
  for (ec in names(ext_data)) {
    d <- ext_data[[ec]]
    EF_sub <- fit$EF[d$train_idx, , drop = FALSE]
    beta_to_use <- if (use_pooled_beta) fit$EBeta_pooled else fit$EBeta
    pred <- predict_cox_on_yf(d$Y_ext, EF_sub, beta_to_use, EF_norms = fit$EF_norms)
    out[[ec]] <- list(risk = pred$risk_scores, time = d$time, status = d$status)
  }
  out
}

arm_risk <- list(
  joint_yfb           = score_arm(fits$joint_yfb, use_pooled_beta = FALSE),
  joint_yfb_cohort_L  = score_arm(fits$joint_yfb_cohort_L, use_pooled_beta = FALSE),
  joint_yfb_beta_c    = score_arm(fits$joint_yfb_beta_c, use_pooled_beta = TRUE),
  joint_yfb_all_c     = score_arm(fits$joint_yfb_all_c, use_pooled_beta = TRUE)
)

# Two-step baseline: reuse cached risk scores (K=40, LASSO stage 2), NOT refit.
two_step_rds <- "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_riskscores.rds"
two_step <- if (file.exists(two_step_rds)) readRDS(two_step_rds)$k40 else NULL
if (!is.null(two_step)) arm_risk[["two_step_ebmf_cox"]] <- two_step

# --------------------------------------------------------------------------
# Pairwise bootstrap CIs: every arm vs. joint_yfb (the baseline this
# comparison is about), per cohort and pooled.
# --------------------------------------------------------------------------
BASELINE <- "joint_yfb"
comparison_arms <- setdiff(names(arm_risk), BASELINE)

rows <- list()
for (arm in comparison_arms) {
  pooled_a <- pooled_b <- pooled_time <- pooled_status <- c()
  for (ec in EXTERNAL_COHORTS) {
    a <- arm_risk[[BASELINE]][[ec]]; b <- arm_risk[[arm]][[ec]]
    if (is.null(a) || is.null(b)) next
    stopifnot(identical(a$time, b$time), identical(a$status, b$status))

    ci <- bootstrap_concordance_diff_ci(a$risk, b$risk, a$time, a$status, B = 2000, seed = 1)
    rows[[length(rows) + 1]] <- data.frame(
      comparison = sprintf("%s_minus_%s", BASELINE, arm), cohort = ec, n = length(a$time),
      diff_estimate = round(ci$estimate, 4), diff_lower = round(ci$lower, 4),
      diff_upper = round(ci$upper, 4), significant = ci$significant, stringsAsFactors = FALSE
    )
    pooled_a <- c(pooled_a, a$risk); pooled_b <- c(pooled_b, b$risk)
    pooled_time <- c(pooled_time, a$time); pooled_status <- c(pooled_status, a$status)
  }
  if (length(pooled_time) > 0) {
    ci_p <- bootstrap_concordance_diff_ci(pooled_a, pooled_b, pooled_time, pooled_status, B = 2000, seed = 1)
    rows[[length(rows) + 1]] <- data.frame(
      comparison = sprintf("%s_minus_%s", BASELINE, arm), cohort = "POOLED", n = length(pooled_time),
      diff_estimate = round(ci_p$estimate, 4), diff_lower = round(ci_p$lower, 4),
      diff_upper = round(ci_p$upper, 4), significant = ci_p$significant, stringsAsFactors = FALSE
    )
    cat(sprintf("%-35s POOLED (n=%d): diff=%.4f [%.4f, %.4f] sig=%s\n",
                sprintf("%s - %s", BASELINE, arm), length(pooled_time),
                ci_p$estimate, ci_p$lower, ci_p$upper, ci_p$significant))
  }
}

results <- do.call(rbind, rows)
out_csv <- file.path(OUT_DIR, "cohort_beta_bootstrap_ci.csv")
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nResults saved: %s\n", out_csv))
