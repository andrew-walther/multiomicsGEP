# ============================================================
# Script:  results/benchmark_sim/run_ebmf_cox_regularized.R
# Purpose: A fairer two-step (EBMF -> Cox) baseline: a large K NOT informed
#          by YFB's answer (reusing the already-fitted K=20 and K=40 EBMF
#          factors from run_ebmf_cox_external.R --k 20 / --k 40, no
#          re-fitting of the unsupervised step), but with stage 2 replaced
#          by a REGULARIZED (LASSO) Cox model (glmnet, family="cox",
#          lambda chosen by cross-validation) instead of a plain
#          unregularized coxph() on all K factor scores.
#
#          Motivation (DECISIONS.md 2026-08-20 discussion): giving EBMF the
#          same K=7 as YFB risks "leaking" the joint model's answer about
#          model complexity. But giving EBMF a large, uninformed K (20 or
#          40) with a plain unregularized coxph() at stage 2 confounds two
#          separate questions -- "does the unsupervised step need a large K"
#          and "does an unregularized multi-covariate Cox model overfit
#          with that many covariates relative to events" -- and the K=40
#          run's training-C-up/external-C-down pattern (0.731 / 0.578) is
#          a textbook overfitting signature. A LASSO Cox is also a more
#          faithful implementation of the two-step method's actual intent
#          (identify WHICH factors are prognostic), which plain coxph
#          never does (every factor gets a nonzero coefficient regardless).
#
#   Output: results/benchmark_sim/outputs/ebmf_cox_external/
#             ebmf_cox_regularized_results.csv (one row per K x cohort)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   PDAC_DATA_ROOT=<path> Rscript results/benchmark_sim/run_ebmf_cox_regularized.R
# Requires: ebmf_cox_external_fit_k20.rds / _k40.rds already exist (run
#           run_ebmf_cox_external.R --k 20 and --k 40 first if missing;
#           NOTE the K=20 run's default filenames are unsuffixed -- see
#           that script's k_suffix logic).
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival); library(glmnet) })
cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")

TOP_N_DESURV     <- cfg$preprocessing$top_n_genes_desurv
TRAIN_COHORTS    <- cfg$pdac$training_cohorts
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

OUT_DIR <- "results/benchmark_sim/outputs/ebmf_cox_external"

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}
project_scores <- function(Ymat, Fmat, F_norms) Ymat %*% sweep(Fmat, 2, F_norms, "/")

# --------------------------------------------------------------------------
# 1. Load + preprocess training data (identical to run_ebmf_cox_external.R)
# --------------------------------------------------------------------------
cat("--- Loading TCGA_PAAD + CPTAC ---\n")
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))

pp <- preprocess_merged_cohorts(
  cohort_raw_list = train_raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = TOP_N_DESURV, rank_transform = FALSE, per_platform_standardize = TRUE,
  normalize_method = "none", selection_per_cohort = TRUE, selection_method = "combined_rank"
)
Y_train <- pp$Y

# --------------------------------------------------------------------------
# 2. For each cached EBMF fit (K=20, K=40): refit stage 2 with LASSO Cox
# --------------------------------------------------------------------------
FITS <- list(
  k20 = "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit.rds",     # unsuffixed = default K=20
  k40 = "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit_k40.rds"
)

results_rows <- list()
risk_cache   <- list()   # [[tag]][[cohort]] = {risk, time, status} for post-hoc bootstrap CIs

for (tag in names(FITS)) {
  path <- FITS[[tag]]
  if (!file.exists(path)) { message(sprintf("Skipping %s: %s not found", tag, path)); next }
  f <- readRDS(path)
  K_ebmf <- f$K_ebmf
  cat(sprintf("\n=== %s (K=%d) ===\n", tag, K_ebmf))

  # f$train_genes is pp$gene_names from the ORIGINAL fitting run (same preprocessing
  # call, deterministic) -- F_ebmf's rows are already in that same gene order, so
  # S_train can be recomputed directly against the freshly-loaded Y_train (same order).
  stopifnot(identical(f$train_genes, pp$gene_names))
  S_train <- project_scores(Y_train, f$F_ebmf, f$F_norms)
  colnames(S_train) <- paste0("EBMF", seq_len(K_ebmf))

  cv_fit <- cv.glmnet(S_train, Surv(time_train, status_train), family = "cox", alpha = 1)
  beta_lasso <- as.numeric(coef(cv_fit, s = "lambda.min"))
  names(beta_lasso) <- colnames(S_train)
  k_eff <- sum(abs(beta_lasso) > 1e-8)
  cat(sprintf("  LASSO (lambda.min): %d / %d factors with nonzero coefficient\n", k_eff, K_ebmf))

  risk_train <- as.numeric(S_train %*% beta_lasso)
  c_train <- oriented_cindex(risk_train, time_train, status_train)
  cat(sprintf("  Training oriented C-index = %.4f\n", c_train))

  for (ext_cohort in EXTERNAL_COHORTS) {
    raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
    pre_ext <- preprocess_desurv_cohort(
      Y = raw_ext$Y, gene_names = raw_ext$gene_names, top_n = NULL,
      log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]], cohort_name = ext_cohort,
      rank_transform = FALSE, per_platform_standardize = TRUE
    )
    common <- intersect(f$train_genes, pre_ext$gene_names)
    if (length(common) < 100) { cat(sprintf("  Skipping %s: only %d common genes\n", ext_cohort, length(common))); next }

    Y_ext  <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    tr_idx <- match(common, f$train_genes)
    F_sub  <- f$F_ebmf[tr_idx, , drop = FALSE]
    S_ext  <- project_scores(Y_ext, F_sub, f$F_norms)
    risk   <- as.numeric(S_ext %*% beta_lasso)
    c_val  <- oriented_cindex(risk, raw_ext$time, raw_ext$status)
    cat(sprintf("  %-20s C=%.4f (%d common genes)\n", ext_cohort, c_val, length(common)))

    if (is.null(risk_cache[[tag]])) risk_cache[[tag]] <- list()
    risk_cache[[tag]][[ext_cohort]] <- list(risk = risk, time = raw_ext$time, status = raw_ext$status)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      k_tag = tag, K = K_ebmf, cohort = ext_cohort, c_index = round(c_val, 4),
      k_eff_lasso = k_eff, training_c = round(c_train, 4), n_common_genes = length(common),
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, results_rows)
out_csv <- file.path(OUT_DIR, "ebmf_cox_regularized_results.csv")
write.csv(results, out_csv, row.names = FALSE)
saveRDS(risk_cache, file.path(OUT_DIR, "ebmf_cox_regularized_riskscores.rds"))

cat("\n============================================================\n")
cat(" Mean external C by K (LASSO stage 2)\n")
cat("============================================================\n")
agg <- aggregate(c_index ~ k_tag + K, data = results, FUN = mean)
for (i in seq_len(nrow(agg)))
  cat(sprintf("  %s (K=%d): mean external C=%.4f\n", agg$k_tag[i], agg$K[i], agg$c_index[i]))
cat(sprintf("\nResults saved: %s\n", out_csv))
