# ============================================================
# Script:  results/benchmark_sim/run_k_init_sweep.R
# Purpose: Analysis A (docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md).
#          Fits YFB on the real TCGA+CPTAC training data at several K_init
#          values spanning below, at, and above the current CV-selected K=7,
#          and classifies factors at each fit with classify_factors(), to
#          check whether K_eff_survival stays at 2 regardless of the
#          starting K — i.e. whether ARD pruning from an over-specified K is
#          a stable alternative to CV-selecting K directly from held-out
#          C-index.
#
#          Also records each fit's full ELBO (elbo_full, alpha-weighted
#          genomics + survival + KL) and final RMSE, so K_init itself can be
#          selected by the ELBO criterion this project's own K-selection
#          policy prefers over CV (DECISIONS.md 2026-04-24: "ARD preferred
#          over ELBO grid search for K selection" — ARD determines K_eff
#          within one large-K fit, but comparing ELBO *across* K_init values
#          checks whether that single large fit was actually the best one,
#          since ELBO is not guaranteed to be monotone non-decreasing in K
#          once the prior's KL cost is counted).
#
#          Preprocessing matches the D4 configuration in
#          run_desurv_comparison.R exactly: YFB + per-platform z-std +
#          combined_rank gene selection (top-3000 per cohort, before
#          normalization) + no cohort_id.
#
#   Output: results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-19
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R --quick
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

# Navigate to project root if invoked from a subdirectory
if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")
source("code/select_K.R")

YML_PATH <- "config/globals.yml"
cfg      <- yaml::read_yaml(YML_PATH)
b        <- cfg$benchmark
p        <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
PVE_THRESH   <- cfg$k_selection$pve_threshold
TOP_N_DESURV <- p$top_n_genes_desurv  # 3000 — DeSurv-aligned, D4 config

# K_init values to sweep: below, at, and above the current CV-selected K=7,
# to test both ARD pruning stability and whether K=7 is actually ELBO-optimal
# among nearby candidates (not just among 7/10/15/20).
K_INIT_VALUES <- c(5L, 6L, 7L, 8L, 9L, 10L, 15L, 20L)

OUT_DIR <- "results/benchmark_sim/outputs/k_init_sweep"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | K_init values: %s\n",
            QUICK_MODE, MAX_ITER, paste(K_INIT_VALUES, collapse = ", ")))

# --------------------------------------------------------------------------
# 1. Load training data (TCGA_PAAD + CPTAC)
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Preprocess — D4 config: YFB + per-platform z-std + combined_rank +
#    top-3000 per-cohort selection (before normalization), no cohort_id.
# --------------------------------------------------------------------------

cat("--- Preprocessing training data (D4 config) ---\n")
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
Y_train    <- pp$Y
train_genes <- pp$gene_names
cat(sprintf("  n=%d, p=%d\n\n", nrow(Y_train), ncol(Y_train)))

# --------------------------------------------------------------------------
# 3. Fit YFB at each K_init, classify factors
# --------------------------------------------------------------------------

cat("=== Fitting YFB at each K_init ===\n\n")
fits <- list()

for (K_init in K_INIT_VALUES) {
  cat(sprintf("--- K_init=%d ---\n", K_init))
  set.seed(42L)
  fit <- suppressMessages(
    fit_cox_on_yf(
      Y_train, time_train, status_train,
      K = K_init, max_iter = MAX_ITER, alpha = ALPHA,
      prior_beta = PRIOR_BETA, verbose = TRUE
    )
  )
  fits[[as.character(K_init)]] <- fit
  n_iter <- fit$history$n_iter
  cat(sprintf("  |beta|: [%s]\n", paste(sprintf("%.4f", abs(fit$EBeta)), collapse = ", ")))
  cat(sprintf("  final elbo_full=%.4f | final rmse=%.6f\n\n",
              fit$history$elbo_full[n_iter], fit$history$rmse[n_iter]))
}

# --------------------------------------------------------------------------
# 4. External validation on 5 held-out cohorts, per K_init fit
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

ext_cindex <- vector("list", length(K_INIT_VALUES))
names(ext_cindex) <- as.character(K_INIT_VALUES)

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  # top_n=NULL: keep all external genes; intersection with train_genes controls
  # the final gene set — matches run_desurv_comparison.R's D4 external eval.
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE
  )

  common    <- intersect(train_genes, pre_ext$gene_names)
  if (length(common) < 100) {
    cat(sprintf("    Skipping %s: only %d common genes\n", ext_cohort, length(common)))
    next
  }
  Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
  train_idx <- match(common, train_genes)

  for (K_init in K_INIT_VALUES) {
    fit    <- fits[[as.character(K_init)]]
    EF_sub <- fit$EF[train_idx, , drop = FALSE]
    pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
    c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
    ext_cindex[[as.character(K_init)]][[ext_cohort]] <- c_val
  }
}

# --------------------------------------------------------------------------
# 5. Classify factors at each K_init and assemble results table
# --------------------------------------------------------------------------

cat("\n=== Factor classification per K_init ===\n\n")
results_rows <- list()

for (K_init in K_INIT_VALUES) {
  fit    <- fits[[as.character(K_init)]]
  n_iter <- fit$history$n_iter
  cls    <- classify_factors(fit, Y_train, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)

  K_survival_active <- sum(cls$category == "survival_active")
  K_genomics_only   <- sum(cls$category == "genomics_only")
  K_dead            <- sum(cls$category == "dead")
  K_eff_total       <- K_survival_active + K_genomics_only

  cohort_c <- ext_cindex[[as.character(K_init)]]
  mean_c   <- if (length(cohort_c) > 0) mean(unlist(cohort_c)) else NA_real_

  cat(sprintf("K_init=%2d: K_survival_active=%d, K_genomics_only=%d, K_dead=%d, K_eff_total=%d | mean external C=%.4f | elbo_full=%.4f | rmse=%.6f\n",
              K_init, K_survival_active, K_genomics_only, K_dead, K_eff_total, mean_c,
              fit$history$elbo_full[n_iter], fit$history$rmse[n_iter]))

  results_rows[[length(results_rows) + 1]] <- data.frame(
    K_init            = K_init,
    K_total           = K_init,
    K_survival_active = K_survival_active,
    K_genomics_only   = K_genomics_only,
    K_dead            = K_dead,
    K_eff_total       = K_eff_total,
    elbo_full         = round(fit$history$elbo_full[n_iter], 4),
    rmse              = round(fit$history$rmse[n_iter], 6),
    n_iter            = n_iter,
    mean_external_c   = round(mean_c, 4),
    stringsAsFactors  = FALSE
  )
}

results <- do.call(rbind, results_rows)

# ELBO-preferred K_init: this project's stated K-selection policy
# (DECISIONS.md 2026-04-24) prefers ELBO over CV/C-index for model comparison.
elbo_best_idx <- which.max(results$elbo_full)
cat(sprintf("\nELBO-preferred K_init = %d (elbo_full=%.4f, K_eff_total=%d)\n",
            results$K_init[elbo_best_idx], results$elbo_full[elbo_best_idx],
            results$K_eff_total[elbo_best_idx]))

# Attach per-cohort C-index as separate columns for full traceability.
for (ext_cohort in EXTERNAL_COHORTS) {
  results[[paste0("c_", ext_cohort)]] <- sapply(results$K_init, function(K_init) {
    v <- ext_cindex[[as.character(K_init)]][[ext_cohort]]
    if (is.null(v)) NA_real_ else round(v, 4)
  })
}

out_csv <- file.path(OUT_DIR, "k_init_sweep_results.csv")
write.csv(results, out_csv, row.names = FALSE)
saveRDS(fits, file.path(OUT_DIR, "k_init_sweep_fits.rds"))

cat(sprintf("\n=== Results saved: %s ===\n", out_csv))
