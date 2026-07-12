# ============================================================
# Script:  results/benchmark_sim/run_merged_benchmark.R
# Purpose: Comprehensive merged-cohort benchmark — all preprocessing options.
#          Fits 18 model configurations at CV-selected K (biological floor K>=3)
#          and evaluates external C-index on 5 held-out PDAC cohorts.
#
#          Model IDs:
#            M1–M6:   existing (joint_quantile_rank, perplatform_zstd) x LB/YFB x +-cohort
#            M7–M18:  new (joint_quantile_norank, joint_zstd, log_only) x LB/YFB x +-cohort
#
#          K values read from globals.yml. Run run_merged_kcv.R first.
#          YFB x joint QN excluded (structural beta->0, DECISIONS.md 2026-05-22).
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_benchmark.R [--quick]
#          --quick: max_iter=30, skip top-gene table (smoke test)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")
b   <- cfg$benchmark

# --------------------------------------------------------------------------
# Guard: all K values must be filled
# --------------------------------------------------------------------------

REQUIRED_KEYS <- c("k_merged_lb_joint", "k_merged_lb_perplatform",
                   "k_merged_yfb_perplatform", "k_merged_lb_joint_norank",
                   "k_merged_yfb_joint_norank", "k_merged_lb_zstd",
                   "k_merged_yfb_zstd", "k_merged_lb_logonly", "k_merged_yfb_logonly")
missing_keys <- REQUIRED_KEYS[sapply(REQUIRED_KEYS, function(k) is.null(b[[k]]))]
if (length(missing_keys) > 0) {
  stop(sprintf("K values not set: %s\nRun: Rscript results/benchmark_sim/run_merged_kcv.R",
               paste(missing_keys, collapse = ", ")))
}

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

ALPHA       <- b$alpha
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA  <- "normal"
SIGMA_COH   <- 1.0
BETA_THRESH <- cfg$k_selection$beta_threshold

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load training data
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
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
# 2. Model configuration table
#    Each row: model ID, model type, preprocessing params, cohort_id flag, K key
# --------------------------------------------------------------------------

MODEL_CONFIGS <- list(
  # --- existing (joint_quantile_rank) ---
  list(id="M1",  model="LB",  per_plat=FALSE, norm="quantile", rank=TRUE,  cohort=FALSE, k_key="k_merged_lb_joint"),
  list(id="M2",  model="LB",  per_plat=FALSE, norm="quantile", rank=TRUE,  cohort=TRUE,  k_key="k_merged_lb_joint"),
  # --- existing (perplatform_zstd); M5/M6 re-run at K=3 floor ---
  list(id="M3",  model="LB",  per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_lb_perplatform"),
  list(id="M4",  model="LB",  per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_perplatform"),
  list(id="M5",  model="YFB", per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_perplatform"),
  list(id="M6",  model="YFB", per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_perplatform"),
  # --- new (joint_quantile_norank) ---
  list(id="M7",  model="LB",  per_plat=FALSE, norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_lb_joint_norank"),
  list(id="M8",  model="LB",  per_plat=FALSE, norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_joint_norank"),
  list(id="M13", model="YFB", per_plat=FALSE, norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_joint_norank"),
  list(id="M14", model="YFB", per_plat=FALSE, norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_joint_norank"),
  # --- new (joint_zstd) ---
  list(id="M9",  model="LB",  per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=FALSE, k_key="k_merged_lb_zstd"),
  list(id="M10", model="LB",  per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_zstd"),
  list(id="M15", model="YFB", per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_zstd"),
  list(id="M16", model="YFB", per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_zstd"),
  # --- new (log_only) ---
  list(id="M11", model="LB",  per_plat=FALSE, norm="none",     rank=FALSE, cohort=FALSE, k_key="k_merged_lb_logonly"),
  list(id="M12", model="LB",  per_plat=FALSE, norm="none",     rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_logonly"),
  list(id="M17", model="YFB", per_plat=FALSE, norm="none",     rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_logonly"),
  list(id="M18", model="YFB", per_plat=FALSE, norm="none",     rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_logonly")
)

# --------------------------------------------------------------------------
# Helper: oriented C-index
# --------------------------------------------------------------------------

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 3. Preprocess training data — cache by (per_plat, norm, rank) combo
# --------------------------------------------------------------------------

cat("--- Preprocessing training data (all modes) ---\n")
preproc_cache  <- list()
gene_set_cache <- list()

for (mcfg in MODEL_CONFIGS) {
  cache_key <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
  if (!cache_key %in% names(preproc_cache)) {
    cat(sprintf("  Preprocessing: per_plat=%s, norm=%s, rank=%s ...\n",
                mcfg$per_plat, mcfg$norm, mcfg$rank))
    pp <- preprocess_merged_cohorts(
      train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
      top_n = cfg$preprocessing$top_n_genes,
      rank_transform = mcfg$rank, per_platform_standardize = mcfg$per_plat,
      normalize_method = mcfg$norm
    )
    preproc_cache[[cache_key]]  <- pp$Y
    gene_set_cache[[cache_key]] <- pp$gene_names
    cat(sprintf("    n=%d, p=%d\n", nrow(pp$Y), ncol(pp$Y)))
  }
}

# --------------------------------------------------------------------------
# 4. Fit all configurations
# --------------------------------------------------------------------------

cat("\n=== Fitting all 18 configurations ===\n\n")
fits <- list()

for (mcfg in MODEL_CONFIGS) {
  cache_key <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
  Y_train   <- preproc_cache[[cache_key]]
  K         <- b[[mcfg$k_key]]
  cohort_id <- if (mcfg$cohort) cohort_labels else NULL

  cat(sprintf("--- Fitting %s [%s, per_plat=%s, norm=%s, rank=%s, cohort=%s, K=%d] ---\n",
              mcfg$id, mcfg$model, mcfg$per_plat, mcfg$norm, mcfg$rank, mcfg$cohort, K))
  set.seed(42L)

  fit <- suppressMessages(
    if (mcfg$model == "LB")
      fit_supervised_mf_modular(Y_train, time_train, status_train,
                                K = K, max_iter = MAX_ITER, alpha = ALPHA,
                                prior_beta = PRIOR_BETA,
                                verbose = TRUE, cohort_id = cohort_id,
                                sigma_F_cohort = SIGMA_COH)
    else
      fit_cox_on_yf(Y_train, time_train, status_train,
                   K = K, max_iter = MAX_ITER, alpha = ALPHA,
                   prior_beta = PRIOR_BETA,
                   verbose = TRUE, cohort_id = cohort_id,
                   sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fits[[mcfg$id]] <- fit
}

# --------------------------------------------------------------------------
# 5. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y = raw_ext$Y, gene_names = raw_ext$gene_names,
    top_n = cfg$preprocessing$top_n_genes,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name = ext_cohort
  )

  for (mcfg in MODEL_CONFIGS) {
    cache_key   <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
    train_genes <- gene_set_cache[[cache_key]]
    common      <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) {
      cat(sprintf("    Skipping %s x %s: only %d common genes\n",
                  mcfg$id, ext_cohort, length(common)))
      next
    }

    Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    train_idx <- match(common, train_genes)
    fit       <- fits[[mcfg$id]]
    EF_sub    <- fit$EF[train_idx, , drop = FALSE]

    pred  <- if (mcfg$model == "LB")
      predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)
    else
      predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)

    c_val <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      model      = mcfg$id,
      model_type = mcfg$model,
      preprocess = paste(mcfg$norm, if(mcfg$per_plat) "perplat" else "joint",
                         if(mcfg$rank) "rank" else "norank", sep = "_"),
      has_cohort = mcfg$cohort,
      K          = b[[mcfg$k_key]],
      cohort     = ext_cohort,
      c_index    = round(c_val, 4),
      k_eff      = sum(abs(fit$EBeta) > BETA_THRESH),
      beta_max   = round(max(abs(fit$EBeta)), 4),
      n_iters    = fit$history$n_iter,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, results_rows)

# --------------------------------------------------------------------------
# 6. Print summary
# --------------------------------------------------------------------------

cat("\n============================================================\n")
cat(" Merged Benchmark Results - Mean External C-index by Model\n")
cat("============================================================\n")
model_order <- c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10",
                 "M11","M12","M13","M14","M15","M16","M17","M18")
for (mid in model_order) {
  sub <- results_df[results_df$model == mid, ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %s [%s, %s, cohort=%s, K=%d]: mean C=%.3f | K_eff=%d | beta_max=%.4f\n",
              mid, sub$model_type[1], sub$preprocess[1], sub$has_cohort[1],
              sub$K[1], mean(sub$c_index), sub$k_eff[1], sub$beta_max[1]))
}

# --------------------------------------------------------------------------
# 7. Top-20 gene table (skipped in quick mode)
# --------------------------------------------------------------------------

if (!QUICK_MODE) {
  cat("\n--- Factor top-20 genes ---\n")
  top_genes_rows <- list()
  for (mcfg in MODEL_CONFIGS) {
    cache_key <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
    genes     <- gene_set_cache[[cache_key]]
    fit       <- fits[[mcfg$id]]
    for (k in seq_len(ncol(fit$EF))) {
      idx <- order(abs(fit$EF[, k]), decreasing = TRUE)[1:min(20, nrow(fit$EF))]
      top_genes_rows[[length(top_genes_rows) + 1]] <- data.frame(
        model = mcfg$id, factor = k, gene = genes[idx],
        loading = round(fit$EF[idx, k], 4), stringsAsFactors = FALSE
      )
    }
  }
  top_genes_df <- do.call(rbind, top_genes_rows)
  write.csv(top_genes_df,
            file.path(OUT_DIR, "merged_benchmark_top_genes_extended.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved: %s\n",
              file.path(OUT_DIR, "merged_benchmark_top_genes_extended.csv")))
}

# --------------------------------------------------------------------------
# 8. Save results
# --------------------------------------------------------------------------

write.csv(results_df,
          file.path(OUT_DIR, "merged_benchmark_results_extended.csv"),
          row.names = FALSE)

compact_fits <- lapply(fits, function(f) list(
  EBeta      = f$EBeta,
  EBeta2     = f$EBeta2,
  EF         = f$EF,
  EF_cohort  = f$EF_cohort,
  EF2_cohort = f$EF2_cohort,
  EF_norms   = if (!is.null(f$EF_norms)) f$EF_norms else NULL,
  history    = f$history
))
saveRDS(list(fits      = compact_fits,
             results   = results_df,
             gene_sets = gene_set_cache,
             params    = list(ALPHA = ALPHA, PRIOR_BETA = PRIOR_BETA, K_floor = 3L),
             date      = Sys.time()),
        file.path(OUT_DIR, "merged_benchmark_fits_extended.rds"))

cat(sprintf("\nResults saved to: %s\n", OUT_DIR))
cat("============================================================\n")
