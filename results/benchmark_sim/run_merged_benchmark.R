# ============================================================
# Script:  results/benchmark_sim/run_merged_benchmark.R
# Purpose: Comprehensive merged-cohort benchmark (TCGA_PAAD + CPTAC).
#          Fits 6 model configurations at CV-selected K and evaluates
#          external C-index on 5 held-out PDAC cohorts.
#
#          Configurations (all prior_beta="normal"):
#            M1: LB  x joint quantile+rank  x no cohort_id
#            M2: LB  x joint quantile+rank  x cohort_id
#            M3: LB  x per-platform z-std   x no cohort_id
#            M4: LB  x per-platform z-std   x cohort_id
#            M5: YFB x per-platform z-std   x no cohort_id
#            M6: YFB x per-platform z-std   x cohort_id
#
#          K for each configuration is read from globals.yml
#          (benchmark.k_merged_lb_joint, k_merged_lb_perplatform,
#           k_merged_yfb_perplatform). Run run_merged_kcv.R first.
#
#          Excluded (documented beta->0 structural failure, all V0-V11 exhausted):
#            YFB x joint quantile+rank x No/Yes
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_benchmark.R [--quick]
#          --quick: max_iter=30, skips interpretability output (smoke test)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../../")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

cfg <- yaml::read_yaml("config/globals.yml")

# Guard: K values must be filled by run_merged_kcv.R first
k_lb_joint        <- cfg$benchmark$k_merged_lb_joint
k_lb_perplatform  <- cfg$benchmark$k_merged_lb_perplatform
k_yfb_perplatform <- cfg$benchmark$k_merged_yfb_perplatform

if (is.null(k_lb_joint) || is.null(k_lb_perplatform) || is.null(k_yfb_perplatform)) {
  missing <- c(
    if (is.null(k_lb_joint))        "k_merged_lb_joint",
    if (is.null(k_lb_perplatform))  "k_merged_lb_perplatform",
    if (is.null(k_yfb_perplatform)) "k_merged_yfb_perplatform"
  )
  stop(sprintf(
    "K values not yet set in globals.yml: %s\nRun: Rscript results/benchmark_sim/run_merged_kcv.R",
    paste(missing, collapse = ", ")
  ))
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

ALPHA      <- cfg$benchmark$alpha
LAMBDA     <- cfg$benchmark$lambda
MAX_ITER   <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA <- "normal"
SIGMA_COH  <- 1.0
BETA_THRESH <- cfg$k_selection$beta_threshold

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Merged Cohort Benchmark - 6-way comparison ===\n")
cat(sprintf("    alpha=%.2f | prior_beta=%s | max_iter=%d | QUICK=%s\n",
            ALPHA, PRIOR_BETA, MAX_ITER, QUICK_MODE))
cat(sprintf("    K: LB_joint=%d | LB_perplatform=%d | YFB_perplatform=%d\n\n",
            k_lb_joint, k_lb_perplatform, k_yfb_perplatform))

# --------------------------------------------------------------------------
# 1. Load and preprocess training data -- both versions
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

merged_joint <- preprocess_merged_cohorts(
  train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes,
  rank_transform = TRUE, per_platform_standardize = FALSE
)

merged_perplatform <- preprocess_merged_cohorts(
  train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes,
  rank_transform = FALSE, per_platform_standardize = TRUE
)

cat(sprintf("  joint: n=%d, p=%d\n", nrow(merged_joint$Y), ncol(merged_joint$Y)))
cat(sprintf("  perplatform: n=%d, p=%d\n\n", nrow(merged_perplatform$Y),
            ncol(merged_perplatform$Y)))

# --------------------------------------------------------------------------
# Helper: oriented C-index (always >= 0.5 by flipping sign if needed)
# --------------------------------------------------------------------------
oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. Fit all 6 configurations
# --------------------------------------------------------------------------

fits      <- list()
gene_sets <- list()   # store gene_names per preprocessing type for external validation

fit_lb <- function(Y, time, status, K, cohort_id = NULL, label = "") {
  cat(sprintf("--- Fitting %s (K=%d) ---\n", label, K))
  set.seed(42L)
  fit <- suppressMessages(
    fit_supervised_mf_modular(Y, time, status,
                              K          = K,
                              max_iter   = MAX_ITER,
                              alpha      = ALPHA,
                              lambda     = LAMBDA,
                              prior_beta = PRIOR_BETA,
                              verbose    = TRUE,
                              cohort_id  = cohort_id,
                              sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fit
}

fit_yfb <- function(Y, time, status, K, cohort_id = NULL, label = "") {
  cat(sprintf("--- Fitting %s (K=%d) ---\n", label, K))
  set.seed(42L)
  fit <- suppressMessages(
    fit_cox_on_yf(Y, time, status,
                  K          = K,
                  max_iter   = MAX_ITER,
                  alpha      = ALPHA,
                  lambda     = LAMBDA,
                  prior_beta = PRIOR_BETA,
                  verbose    = TRUE,
                  cohort_id  = cohort_id,
                  sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fit
}

fits$M1 <- fit_lb(merged_joint$Y,       time_train, status_train, k_lb_joint,
                   cohort_id = NULL,         label = "M1 LB joint no-cohort")
fits$M2 <- fit_lb(merged_joint$Y,       time_train, status_train, k_lb_joint,
                   cohort_id = cohort_labels, label = "M2 LB joint cohort_id")
fits$M3 <- fit_lb(merged_perplatform$Y, time_train, status_train, k_lb_perplatform,
                   cohort_id = NULL,         label = "M3 LB perplatform no-cohort")
fits$M4 <- fit_lb(merged_perplatform$Y, time_train, status_train, k_lb_perplatform,
                   cohort_id = cohort_labels, label = "M4 LB perplatform cohort_id")
fits$M5 <- fit_yfb(merged_perplatform$Y, time_train, status_train, k_yfb_perplatform,
                    cohort_id = NULL,         label = "M5 YFB perplatform no-cohort")
fits$M6 <- fit_yfb(merged_perplatform$Y, time_train, status_train, k_yfb_perplatform,
                    cohort_id = cohort_labels, label = "M6 YFB perplatform cohort_id")

gene_sets$joint       <- merged_joint$gene_names
gene_sets$perplatform <- merged_perplatform$gene_names

# --------------------------------------------------------------------------
# 3. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")

EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = cfg$preprocessing$top_n_genes,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort
  )

  # For each training gene set, intersect and evaluate all models using that set
  for (prep_type in c("joint", "perplatform")) {
    train_genes <- gene_sets[[prep_type]]
    common      <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) next
    Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    train_idx <- match(common, train_genes)

    models_for_prep <- if (prep_type == "joint") c("M1", "M2") else c("M3", "M4", "M5", "M6")
    for (mid in models_for_prep) {
      fit    <- fits[[mid]]
      is_yfb <- mid %in% c("M5", "M6")
      if (!is_yfb) {
        EF_sub <- fit$EF[train_idx, , drop = FALSE]
        pred   <- predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)
        c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
      } else {
        EF_sub <- fit$EF[train_idx, , drop = FALSE]
        pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
        c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
      }
      results_rows[[length(results_rows) + 1]] <- data.frame(
        model      = mid,
        cohort     = ext_cohort,
        c_index    = round(c_val, 4),
        k_eff      = sum(abs(fit$EBeta) > BETA_THRESH),
        beta_max   = round(max(abs(fit$EBeta)), 4),
        n_iters    = fit$history$n_iter,
        preprocess = prep_type,
        has_cohort = mid %in% c("M2", "M4", "M6"),
        stringsAsFactors = FALSE
      )
    }
  }
}

results_df <- do.call(rbind, results_rows)

# --------------------------------------------------------------------------
# 4. Summary table
# --------------------------------------------------------------------------

cat("\n============================================================\n")
cat(" Merged Benchmark Results - External C-index\n")
cat("============================================================\n")

model_ids    <- c("M1","M2","M3","M4","M5","M6")
model_labels <- c("LB_joint","LB_joint_coh","LB_perplat","LB_perplat_coh",
                  "YFB_perplat","YFB_perplat_coh")

for (i in seq_along(model_ids)) {
  mid <- model_ids[i]
  sub <- results_df[results_df$model == mid, ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %s (%s):\n", mid, model_labels[i]))
  for (j in seq_len(nrow(sub))) {
    cat(sprintf("    %-22s C=%.3f\n", sub$cohort[j], sub$c_index[j]))
  }
  cat(sprintf("    %-22s C=%.3f  (K_eff=%d, beta_max=%.4f)\n\n",
              "MEAN", mean(sub$c_index), sub$k_eff[1], sub$beta_max[1]))
}

# --------------------------------------------------------------------------
# 5. Factor top-20 gene table (interpretability)
# --------------------------------------------------------------------------

if (!QUICK_MODE) {
  cat("--- Factor top-20 genes per configuration ---\n")
  top_genes_list <- list()
  for (mid in model_ids) {
    fit <- fits[[mid]]
    K_fit <- ncol(fit$EF)
    prep_type <- if (mid %in% c("M1","M2")) "joint" else "perplatform"
    genes <- gene_sets[[prep_type]]
    top_genes_list[[mid]] <- lapply(seq_len(K_fit), function(k) {
      idx <- order(abs(fit$EF[, k]), decreasing = TRUE)[1:min(20, nrow(fit$EF))]
      data.frame(model = mid, factor = k, gene = genes[idx],
                 loading = round(fit$EF[idx, k], 4))
    })
  }
  top_genes_df <- do.call(rbind, lapply(top_genes_list, function(x) do.call(rbind, x)))
  write.csv(top_genes_df,
            file.path(OUT_DIR, "merged_benchmark_top_genes.csv"),
            row.names = FALSE)
  cat(sprintf("  Top genes saved: %s\n\n",
              file.path(OUT_DIR, "merged_benchmark_top_genes.csv")))
}

# --------------------------------------------------------------------------
# 6. Save results
# --------------------------------------------------------------------------

write.csv(results_df, file.path(OUT_DIR, "merged_benchmark_results.csv"), row.names = FALSE)

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
             gene_sets = gene_sets,
             params    = list(k_lb_joint        = k_lb_joint,
                              k_lb_perplatform  = k_lb_perplatform,
                              k_yfb_perplatform = k_yfb_perplatform,
                              ALPHA      = ALPHA,
                              PRIOR_BETA = PRIOR_BETA),
             date = Sys.time()),
        file.path(OUT_DIR, "merged_benchmark_fits.rds"))

cat(sprintf("Results saved to: %s\n", OUT_DIR))
cat("============================================================\n")
