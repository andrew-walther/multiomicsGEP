# ============================================================
# Script:  results/benchmark_sim/run_cohort_lmm_benchmark.R
# Purpose: Stage 2 real-data benchmark for the cohort indicator extension.
#          Fits four models on merged TCGA_PAAD + CPTAC training data and
#          evaluates external C-index on 5 held-out PDAC cohorts.
#
#          Models compared (all use prior_beta="normal"):
#            LB_base    — LB (eta = L·beta), no cohort adjustment
#            LB_cohort  — LB + cohort_id (TCGA vs CPTAC corner-point column)
#            YFB_base   — YFB (eta = (Y·F)·beta), no cohort adjustment
#            YFB_cohort — YFB + cohort_id
#
#          Baseline reference (from main-branch benchmark runs):
#            LB merged:  K_eff=3, beta_max=0.021, C-ext: 0.51–0.67
#            YFB merged: K_eff=0 (beta->0), C-ext: 0.50 (point_normal)
#                        / 0.54–0.64 (normal)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-22
# Usage:   Rscript results/benchmark_sim/run_cohort_lmm_benchmark.R [--quick]
#          --quick: K_LB=5, K_YFB=2, max_iter=30 (smoke test)
# ============================================================

# --------------------------------------------------------------------------
# 0. Setup
# --------------------------------------------------------------------------

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args
LOW_K_MODE <- "--low-k" %in% args

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R")
source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

K_LB       <- if (QUICK_MODE) 5L  else if (LOW_K_MODE) 5L  else cfg$benchmark$k_pdac       # 20
K_YFB      <- if (QUICK_MODE) 2L  else if (LOW_K_MODE) 3L  else cfg$benchmark$k_pdac_yfb_merged  # 3
ALPHA      <- cfg$benchmark$alpha        # 0.5
LAMBDA     <- cfg$benchmark$lambda
MAX_ITER   <- if (QUICK_MODE) 30L else cfg$cavi$max_iter           # 300
PRIOR_BETA <- "normal"   # avoids beta->0 collapse; consistent with best main-branch result
SIGMA_COH  <- 1.0        # sigma_F_cohort: prior SD for cohort F rows
BETA_THRESH <- cfg$k_selection$beta_threshold

cat("=== Cohort LMM Benchmark — 4-way comparison ===\n")
cat(sprintf("    K_LB=%d | K_YFB=%d | alpha=%.2f | prior_beta=%s\n",
            K_LB, K_YFB, ALPHA, PRIOR_BETA))
cat(sprintf("    max_iter=%d | sigma_F_cohort=%.1f | QUICK=%s | LOW_K=%s\n\n",
            MAX_ITER, SIGMA_COH, QUICK_MODE, LOW_K_MODE))

OUT_DIR <- if (LOW_K_MODE) "results/benchmark_sim/outputs/cohort_lmm_benchmark_low_k" else
                            "results/benchmark_sim/outputs/cohort_lmm_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load and preprocess merged TCGA_PAAD + CPTAC training data
# --------------------------------------------------------------------------

if (!dir.exists(PDAC_DATA_ROOT)) {
  stop(sprintf(
    "PDAC data not found at %s.\nSet PDAC_DATA_ROOT env var or check OneDrive sync.",
    PDAC_DATA_ROOT
  ))
}

cat("--- Loading training cohorts (TCGA_PAAD + CPTAC) ---\n")
TRAIN_COHORTS <- c("TCGA_PAAD", "CPTAC")
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  Loading %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})

cat("  Preprocessing (intersect-first + rank-transform) ...\n")
merged <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n               = cfg$preprocessing$top_n_genes,
  rank_transform      = TRUE
)
Y_train      <- merged$Y
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
train_genes  <- merged$gene_names

# Build cohort indicator: one label per patient, matching row order in Y_train
n_tcga <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))

cat(sprintf("  Training matrix: n=%d (TCGA=%d, CPTAC=%d), p=%d\n",
            nrow(Y_train), n_tcga, n_cptac, ncol(Y_train)))
cat(sprintf("  Events: %d / %d (%.1f%% censored)\n\n",
            sum(status_train), nrow(Y_train),
            100 * mean(status_train == 0)))

# --------------------------------------------------------------------------
# Helper: compute oriented C-index (max of c, 1-c)
# --------------------------------------------------------------------------
oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. Fit all four models
# --------------------------------------------------------------------------

fits <- list()

# ------- 2a. LB_base -------
cat("--- Fitting LB_base ---\n")
set.seed(42L)
fits$LB_base <- suppressMessages(
  fit_supervised_mf_modular(Y_train, time_train, status_train,
                            K        = K_LB,
                            max_iter = MAX_ITER,
                            alpha    = ALPHA,
                            lambda   = LAMBDA,
                            prior_beta = PRIOR_BETA,
                            verbose  = TRUE,
                            cohort_id = NULL)
)
cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
            sum(abs(fits$LB_base$EBeta) > BETA_THRESH),
            max(abs(fits$LB_base$EBeta)),
            fits$LB_base$history$n_iter))

# ------- 2b. LB_cohort -------
cat("--- Fitting LB_cohort ---\n")
set.seed(42L)
fits$LB_cohort <- suppressMessages(
  fit_supervised_mf_modular(Y_train, time_train, status_train,
                            K        = K_LB,
                            max_iter = MAX_ITER,
                            alpha    = ALPHA,
                            lambda   = LAMBDA,
                            prior_beta     = PRIOR_BETA,
                            verbose        = TRUE,
                            cohort_id      = cohort_labels,
                            sigma_F_cohort = SIGMA_COH)
)
cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n",
            sum(abs(fits$LB_cohort$EBeta) > BETA_THRESH),
            max(abs(fits$LB_cohort$EBeta)),
            fits$LB_cohort$history$n_iter))
if (!is.null(fits$LB_cohort$EF_cohort)) {
  cohort_norm <- sqrt(sum(fits$LB_cohort$EF_cohort^2))
  cat(sprintf("  ||EF_cohort||_F = %.4f  (L2 norm of platform offset row)\n\n",
              cohort_norm))
}

# ------- 2c. YFB_base -------
cat("--- Fitting YFB_base ---\n")
set.seed(42L)
fits$YFB_base <- suppressMessages(
  fit_cox_on_yf(Y_train, time_train, status_train,
                K          = K_YFB,
                max_iter   = MAX_ITER,
                alpha      = ALPHA,
                lambda     = LAMBDA,
                prior_beta = PRIOR_BETA,
                verbose    = TRUE,
                cohort_id  = NULL)
)
cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
            sum(abs(fits$YFB_base$EBeta) > BETA_THRESH),
            max(abs(fits$YFB_base$EBeta)),
            fits$YFB_base$history$n_iter))

# ------- 2d. YFB_cohort -------
cat("--- Fitting YFB_cohort ---\n")
set.seed(42L)
fits$YFB_cohort <- suppressMessages(
  fit_cox_on_yf(Y_train, time_train, status_train,
                K              = K_YFB,
                max_iter       = MAX_ITER,
                alpha          = ALPHA,
                lambda         = LAMBDA,
                prior_beta     = PRIOR_BETA,
                verbose        = TRUE,
                cohort_id      = cohort_labels,
                sigma_F_cohort = SIGMA_COH)
)
cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n",
            sum(abs(fits$YFB_cohort$EBeta) > BETA_THRESH),
            max(abs(fits$YFB_cohort$EBeta)),
            fits$YFB_cohort$history$n_iter))
if (!is.null(fits$YFB_cohort$EF_cohort)) {
  cohort_norm_yfb <- sqrt(sum(fits$YFB_cohort$EF_cohort^2))
  cat(sprintf("  ||EF_cohort||_F = %.4f\n\n", cohort_norm_yfb))
}

# --------------------------------------------------------------------------
# 3. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation ---\n")

results_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = cfg$preprocessing$top_n_genes,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort
  )

  # Intersect to training gene set
  common <- intersect(train_genes, pre_ext$gene_names)
  if (length(common) < 100) {
    cat(sprintf("    [SKIP] Only %d common genes — too few for reliable prediction.\n",
                length(common)))
    next
  }
  Y_ext <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]

  # Subset EF matrices to common genes
  train_idx_genes <- match(common, train_genes)

  # --- LB_base ---
  EF_lb_sub  <- fits$LB_base$EF[train_idx_genes, , drop = FALSE]
  pred_lb    <- predict_supervised_mf(Y_ext, EF_lb_sub, fits$LB_base$EBeta)
  c_lb       <- oriented_cindex(pred_lb$risk_scores, raw_ext$time, raw_ext$status)

  # --- LB_cohort ---
  EF_lbc_sub <- fits$LB_cohort$EF[train_idx_genes, , drop = FALSE]
  pred_lbc   <- predict_supervised_mf(Y_ext, EF_lbc_sub, fits$LB_cohort$EBeta)
  c_lbc      <- oriented_cindex(pred_lbc$risk_scores, raw_ext$time, raw_ext$status)

  # --- YFB_base ---
  EF_yfb_sub  <- fits$YFB_base$EF[train_idx_genes, , drop = FALSE]
  pred_yfb    <- predict_cox_on_yf(Y_ext, EF_yfb_sub, fits$YFB_base$EBeta,
                                    EF_norms = fits$YFB_base$EF_norms)
  c_yfb       <- oriented_cindex(pred_yfb$risk_scores, raw_ext$time, raw_ext$status)

  # --- YFB_cohort ---
  EF_yfbc_sub <- fits$YFB_cohort$EF[train_idx_genes, , drop = FALSE]
  pred_yfbc   <- predict_cox_on_yf(Y_ext, EF_yfbc_sub, fits$YFB_cohort$EBeta,
                                    EF_norms = fits$YFB_cohort$EF_norms)
  c_yfbc      <- oriented_cindex(pred_yfbc$risk_scores, raw_ext$time, raw_ext$status)

  cat(sprintf(
    "    %-22s | LB_base=%.3f  LB_cohort=%.3f  YFB_base=%.3f  YFB_cohort=%.3f\n",
    ext_cohort, c_lb, c_lbc, c_yfb, c_yfbc
  ))

  for (nm in c("LB_base","LB_cohort","YFB_base","YFB_cohort")) {
    c_val <- c(LB_base=c_lb, LB_cohort=c_lbc, YFB_base=c_yfb, YFB_cohort=c_yfbc)[nm]
    model_type <- sub("_.*", "", nm)
    has_cohort <- grepl("cohort", nm)
    fit_obj    <- fits[[nm]]
    results_rows[[length(results_rows) + 1]] <- data.frame(
      model      = nm,
      model_type = model_type,
      has_cohort = has_cohort,
      cohort     = ext_cohort,
      c_index    = round(c_val, 4),
      k_eff      = sum(abs(fit_obj$EBeta) > BETA_THRESH),
      beta_max   = round(max(abs(fit_obj$EBeta)), 4),
      n_iters    = fit_obj$history$n_iter,
      stringsAsFactors = FALSE
    )
  }
}

# --------------------------------------------------------------------------
# 4. Summary table
# --------------------------------------------------------------------------

results_df <- do.call(rbind, results_rows)

cat("\n============================================================\n")
cat(" Stage 2 Results — External C-index (oriented)\n")
cat("============================================================\n")

# Wide format: one row per cohort, 4 model columns
cohorts_eval <- unique(results_df$cohort)
wide_rows <- lapply(cohorts_eval, function(coh) {
  sub_df <- results_df[results_df$cohort == coh, ]
  row_vals <- setNames(sub_df$c_index, sub_df$model)
  c(cohort = coh, as.list(row_vals[c("LB_base","LB_cohort","YFB_base","YFB_cohort")]))
})
wide_df <- do.call(rbind, lapply(wide_rows, as.data.frame, stringsAsFactors=FALSE))
for (col in c("LB_base","LB_cohort","YFB_base","YFB_cohort")) {
  wide_df[[col]] <- as.numeric(wide_df[[col]])
}

cat(sprintf("%-22s  %8s  %10s  %8s  %11s\n",
            "Cohort", "LB_base", "LB_cohort", "YFB_base", "YFB_cohort"))
cat(sprintf("%-22s  %8s  %10s  %8s  %11s\n",
            "------", "-------", "---------", "--------", "----------"))
for (i in seq_len(nrow(wide_df))) {
  cat(sprintf("%-22s  %8.3f  %10.3f  %8.3f  %11.3f\n",
              wide_df$cohort[i],
              wide_df$LB_base[i], wide_df$LB_cohort[i],
              wide_df$YFB_base[i], wide_df$YFB_cohort[i]))
}

# Mean C-index across external cohorts
means <- colMeans(wide_df[, c("LB_base","LB_cohort","YFB_base","YFB_cohort")],
                  na.rm = TRUE)
cat(sprintf("%-22s  %8.3f  %10.3f  %8.3f  %11.3f\n",
            "MEAN", means["LB_base"], means["LB_cohort"],
            means["YFB_base"], means["YFB_cohort"]))

# Cohort ||EF_cohort|| norms for offset absorption
cat("\n--- Cohort offset absorption ---\n")
for (nm in c("LB_cohort","YFB_cohort")) {
  ef_c <- fits[[nm]]$EF_cohort
  if (!is.null(ef_c)) {
    cat(sprintf("  %s: ||EF_cohort||_F = %.4f | EF_cohort range [%.4f, %.4f]\n",
                nm, sqrt(sum(ef_c^2)), min(ef_c), max(ef_c)))
  }
}

cat("\n--- Model summary (K_eff / beta_max) ---\n")
for (nm in names(fits)) {
  fit_obj <- fits[[nm]]
  cat(sprintf("  %-12s: K_eff=%d | beta_max=%.4f | iters=%d\n",
              nm,
              sum(abs(fit_obj$EBeta) > BETA_THRESH),
              max(abs(fit_obj$EBeta)),
              fit_obj$history$n_iter))
}

cat("\n--- Baseline reference (main-branch merged, prior_beta=normal) ---\n")
cat("  LB_base (main):  C-ext 0.571 / 0.513 / 0.648 / 0.671 / 0.644 (mean~0.609)\n")
cat("  YFB_base (main): C-ext 0.573 / 0.537 / 0.607 / 0.637 / 0.562 (mean~0.583)\n")

# --------------------------------------------------------------------------
# 5. Save results
# --------------------------------------------------------------------------

out_csv <- file.path(OUT_DIR, "cohort_lmm_benchmark_results.csv")
write.csv(results_df, out_csv, row.names = FALSE)

out_rds <- file.path(OUT_DIR, "cohort_lmm_benchmark_fits.rds")
# Save EBeta, EF_cohort, EF norms — not full EL/EF matrices (large)
compact_fits <- lapply(fits, function(f) {
  list(EBeta    = f$EBeta,
       EBeta2   = f$EBeta2,
       EF_cohort  = f$EF_cohort,
       EF2_cohort = f$EF2_cohort,
       EF_norms   = if (!is.null(f$EF_norms)) f$EF_norms else NULL,
       history    = f$history,
       n_iter     = f$history$n_iter)
})
saveRDS(list(fits = compact_fits, results = results_df, wide = wide_df,
             params = list(K_LB=K_LB, K_YFB=K_YFB, ALPHA=ALPHA,
                           PRIOR_BETA=PRIOR_BETA, SIGMA_COH=SIGMA_COH),
             date = Sys.time()),
        out_rds)

cat(sprintf("\nResults saved:\n  %s\n  %s\n", out_csv, out_rds))
cat("============================================================\n")
