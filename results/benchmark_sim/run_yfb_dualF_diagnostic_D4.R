# ============================================================
# Script:  results/benchmark_sim/run_yfb_dualF_diagnostic_D4.R
# Purpose: Step 0 of docs/plans/yfb_dual_source_F_experiments_08_20_2026.md.
#          Reconfirm the alpha_F>0 dual-source-F failure (DECISIONS.md
#          2026-05-22) under the CURRENT recommended D4 preprocessing
#          (per-platform z-std + DeSurv combined_rank top-3000 per-cohort
#          gene selection, K=7 fixed) instead of the original diagnostic's
#          merged-rank-transform preprocessing, which CLAUDE.md documents as
#          collapsing beta->0 regardless of model. That confound made the
#          original 2026-05-22 alpha_F sweep (V6-V8 in
#          run_yfb_beta_fix_diagnostic.R) uninformative about whether
#          alpha_F specifically is the problem under the actually-recommended
#          config.
#
#          Sweeps alpha_F in {0 (baseline), 0.1, 0.3, 0.5} with everything
#          else fixed at the D4 recommended config (see
#          run_desurv_comparison.R's D4 entry). Full 5-cohort external
#          validation, matching the exact D4 protocol, for a like-for-like
#          comparison against the recommended config's 0.627 mean C-index.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   Rscript results/benchmark_sim/run_yfb_dualF_diagnostic_D4.R [--quick]
#          --quick: max_iter=30 (fast sanity check)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/predict_cox_on_yf.R")
source("code/preprocess_desurv.R")

b <- cfg$benchmark
p <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
TOP_N_DESURV <- p$top_n_genes_desurv
K_D4         <- b$k_merged_yfb_desurv
ALPHA_F_GRID <- c(0, 0.1, 0.3, 0.5)

if (is.null(K_D4))
  stop("config/globals.yml benchmark$k_merged_yfb_desurv is NULL -- D4's K has not been set.")

OUT_DIR <- "results/benchmark_sim/outputs/yfb_dualF_diagnostic_D4"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== YFB Dual-Source F Diagnostic (D4 preprocessing, alpha_F sweep) ===\n")
cat(sprintf("    K=%d (fixed, from globals.yml k_merged_yfb_desurv) | alpha=%.2f | max_iter=%d\n",
            K_D4, ALPHA, MAX_ITER))
cat(sprintf("    alpha_F grid: %s\n\n", paste(ALPHA_F_GRID, collapse = ", ")))

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 1. Load + preprocess training data (D4: per-platform z-std, combined_rank,
#    top-3000 per-cohort, no rank transform, no quantile normalization)
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d, events=%d\n\n", length(time_train), sum(status_train)))

cat("--- Preprocessing (D4: per-platform z-std + combined_rank top-3000 per-cohort) ---\n")
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
Y_train     <- pp$Y
train_genes <- pp$gene_names
cat(sprintf("  Training matrix: n=%d, p=%d\n\n", nrow(Y_train), ncol(Y_train)))

# --------------------------------------------------------------------------
# 2. Fit alpha_F sweep at fixed K=7
# --------------------------------------------------------------------------

cat(sprintf("--- Fitting alpha_F sweep at K=%d ---\n\n", K_D4))

fits    <- list()
results <- list()

for (alpha_f in ALPHA_F_GRID) {
  vname <- sprintf("alphaF_%.2f", alpha_f)
  cat(sprintf("--- %s ---\n", vname))
  set.seed(42L)
  fit <- suppressMessages(
    fit_cox_on_yf(
      Y_train, time_train, status_train,
      K = K_D4, max_iter = MAX_ITER, tol = cfg$cavi$tol,
      prior_LF = "point_exponential", prior_beta = PRIOR_BETA,
      alpha = ALPHA, alpha_F = alpha_f, verbose = TRUE
    )
  )
  fits[[vname]] <- fit

  k_eff    <- sum(abs(fit$EBeta) > BETA_THRESH)
  beta_max <- max(abs(fit$EBeta))
  n_iters  <- fit$history$n_iter
  cat(sprintf("  -> K_eff=%d | beta_max=%.4f | iters=%d\n\n", k_eff, beta_max, n_iters))

  results[[vname]] <- data.frame(
    variant  = vname, alpha_F = alpha_f, K = K_D4,
    k_eff = k_eff, beta_max = round(beta_max, 4), n_iters = n_iters,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# 3. External validation on all 5 held-out cohorts (matches D4 protocol
#    exactly: rank_transform=FALSE, per_platform_standardize=TRUE, top_n=NULL)
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
ext_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE
  )
  common <- intersect(train_genes, pre_ext$gene_names)
  if (length(common) < 100) {
    cat(sprintf("    Skipping %s: only %d common genes\n", ext_cohort, length(common)))
    next
  }
  Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
  train_idx <- match(common, train_genes)

  for (vname in names(fits)) {
    fit    <- fits[[vname]]
    EF_sub <- fit$EF[train_idx, , drop = FALSE]
    pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
    c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

    ext_rows[[length(ext_rows) + 1]] <- data.frame(
      variant = vname, cohort = ext_cohort, c_index = round(c_val, 4),
      n_common_genes = length(common),
      stringsAsFactors = FALSE
    )
  }
}

ext_tbl <- do.call(rbind, ext_rows)
# alpha_F wasn't stored on the fit object; recover it from the variant name.
ext_tbl$alpha_F <- as.numeric(sub("alphaF_", "", ext_tbl$variant))

# --------------------------------------------------------------------------
# 4. Summary
# --------------------------------------------------------------------------

mean_c <- aggregate(c_index ~ variant + alpha_F, data = ext_tbl, FUN = mean)
mean_c <- mean_c[order(mean_c$alpha_F), ]

cat("\n============================================================\n")
cat(" Diagnostic Summary — D4 preprocessing, alpha_F sweep\n")
cat(sprintf(" K=%d (fixed) | Training n=%d, p=%d\n", K_D4, nrow(Y_train), ncol(Y_train)))
cat(" Reference: current recommended config (alpha_F=0), mean external C = 0.627\n")
cat("============================================================\n")
print(mean_c, row.names = FALSE)
cat("============================================================\n")

fit_summary <- do.call(rbind, results)
write.csv(fit_summary, file.path(OUT_DIR, "fit_summary.csv"), row.names = FALSE)
write.csv(ext_tbl,     file.path(OUT_DIR, "external_validation.csv"), row.names = FALSE)
write.csv(mean_c,      file.path(OUT_DIR, "mean_c_by_alphaF.csv"), row.names = FALSE)

cat("\nResults saved to:\n")
cat(sprintf("  %s/fit_summary.csv\n", OUT_DIR))
cat(sprintf("  %s/external_validation.csv\n", OUT_DIR))
cat(sprintf("  %s/mean_c_by_alphaF.csv\n", OUT_DIR))
cat("============================================================\n")
