# ============================================================
# Script:  results/benchmark_sim/run_yfb_dualF_frozen_grid_D4.R
# Purpose: Step 2 (Experiment B) of
#          docs/plans/yfb_dual_source_F_experiments_08_20_2026.md.
#
#          Tests whether letting beta warm up (via the existing N_frozen
#          frozen-F preconditioning, DECISIONS.md 2026-05-22) BEFORE
#          activating dual-source F (alpha_F > 0) lets the survival gradient
#          actually improve on genomics-only F -- rather than the cold-start
#          alpha_F sweep in run_yfb_dualF_diagnostic_D4.R (Step 0), which
#          gave the survival term a non-zero weight from iteration 1 and
#          found it statistically flat vs. the alpha_F=0 baseline.
#
#          N_frozen and alpha_F already compose correctly in fit_cox_on_yf.R
#          with NO code changes needed: during iter <= N_frozen, F is held
#          fixed at its SVD init (beta and L update freely); the moment F
#          unfreezes at iter = N_frozen+1, it immediately uses whatever
#          alpha_F was passed. This combination (N_frozen>0 AND alpha_F>0
#          together) was never tested in the original 2026-05-22 diagnostic
#          (which only varied them one at a time).
#
#          Grid: N_frozen in {0, 10, 20, 30} x alpha_F in {0, 0.1, 0.3, 0.5},
#          under the D4 preprocessing (per-platform z-std, combined_rank
#          top-3000 per-cohort, K=7 fixed), full 5-cohort external validation.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   Rscript results/benchmark_sim/run_yfb_dualF_frozen_grid_D4.R [--quick]
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

N_FROZEN_GRID <- c(0L, 10L, 20L, 30L)
ALPHA_F_GRID  <- c(0, 0.1, 0.3, 0.5)

if (is.null(K_D4))
  stop("config/globals.yml benchmark$k_merged_yfb_desurv is NULL -- D4's K has not been set.")

OUT_DIR <- "results/benchmark_sim/outputs/yfb_dualF_frozen_grid_D4"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== YFB Dual-Source F Diagnostic (D4 preprocessing, N_frozen x alpha_F grid) ===\n")
cat(sprintf("    K=%d (fixed) | alpha=%.2f | max_iter=%d\n", K_D4, ALPHA, MAX_ITER))
cat(sprintf("    N_frozen grid: %s | alpha_F grid: %s\n\n",
            paste(N_FROZEN_GRID, collapse = ", "), paste(ALPHA_F_GRID, collapse = ", ")))

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 1. Load + preprocess training data (D4 preprocessing, same as Step 0)
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
# 2. Fit the N_frozen x alpha_F grid at fixed K=7
# --------------------------------------------------------------------------

grid <- expand.grid(N_frozen = N_FROZEN_GRID, alpha_F = ALPHA_F_GRID)
cat(sprintf("--- Fitting %d (N_frozen, alpha_F) combinations at K=%d ---\n\n",
            nrow(grid), K_D4))

fits    <- list()
results <- list()

for (gi in seq_len(nrow(grid))) {
  nf <- grid$N_frozen[gi]; af <- grid$alpha_F[gi]
  vname <- sprintf("Nf%02d_aF%.2f", nf, af)
  cat(sprintf("--- %s (N_frozen=%d, alpha_F=%.2f) ---\n", vname, nf, af))
  set.seed(42L)
  fit <- suppressMessages(
    fit_cox_on_yf(
      Y_train, time_train, status_train,
      K = K_D4, max_iter = MAX_ITER, tol = cfg$cavi$tol,
      prior_LF = "point_exponential", prior_beta = PRIOR_BETA,
      alpha = ALPHA, alpha_F = af, N_frozen = nf, verbose = TRUE
    )
  )
  fits[[vname]] <- fit

  k_eff    <- sum(abs(fit$EBeta) > BETA_THRESH)
  beta_max <- max(abs(fit$EBeta))
  n_iters  <- fit$history$n_iter
  cat(sprintf("  -> K_eff=%d | beta_max=%.4f | iters=%d\n\n", k_eff, beta_max, n_iters))

  results[[vname]] <- data.frame(
    variant = vname, N_frozen = nf, alpha_F = af, K = K_D4,
    k_eff = k_eff, beta_max = round(beta_max, 4), n_iters = n_iters,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# 3. External validation on all 5 held-out cohorts (matches D4 protocol)
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
ext_tbl <- merge(ext_tbl, grid_labeled <- data.frame(
  variant  = names(fits),
  N_frozen = grid$N_frozen,
  alpha_F  = grid$alpha_F,
  stringsAsFactors = FALSE
), by = "variant")

# --------------------------------------------------------------------------
# 4. Summary
# --------------------------------------------------------------------------

mean_c <- aggregate(c_index ~ variant + N_frozen + alpha_F, data = ext_tbl, FUN = mean)
mean_c <- mean_c[order(mean_c$N_frozen, mean_c$alpha_F), ]

cat("\n============================================================\n")
cat(" Diagnostic Summary — D4 preprocessing, N_frozen x alpha_F grid\n")
cat(sprintf(" K=%d (fixed) | Training n=%d, p=%d\n", K_D4, nrow(Y_train), ncol(Y_train)))
cat(" Reference: current recommended config (N_frozen=0, alpha_F=0), mean external C = 0.627\n")
cat("============================================================\n")
print(mean_c, row.names = FALSE)
cat("============================================================\n")

fit_summary <- do.call(rbind, results)
write.csv(fit_summary, file.path(OUT_DIR, "fit_summary.csv"), row.names = FALSE)
write.csv(ext_tbl,     file.path(OUT_DIR, "external_validation.csv"), row.names = FALSE)
write.csv(mean_c,      file.path(OUT_DIR, "mean_c_by_config.csv"), row.names = FALSE)

cat("\nResults saved to:\n")
cat(sprintf("  %s/fit_summary.csv\n", OUT_DIR))
cat(sprintf("  %s/external_validation.csv\n", OUT_DIR))
cat(sprintf("  %s/mean_c_by_config.csv\n", OUT_DIR))
cat("============================================================\n")
