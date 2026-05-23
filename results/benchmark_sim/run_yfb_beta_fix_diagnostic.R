# ============================================================
# Script:  results/benchmark_sim/run_yfb_beta_fix_diagnostic.R
# Purpose: Diagnose and fix YFB β→0 collapse on merged TCGA+CPTAC.
#
#          Step 1. K-CV on merged training data to select K for YFB
#                  (replaces the previous K=20 + ARD approach).
#                  K-CV uses C-index as selection criterion (not ELBO:
#                  ELBO always favours larger K regardless of survival
#                  signal, making it unsuitable for K selection here).
#                  The reference run uses V1 (cox_warmstart=TRUE) since
#                  V0 collapses to β=0 in all folds (C=0.5 everywhere),
#                  which gives an uninformative K-CV curve.
#
#          Step 2. Fit all variants at K_cv and diagnose β→0:
#            V0: Baseline — cox_warmstart=FALSE, N_burnin=0, α=0.50
#            V1: Cox warm-start — cox_warmstart=TRUE, N_burnin=0, α=0.50
#            V2: β burn-in — cox_warmstart=FALSE, N_burnin=10, α=0.50
#            V3: Combined — cox_warmstart=TRUE, N_burnin=10, α=0.50
#            V4: High α — cox_warmstart=FALSE, N_burnin=0, α=0.75
#            V5: High α + warm-start — cox_warmstart=TRUE, N_burnin=0, α=0.75
#
#          Step 3. Quick external validation on Dijk (proxy for all 5
#                  cohorts — fastest external cohort to load).
#
# Key diagnostic printed at iter 1:
#   A_surv/A_gen ratio — the scale imbalance that drives β→0.
#   Target: ratio > 0.01 at iter 1 (comparable to single-cohort TCGA).
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-22
# Usage:   Rscript results/benchmark_sim/run_yfb_beta_fix_diagnostic.R [--quick]
#          --quick: K_grid=2:5, max_iter=50, 3 folds (fast sanity check ~5 min)
# ============================================================

# --------------------------------------------------------------------------
# 0. Setup
# --------------------------------------------------------------------------

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
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
source("code/predict.R"); source("code/predict_cox_on_yf.R")
source("code/select_K.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

K_GRID    <- if (QUICK_MODE) 2L:5L else 2L:10L
N_FOLDS   <- if (QUICK_MODE) 3L else 5L
MAX_ITER  <- if (QUICK_MODE) 50L else 150L
PRIOR_BETA <- "normal"
ALPHA_BASE <- 0.50
ALPHA_HIGH <- 0.75   # p/(n+p) ≈ 2000/2273 ≈ 0.88 would be extreme; test 0.75 first

cat("=== YFB β→0 Diagnostic ===\n")
cat(sprintf("    K_grid=%s | n_folds=%d | max_iter=%d | QUICK=%s\n\n",
            paste(range(K_GRID), collapse = ":"), N_FOLDS, MAX_ITER, QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/yfb_beta_fix"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load training data
# --------------------------------------------------------------------------

cat("--- Loading merged TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- c("TCGA_PAAD", "CPTAC")
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  Loading %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})

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
n_tcga       <- train_raw$TCGA_PAAD$n
n_cptac      <- train_raw$CPTAC$n

cat(sprintf("  Training: n=%d (TCGA=%d, CPTAC=%d), p=%d, events=%d (%.0f%% censored)\n\n",
            nrow(Y_train), n_tcga, n_cptac, ncol(Y_train),
            sum(status_train), 100 * mean(status_train == 0)))

# --------------------------------------------------------------------------
# Helper: oriented C-index
# --------------------------------------------------------------------------

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. K-CV on merged data (reference variant: V1, cox_warmstart=TRUE)
# --------------------------------------------------------------------------
# V0 (baseline) collapses β→0 in every fold, giving C=0.5 for all K.
# K-CV on V0 is uninformative. We use V1 (cox_warmstart) as the reference
# because it is the minimal change needed to escape β→0 and thus produces
# a meaningful C-index curve for K selection.

cat("--- Step 1: K-CV on merged data (reference: cox_warmstart=TRUE) ---\n")
cat(sprintf("    K_grid=%d:%d, %d folds, max_iter=%d\n\n",
            min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER))

set.seed(42L)
kcv_ref <- select_K_cv(
  Y_train, time_train, status_train,
  K_grid   = K_GRID,
  n_folds  = N_FOLDS,
  model    = "YFB",
  seed     = 42L,
  verbose  = TRUE,
  # variant V1 parameters (the reference):
  max_iter      = MAX_ITER,
  prior_beta    = PRIOR_BETA,
  alpha         = ALPHA_BASE,
  cox_warmstart = TRUE,
  N_burnin      = 0L
)

K_cv <- kcv_ref$K_opt
cat(sprintf("\n  K-CV selected K = %d (1-SE rule)\n", K_cv))
cat("\n  K-CV table (V1 reference):\n")
print(kcv_ref$cv_table[, c("K", "mean_cindex", "se_cindex")], row.names = FALSE)
write.csv(kcv_ref$cv_table,
          file.path(OUT_DIR, "kcv_merged_V1_reference.csv"),
          row.names = FALSE)

# --------------------------------------------------------------------------
# 3. Define variants
# --------------------------------------------------------------------------

variants <- list(
  V0 = list(label = "Baseline",
            cox_warmstart = FALSE, N_burnin = 0L,  alpha = ALPHA_BASE, alpha_F = 0),
  V1 = list(label = "Cox warm-start",
            cox_warmstart = TRUE,  N_burnin = 0L,  alpha = ALPHA_BASE, alpha_F = 0),
  V2 = list(label = "β burn-in (N=10)",
            cox_warmstart = FALSE, N_burnin = 10L, alpha = ALPHA_BASE, alpha_F = 0),
  V3 = list(label = "Warm-start + burn-in",
            cox_warmstart = TRUE,  N_burnin = 10L, alpha = ALPHA_BASE, alpha_F = 0),
  V4 = list(label = "High α=0.75",
            cox_warmstart = FALSE, N_burnin = 0L,  alpha = ALPHA_HIGH, alpha_F = 0),
  V5 = list(label = "High α=0.75 + warm-start",
            cox_warmstart = TRUE,  N_burnin = 0L,  alpha = ALPHA_HIGH, alpha_F = 0),
  # Dual-source F variants (alpha_F > 0 activates survival gradient in F update)
  V6 = list(label = "Dual-F α_F=0.1 + warm-start",
            cox_warmstart = TRUE,  N_burnin = 0L,  alpha = ALPHA_BASE, alpha_F = 0.1),
  V7 = list(label = "Dual-F α_F=0.3 + warm-start",
            cox_warmstart = TRUE,  N_burnin = 0L,  alpha = ALPHA_BASE, alpha_F = 0.3),
  V8 = list(label = "Dual-F α_F=0.5 + warm-start",
            cox_warmstart = TRUE,  N_burnin = 0L,  alpha = ALPHA_BASE, alpha_F = 0.5)
)

# --------------------------------------------------------------------------
# 4. Load Dijk for quick external validation
# --------------------------------------------------------------------------

cat("\n--- Loading Dijk (external validation proxy) ---\n")
dijk_raw   <- load_pdac_raw("Dijk", PDAC_DATA_ROOT)
dijk_genes <- dijk_raw$gene_names
common     <- intersect(train_genes, dijk_genes)
cat(sprintf("  Dijk: n=%d, common genes with training=%d\n\n",
            dijk_raw$n, length(common)))

# --------------------------------------------------------------------------
# 5. Fit all variants at K_cv and evaluate
# --------------------------------------------------------------------------

cat(sprintf("--- Step 2: Fit all variants at K=%d ---\n\n", K_cv))

results <- vector("list", length(variants))
fits    <- list()

for (vi in seq_along(variants)) {
  v    <- variants[[vi]]
  vname <- names(variants)[vi]
  cat(sprintf("--- %s: %s ---\n", vname, v$label))
  cat(sprintf("    cox_warmstart=%s | N_burnin=%d | alpha=%.2f | alpha_F=%.2f\n",
              v$cox_warmstart, v$N_burnin, v$alpha, v$alpha_F))

  set.seed(42L)
  fit <- suppressMessages(
    fit_cox_on_yf(Y_train, time_train, status_train,
                  K             = K_cv,
                  max_iter      = MAX_ITER,
                  prior_LF      = "point_exponential",
                  prior_beta    = PRIOR_BETA,
                  alpha         = v$alpha,
                  cox_warmstart = v$cox_warmstart,
                  N_burnin      = v$N_burnin,
                  alpha_F       = v$alpha_F,
                  verbose       = TRUE)
  )
  fits[[vname]] <- fit

  # --- training concordance ---
  zf_train  <- Y_train %*% (fit$EF / pmax(fit$EF_norms, 1e-10))
  eta_train <- as.vector(zf_train %*% fit$EBeta)
  c_train   <- oriented_cindex(eta_train, time_train, status_train)

  # --- external validation on Dijk ---
  # Subset EF to common genes (same gene order as dijk raw)
  ef_sub <- fit$EF[match(common, train_genes), , drop = FALSE]
  Y_dijk <- dijk_raw$Y[, match(common, dijk_genes), drop = FALSE]

  # Column-center Dijk using training gene means over common genes
  train_means_common <- colMeans(Y_train[, match(common, train_genes), drop = FALSE])
  Y_dijk_c <- sweep(Y_dijk, 2, train_means_common, "-")

  pred_dijk <- predict_cox_on_yf(Y_dijk_c, ef_sub, fit$EBeta,
                                  EF_norms = fit$EF_norms)
  c_dijk <- oriented_cindex(pred_dijk$risk_scores,
                             dijk_raw$time, dijk_raw$status)

  # fit_cox_on_yf does not return K_eff; derive from |EBeta| > 0.01 threshold
  k_eff_val <- sum(abs(fit$EBeta) > 0.01)
  beta_max_val <- if (length(fit$EBeta) > 0) max(abs(fit$EBeta)) else 0.0

  n_iters <- sum(fit$history$elbo_full != 0)
  results[[vi]] <- data.frame(
    variant       = vname,
    label         = v$label,
    K_cv          = K_cv,
    K_eff         = k_eff_val,
    beta_max      = round(beta_max_val, 4),
    alpha         = v$alpha,
    alpha_F       = v$alpha_F,
    cox_warmstart = v$cox_warmstart,
    N_burnin      = v$N_burnin,
    iters         = n_iters,
    C_train       = round(c_train, 3),
    C_dijk        = round(c_dijk, 3),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  → K_eff=%d | beta_max=%.4f | iters=%d | C_train=%.3f | C_dijk=%.3f\n\n",
              k_eff_val, beta_max_val, n_iters, c_train, c_dijk))
}

# --------------------------------------------------------------------------
# 6. Summary table
# --------------------------------------------------------------------------

res_tbl <- do.call(rbind, results)

cat("============================================================\n")
cat(" Diagnostic Summary — YFB β→0 Fix Strategies\n")
cat(sprintf(" K_cv (1-SE rule on merged data) = %d\n", K_cv))
cat(sprintf(" Training: n=%d (TCGA=%d, CPTAC=%d), p=%d\n",
            nrow(Y_train), n_tcga, n_cptac, ncol(Y_train)))
cat("============================================================\n")
cat(sprintf("%-6s %-32s %5s %5s %8s %7s %5s %5s %5s\n",
            "Var", "Label", "K_eff", "β_max", "C_train", "C_Dijk", "α", "α_F", "Brnin"))
cat(strrep("-", 78), "\n")
for (i in seq_len(nrow(res_tbl))) {
  r <- res_tbl[i, ]
  cat(sprintf("%-6s %-32s %5d %5.4f %7.3f %7.3f %5.2f %5.2f %5d\n",
              r$variant, r$label, r$K_eff, r$beta_max,
              r$C_train, r$C_dijk, r$alpha, r$alpha_F, r$N_burnin))
}
cat("============================================================\n")

# Reference
cat("\nBaseline reference (main-branch YFB merged, prior_β=normal):\n")
cat("  YFB_base: K_eff=0, β→0, C-ext ≈ 0.54 (Dijk: 0.506)\n")
cat("  YFB_perplatform (Phase 1 fix): C-ext Dijk ≈ 0.573 (target to beat)\n\n")

# Save
write.csv(res_tbl,
          file.path(OUT_DIR, "yfb_beta_fix_results.csv"),
          row.names = FALSE)

# Compact RDS: save EBeta and convergence history only
# K_eff derived from |EBeta| threshold (fit_cox_on_yf does not return K_eff)
fits_compact <- lapply(fits, function(f)
  list(EBeta = f$EBeta, EBeta2 = f$EBeta2, EF_norms = f$EF_norms,
       K_eff = sum(abs(f$EBeta) > 0.01), history = f$history))
saveRDS(fits_compact,
        file.path(OUT_DIR, "yfb_beta_fix_fits_compact.rds"))

cat("Results saved to:\n")
cat(sprintf("  %s/yfb_beta_fix_results.csv\n", OUT_DIR))
cat(sprintf("  %s/kcv_merged_V1_reference.csv\n", OUT_DIR))
cat("============================================================\n")
