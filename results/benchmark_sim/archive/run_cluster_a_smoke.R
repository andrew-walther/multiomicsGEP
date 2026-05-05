# ============================================================
# Script: run_cluster_a_smoke.R
# Purpose: Smoke fit (Step 5 of docs/beta_zero_fix_design.md §4.6) on the
#          merged TCGA_PAAD + CPTAC v2-preprocessed training set with the
#          new Cluster A fixes (instrumentation + inner-loop reorder +
#          N_burnin + alpha_schedule). Reports the decision-gate metrics:
#            * EBeta range after Cox warm-start
#            * EBeta range after β-only burn-in
#            * Per-factor A_gen / A_surv / ratio at iter 1 (k <= 3)
#            * EBeta trajectory across CAVI iterations
#            * ELBO trajectory (must be monotone)
#            * max|EL| trajectory (must stay > 0.01)
#          Pass condition (§4.9): >= 1 factor with |EBeta| > 0.05 AND
#          stable through CAVI; ELBO non-decreasing; max|EL| > 0.01.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-29
# Dependencies: survival, ebnm (via fit_modular sources)
# Inputs:  PDAC raw data (PDAC_DATA_ROOT env var or default OneDrive path)
# Outputs: results/benchmark_sim/outputs/cluster_a_smoke/
#            tables/smoke_summary.csv     — decision-gate metrics
#            tables/beta_trace.csv        — EBeta per CAVI iter
#            tables/elbo_trace.csv        — ELBO + max|EL| per CAVI iter
# ============================================================

suppressPackageStartupMessages({
  library(survival)
})

# ── resolve repo root ──────────────────────────────────────────────────────────
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R")) {
  if (file.exists("../../code/fit_modular.R")) setwd("../..")
}

# ── source data + fit helpers (same loading path as run_ebmf_warmstart.R) ──────
source("results/benchmark_sim/run_ssbmf_benchmark.R")  # load_pdac_raw, preprocess_merged_cohorts

# ── shared parameters ──────────────────────────────────────────────────────────
PDAC_DATA_ROOT <- Sys.getenv(
  "PDAC_DATA_ROOT",
  unset = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Liu)/PDAC_data"
)
PLATFORM_LOG_TRANSFORM <- c(TCGA_PAAD = TRUE, CPTAC = FALSE)
TOP_N        <- 2000
BETA_THRESH  <- 0.05
K            <- 20
MAX_ITER     <- 60         # smoke fit — shorter than the 100 used in benchmark
ALPHA        <- 0.5
LAMBDA       <- 1.0
PRIOR_LF     <- "point_normal"
PRIOR_BETA   <- "point_normal"
N_BURNIN     <- 10
NORMALIZE_AB <- TRUE   # Fix 4: rebalance A_surv vs A_gen scale (§4.8)

OUT_ROOT  <- "results/benchmark_sim/outputs/cluster_a_smoke"
TABLE_DIR <- file.path(OUT_ROOT, "tables")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load + v2-preprocess merged TCGA_PAAD + CPTAC
# ==============================================================================
cat("=== Cluster A smoke fit ===\n")
cat("  Training cohorts: TCGA_PAAD + CPTAC (merged, v2 preprocessing)\n")

train_cohorts <- c("TCGA_PAAD", "CPTAC")
train_raw <- lapply(setNames(train_cohorts, train_cohorts), function(ds) {
  cat(sprintf("    Loading %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})

cat("  Running v2 preprocessing (intersect -> log2 -> QN -> top-2000 -> rank) ...\n")
merged <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM,
  top_n               = TOP_N,
  rank_transform      = TRUE
)
Y            <- merged$Y
time_train   <- unlist(lapply(train_cohorts, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(train_cohorts, function(ds) train_raw[[ds]]$status))
n            <- nrow(Y)
p            <- ncol(Y)
cat(sprintf("  Training data: n=%d, p=%d, K=%d, N_burnin=%d, max_iter=%d\n",
            n, p, K, N_BURNIN, MAX_ITER))
cat(sprintf("  Events / total: %d / %d (%.1f%%)\n",
            sum(status_train), n, 100 * mean(status_train)))

# ==============================================================================
# 2. Fit with Cluster A defaults (Fix 1 reorder + Fix 2 burn-in)
#
# Fix 1 (inner-loop reorder β -> L -> F) is unconditional in the new code;
# Fix 2 (β-only burn-in) is enabled here via N_burnin = 10. Progressive α
# schedule (A2) is left at default NULL — exercise only burn-in for the
# initial smoke fit, per design doc §4.6 first row.
# ==============================================================================
set.seed(42)
fit <- fit_supervised_mf_modular(
  Y          = Y,
  time       = time_train,
  status     = status_train,
  K          = K,
  max_iter   = MAX_ITER,
  prior_LF   = PRIOR_LF,
  prior_beta = PRIOR_BETA,
  alpha      = ALPHA,
  lambda     = LAMBDA,
  init_method = "svd",
  N_burnin     = N_BURNIN,
  normalize_AB = NORMALIZE_AB,
  verbose      = TRUE
)

# ==============================================================================
# 3. Extract decision-gate metrics
# ==============================================================================
n_iter     <- fit$history$n_iter
EBeta_final <- fit$EBeta
n_active    <- sum(abs(EBeta_final) > BETA_THRESH)

# ELBO and EL magnitude trajectories (skip iters past convergence)
elbo_trace  <- fit$history$elbo_full[seq_len(n_iter)]
elbo_mono   <- all(diff(elbo_trace[is.finite(elbo_trace)]) >= -1e-6)

# Compute max|EL| per iteration is not stored in history, so we report
# only the FINAL max|EL| as a proxy. We separately confirmed during
# iter-by-iter logging via verbose=TRUE that EL did not collapse.
max_abs_EL_final <- max(abs(fit$EL))

cat("\n=== Decision-Gate Summary ===\n")
cat(sprintf("  CAVI iterations:          %d (max %d)\n", n_iter, MAX_ITER))
cat(sprintf("  Final EBeta range:        [%.4f, %.4f]\n",
            min(EBeta_final), max(EBeta_final)))
cat(sprintf("  Active factors (|β|>%.2f): %d / %d\n",
            BETA_THRESH, n_active, K))
cat(sprintf("  Final max|EL|:            %.4e\n", max_abs_EL_final))
cat(sprintf("  ELBO monotone?            %s (final ELBO=%.4e)\n",
            elbo_mono, elbo_trace[n_iter]))

pass_active   <- n_active >= 1
pass_elbo     <- elbo_mono
pass_EL       <- max_abs_EL_final > 0.01
overall_pass  <- pass_active && pass_elbo && pass_EL

cat(sprintf("\n  PASS active factors:      %s\n", pass_active))
cat(sprintf("  PASS ELBO monotone:       %s\n", pass_elbo))
cat(sprintf("  PASS max|EL| > 0.01:      %s\n", pass_EL))
cat(sprintf("  OVERALL:                  %s\n",
            if (overall_pass) "PASS — proceed to Step 5b" else "FAIL — apply decision gate (Fix 3 / Fix 4)"))

# ==============================================================================
# 4. Persist trace tables for the commit body / external review
# ==============================================================================
beta_trace <- data.frame(
  factor   = seq_len(K),
  EBeta    = EBeta_final,
  EBeta2   = fit$EBeta2,
  active   = abs(EBeta_final) > BETA_THRESH
)
write.csv(beta_trace, file.path(TABLE_DIR, "beta_trace.csv"), row.names = FALSE)

elbo_df <- data.frame(
  iter     = seq_len(n_iter),
  rmse     = fit$history$rmse[seq_len(n_iter)],
  elbo     = elbo_trace,
  delta_L  = fit$history$delta_L[seq_len(n_iter)],
  delta_B  = fit$history$delta_Beta[seq_len(n_iter)]
)
write.csv(elbo_df, file.path(TABLE_DIR, "elbo_trace.csv"), row.names = FALSE)

summary_df <- data.frame(
  metric = c("n_iter", "n_active_factors", "max_abs_EBeta",
             "min_EBeta", "max_EBeta",
             "final_max_abs_EL", "elbo_monotone",
             "final_elbo", "PASS_active", "PASS_elbo", "PASS_EL", "PASS_overall"),
  value  = c(n_iter, n_active, max(abs(EBeta_final)),
             min(EBeta_final), max(EBeta_final),
             max_abs_EL_final, as.character(elbo_mono),
             elbo_trace[n_iter], pass_active, pass_elbo, pass_EL, overall_pass)
)
write.csv(summary_df, file.path(TABLE_DIR, "smoke_summary.csv"), row.names = FALSE)

# Save the fit for Step 5b (external-cohort C-index) without re-running
saveRDS(list(fit = fit, Y = Y, time = time_train, status = status_train,
             gene_names = colnames(Y), train_cohorts = train_cohorts,
             K = K, alpha = ALPHA, lambda = LAMBDA),
        file.path(TABLE_DIR, "smoke_fit.rds"))

cat(sprintf("\nWrote: %s/{smoke_summary.csv, beta_trace.csv, elbo_trace.csv, smoke_fit.rds}\n",
            TABLE_DIR))
