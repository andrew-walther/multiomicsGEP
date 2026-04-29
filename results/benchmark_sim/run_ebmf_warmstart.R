# ============================================================
# Script: run_ebmf_warmstart.R
# Purpose: EBMF warm-start diagnostic — two experiments to distinguish
#          whether SSBMF's β=0 failure is due to (a) broken β update or
#          (b) poor initialization / early CAVI dynamics.
#
#          Background: SSBMF trained on merged TCGA_PAAD + CPTAC gives
#          EBeta = 0 for all K factors. The EBMF diagnostic (run_ebmf_diagnostic.R)
#          confirmed that 5/20 unsupervised EBMF factors are Cox-significant
#          (C-index up to 0.629). So the data HAS survival signal; SSBMF
#          cannot find it. Two hypotheses remain:
#
#          H1 — Initialization/dynamics: CAVI starts from SVD and the survival
#               gradient in the L update is too weak to pull factors toward
#               the biologically-informative EBMF directions. The β update
#               itself is fine, but it never gets factors worth selecting.
#
#          H2 — Broken β update: even given factors known to be survival-
#               associated, the β CAVI update cannot assign non-zero
#               coefficients (prior too aggressive, gradient wrong, or bug).
#
# Experiments:
#   1. β-ONLY UPDATE — Hold EL fixed at the EBMF loading matrix. Run only
#      the β update (update_beta_k) iteratively for 30 iterations. If β
#      becomes non-zero: H1 is correct (β update works, init is the problem).
#      If β stays zero: H2 is correct (β update is broken).
#
#   2. FULL CAVI WARM-START — Initialize EL and EF from the EBMF solution,
#      then run the complete CAVI loop normally. If β is non-zero and CAVI
#      maintains the EBMF structure: warm-starting is a viable fix. If β
#      collapses back to zero: the L update is washing out the survival signal.
#
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-29
# Dependencies: flashier (already fitted), survival
# Inputs:  results/benchmark_sim/outputs/ebmf_diagnostic/tables/ebmf_fit.rds
#          code/fit_modular.R, code/update_beta.R
# Outputs: results/benchmark_sim/outputs/ebmf_warmstart/
#            tables/exp1_beta_trace.csv       — EBeta per iteration (β-only)
#            tables/exp1_beta_final.csv       — final EBeta + Cox c-index (β-only)
#            tables/exp2_warmstart_beta.csv   — final EBeta from full CAVI
#            tables/exp2_warmstart_elbo.csv   — ELBO trace from full CAVI
#            tables/exp2_warmstart_beta_summary.csv — per-factor β + active status
# ============================================================

suppressPackageStartupMessages({
  library(flashier)
  library(survival)
  library(pheatmap)
})

# ── resolve repo root ──────────────────────────────────────────────────────────
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R")) {
  if (file.exists("../../code/fit_modular.R")) setwd("../..")
}

# ── source all code dependencies ───────────────────────────────────────────────
# run_ssbmf_benchmark.R sources fit_modular.R (via tryCatch to suppress the
# entry-point block's stopifnot), update_beta.R, and all other dependencies.
# The if (sys.nframe() == 0) guard prevents the benchmark entry-point from
# running when sourced.  We source it unconditionally here so that all helper
# functions (load_pdac_raw, calc_cox_taylor, compute_z_no_k, update_beta_k,
# fit_supervised_mf_modular, preprocess_merged_cohorts) are in scope.
source("results/benchmark_sim/run_ssbmf_benchmark.R")
source("results/benchmark_sim/run_phase1_diagnostics.R")  # plot_cohort_loading_heatmap

# ── output directories ─────────────────────────────────────────────────────────
OUT_ROOT  <- "results/benchmark_sim/outputs/ebmf_warmstart"
TABLE_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR   <- file.path(OUT_ROOT, "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

# ── shared parameters (must match the EBMF diagnostic run) ─────────────────────
PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT",
                             unset = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/PDAC_data")
PLATFORM_LOG_TRANSFORM <- c(TCGA_PAAD = TRUE, CPTAC = FALSE)
TOP_N       <- 2000
BETA_THRESH <- 0.05   # |β| > BETA_THRESH counts as "active"

# ==============================================================================
# 1. Load saved EBMF results and re-derive training data
# ==============================================================================
cat("=== EBMF Warm-Start: loading data ===\n")

ebmf_rds <- "results/benchmark_sim/outputs/ebmf_diagnostic/tables/ebmf_fit.rds"
if (!file.exists(ebmf_rds))
  stop(sprintf("EBMF fit not found at %s — run run_ebmf_diagnostic.R first.", ebmf_rds))

ebmf_obj    <- readRDS(ebmf_rds)
flash_fit   <- ebmf_obj$flash_fit   # full flashier fit object
L_ebmf_unit <- ebmf_obj$L          # n × K, unit-L2-norm columns (from ldf type="2")
D_ebmf      <- ebmf_obj$D          # K-vector: scale per factor
time_train  <- ebmf_obj$time
status_train <- ebmf_obj$status

# Re-derive Y via v2 preprocessing (fast — no modelling, just data transforms).
# Must match the preprocessing used for run_ebmf_diagnostic.R exactly.
train_cohorts <- c("TCGA_PAAD", "CPTAC")
train_raw <- lapply(setNames(train_cohorts, train_cohorts), function(ds) {
  cat(sprintf("  Loading %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})

cat("  Running v2 preprocessing (intersect → log2 → QN → top-2000 → rank) ...\n")
merged <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM,
  top_n               = TOP_N,
  rank_transform      = TRUE
)
Y             <- merged$Y
cohort_labels <- merged$dataset_labels
n             <- nrow(Y)
p             <- ncol(Y)

# Extract gene-loading matrix F from flash_fit.
# ldf(type="2") gives unit-L2-norm L and F columns scaled by D so that
# L_unit %*% diag(D) %*% t(F_unit) ≈ Y.
ldf_res  <- ldf(flash_fit, type = "2")
F_ebmf_unit <- ldf_res$F   # p × K, unit-L2-norm columns

K <- ncol(L_ebmf_unit)   # number of EBMF factors retained

cat(sprintf("  Training data: n=%d, p=%d | EBMF factors: K=%d\n", n, p, K))
cat(sprintf("  Event rate: %.1f%%\n", 100 * mean(status_train)))

# ==============================================================================
# Experiment 1: β-ONLY UPDATE
#
# Fix EL = L_ebmf and run only the β CAVI update for N_ITER_BETA iterations.
# EL2 is set to EL^2 (treating the EBMF posterior means as point estimates —
# zero posterior variance). This isolates the β update from any L/F dynamics.
#
# Interpretation:
#   EBeta moves away from 0 → β update works; initialization is the culprit
#   EBeta stays at 0      → β update is broken regardless of the factors
# ==============================================================================
cat("\n=================================================================\n")
cat(" EXPERIMENT 1: β-ONLY UPDATE (EL fixed at EBMF loadings)\n")
cat("=================================================================\n")

N_ITER_BETA <- 30
PRIOR_BETA  <- "point_normal"
ALPHA       <- 0.5

#' Run the β CAVI update in isolation with EL fixed.
#'
#' EL is held constant at the EBMF loading matrix throughout. Only
#' update_beta_k() is called at each iteration; L, F, and τ are never updated.
#'
#' @param EL         n × K matrix: fixed subject loadings (EBMF L)
#' @param time       numeric n-vector: follow-up times
#' @param status     integer n-vector: event indicators (1 = event)
#' @param prior_beta character: EBNM prior for β (default "point_normal")
#' @param alpha      numeric in [0,1]: survival mixing weight (default 0.5)
#' @param n_iter     number of β-update iterations (default 30)
#' @param verbose    print EBeta each iteration? (default TRUE)
#' @return list(EBeta, EBeta2, beta_trace n_iter × K matrix)
run_beta_only <- function(EL, time, status,
                          prior_beta = "point_normal",
                          alpha      = 0.5,
                          n_iter     = 30,
                          verbose    = TRUE) {
  K_loc  <- ncol(EL)
  # EL2 = EL^2: zero posterior variance — treats EBMF means as fixed constants.
  EL2    <- EL^2
  EBeta  <- rep(0, K_loc)
  EBeta2 <- rep(0, K_loc)

  beta_trace <- matrix(NA_real_, nrow = n_iter, ncol = K_loc,
                       dimnames = list(NULL, paste0("EBMF", seq_len(K_loc))))

  for (iter in seq_len(n_iter)) {
    # Cox Taylor expansion at current linear predictor η = EL β
    eta    <- as.vector(EL %*% EBeta)
    taylor <- calc_cox_taylor(eta, time, status)
    z      <- eta + taylor$u / taylor$w   # working response
    w      <- taylor$w                    # diagonal weights W_{ii}

    for (k in seq_len(K_loc)) {
      z_no_k   <- compute_z_no_k(z, EL, EBeta, k)
      res      <- update_beta_k(w, z_no_k, EL[, k], EL2[, k],
                                prior_family = prior_beta, alpha = alpha)
      EBeta[k]  <- res$mean
      EBeta2[k] <- res$second
    }

    beta_trace[iter, ] <- EBeta

    if (verbose) {
      n_active <- sum(abs(EBeta) > BETA_THRESH)
      cat(sprintf("  [iter %2d] active=%d  max|β|=%.4f  β=[%s]\n",
                  iter, n_active, max(abs(EBeta)),
                  paste(sprintf("%+.3f", EBeta), collapse = ", ")))
    }
  }

  list(EBeta = EBeta, EBeta2 = EBeta2, beta_trace = beta_trace)
}

exp1 <- run_beta_only(
  EL         = L_ebmf_unit,
  time       = time_train,
  status     = status_train,
  prior_beta = PRIOR_BETA,
  alpha      = ALPHA,
  n_iter     = N_ITER_BETA,
  verbose    = TRUE
)

# Compute C-index for the final β-only risk score
risk_exp1 <- as.vector(L_ebmf_unit %*% exp1$EBeta)
cindex_exp1 <- tryCatch({
  surv_obj <- Surv(time_train, status_train)
  conc     <- concordance(surv_obj ~ risk_exp1)
  conc$concordance
}, error = function(e) NA_real_)

n_active_exp1 <- sum(abs(exp1$EBeta) > BETA_THRESH)
cat(sprintf("\n  Experiment 1 summary:\n"))
cat(sprintf("    Active factors (|β| > %.2f): %d / %d\n", BETA_THRESH, n_active_exp1, K))
cat(sprintf("    Max |β|: %.6f\n", max(abs(exp1$EBeta))))
cat(sprintf("    C-index (training): %.4f\n", cindex_exp1))

# Save results
beta_trace_df <- as.data.frame(exp1$beta_trace)
beta_trace_df$iter <- seq_len(nrow(beta_trace_df))
write.csv(beta_trace_df,
          file.path(TABLE_DIR, "exp1_beta_trace.csv"), row.names = FALSE)

beta_final_df <- data.frame(
  Factor      = paste0("EBMF", seq_len(K)),
  EBeta       = round(exp1$EBeta,  6),
  Abs_EBeta   = round(abs(exp1$EBeta), 6),
  Active      = abs(exp1$EBeta) > BETA_THRESH,
  Cox_p_value = NA_real_   # placeholder; filled below
)
# Annotate with Cox p-values from the original EBMF diagnostic
ebmf_cox <- read.csv("results/benchmark_sim/outputs/ebmf_diagnostic/tables/ebmf_cox_summary.csv",
                     stringsAsFactors = FALSE)
beta_final_df$Cox_p_value <- ebmf_cox$p_value[match(beta_final_df$Factor, ebmf_cox$Factor)]
write.csv(beta_final_df,
          file.path(TABLE_DIR, "exp1_beta_final.csv"), row.names = FALSE)

# ==============================================================================
# Experiment 2: FULL CAVI WARM-START
#
# Initialize SSBMF EL and EF from the EBMF posterior means (flash_fit$L_pm,
# flash_fit$F_pm), then run the full CAVI loop normally — all updates active.
#
# Scaling: flash_fit$L_pm %*% t(flash_fit$F_pm) ≈ Y by construction, so this
# is a properly scaled warm start (EL %*% t(EF) ≈ Y from the first iteration).
#
# prior_LF is set to "point_normal" (allows negative values) rather than the
# default "point_exponential" (non-negative only), because EBMF factors are
# signed and pmax-clipping them to 0 would destroy the initialisation.
#
# Interpretation:
#   β non-zero + C-index improves → warm-start works; poor init was the culprit
#   β collapses back to 0          → L update washes out the survival signal
#   L matrix drifts far from EBMF  → CAVI abandons the warm-start direction
# ==============================================================================
cat("\n=================================================================\n")
cat(" EXPERIMENT 2: FULL CAVI WARM-START (EL, EF from EBMF)\n")
cat("=================================================================\n")

# flash_fit$L_pm: n × K posterior means of L (columns scale with factor magnitude)
# flash_fit$F_pm: p × K posterior means of F
EL_init_exp2 <- flash_fit$L_pm   # n × K
EF_init_exp2 <- flash_fit$F_pm   # p × K

cat(sprintf("  EL_init: %d x %d  (range [%.3f, %.3f])\n",
            nrow(EL_init_exp2), ncol(EL_init_exp2),
            min(EL_init_exp2), max(EL_init_exp2)))
cat(sprintf("  EF_init: %d x %d  (range [%.3f, %.3f])\n",
            nrow(EF_init_exp2), ncol(EF_init_exp2),
            min(EF_init_exp2), max(EF_init_exp2)))

fit_ws <- fit_supervised_mf_modular(
  Y          = Y,
  time       = time_train,
  status     = status_train,
  K          = K,
  max_iter   = 100,
  tol        = 1e-5,
  prior_LF   = "point_normal",    # signed factors; matches EBMF assumption
  prior_beta = PRIOR_BETA,
  alpha      = ALPHA,
  lambda     = 1.0,
  EL_init    = EL_init_exp2,
  EF_init    = EF_init_exp2,
  verbose    = TRUE
)

n_active_ws <- sum(abs(fit_ws$EBeta) > BETA_THRESH)
risk_ws     <- as.vector(fit_ws$EL %*% fit_ws$EBeta)
cindex_ws   <- tryCatch({
  conc <- concordance(Surv(time_train, status_train) ~ risk_ws)
  conc$concordance
}, error = function(e) NA_real_)

cat(sprintf("\n  Experiment 2 summary:\n"))
cat(sprintf("    Converged: %s (%d iters)\n",
            fit_ws$history$converged, fit_ws$history$n_iter))
cat(sprintf("    Active factors (|β| > %.2f): %d / %d\n", BETA_THRESH, n_active_ws, K))
cat(sprintf("    Max |β|: %.6f\n", max(abs(fit_ws$EBeta))))
cat(sprintf("    C-index (training): %.4f\n", cindex_ws))
cat(sprintf("    Final ELBO: %.2f\n", tail(fit_ws$history$elbo_full, 1)))

# Save results
ws_elbo_df <- data.frame(
  iter       = seq_along(fit_ws$history$elbo_full),
  elbo_full  = fit_ws$history$elbo_full,
  elbo_proxy = fit_ws$history$elbo_proxy,
  delta_L    = fit_ws$history$delta_L,
  delta_Beta = fit_ws$history$delta_Beta
)
write.csv(ws_elbo_df,
          file.path(TABLE_DIR, "exp2_warmstart_elbo.csv"), row.names = FALSE)

ws_beta_df <- data.frame(
  Factor     = paste0("EBMF", seq_len(K)),
  EBeta      = round(fit_ws$EBeta,  6),
  Abs_EBeta  = round(abs(fit_ws$EBeta), 6),
  Active     = abs(fit_ws$EBeta) > BETA_THRESH
)
ws_beta_df$Cox_p_value <- ebmf_cox$p_value[match(ws_beta_df$Factor, ebmf_cox$Factor)]
write.csv(ws_beta_df,
          file.path(TABLE_DIR, "exp2_warmstart_beta.csv"), row.names = FALSE)

# Loading heatmap: warm-start final L with Cox p-value annotation
plot_cohort_loading_heatmap(
  EL            = fit_ws$EL,
  EBeta         = fit_ws$EBeta,
  cohort_labels = factor(cohort_labels),
  out_stub      = file.path(FIG_DIR, "exp2_warmstart_loading_heatmap"),
  beta_thresh   = BETA_THRESH
)

# ==============================================================================
# Summary comparison
# ==============================================================================
cat("\n=================================================================\n")
cat(" WARM-START DIAGNOSTIC SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("  Exp 1 (β-only, EL fixed): %d / %d active | max|β|=%.4f | C-index=%.4f\n",
            n_active_exp1, K, max(abs(exp1$EBeta)), cindex_exp1))
cat(sprintf("  Exp 2 (full CAVI ws):     %d / %d active | max|β|=%.4f | C-index=%.4f\n",
            n_active_ws,   K, max(abs(fit_ws$EBeta)), cindex_ws))

cat("\n  Interpretation:\n")
if (n_active_exp1 > 0) {
  cat("  ✓ Exp 1 β ≠ 0 → β update CAN select survival factors given EBMF loadings.\n")
  cat("    Root cause of SSBMF failure: initialization / early CAVI dynamics.\n")
  cat("    Fix direction: EBMF pre-warming or improved initialization strategy.\n")
} else {
  cat("  ✗ Exp 1 β = 0 even with EBMF loadings fixed.\n")
  cat("    Root cause: β CAVI update is broken (prior too aggressive, gradient\n")
  cat("    wrong, or code bug). Debug update_beta_k() directly.\n")
}
if (n_active_ws > 0) {
  cat("  ✓ Exp 2 β ≠ 0 → full CAVI warm-start works.\n")
  cat("    EBMF initialization is a viable production strategy.\n")
} else {
  cat("  ✗ Exp 2 β = 0 after full CAVI. L/F updates are washing out survival signal.\n")
  cat("    Investigate update_L_k() — survival gradient in A_surv may be dominated\n")
  cat("    by the genomics reconstruction gradient.\n")
}

summary_df <- data.frame(
  Experiment      = c("Exp1_beta_only", "Exp2_full_warmstart"),
  n_active_beta   = c(n_active_exp1, n_active_ws),
  max_abs_beta    = round(c(max(abs(exp1$EBeta)), max(abs(fit_ws$EBeta))), 6),
  cindex_training = round(c(cindex_exp1, cindex_ws), 4),
  stringsAsFactors = FALSE
)
write.csv(summary_df,
          file.path(TABLE_DIR, "warmstart_summary.csv"), row.names = FALSE)

cat(sprintf("\n  Results written to %s\n", OUT_ROOT))
