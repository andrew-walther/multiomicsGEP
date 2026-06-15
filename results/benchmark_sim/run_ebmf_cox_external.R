# ============================================================
# Script:  results/benchmark_sim/run_ebmf_cox_external.R
# Purpose: Real-data UNSUPERVISED 2-step baseline (EBMF -> Cox) evaluated on the
#          same external-validation protocol as the supervised LB / YFB models in
#          run_desurv_comparison.R.
#
#          Pipeline (identical preprocessing, split, and scoring to YFB):
#            1. Load merged TCGA_PAAD + CPTAC training data.
#            2. DeSurv-aligned preprocessing (per-platform z-std + combined-rank
#               top-3000-per-cohort gene selection); the SAME settings used for
#               the recommended YFB model (D4) in run_desurv_comparison.R.
#            3. Fit EBMF (flashier, var_type=2, greedy + backfit) -> gene weights
#               F_ebmf and factor scores L_ebmf. This is UNSUPERVISED: survival
#               plays no role in learning F.
#            4. Fit survival::coxph on the EBMF factor scores -> beta_ebmf (2-step).
#            5. Score each of the 5 held-out cohorts by the EXACT same projection
#               used for YFB: eta_new = (Y_new F_norm) beta. The only difference
#               from YFB is that F here is learned unsupervised.
#            6. Per-cohort external C-index via the oriented concordance used in
#               run_desurv_comparison.R (training-concordance sign correction).
#
#          FAIL LOUD: if the Cox arm collapses (all |beta| ~ 0) or the training
#          oriented concordance is at chance, the script stops with an explicit
#          error rather than writing a degenerate result silently.
#
#   Output: results/benchmark_sim/outputs/ebmf_cox_external/
#     ebmf_cox_external_results.csv   (tidy, columns aligned with desurv CSV)
#     ebmf_cox_external_fit.rds       (flashier fit, F, beta, for post-hoc use)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-15
# Dependencies: flashier (>= 1.0.0), survival, yaml
# Usage:   PDAC_DATA_ROOT=<path> caffeinate -i Rscript \
#            results/benchmark_sim/run_ebmf_cox_external.R
#          PDAC_DATA_ROOT=<path> Rscript \
#            results/benchmark_sim/run_ebmf_cox_external.R --quick
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

# Navigate to project root if invoked from a subdirectory
if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
  library(flashier)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")  # load_pdac_raw, PLATFORM_LOG_TRANSFORM
source("code/preprocess_desurv.R")                   # preprocess_merged_cohorts, preprocess_desurv_cohort

b <- cfg$benchmark
p <- cfg$preprocessing

# DeSurv-aligned settings — must match the recommended YFB config (D4) so the
# baseline is evaluated under identical preprocessing and gene selection.
TOP_N_DESURV     <- p$top_n_genes_desurv   # 3000 per cohort
BETA_THRESH      <- cfg$k_selection$beta_threshold
K_EBMF_MAX       <- if (QUICK_MODE) 5L else b$k_pdac   # greedy ceiling (flashier may use fewer)
TRAIN_COHORTS    <- cfg$pdac$training_cohorts
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

OUT_DIR <- "results/benchmark_sim/outputs/ebmf_cox_external"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | K_EBMF_MAX: %d | top_n: %d\n",
            QUICK_MODE, K_EBMF_MAX, TOP_N_DESURV))

# --------------------------------------------------------------------------
# Helper: oriented C-index (same definition as run_desurv_comparison.R)
# --------------------------------------------------------------------------
#' Oriented concordance: orient the risk score so higher = worse prognosis.
#'
#' Cox maximises the partial likelihood, so the training concordance is > 0.5
#' by construction; this max(c, 1-c) form applies the training-concordance sign
#' correction at scoring time, matching the supervised LB / YFB external scoring.
#'
#' @param risk    numeric vector of linear-predictor risk scores.
#' @param time    numeric vector of follow-up times.
#' @param status  integer event indicator (1 = event).
#' @return scalar oriented C-index in [0.5, 1].
oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# Helper: YFB-style projection scoring  eta = (Y F_norm) beta
# --------------------------------------------------------------------------
#' Project an expression matrix onto unit-normalised EBMF gene weights.
#'
#' Mirrors predict_cox_on_yf(): EF columns are scaled to unit L2 norm using the
#' TRAINING column norms (F_norms), so a row-subset of F for external validation
#' stays on the training scale. Used for BOTH the training Cox fit and external
#' scoring, guaranteeing train/test consistency.
#'
#' @param Ymat    n x m expression matrix (m = number of genes used).
#' @param Fmat    m x K gene-weight matrix (row-subset of the full F is allowed).
#' @param F_norms length-K training column norms of the full F.
#' @return n x K matrix of projection scores.
project_scores <- function(Ymat, Fmat, F_norms) {
  Ymat %*% sweep(Fmat, 2, F_norms, "/")
}

# --------------------------------------------------------------------------
# 1. Load + preprocess merged TCGA_PAAD + CPTAC training data (DeSurv-aligned)
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (%s), events=%d\n",
            length(time_train),
            paste(sprintf("%s=%d", TRAIN_COHORTS,
                          vapply(train_raw, function(x) x$n, integer(1))), collapse = ", "),
            sum(status_train)))

cat("\n--- DeSurv-aligned preprocessing (per-platform z-std + combined-rank top-3000/cohort) ---\n")
pp <- preprocess_merged_cohorts(
  cohort_raw_list          = train_raw,
  log_transform_flags      = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n                    = TOP_N_DESURV,
  rank_transform           = FALSE,         # matches D4 training preprocessing
  per_platform_standardize = TRUE,
  normalize_method         = "none",
  selection_per_cohort     = TRUE,
  selection_method         = "combined_rank"
)
Y_train     <- pp$Y
train_genes <- pp$gene_names
cat(sprintf("  Training matrix: n=%d, p=%d genes\n", nrow(Y_train), ncol(Y_train)))

# --------------------------------------------------------------------------
# 2. Fit EBMF (flashier) — unsupervised; survival not used here
# --------------------------------------------------------------------------

cat(sprintf("\n--- Fitting EBMF (flashier, var_type=2, greedy_Kmax=%d, backfit) ---\n",
            K_EBMF_MAX))
set.seed(42L)
# var_type = 2: per-gene (per-column) residual variance — each gene has its own
# noise level, matching the genomics heteroscedasticity. greedy + backfit lets
# flashier choose how many factors to keep (>= 1; may be < K_EBMF_MAX).
flash_fit <- flash(
  Y_train,
  var_type    = 2,
  greedy_Kmax = K_EBMF_MAX,
  backfit     = TRUE,
  verbose     = 1
)
K_ebmf <- flash_fit$n_factors
if (is.null(K_ebmf) || K_ebmf < 1)
  stop("EBMF returned no factors — cannot build a 2-step baseline.")
cat(sprintf("  EBMF converged with K=%d factors (greedy ceiling %d).\n",
            K_ebmf, K_EBMF_MAX))

# LDF decomposition, type "2": unit-L2-norm L and F columns, scale absorbed by D.
ldf_res <- ldf(flash_fit, type = "2")
F_ebmf  <- ldf_res$F   # p x K gene weights (unit-norm columns)
if (nrow(F_ebmf) != ncol(Y_train))
  stop(sprintf("F_ebmf has %d rows but training data has %d genes.",
               nrow(F_ebmf), ncol(Y_train)))
F_norms <- sqrt(colSums(F_ebmf^2))          # ~1 by construction; kept for the
F_norms <- pmax(F_norms, 1e-12)             # exact YFB-style scoring contract

# --------------------------------------------------------------------------
# 3. Fit Cox on EBMF factor scores (2-step), on the SAME projection used at test
# --------------------------------------------------------------------------

cat("\n--- Cox on EBMF factor scores (2-step) ---\n")
S_train <- project_scores(Y_train, F_ebmf, F_norms)   # n x K
colnames(S_train) <- paste0("EBMF", seq_len(K_ebmf))

cox_fit <- coxph(Surv(time_train, status_train) ~ S_train)
beta_ebmf <- as.numeric(coef(cox_fit))
beta_ebmf[is.na(beta_ebmf)] <- 0           # rank-deficient columns -> 0 (loud below)
names(beta_ebmf) <- colnames(S_train)

k_eff_cox <- sum(abs(beta_ebmf) > BETA_THRESH)
cat(sprintf("  Cox coefficients: %d / %d with |beta| > %.3g | max|beta|=%.4f\n",
            k_eff_cox, K_ebmf, BETA_THRESH, max(abs(beta_ebmf))))

# FAIL LOUD: a fully collapsed Cox arm is not a usable baseline.
if (max(abs(beta_ebmf)) < 1e-6)
  stop("EBMF->Cox COLLAPSED: all Cox coefficients ~ 0. ",
       "The 2-step baseline produced no survival signal on the training data.")

# Training-concordance sanity / sign-correction check.
risk_train  <- as.numeric(S_train %*% beta_ebmf)
c_train     <- oriented_cindex(risk_train, time_train, status_train)
cat(sprintf("  Training oriented C-index = %.4f\n", c_train))
if (c_train < 0.55)
  stop(sprintf("EBMF->Cox training concordance at chance (C=%.3f). ", c_train),
       "The fitted model has no usable training survival signal — reporting this ",
       "would be a degenerate baseline; investigate before scoring externally.")

# --------------------------------------------------------------------------
# 4. External validation on the 5 held-out cohorts (identical protocol to YFB)
# --------------------------------------------------------------------------

cat("\n--- External validation (5 held-out cohorts) ---\n")
results_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)

  # top_n = NULL: keep all external genes; the intersection with train_genes
  # controls the final gene set (mirrors run_desurv_comparison.R).
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort
  )

  common <- intersect(train_genes, pre_ext$gene_names)
  if (length(common) < 100) {
    cat(sprintf("    Skipping %s: only %d common genes\n", ext_cohort, length(common)))
    next
  }

  Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
  train_idx <- match(common, train_genes)
  F_sub     <- F_ebmf[train_idx, , drop = FALSE]

  S_ext <- project_scores(Y_ext, F_sub, F_norms)
  risk  <- as.numeric(S_ext %*% beta_ebmf)
  c_val <- oriented_cindex(risk, raw_ext$time, raw_ext$status)

  cat(sprintf("    common genes=%d | external C=%.4f\n", length(common), c_val))

  results_rows[[length(results_rows) + 1]] <- data.frame(
    model          = "EBMF_COX",
    label          = "EBMF->Cox (unsupervised)",
    model_type     = "EBMF",
    cohort         = ext_cohort,
    c_index        = round(c_val, 4),
    K              = K_ebmf,
    k_eff          = k_eff_cox,
    beta_max       = round(max(abs(beta_ebmf)), 4),
    top_n          = TOP_N_DESURV,
    sel_method     = "combined_rank",
    per_cohort     = TRUE,
    n_common_genes = length(common),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# 5. Save and report
# --------------------------------------------------------------------------

if (length(results_rows) == 0)
  stop("No external validation rows: all cohorts had < 100 common genes.")

results <- do.call(rbind, results_rows)
out_csv <- file.path(OUT_DIR, "ebmf_cox_external_results.csv")
write.csv(results, out_csv, row.names = FALSE)

saveRDS(
  list(flash_fit = flash_fit, F_ebmf = F_ebmf, F_norms = F_norms,
       beta_ebmf = beta_ebmf, cox_fit = cox_fit, K_ebmf = K_ebmf,
       train_genes = train_genes, c_train = c_train),
  file.path(OUT_DIR, "ebmf_cox_external_fit.rds")
)

mean_c <- mean(results$c_index)

# FAIL LOUD: a mean external C at chance means the baseline is degenerate.
if (mean_c < 0.52)
  stop(sprintf("EBMF->Cox external mean C=%.3f is at chance — degenerate baseline. ",
               mean_c),
       "Results were still written for inspection, but treat as a COLLAPSE.")

cat(sprintf("\n=== Results saved: %s ===\n", out_csv))
cat(sprintf("EBMF->Cox: K=%d | Cox-active=%d | training C=%.3f | mean external C=%.3f\n",
            K_ebmf, k_eff_cox, c_train, mean_c))
cat("Per-cohort external C-index:\n")
for (i in seq_len(nrow(results)))
  cat(sprintf("  %-20s C=%.3f  (%d common genes)\n",
              results$cohort[i], results$c_index[i], results$n_common_genes[i]))
cat("\nDone.\n")
