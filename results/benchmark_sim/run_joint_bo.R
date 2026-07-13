# ============================================================
# Script:  results/benchmark_sim/run_joint_bo.R
# Purpose: Step 3 of the K-parsimony follow-up plan
#          (docs/plans/ssbmf_k_parsimony_followup_plan_07_13_2026.md).
#          Jointly tunes (K, alpha) via Bayesian optimization
#          (code/select_k_alpha_bo.R) on the D4 recommended config, testing
#          whether DeSurv's own approach -- joint (k, alpha, penalty) BO
#          rather than our fixed-alpha, K-only CV -- finds a smaller K with
#          comparable external performance. Motivated by Step 1/2's negative
#          results: neither multistart nor deflation-init (both fixed at
#          alpha=0.5) rescue K=4/K=5 on their own; this asks whether varying
#          alpha itself reshapes the fitting landscape rather than just the
#          starting point within it.
#
#          Sanity check (plan Step 3 item 4): the BO result should be at
#          least as good, on the internal CV objective, as the existing
#          fixed-alpha K-CV grid's K=7 internal score (DECISIONS.md
#          2026-07-12 "Fresh K-CV" entry: mean C=0.633, SE=0.029 at alpha=0.5).
#
#          Output: results/benchmark_sim/outputs/joint_bo/
#            joint_bo_history.csv        (every (K, alpha) point evaluated)
#            joint_bo_external_val.csv   (external validation of the BO winner)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_joint_bo.R
#          caffeinate -i Rscript results/benchmark_sim/run_joint_bo.R --quick
# Requires: PDAC_DATA_ROOT set (real data not in git); rBayesianOptimization
#           package (install.packages("rBayesianOptimization")).
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({
  library(yaml); library(survival); library(rBayesianOptimization)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"), error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/select_K.R")
source("code/select_k_alpha_bo.R")
source("code/preprocess_desurv.R")

b <- cfg$benchmark
p <- cfg$preprocessing

MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
N_FOLDS      <- if (QUICK_MODE) 3L  else 5L
INIT_POINTS  <- if (QUICK_MODE) 3L  else 8L
N_ITER       <- if (QUICK_MODE) 2L  else 15L
K_BOUNDS     <- c(2L, 10L)
ALPHA_BOUNDS <- c(0, 1)
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
TOP_N_DESURV <- p$top_n_genes_desurv

K7_INTERNAL_CV_REF <- 0.633   # DECISIONS.md 2026-07-12, alpha=0.5 fixed K-CV at K=7

OUT_DIR <- "results/benchmark_sim/outputs/joint_bo"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | N_FOLDS: %d | init_points: %d | n_iter: %d\n",
            QUICK_MODE, MAX_ITER, N_FOLDS, INIT_POINTS, N_ITER))

# --------------------------------------------------------------------------
# 1. Load and preprocess training data (D4 config)
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d, events=%d\n\n", length(time_train), sum(status_train)))

cat("--- Preprocessing (D4 config) ---\n")
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
cat(sprintf("  n=%d, p=%d\n\n", nrow(Y_train), ncol(Y_train)))

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. Joint (K, alpha) Bayesian optimization
# --------------------------------------------------------------------------

cat("=== Joint (K, alpha) Bayesian optimization ===\n\n")
bo_res <- select_k_alpha_bayesopt(
  Y_train, time_train, status_train,
  K_bounds = K_BOUNDS, alpha_bounds = ALPHA_BOUNDS,
  n_folds = N_FOLDS, seed = 42L, model = "YFB",
  init_points = INIT_POINTS, n_iter = N_ITER, bo_seed = 42L,
  prior_beta = PRIOR_BETA, max_iter = MAX_ITER, verbose = TRUE)

write.csv(bo_res$History, file.path(OUT_DIR, "joint_bo_history.csv"), row.names = FALSE)

cat(sprintf("\nBO raw Best_Par: K=%d, alpha=%.4f, internal CV mean C=%.4f\n",
            round(bo_res$Best_Par[["K"]]), bo_res$Best_Par[["alpha"]], bo_res$Best_Value))

# --------------------------------------------------------------------------
# 2b. Validity gate: BO's raw Best_Par can be a degenerate alpha extreme that
#     scores well by incidental unsupervised-reconstruction alignment with
#     the outcome (K_eff=0), not genuine survival modeling -- a documented
#     failure mode in this project (DECISIONS.md 2026-05-05). Re-check the
#     top 5 candidates and require K_eff > 0 before trusting a winner.
# --------------------------------------------------------------------------

cat("\n--- Validity gate: re-checking top 5 candidates for K_eff > 0 ---\n")
winner <- pick_trustworthy_bo_winner(
  bo_res$History, Y_train, time_train, status_train,
  n_candidates = 5, beta_threshold = BETA_THRESH,
  prior_beta = PRIOR_BETA, max_iter = MAX_ITER)

print(winner$candidates_checked)

K_bo     <- winner$K
alpha_bo <- winner$alpha
cv_bo    <- winner$cv_value
fit_bo   <- winner$fit   # already fit during the validity check -- reused below, no need to refit

cat(sprintf("\nTrustworthy BO winner: K=%d, alpha=%.4f, internal CV mean C=%.4f, K_eff=%d\n",
            K_bo, alpha_bo, cv_bo, winner$k_eff))

# --------------------------------------------------------------------------
# 3. Sanity check vs. the existing fixed-alpha K-CV grid (plan Step 3 item 4)
# --------------------------------------------------------------------------

cat(sprintf("\n--- Sanity check: BO internal CV (%.4f) vs. fixed-alpha K=7 K-CV reference (%.4f) ---\n",
            cv_bo, K7_INTERNAL_CV_REF))
if (cv_bo < K7_INTERNAL_CV_REF - 0.02) {
  warning(sprintf(
    "BO's internal CV score (%.4f) is more than 0.02 below the existing K=7/alpha=0.5 K-CV
reference (%.4f) -- the search may not be finding a competitive region; inspect joint_bo_history.csv
before trusting K_bo/alpha_bo.", cv_bo, K7_INTERNAL_CV_REF))
} else {
  cat("BO's internal CV score is competitive with (not worse than) the existing marginal optimum.\n")
}

# --------------------------------------------------------------------------
# 4. External validation, using the full-training-data fit pick_trustworthy_bo_winner()
#    already produced (fit_bo <- winner$fit above) -- no need to refit.
# --------------------------------------------------------------------------

cat(sprintf("\n=== BO winner already fit on full training data (K=%d, alpha=%.4f) ===\n",
            K_bo, alpha_bo))
cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
            sum(abs(fit_bo$EBeta) > BETA_THRESH), max(abs(fit_bo$EBeta)),
            fit_bo$history$n_iter))

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows     <- list()

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

  EF_sub <- fit_bo$EF[train_idx, , drop = FALSE]
  pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit_bo$EBeta, EF_norms = fit_bo$EF_norms)
  c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

  results_rows[[length(results_rows) + 1]] <- data.frame(
    K = K_bo, alpha = round(alpha_bo, 4), cohort = ext_cohort, c_index = round(c_val, 4),
    k_eff = sum(abs(fit_bo$EBeta) > BETA_THRESH), beta_max = round(max(abs(fit_bo$EBeta)), 4),
    n_iter = fit_bo$history$n_iter, n_common_genes = length(common),
    stringsAsFactors = FALSE
  )
}

if (length(results_rows) == 0)
  stop("No external validation rows: all cohorts had < 100 common genes.")

results <- do.call(rbind, results_rows)
out_csv <- file.path(OUT_DIR, "joint_bo_external_val.csv")
write.csv(results, out_csv, row.names = FALSE)

cat("\n============================================================\n")
cat(" Joint (K, alpha) BO winner: external validation\n")
cat("============================================================\n")
cat(sprintf("  K=%d, alpha=%.4f: mean C=%.4f (SE=%.4f), K_eff=%d\n",
            K_bo, alpha_bo, mean(results$c_index), sd(results$c_index) / sqrt(nrow(results)),
            unique(results$k_eff)))
cat(sprintf("\nFor comparison (Step 1, DECISIONS.md 2026-07-13):\n"))
cat("  K=7 (alpha=0.5, fresh SVD):        mean C=0.6267 (SE=0.0199), K_eff=2\n")
cat("  K=4 (alpha=0.5, warm-start-from-7): mean C=0.6270 (SE=0.0198), K_eff=2\n")
cat("  K=5 (alpha=0.5, warm-start-from-7): mean C=0.6270 (SE=0.0198), K_eff=2\n")
cat(sprintf("\nResults: %s\n", out_csv))
cat("============================================================\n")
