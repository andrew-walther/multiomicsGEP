# ============================================================
# Script:  results/benchmark_sim/run_merged_kcv.R
# Purpose: Select K via cross-validated C-index (1-SE rule) on merged
#          TCGA_PAAD + CPTAC training data for three preprocessing x
#          model configurations:
#            (a) LB    x joint quantile+rank   (for M1, M2)
#            (b) LB    x per-platform z-std    (for M3, M4)
#            (c) YFB   x per-platform z-std    (for M5, M6)
#
#          Results are written into config/globals.yml under
#          benchmark.k_merged_lb_joint, benchmark.k_merged_lb_perplatform,
#          benchmark.k_merged_yfb_perplatform. run_merged_benchmark.R
#          reads these values.
#
#          NOTE on preprocessing leakage: preprocessing is applied to the
#          full merged matrix before CV folds split. Held-out fold data
#          therefore influences quantile normalization / z-std parameters.
#          This is standard practice for K selection (not final model
#          parameters) and the effect is small, but it is not leakage-free.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_kcv.R [--quick]
#          --quick: K_grid=2:5, n_folds=3, max_iter=50 (~5 min)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../../")
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
N_FOLDS   <- if (QUICK_MODE) 3L else cfg$cavi$n_cv_folds    # 5
MAX_ITER  <- if (QUICK_MODE) 50L else cfg$cavi$max_iter      # 300
ALPHA     <- cfg$benchmark$alpha                              # 0.5
LAMBDA    <- cfg$benchmark$lambda                             # 1.0
PRIOR_BETA <- "normal"

cat("=== Merged-Cohort K-CV ===\n")
cat(sprintf("    K_grid=%d:%d | n_folds=%d | max_iter=%d | QUICK=%s\n\n",
            min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER, QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load raw training data
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))

time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))

cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n", n_tcga + n_cptac,
            n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Preprocess — two versions
# --------------------------------------------------------------------------

cat("--- Preprocessing (joint quantile+rank) ---\n")
merged_joint <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n               = cfg$preprocessing$top_n_genes,
  rank_transform      = TRUE,
  per_platform_standardize = FALSE
)
Y_joint <- merged_joint$Y
cat(sprintf("  Y_joint: n=%d, p=%d\n\n", nrow(Y_joint), ncol(Y_joint)))

cat("--- Preprocessing (per-platform z-std) ---\n")
merged_perplatform <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n               = cfg$preprocessing$top_n_genes,
  rank_transform      = FALSE,
  per_platform_standardize = TRUE
)
Y_perplatform <- merged_perplatform$Y
cat(sprintf("  Y_perplatform: n=%d, p=%d\n\n", nrow(Y_perplatform), ncol(Y_perplatform)))

# --------------------------------------------------------------------------
# 3. K-CV for each configuration
# --------------------------------------------------------------------------

run_kcv <- function(Y, time, status, model, label) {
  cat(sprintf("--- K-CV: %s ---\n", label))
  cat(sprintf("    model=%s | K_grid=%d:%d | n_folds=%d | max_iter=%d\n\n",
              model, min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER))
  set.seed(42L)
  result <- select_K_cv(
    Y, time, status,
    K_grid     = K_GRID,
    n_folds    = N_FOLDS,
    model      = model,
    seed       = 42L,
    verbose    = TRUE,
    max_iter   = MAX_ITER,
    prior_beta = PRIOR_BETA,
    alpha      = ALPHA,
    lambda     = LAMBDA,
    sign_correction = FALSE    # must be FALSE in CV folds
  )
  cat(sprintf("\n  -> K_opt=%d (1-SE rule)\n", result$K_opt))
  cat("\n  CV table:\n")
  print(result$cv_table[, c("K", "mean_cindex", "se_cindex")], row.names = FALSE)
  out_csv <- file.path(OUT_DIR, sprintf("kcv_%s.csv", gsub(" ", "_", tolower(label))))
  write.csv(result$cv_table, out_csv, row.names = FALSE)
  cat(sprintf("  Saved: %s\n\n", out_csv))
  result
}

kcv_lb_joint        <- run_kcv(Y_joint,        time_train, status_train, "LB",  "LB joint")
kcv_lb_perplatform  <- run_kcv(Y_perplatform,  time_train, status_train, "LB",  "LB per-platform")
kcv_yfb_perplatform <- run_kcv(Y_perplatform,  time_train, status_train, "YFB", "YFB per-platform")

# --------------------------------------------------------------------------
# 4. Write K values back into globals.yml
# --------------------------------------------------------------------------

cat("--- Updating config/globals.yml with CV-selected K values ---\n")

globals_text <- readLines("config/globals.yml")

replace_null_key <- function(lines, key, value) {
  pattern <- paste0("^(\\s*", key, ":\\s*)null")
  idx <- grep(pattern, lines)
  if (length(idx) == 0) stop(sprintf("Key '%s' not found in globals.yml", key))
  lines[idx] <- sub("null", as.character(value), lines[idx])
  lines
}

globals_text <- replace_null_key(globals_text, "k_merged_lb_joint",
                                  kcv_lb_joint$K_opt)
globals_text <- replace_null_key(globals_text, "k_merged_lb_perplatform",
                                  kcv_lb_perplatform$K_opt)
globals_text <- replace_null_key(globals_text, "k_merged_yfb_perplatform",
                                  kcv_yfb_perplatform$K_opt)

writeLines(globals_text, "config/globals.yml")
cat(sprintf("  k_merged_lb_joint        = %d\n", kcv_lb_joint$K_opt))
cat(sprintf("  k_merged_lb_perplatform  = %d\n", kcv_lb_perplatform$K_opt))
cat(sprintf("  k_merged_yfb_perplatform = %d\n", kcv_yfb_perplatform$K_opt))
cat("  globals.yml updated.\n\n")

cat("============================================================\n")
cat(" K-CV complete. Run run_merged_benchmark.R to fit all models.\n")
cat("============================================================\n")
