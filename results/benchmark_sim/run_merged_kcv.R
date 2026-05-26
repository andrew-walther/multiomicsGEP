# ============================================================
# Script:  results/benchmark_sim/run_merged_kcv.R
# Purpose: Select K via cross-validated C-index (1-SE rule, biological floor)
#          on merged TCGA_PAAD + CPTAC training data for all preprocessing x
#          model configurations used in the merged benchmark.
#
#          Biological K floor: K_final = max(K_1se, K_MIN_BIOLOGICAL) where
#          K_MIN_BIOLOGICAL = 3. Motivation: DeSurv (Young et al. 2025) uses
#          K=3-4; the model identifies biological programs (plural requires K>=3);
#          for YFB per-platform, 1-SE selects K=2 but K=3 has higher CV mean C
#          (0.634 vs 0.625) and lower fold SE (0.022 vs 0.030).
#
#          Preprocessing configurations:
#            joint_quantile_rank    - joint QN + rank transform (M1, M2)
#            perplatform_zstd       - per-platform z-std + joint QN (M3-M6)
#            joint_quantile_norank  - joint QN, no rank transform (M7, M8, M13, M14)
#            joint_zstd             - joint z-standardization (M9, M10, M15, M16)
#            log_only               - log transform only (M11, M12, M17, M18)
#
#          Results written to config/globals.yml under benchmark.k_merged_*.
#          Previously computed values are overwritten (not just nulls).
#
#          NOTE on leakage: preprocessing applied to full training matrix before
#          CV fold split. This is standard for K selection, not final inference.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_kcv.R [--quick]
#          --quick: K_grid=2:4, n_folds=3, max_iter=30 (~3 min)
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

K_GRID            <- if (QUICK_MODE) 2L:4L  else 2L:10L
N_FOLDS           <- if (QUICK_MODE) 3L      else cfg$cavi$n_cv_folds    # 5
MAX_ITER          <- if (QUICK_MODE) 30L     else cfg$cavi$max_iter      # 300
ALPHA             <- cfg$benchmark$alpha                                  # 0.5
LAMBDA            <- cfg$benchmark$lambda                                 # 1.0
PRIOR_BETA        <- "normal"
K_MIN_BIOLOGICAL  <- 3L   # biological floor: K_final = max(K_1se, K_MIN_BIOLOGICAL)

cat("=== Merged-Cohort K-CV (all preprocessing options) ===\n")
cat(sprintf("    K_grid=%d:%d | n_folds=%d | max_iter=%d | K_floor=%d | QUICK=%s\n\n",
            min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER, K_MIN_BIOLOGICAL, QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load raw training data (done once; preprocessing varies below)
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga        <- train_raw$TCGA_PAAD$n
n_cptac       <- train_raw$CPTAC$n
time_train    <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train  <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Preprocessing configuration table
#    Each row: label, globals_key_lb, globals_key_yfb, preprocess_args
# --------------------------------------------------------------------------

PREPROC_CONFIGS <- list(
  list(
    label       = "joint_quantile_rank",
    key_lb      = "k_merged_lb_joint",
    key_yfb     = NULL,  # YFB x joint QN excluded (beta->0 structural, see DECISIONS.md 2026-05-22)
    per_plat    = FALSE,
    norm_method = "quantile",
    rank        = TRUE
  ),
  list(
    label       = "perplatform_zstd",
    key_lb      = "k_merged_lb_perplatform",
    key_yfb     = "k_merged_yfb_perplatform",
    per_plat    = TRUE,
    norm_method = "quantile",
    rank        = FALSE
  ),
  list(
    label       = "joint_quantile_norank",
    key_lb      = "k_merged_lb_joint_norank",
    key_yfb     = "k_merged_yfb_joint_norank",
    per_plat    = FALSE,
    norm_method = "quantile",
    rank        = FALSE
  ),
  list(
    label       = "joint_zstd",
    key_lb      = "k_merged_lb_zstd",
    key_yfb     = "k_merged_yfb_zstd",
    per_plat    = FALSE,
    norm_method = "z_score",
    rank        = FALSE
  ),
  list(
    label       = "log_only",
    key_lb      = "k_merged_lb_logonly",
    key_yfb     = "k_merged_yfb_logonly",
    per_plat    = FALSE,
    norm_method = "none",
    rank        = FALSE
  )
)

# --------------------------------------------------------------------------
# Helper: write any value (null or integer) back to globals.yml key
# --------------------------------------------------------------------------

set_key <- function(lines, key, value) {
  pattern <- paste0("^\\s*", key, ":")
  idx     <- grep(pattern, lines)
  if (length(idx) == 0) stop(sprintf("Key '%s' not found in globals.yml", key))
  # Replace the value token (first non-whitespace after the colon) with the new value,
  # preserving any trailing comment.
  lines[idx] <- sub(
    paste0("(^\\s*", key, ":\\s*)\\S+"),
    paste0("\\1", as.character(value)),
    lines[idx]
  )
  lines
}

# --------------------------------------------------------------------------
# Helper: run K-CV for one (preprocessing, model) combination
# --------------------------------------------------------------------------

run_kcv <- function(Y, time, status, model, label) {
  cat(sprintf("--- K-CV: %s x %s ---\n", label, model))
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
    sign_correction = FALSE
  )
  K_1se   <- result$K_opt
  K_final <- max(K_1se, K_MIN_BIOLOGICAL)
  if (K_final > K_1se) {
    cat(sprintf("\n  NOTE: 1-SE selected K=%d; biological floor (K>=%d) applied -> K_final=%d\n",
                K_1se, K_MIN_BIOLOGICAL, K_final))
  } else {
    cat(sprintf("\n  -> K_opt=%d (1-SE rule; no floor applied)\n", K_1se))
  }
  cat("\n  CV table:\n")
  print(result$cv_table[, c("K", "mean_cindex", "se_cindex")], row.names = FALSE)
  out_csv <- file.path(OUT_DIR, sprintf("kcv_%s_%s.csv",
                                        gsub(" ", "_", tolower(label)),
                                        tolower(model)))
  write.csv(result$cv_table, out_csv, row.names = FALSE)
  cat(sprintf("  Saved: %s\n\n", out_csv))
  list(K_opt = K_final, K_1se = K_1se, cv_table = result$cv_table)
}

# --------------------------------------------------------------------------
# 3. Run K-CV for all configurations, write results back to globals.yml
# --------------------------------------------------------------------------

globals_text <- readLines("config/globals.yml")
k_results    <- list()

for (pcfg in PREPROC_CONFIGS) {

  cat(sprintf("\n========== Preprocessing: %s ==========\n\n", pcfg$label))

  Y_prep <- preprocess_merged_cohorts(
    cohort_raw_list          = train_raw,
    log_transform_flags      = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
    top_n                    = cfg$preprocessing$top_n_genes,
    rank_transform           = pcfg$rank,
    per_platform_standardize = pcfg$per_plat,
    normalize_method         = pcfg$norm_method
  )$Y
  cat(sprintf("  Y: n=%d, p=%d\n\n", nrow(Y_prep), ncol(Y_prep)))

  if (!is.null(pcfg$key_lb)) {
    r_lb <- run_kcv(Y_prep, time_train, status_train, "LB", pcfg$label)
    k_results[[pcfg$key_lb]] <- r_lb$K_opt
    globals_text <- set_key(globals_text, pcfg$key_lb, r_lb$K_opt)
  }
  if (!is.null(pcfg$key_yfb)) {
    r_yfb <- run_kcv(Y_prep, time_train, status_train, "YFB", pcfg$label)
    k_results[[pcfg$key_yfb]] <- r_yfb$K_opt
    globals_text <- set_key(globals_text, pcfg$key_yfb, r_yfb$K_opt)
  }
}

# --------------------------------------------------------------------------
# 4. Write all K values back to globals.yml in one atomic write
# --------------------------------------------------------------------------

writeLines(globals_text, "config/globals.yml")
cat("\n--- K values written to config/globals.yml ---\n")
for (nm in names(k_results)) cat(sprintf("  %-35s = %d\n", nm, k_results[[nm]]))

cat("\n============================================================\n")
cat(" K-CV complete. Run run_merged_benchmark.R to fit all models.\n")
cat("============================================================\n")
