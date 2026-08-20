# ============================================================
# Script:  results/benchmark_sim/run_k_init_multistart_check.R
# Purpose: Follow-up to run_k_init_sweep.R (Analysis A). The single-init
#          sweep found K_init=5/6 to have the BEST training ELBO of the
#          whole grid, yet a qualitatively different factor split (3
#          survival-active factors instead of 2) and far worse external
#          C-index (~0.596-0.597 vs 0.626-0.628 for K>=7) -- a genuine
#          ELBO-vs-generalization disagreement. Before treating that as a
#          real finding, this script checks whether it survives multistart:
#          a single SVD init at any given K is not guaranteed to reach that
#          K's true ELBO optimum (CAVI is non-convex), so K=5/6's apparent
#          ELBO lead could simply be a lucky single-init artifact rather
#          than evidence that a smaller, differently-structured model is
#          genuinely preferred.
#
#          Method: fit_cox_on_yf_multistart() (already used for exactly
#          this purpose in run_k_parsimony_followup.R, 2026-07-13) --
#          n_init=15 restarts per K (1 SVD + 14 random inits), best-ELBO
#          selection -- at K in {5,6,7,8,9,10}. K=15/20 are excluded: they
#          already lose to K=7-10 on both ELBO and external C-index in the
#          single-init sweep, so multistart there answers no open question.
#
#          Same D4 preprocessing as run_k_init_sweep.R (per-platform z-std,
#          combined_rank gene selection, top-3000 per cohort, no cohort_id).
#
#   Output: results/benchmark_sim/outputs/k_init_sweep/
#             k_init_multistart_results.csv   (one row per K, best-of-multistart)
#             k_init_multistart_restarts.csv  (one row per K x restart, full diagnostics)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-19
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_k_init_multistart_check.R
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_multistart_check.R --quick
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
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/fit_modular_multistart.R")   # fit_cox_on_yf_multistart()
source("code/preprocess_desurv.R")
source("code/select_K.R")                 # classify_factors()

YML_PATH <- "config/globals.yml"
cfg      <- yaml::read_yaml(YML_PATH)
b        <- cfg$benchmark
p        <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
N_INIT       <- if (QUICK_MODE) 3L  else 15L
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
PVE_THRESH   <- cfg$k_selection$pve_threshold
TOP_N_DESURV <- p$top_n_genes_desurv

K_VALUES <- c(5L, 6L, 7L, 8L, 9L, 10L)

OUT_DIR <- "results/benchmark_sim/outputs/k_init_sweep"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | N_INIT: %d | K values: %s\n",
            QUICK_MODE, MAX_ITER, N_INIT, paste(K_VALUES, collapse = ", ")))

# --------------------------------------------------------------------------
# 1. Load + preprocess training data (identical to run_k_init_sweep.R)
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

cat("--- Preprocessing training data (D4 config) ---\n")
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

# --------------------------------------------------------------------------
# 2. Multistart fit at each K
# --------------------------------------------------------------------------

cat("=== Multistart fitting (best-ELBO) at each K ===\n\n")
best_fits    <- list()
restart_rows <- list()

for (K in K_VALUES) {
  cat(sprintf("--- K=%d (n_init=%d) ---\n", K, N_INIT))
  ms <- fit_cox_on_yf_multistart(
    Y_train, time_train, status_train,
    K = K, max_iter = MAX_ITER, alpha = ALPHA, prior_beta = PRIOR_BETA,
    n_init = N_INIT, init_seed_base = 42L
  )
  best_fits[[as.character(K)]] <- ms$best
  cat(sprintf("  best_idx=%d/%d | final_elbo=%.4f | |beta|: [%s]\n\n",
              ms$best_idx, N_INIT, tail(ms$best$history$elbo_full, 1),
              paste(sprintf("%.4f", abs(ms$best$EBeta)), collapse = ", ")))
  ms$restarts$K <- K
  restart_rows[[length(restart_rows) + 1]] <- ms$restarts
}

restarts_all <- do.call(rbind, restart_rows)
write.csv(restarts_all, file.path(OUT_DIR, "k_init_multistart_restarts.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# 3. External validation on 5 held-out cohorts, for the best-of-multistart
#    fit at each K
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts), best-of-multistart fits ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

ext_cindex <- vector("list", length(K_VALUES))
names(ext_cindex) <- as.character(K_VALUES)

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

  for (K in K_VALUES) {
    fit    <- best_fits[[as.character(K)]]
    EF_sub <- fit$EF[train_idx, , drop = FALSE]
    pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
    c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
    ext_cindex[[as.character(K)]][[ext_cohort]] <- c_val
  }
}

# --------------------------------------------------------------------------
# 4. Classify factors for the best-of-multistart fit at each K, assemble table
# --------------------------------------------------------------------------

cat("\n=== Best-of-multistart factor classification per K ===\n\n")
results_rows <- list()

for (K in K_VALUES) {
  fit    <- best_fits[[as.character(K)]]
  n_iter <- fit$history$n_iter
  cls    <- classify_factors(fit, Y_train, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)

  K_survival_active <- sum(cls$category == "survival_active")
  K_genomics_only   <- sum(cls$category == "genomics_only")
  K_dead            <- sum(cls$category == "dead")
  K_eff_total       <- K_survival_active + K_genomics_only

  cohort_c <- ext_cindex[[as.character(K)]]
  mean_c   <- if (length(cohort_c) > 0) mean(unlist(cohort_c)) else NA_real_

  cat(sprintf("K=%2d: K_survival_active=%d, K_genomics_only=%d, K_dead=%d, K_eff_total=%d | mean external C=%.4f | best_elbo=%.4f\n",
              K, K_survival_active, K_genomics_only, K_dead, K_eff_total, mean_c,
              tail(fit$history$elbo_full, 1)))

  results_rows[[length(results_rows) + 1]] <- data.frame(
    K                 = K,
    K_survival_active = K_survival_active,
    K_genomics_only   = K_genomics_only,
    K_dead            = K_dead,
    K_eff_total       = K_eff_total,
    best_elbo_full    = round(tail(fit$history$elbo_full, 1), 4),
    n_iter            = n_iter,
    mean_external_c   = round(mean_c, 4),
    stringsAsFactors  = FALSE
  )
}

results <- do.call(rbind, results_rows)

for (ext_cohort in EXTERNAL_COHORTS) {
  results[[paste0("c_", ext_cohort)]] <- sapply(results$K, function(K) {
    v <- ext_cindex[[as.character(K)]][[ext_cohort]]
    if (is.null(v)) NA_real_ else round(v, 4)
  })
}

out_csv <- file.path(OUT_DIR, "k_init_multistart_results.csv")
write.csv(results, out_csv, row.names = FALSE)
saveRDS(best_fits, file.path(OUT_DIR, "k_init_multistart_best_fits.rds"))

elbo_best_idx <- which.max(results$best_elbo_full)
cat(sprintf("\nBest-of-multistart ELBO-preferred K = %d (best_elbo_full=%.4f, K_eff_total=%d)\n",
            results$K[elbo_best_idx], results$best_elbo_full[elbo_best_idx],
            results$K_eff_total[elbo_best_idx]))
cat(sprintf("\n=== Results saved: %s ===\n", out_csv))
