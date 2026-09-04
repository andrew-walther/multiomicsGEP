# ============================================================
# Script:  results/benchmark_sim/run_cohort_beta_comparison.R
# Purpose: A focused preview of the 9/4 plan's Stage 3a performance
#          comparison, run immediately after Stage 2 (beta_cohort_id)
#          lands, to answer directly: does letting cohort membership reach
#          the SURVIVAL COEFFICIENTS add value over the current recommended
#          model, and how does it compare to the pre-existing genomics-only
#          cohort_id extension and to the two-step EBMF->Cox baseline?
#
#          This is NOT the full Stage 3 analysis (no bootstrap CIs,
#          no leave-one-study-out, no factor-comparison/overlap analysis --
#          those remain Stage 3 proper). It is a single external-C-index
#          table across five arms, all at the recommended K_init=7, D4
#          preprocessing (matches run_k_init_sweep.R and
#          run_desurv_comparison.R exactly):
#
#            1. joint_yfb           -- current recommended model (no cohort
#                                       information at all)
#            2. joint_yfb_cohort_L  -- + cohort_id (2026-05-22: genomics
#                                       offsets in L/F, beta_cohort=0 by
#                                       construction). Never benchmarked
#                                       under the CURRENT D4/K=7 config --
#                                       last benchmarked 2026-05-22 at K=20
#                                       under older preprocessing, where it
#                                       lost (mean external C 0.614 vs 0.627).
#            3. joint_yfb_beta_c    -- + beta_cohort_id (Stage 2, this
#                                       session): cohort-specific SURVIVAL
#                                       coefficients, the new contribution.
#            4. joint_yfb_all_c     -- cohort_id + strata_id + beta_cohort_id
#                                       together (all three cohort-aware
#                                       extensions at once).
#            5. two_step_ebmf_cox   -- reused from
#                                       ebmf_cox_regularized_results.csv
#                                       (K=40, LASSO stage 2) -- NOT refit,
#                                       so the number stays comparable to
#                                       the deck's own +0.026 headline.
#
#          Unsupervised EBMF/NMF alone is NOT included: it has no survival
#          coefficient, so there is no risk score to report a C-index for.
#
#          The grouping vector for cohort_id/strata_id/beta_cohort_id is the
#          same dataset_labels vector (TCGA_PAAD vs CPTAC) in every arm that
#          uses it -- C=2, so update_F_cohort_all() is exact, not the
#          Gauss-Seidel approximation it falls back to for C>2.
#
#   Sanity gate (must hold before interpreting anything else): arm 1's mean
#   external C and K_survival_active must match
#   desurv_comparison_results.csv (model D4) -- 0.627, K_eff=2.
#
#   Output: results/benchmark_sim/outputs/cohort_beta_comparison/
#             cohort_beta_comparison_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_cohort_beta_comparison.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_beta_cohort.R")
source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/compute_bic.R")
source("code/preprocess_desurv.R")
source("code/select_K.R")
source("code/concordance_ci.R")  # frozen_reverse_cindex()

b   <- cfg$benchmark
pcf <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- cfg$cavi$max_iter
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
PVE_THRESH   <- cfg$k_selection$pve_threshold
TOP_N_DESURV <- pcf$top_n_genes_desurv
K_REC        <- 7L

OUT_DIR <- "results/benchmark_sim/outputs/cohort_beta_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load + preprocess training data -- identical D4 config to
#    run_k_init_sweep.R / run_desurv_comparison.R.
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
dataset_labels <- c(rep(TRAIN_COHORTS[1], n_tcga), rep(TRAIN_COHORTS[2], n_cptac))
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
# 2. Load the 5 external cohorts once.
# --------------------------------------------------------------------------

cat("--- Loading external cohorts (5) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
ext_data <- list()
for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y = raw_ext$Y, gene_names = raw_ext$gene_names, top_n = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]], cohort_name = ext_cohort,
    rank_transform = FALSE, per_platform_standardize = TRUE
  )
  common <- intersect(train_genes, pre_ext$gene_names)
  if (length(common) < 100) { cat(sprintf("    Skipping %s\n", ext_cohort)); next }
  ext_data[[ext_cohort]] <- list(
    Y_ext = pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE],
    train_idx = match(common, train_genes),
    time = raw_ext$time, status = raw_ext$status
  )
}
cat(sprintf("  %d/%d external cohorts usable\n\n", length(ext_data), length(EXTERNAL_COHORTS)))

# frozen_cindex(): the orientation is decided ONCE, at fit time, by
# fit_cox_on_yf()'s Phase C (fixed 2026-09-04, DECISIONS.md -- correct-
# direction concordance, reverse=TRUE). fit$EBeta / fit$EBeta_pooled already
# carry that frozen sign, so external cohorts are scored as-is -- no
# per-cohort max(c, 1-c) re-orientation from the cohort's own outcomes.
frozen_cindex <- function(risk, time, status) frozen_reverse_cindex(risk, time, status)

score_external <- function(fit, use_pooled_beta) {
  cohort_c <- list()
  for (ext_cohort in names(ext_data)) {
    d      <- ext_data[[ext_cohort]]
    EF_sub <- fit$EF[d$train_idx, , drop = FALSE]
    beta_to_use <- if (use_pooled_beta) fit$EBeta_pooled else fit$EBeta
    pred   <- predict_cox_on_yf(d$Y_ext, EF_sub, beta_to_use, EF_norms = fit$EF_norms)
    cohort_c[[ext_cohort]] <- frozen_cindex(pred$risk_scores, d$time, d$status)
  }
  list(cohort_c = cohort_c, mean_c = mean(unlist(cohort_c), na.rm = TRUE))
}

# --------------------------------------------------------------------------
# 3. Fit the four joint-model arms at K_init=7.
# --------------------------------------------------------------------------

run_arm <- function(name, cohort_id = NULL, strata_id = NULL, beta_cohort_id = NULL) {
  cat(sprintf("=== Arm: %s ===\n", name))
  t0 <- proc.time()[["elapsed"]]
  set.seed(42L)
  fit <- fit_cox_on_yf(
    Y_train, time_train, status_train,
    K = K_REC, max_iter = MAX_ITER, alpha = ALPHA, prior_beta = PRIOR_BETA,
    verbose = FALSE, cohort_id = cohort_id, strata_id = strata_id,
    beta_cohort_id = beta_cohort_id
  )
  fit_secs <- proc.time()[["elapsed"]] - t0

  use_pooled <- !is.null(beta_cohort_id)
  sc <- score_external(fit, use_pooled_beta = use_pooled)

  cls <- classify_factors(fit, Y_train, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
  K_survival_active <- sum(cls$category == "survival_active")
  K_genomics_only   <- sum(cls$category == "genomics_only")
  K_eff_total       <- K_survival_active + K_genomics_only

  cat(sprintf("  K_survival_active=%d, K_genomics_only=%d, K_eff_total=%d | mean external C=%.4f | %.1fs\n\n",
              K_survival_active, K_genomics_only, K_eff_total, sc$mean_c, fit_secs))

  row <- data.frame(
    arm = name, K = K_REC,
    K_survival_active = K_survival_active, K_genomics_only = K_genomics_only,
    K_eff_total = K_eff_total, mean_external_c = round(sc$mean_c, 4),
    fit_secs = round(fit_secs, 1), stringsAsFactors = FALSE
  )
  for (ext_cohort in EXTERNAL_COHORTS) {
    v <- sc$cohort_c[[ext_cohort]]
    row[[paste0("c_", ext_cohort)]] <- if (is.null(v)) NA_real_ else round(v, 4)
  }
  list(row = row, fit = fit)
}

arms <- list()
fits <- list()

r <- run_arm("joint_yfb")
arms[["joint_yfb"]] <- r$row; fits[["joint_yfb"]] <- r$fit

r <- run_arm("joint_yfb_cohort_L", cohort_id = dataset_labels)
arms[["joint_yfb_cohort_L"]] <- r$row; fits[["joint_yfb_cohort_L"]] <- r$fit

r <- run_arm("joint_yfb_beta_c", beta_cohort_id = dataset_labels)
arms[["joint_yfb_beta_c"]] <- r$row; fits[["joint_yfb_beta_c"]] <- r$fit

r <- run_arm("joint_yfb_all_c", cohort_id = dataset_labels, strata_id = dataset_labels,
             beta_cohort_id = dataset_labels)
arms[["joint_yfb_all_c"]] <- r$row; fits[["joint_yfb_all_c"]] <- r$fit

# --------------------------------------------------------------------------
# 4. Two-step EBMF->Cox baseline -- REUSED, not refit (K=40, LASSO stage 2,
#    same baseline behind the deck's own +0.026 headline).
# --------------------------------------------------------------------------

two_step_csv <- "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_results.csv"
if (file.exists(two_step_csv)) {
  ts <- read.csv(two_step_csv, stringsAsFactors = FALSE)
  ts40 <- ts[ts$k_tag == "k40", ]
  row <- data.frame(
    arm = "two_step_ebmf_cox", K = 40L,
    K_survival_active = NA_integer_, K_genomics_only = NA_integer_,
    K_eff_total = ts40$k_eff_lasso[1], mean_external_c = round(mean(ts40$c_index), 4),
    fit_secs = NA_real_, stringsAsFactors = FALSE
  )
  for (ext_cohort in EXTERNAL_COHORTS) {
    v <- ts40$c_index[ts40$cohort == ext_cohort]
    row[[paste0("c_", ext_cohort)]] <- if (length(v) == 0) NA_real_ else round(v, 4)
  }
  arms[["two_step_ebmf_cox"]] <- row
  cat(sprintf("=== Arm: two_step_ebmf_cox (reused, not refit) === mean external C=%.4f\n\n",
              row$mean_external_c))
} else {
  cat(sprintf("WARNING: %s not found -- two_step_ebmf_cox arm omitted.\n", two_step_csv))
}

# --------------------------------------------------------------------------
# 5. Sanity gate + save.
# --------------------------------------------------------------------------

results <- do.call(rbind, arms)
rownames(results) <- NULL

desurv_csv <- "results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv"
if (file.exists(desurv_csv)) {
  d4 <- read.csv(desurv_csv, stringsAsFactors = FALSE)
  d4 <- d4[d4$model == "D4", ]
  d4_mean_c <- mean(d4$c_index)
  d4_k_eff  <- d4$k_eff[1]
  gate_c_ok <- abs(results$mean_external_c[results$arm == "joint_yfb"] - d4_mean_c) < 0.005
  gate_k_ok <- results$K_survival_active[results$arm == "joint_yfb"] == d4_k_eff
  cat(sprintf("=== Sanity gate: joint_yfb vs desurv_comparison_results.csv (D4) ===\n"))
  cat(sprintf("  joint_yfb mean C=%.4f vs D4=%.4f -- %s\n",
              results$mean_external_c[results$arm == "joint_yfb"], d4_mean_c,
              if (gate_c_ok) "PASS" else "FAIL"))
  cat(sprintf("  joint_yfb K_survival_active=%d vs D4 k_eff=%d -- %s\n",
              results$K_survival_active[results$arm == "joint_yfb"], d4_k_eff,
              if (gate_k_ok) "PASS" else "FAIL"))
  if (!gate_c_ok || !gate_k_ok) {
    stop("Sanity gate FAILED: joint_yfb has diverged from the D4 reference. ",
         "Fix before interpreting any other arm in this comparison.")
  }
}

cat("\n=== Cohort-beta comparison: mean external C by arm ===\n")
print(results[, c("arm", "K", "K_survival_active", "K_eff_total", "mean_external_c")])

out_csv <- file.path(OUT_DIR, "cohort_beta_comparison_results.csv")
write.csv(results, out_csv, row.names = FALSE)
saveRDS(fits, file.path(OUT_DIR, "cohort_beta_comparison_fits.rds"))
cat(sprintf("\n=== Results saved: %s ===\n", out_csv))
