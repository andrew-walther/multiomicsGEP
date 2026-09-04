# ============================================================
# Script:  results/benchmark_sim/run_cohort_beta_supplementary.R
# Purpose: Two Stage-3 items missed by run_cohort_beta_comparison.R:
#            1. The strata_only decomposition sub-arm (cohort_id and
#               beta_cohort_id already have single-component arms --
#               joint_yfb_cohort_L and joint_yfb_beta_c -- but strata_id
#               alone was never fit).
#            2. Stage 1's held-out survival log-likelihood for every
#               cohort-aware arm, alongside external C-index -- the plan
#               asks for both, only external C was reported originally.
#
#          No re-fitting of the four existing arms' MAIN fits (reused from
#          cohort_beta_comparison_fits.rds); strata_only is fit fresh here
#          since it doesn't exist yet. cv_survival_loglik() itself always
#          refits per fold regardless (that's what "held-out" means), so
#          this script's own cost is 5 arms x 5 folds = 25 fresh fold fits.
#
#   Output: results/benchmark_sim/outputs/cohort_beta_comparison/
#             cohort_beta_supplementary_results.csv  (adds strata_only's
#               external C-index to the existing 5-arm table)
#             cohort_beta_heldout_survival_ll.csv    (held-out survival
#               log-likelihood for all 5 arms)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/benchmark_sim/run_cohort_beta_supplementary.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_beta_cohort.R")
source("code/update_L.R"); source("code/update_F.R"); source("code/update_tau.R")
source("code/compute_elbo.R"); source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/compute_cv_loglik.R")
source("code/preprocess_desurv.R")
source("code/select_K.R")
source("code/concordance_ci.R")  # frozen_reverse_cindex()

OUT_DIR <- "results/benchmark_sim/outputs/cohort_beta_comparison"
CVL <- cfg$cv_loglik
b   <- cfg$benchmark
K_REC <- 7L

# --------------------------------------------------------------------------
# Rebuild D4 training data + external cohorts (identical to
# run_cohort_beta_comparison.R).
# --------------------------------------------------------------------------
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
n_tcga <- train_raw$TCGA_PAAD$n; n_cptac <- train_raw$CPTAC$n
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
dataset_labels <- c(rep(TRAIN_COHORTS[1], n_tcga), rep(TRAIN_COHORTS[2], n_cptac))
pp <- preprocess_merged_cohorts(train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
        top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
        per_platform_standardize = TRUE, normalize_method = "none",
        selection_per_cohort = TRUE, selection_method = "combined_rank")
Y_train <- pp$Y; train_genes <- pp$gene_names

EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
ext_data <- list()
for (ec in EXTERNAL_COHORTS) {
  raw_ext <- load_pdac_raw(ec, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(raw_ext$Y, raw_ext$gene_names, top_n = NULL,
               log_transform = PLATFORM_LOG_TRANSFORM[[ec]], cohort_name = ec,
               rank_transform = FALSE, per_platform_standardize = TRUE)
  common <- intersect(train_genes, pre_ext$gene_names)
  ext_data[[ec]] <- list(Y_ext = pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE],
                          train_idx = match(common, train_genes), time = raw_ext$time, status = raw_ext$status)
}
# frozen_cindex(): orientation frozen at fit time by fit_cox_on_yf()'s Phase C
# (fixed 2026-09-04, DECISIONS.md) -- fit$EBeta/EBeta_pooled already carry the
# correct sign, so no per-cohort max(c,1-c) re-orientation here.
frozen_cindex <- function(risk, time, status) frozen_reverse_cindex(risk, time, status)
score_external <- function(fit, use_pooled_beta) {
  cs <- sapply(names(ext_data), function(ec) {
    d <- ext_data[[ec]]; EF_sub <- fit$EF[d$train_idx, , drop = FALSE]
    beta <- if (use_pooled_beta) fit$EBeta_pooled else fit$EBeta
    pred <- predict_cox_on_yf(d$Y_ext, EF_sub, beta, EF_norms = fit$EF_norms)
    frozen_cindex(pred$risk_scores, d$time, d$status)
  })
  mean(cs, na.rm = TRUE)
}

# --------------------------------------------------------------------------
# 1. strata_only sub-arm (strata_id alone -- fit fresh, not in the cache).
# --------------------------------------------------------------------------
cat("=== Fitting strata_only sub-arm ===\n")
set.seed(42L)
fit_strata_only <- fit_cox_on_yf(Y_train, time_train, status_train, K = K_REC,
                                  max_iter = cfg$cavi$max_iter, alpha = b$alpha, prior_beta = "normal",
                                  verbose = FALSE, strata_id = dataset_labels)
cls_strata <- classify_factors(fit_strata_only, Y_train,
                                beta_thresh = cfg$k_selection$beta_threshold,
                                pve_thresh  = cfg$k_selection$pve_threshold)
mean_c_strata <- score_external(fit_strata_only, use_pooled_beta = FALSE)
cat(sprintf("  strata_only: K_survival_active=%d, mean external C=%.4f\n",
            sum(cls_strata$category == "survival_active"), mean_c_strata))

# --------------------------------------------------------------------------
# 2. Held-out survival log-likelihood for all 5 arms (Stage 1's criterion,
#    requested alongside external C in the Stage 3 evaluation spec).
# --------------------------------------------------------------------------
cat("\n=== Held-out survival log-likelihood, all 5 arms ===\n")
arm_specs <- list(
  joint_yfb           = list(cohort_id = NULL,           strata_id = NULL,           beta_cohort_id = NULL),
  joint_yfb_cohort_L  = list(cohort_id = dataset_labels,  strata_id = NULL,           beta_cohort_id = NULL),
  joint_yfb_beta_c    = list(cohort_id = NULL,            strata_id = NULL,           beta_cohort_id = dataset_labels),
  joint_yfb_all_c     = list(cohort_id = dataset_labels,  strata_id = dataset_labels, beta_cohort_id = dataset_labels),
  strata_only         = list(cohort_id = NULL,            strata_id = dataset_labels, beta_cohort_id = NULL)
)

ll_rows <- list()
for (arm in names(arm_specs)) {
  spec <- arm_specs[[arm]]
  cat(sprintf("  %s ...\n", arm))
  res <- tryCatch(
    # cv_scoring="within_cohort" (explicit, though also the default as of
    # 2026-09-04): this is an ordinary K-fold CV over patients from the same
    # fixed set of training cohorts, so each held-out patient's own cohort
    # label is known -- not a genuinely unseen cohort. See
    # code/compute_cv_loglik.R's cv_scoring roxygen for the full distinction.
    cv_survival_loglik(Y_train, time_train, status_train, K = K_REC,
                        n_folds = CVL$n_folds, seed = CVL$seed,
                        max_iter = cfg$cavi$max_iter, alpha = b$alpha, prior_beta = "normal",
                        cohort_id = spec$cohort_id, strata_id = spec$strata_id,
                        beta_cohort_id = spec$beta_cohort_id, cv_scoring = "within_cohort"),
    error = function(e) { cat(sprintf("    FAILED: %s\n", conditionMessage(e))); NULL }
  )
  ll_rows[[arm]] <- data.frame(
    arm = arm,
    total_logPL = if (is.null(res)) NA_real_ else round(res$total_logPL, 4),
    mean_logPL_per_event = if (is.null(res)) NA_real_ else round(res$mean_logPL_per_event, 6),
    se_logPL = if (is.null(res)) NA_real_ else round(res$se_logPL, 4),
    stringsAsFactors = FALSE
  )
  if (!is.null(res)) cat(sprintf("    total_logPL=%.2f, per-event=%.4f\n", res$total_logPL, res$mean_logPL_per_event))
}
ll_table <- do.call(rbind, ll_rows)
write.csv(ll_table, file.path(OUT_DIR, "cohort_beta_heldout_survival_ll.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# 3. Add strata_only's row to the existing 5-arm comparison table.
# --------------------------------------------------------------------------
existing <- read.csv(file.path(OUT_DIR, "cohort_beta_comparison_results.csv"), stringsAsFactors = FALSE)
strata_row <- data.frame(
  arm = "strata_only", K = K_REC,
  K_survival_active = sum(cls_strata$category == "survival_active"),
  K_genomics_only   = sum(cls_strata$category == "genomics_only"),
  K_eff_total       = sum(cls_strata$category %in% c("survival_active", "genomics_only")),
  mean_external_c   = round(mean_c_strata, 4), fit_secs = NA_real_,
  stringsAsFactors = FALSE
)
for (ec in EXTERNAL_COHORTS) {
  d <- ext_data[[ec]]; EF_sub <- fit_strata_only$EF[d$train_idx, , drop = FALSE]
  pred <- predict_cox_on_yf(d$Y_ext, EF_sub, fit_strata_only$EBeta, EF_norms = fit_strata_only$EF_norms)
  strata_row[[paste0("c_", ec)]] <- round(frozen_cindex(pred$risk_scores, d$time, d$status), 4)
}
common_cols <- intersect(names(existing), names(strata_row))
updated <- rbind(existing[, common_cols], strata_row[, common_cols])
write.csv(updated, file.path(OUT_DIR, "cohort_beta_supplementary_results.csv"), row.names = FALSE)

cat("\n=== Updated arm table (with strata_only) ===\n")
print(updated[, c("arm", "K_survival_active", "K_eff_total", "mean_external_c")])
cat("\n=== Held-out survival log-likelihood by arm ===\n")
print(ll_table)
cat(sprintf("\nOutputs: %s, %s\n",
            file.path(OUT_DIR, "cohort_beta_supplementary_results.csv"),
            file.path(OUT_DIR, "cohort_beta_heldout_survival_ll.csv")))
