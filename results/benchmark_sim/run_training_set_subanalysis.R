# ============================================================
# Script:  results/benchmark_sim/run_training_set_subanalysis.R
# Purpose: Stage 3 leftover -- "does pooling help at all?" The 8/21 notes
#          asserted that using more samples (pooling TCGA_PAAD+CPTAC,
#          n=273) improves factor estimation over training on either
#          cohort alone, but this was never actually tested. Refits the
#          plain joint model (arm 1, K_init=7) on tcga_only (n=144) and
#          cptac_only (n=129), each with its OWN gene selection (not the
#          pooled 2064-gene set) -- the natural consequence of fitting a
#          single cohort, not a choice made for convenience.
#
#          Kept deliberately separate from run_cohort_beta_comparison.R's
#          arm table: this answers "does pooling help," not "does a
#          cohort-aware extension help," and the two questions should not
#          be merged into one table. The pooled-gene-set fit (arm 1 from
#          the main comparison) is the primary result; this is a
#          not-gene-set-matched supplementary check.
#
#   Output: results/benchmark_sim/outputs/cohort_beta_comparison/
#             training_set_subanalysis_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/benchmark_sim/run_training_set_subanalysis.R
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
source("code/preprocess_desurv.R")
source("code/select_K.R")

K_REC <- 7L
b <- cfg$benchmark
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

oriented_cindex <- function(risk, time, status) {
  if (sd(risk) == 0) return(NA_real_)
  c <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance); max(c, 1 - c)
}

run_single_cohort <- function(ds) {
  cat(sprintf("=== Training on %s alone ===\n", ds))
  raw <- load_pdac_raw(ds, PDAC_DATA_ROOT)
  pre <- preprocess_desurv_cohort(
    Y = raw$Y, gene_names = raw$gene_names,
    top_n = cfg$preprocessing$top_n_genes_desurv,  # 3000 -- D4-consistent gene budget, own selection
    log_transform = PLATFORM_LOG_TRANSFORM[[ds]], cohort_name = ds,
    rank_transform = FALSE, per_platform_standardize = TRUE
  )
  Y_train <- pre$Y; train_genes <- pre$gene_names
  cat(sprintf("  n=%d, p=%d (own gene selection, not the pooled 2064)\n", nrow(Y_train), ncol(Y_train)))

  set.seed(42L)
  fit <- fit_cox_on_yf(Y_train, raw$time, raw$status, K = K_REC, max_iter = cfg$cavi$max_iter,
                        alpha = b$alpha, prior_beta = "normal", verbose = FALSE)
  cls <- classify_factors(fit, Y_train, beta_thresh = cfg$k_selection$beta_threshold,
                           pve_thresh = cfg$k_selection$pve_threshold)

  cohort_c <- list()
  for (ec in EXTERNAL_COHORTS) {
    raw_ext <- load_pdac_raw(ec, PDAC_DATA_ROOT)
    pre_ext <- preprocess_desurv_cohort(raw_ext$Y, raw_ext$gene_names, top_n = NULL,
                 log_transform = PLATFORM_LOG_TRANSFORM[[ec]], cohort_name = ec,
                 rank_transform = FALSE, per_platform_standardize = TRUE)
    common <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) { cat(sprintf("    Skipping %s: only %d common genes\n", ec, length(common))); next }
    Y_ext <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    EF_sub <- fit$EF[match(common, train_genes), , drop = FALSE]
    pred <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
    cohort_c[[ec]] <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
  }
  mean_c <- mean(unlist(cohort_c), na.rm = TRUE)
  K_surv <- sum(cls$category == "survival_active")
  cat(sprintf("  K_survival_active=%d, mean external C=%.4f\n\n", K_surv, mean_c))

  row <- data.frame(
    arm = paste0(tolower(ds), "_only"), n = nrow(Y_train), p = ncol(Y_train),
    K_survival_active = K_surv, mean_external_c = round(mean_c, 4), stringsAsFactors = FALSE
  )
  for (ec in EXTERNAL_COHORTS) {
    v <- cohort_c[[ec]]; row[[paste0("c_", ec)]] <- if (is.null(v)) NA_real_ else round(v, 4)
  }
  row
}

results <- rbind(run_single_cohort("TCGA_PAAD"), run_single_cohort("CPTAC"))

pooled <- read.csv("results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_results.csv")
pooled_row <- pooled[pooled$arm == "joint_yfb", ]
cat("=== Pooled (primary, gene-set-matched to itself, n=273, p=2064) vs. single-cohort (own gene selection) ===\n")
cat(sprintf("  pooled (TCGA+CPTAC): n=273, p=2064, K_survival_active=%d, mean external C=%.4f\n",
            pooled_row$K_survival_active, pooled_row$mean_external_c))
print(results[, c("arm", "n", "p", "K_survival_active", "mean_external_c")])

out_csv <- "results/benchmark_sim/outputs/cohort_beta_comparison/training_set_subanalysis_results.csv"
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nOutput: %s\n", out_csv))
cat("NOTE: this comparison is NOT gene-set-matched (each single-cohort fit selects its own genes).\n")
cat("The pooled 2064-gene fit is the primary result; this answers 'does pooling help', a separate question.\n")
