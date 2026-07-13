# ============================================================
# Script:  results/benchmark_sim/run_k_parsimony_curve.R
# Purpose: Phase 3 of the post-lab-meeting action plan -- the internal K-CV
#          curve (ROADMAP.md / DECISIONS.md 2026-07-12: K=2..10 mean CV
#          C-index, peaking at K=8, K=7 selected via 1-SE) tells us nothing
#          about how much of the recommended config's EXTERNAL validation
#          number (D4: mean C=0.627 across 5 held-out PDAC cohorts, K=7)
#          actually requires K=7 vs. would hold at a smaller, more
#          parsimonious K.
#
#          This script refits the YFB D4 configuration (per-platform z-std,
#          DeSurv combined-rank gene selection, top-3000 per cohort, no
#          cohort_id) at K in {2, 3, 4, 5, 7} -- bypassing the single
#          CV-selected K stored in globals.yml -- and re-runs external
#          validation against the same 5 held-out cohorts for each K,
#          producing the actual K-vs-external-performance curve.
#
#          Decision rule (mirrors the CV 1-SE convention already used for
#          internal K selection): report the smallest K whose external
#          mean C-index (computed across the 5 cohorts, one value per
#          cohort) is within 1 SE of K=7's external mean C-index. This is
#          descriptive, not automatically applied -- the result is reported
#          for a human decision, not silently substituted into globals.yml.
#
#          Output: results/benchmark_sim/outputs/k_parsimony_curve/
#            k_parsimony_curve_results.csv   (one row per K x cohort)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_k_parsimony_curve.R
#          caffeinate -i Rscript results/benchmark_sim/run_k_parsimony_curve.R --quick
# Requires: PDAC_DATA_ROOT set (real data not in git) -- see CLAUDE.md /
#           tests/test_real_data_loading.R for the local/Longleaf paths.
#           Fails with an informative error via load_pdac_raw() if unset.
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
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

b <- cfg$benchmark
p <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
K_GRID       <- if (QUICK_MODE) c(2L, 7L) else c(2L, 3L, 4L, 5L, 7L)
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
TOP_N_DESURV <- p$top_n_genes_desurv   # 3000 -- DeSurv-aligned, matches D4

OUT_DIR <- "results/benchmark_sim/outputs/k_parsimony_curve"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | K grid: %s\n",
            QUICK_MODE, MAX_ITER, paste(K_GRID, collapse = ",")))

# --------------------------------------------------------------------------
# 1. Load and preprocess training data (D4 config: per-platform z-std +
#    combined_rank gene selection, per-cohort, top-3000, no cohort_id --
#    identical preprocessing call to D4 in run_desurv_comparison.R)
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
Y_train    <- pp$Y
train_genes <- pp$gene_names
cat(sprintf("  n=%d, p=%d\n\n", nrow(Y_train), ncol(Y_train)))

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. Fit YFB at each K in the grid
# --------------------------------------------------------------------------

cat("=== Fitting YFB D4 config across K grid ===\n\n")
fits <- list()
for (K in K_GRID) {
  cat(sprintf("--- K=%d ---\n", K))
  set.seed(42L)
  fit <- suppressMessages(fit_cox_on_yf(
    Y_train, time_train, status_train,
    K = K, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA, verbose = TRUE))
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fits[[as.character(K)]] <- fit
}

# --------------------------------------------------------------------------
# 3. External validation on the 5 held-out cohorts, for every K
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) x K grid ---\n")
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

  for (K in K_GRID) {
    fit    <- fits[[as.character(K)]]
    EF_sub <- fit$EF[train_idx, , drop = FALSE]
    pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
    c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      K              = K,
      cohort         = ext_cohort,
      c_index        = round(c_val, 4),
      k_eff          = sum(abs(fit$EBeta) > BETA_THRESH),
      beta_max       = round(max(abs(fit$EBeta)), 4),
      n_iter         = fit$history$n_iter,
      n_common_genes = length(common),
      stringsAsFactors = FALSE
    )
  }
}

if (length(results_rows) == 0)
  stop("No external validation rows: all cohorts had < 100 common genes.")

results <- do.call(rbind, results_rows)

# --------------------------------------------------------------------------
# 4. Save and report the K-vs-external-performance curve + 1-SE decision rule
# --------------------------------------------------------------------------

out_csv <- file.path(OUT_DIR, "k_parsimony_curve_results.csv")
write.csv(results, out_csv, row.names = FALSE)

agg <- aggregate(c_index ~ K, data = results,
                 FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
agg <- do.call(data.frame, agg)
names(agg)[2:3] <- c("mean_c", "se_c")
agg$k_eff <- sapply(agg$K, function(k) unique(results$k_eff[results$K == k]))
agg <- agg[order(agg$K), ]

cat("\n============================================================\n")
cat(" External validation: K vs. mean C-index (5 held-out cohorts)\n")
cat("============================================================\n")
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  K=%2d  mean C=%.4f (SE=%.4f)  K_eff=%d\n",
              agg$K[i], agg$mean_c[i], agg$se_c[i], agg$k_eff[i]))
}

# Reference is the empirical best K in this grid (mirrors select_K_cv()'s own
# which.max-based logic), not hardcoded to K=7 -- correct even if the grid or
# the winning K changes on a future rerun.
ref      <- agg[which.max(agg$mean_c), ]
margin   <- ref$mean_c - ref$se_c
eligible <- agg$K[agg$mean_c >= margin]
k_recommended <- min(eligible)
cat(sprintf("\n1-SE decision rule (reference K=%d: mean C=%.4f, SE=%.4f, margin=%.4f):\n",
            ref$K, ref$mean_c, ref$se_c, margin))
cat(sprintf("  Smallest K within 1 SE of best (K=%d): K=%d\n", ref$K, k_recommended))

cat(sprintf("\nResults: %s\n", out_csv))
cat("============================================================\n")
