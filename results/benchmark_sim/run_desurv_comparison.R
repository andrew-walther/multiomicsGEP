# ============================================================
# Script:  results/benchmark_sim/run_desurv_comparison.R
# Purpose: Compare original M4/M5 against DeSurv-aligned preprocessing.
#
#   D1: LB  + per-platform z-std + variance + top-2000 + post-norm selection (= M4)
#   D2: YFB + per-platform z-std + variance + top-2000 + post-norm selection (= M5)
#   D3: LB  + per-platform z-std + combined_rank + top-3000 + per-cohort selection
#   D4: YFB + per-platform z-std + combined_rank + top-3000 + per-cohort selection
#   D5: YFB + per-platform z-std + combined_rank + top-3000 + per-cohort + cohort_id
#
#   D1/D2 reproduce the recommended manuscript configs using the same K already
#   stored in globals.yml.  D3/D4/D5 run K-CV (K_final = max(K_1se, 3)) and store
#   K values in globals.yml before fitting.
#
#   Output: results/benchmark_sim/outputs/desurv_comparison/
#     desurv_comparison_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-27
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_desurv_comparison.R
#          caffeinate -i Rscript results/benchmark_sim/run_desurv_comparison.R --quick
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

# Navigate to project root if invoked from a subdirectory
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
source("code/select_K.R")
source("code/select_alpha_cv.R")

YML_PATH <- "config/globals.yml"
cfg      <- yaml::read_yaml(YML_PATH)
b        <- cfg$benchmark
p        <- cfg$preprocessing

ALPHA            <- b$alpha
MAX_ITER         <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
K_GRID           <- if (QUICK_MODE) 2L:4L else 2L:10L
PRIOR_BETA       <- "normal"
SIGMA_COH        <- 1.0
BETA_THRESH      <- cfg$k_selection$beta_threshold
N_CV_FOLDS       <- cfg$cavi$n_cv_folds
K_MIN_BIOLOGICAL <- 3L
TOP_N_ORIG       <- p$top_n_genes        # 2000 — original setting
TOP_N_DESURV     <- p$top_n_genes_desurv # 3000 — DeSurv-aligned

OUT_DIR <- "results/benchmark_sim/outputs/desurv_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | K_grid: %d:%d\n",
            QUICK_MODE, MAX_ITER, min(K_GRID), max(K_GRID)))

# --------------------------------------------------------------------------
# 1. Load training data
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Configuration table
#    D1/D2 reproduce M4/M5 (for direct comparison, same K from globals.yml).
#    D3/D4 use DeSurv-aligned preprocessing; K-CV runs in section 4.
# --------------------------------------------------------------------------

DESURV_CONFIGS <- list(
  # D1/D2 replicate M4/M5 exactly: per-platform z-std + quantile normalization.
  list(id = "D1", label = "LB orig (M4)",
       model        = "LB",
       top_n        = TOP_N_ORIG,
       sel_method   = "variance",
       per_cohort   = FALSE,
       norm_method  = "quantile",
       cohort_id    = TRUE,
       k_key        = "k_merged_lb_perplatform"),
  list(id = "D2", label = "YFB orig (M5)",
       model        = "YFB",
       top_n        = TOP_N_ORIG,
       sel_method   = "variance",
       per_cohort   = FALSE,
       norm_method  = "quantile",
       cohort_id    = FALSE,
       k_key        = "k_merged_yfb_perplatform"),
  # D3/D4 use per-platform z-std only (no quantile QN on top) — per-cohort
  # selection runs before normalization, so variance signal is preserved.
  list(id = "D3", label = "LB DeSurv-aligned",
       model        = "LB",
       top_n        = TOP_N_DESURV,
       sel_method   = "combined_rank",
       per_cohort   = TRUE,
       norm_method  = "none",
       cohort_id    = TRUE,
       k_key        = "k_merged_lb_desurv"),
  list(id = "D4", label = "YFB DeSurv-aligned",
       model        = "YFB",
       top_n        = TOP_N_DESURV,
       sel_method   = "combined_rank",
       per_cohort   = TRUE,
       norm_method  = "none",
       cohort_id    = FALSE,
       k_key        = "k_merged_yfb_desurv"),
  list(id = "D5", label = "YFB DeSurv-aligned + cohort",
       model        = "YFB",
       top_n        = TOP_N_DESURV,
       sel_method   = "combined_rank",
       per_cohort   = TRUE,
       norm_method  = "none",
       cohort_id    = TRUE,
       k_key        = "k_merged_yfb_desurv_cohort")
)

# --------------------------------------------------------------------------
# Helper: replace any value for a named key in globals.yml.
#
# Reads the file as text lines, finds the line matching the key, and
# substitutes the value token. Preserves all comments and formatting.
# Pattern matches the exact key name at the start of the value token (after
# optional whitespace and colon), so partial key names do not match.
# --------------------------------------------------------------------------

set_key <- function(yml_path, key, value) {
  lines   <- readLines(yml_path)
  pattern <- paste0("^(\\s*", key, ":\\s*)\\S+")
  idx     <- grep(pattern, lines)
  if (length(idx) == 0)
    stop(sprintf("Key '%s' not found in %s", key, yml_path))
  lines[idx[1]] <- sub(pattern, paste0("\\1", value), lines[idx[1]])
  writeLines(lines, yml_path)
  invisible(NULL)
}

# --------------------------------------------------------------------------
# Helper: oriented C-index
# --------------------------------------------------------------------------

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 3. Preprocess training data — one call per unique (top_n, method, per_cohort)
# --------------------------------------------------------------------------

cat("--- Preprocessing training data ---\n")
preproc_cache  <- list()
gene_set_cache <- list()

for (dcfg in DESURV_CONFIGS) {
  ckey <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, dcfg$norm_method, sep = "_")
  if (ckey %in% names(preproc_cache)) next
  cat(sprintf("  top_n=%d, method=%s, per_cohort=%s, norm=%s ...\n",
              dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, dcfg$norm_method))
  pp <- preprocess_merged_cohorts(
    cohort_raw_list          = train_raw,
    log_transform_flags      = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
    top_n                    = dcfg$top_n,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE,
    normalize_method         = dcfg$norm_method,
    selection_per_cohort     = dcfg$per_cohort,
    selection_method         = dcfg$sel_method
  )
  preproc_cache[[ckey]]  <- pp$Y
  gene_set_cache[[ckey]] <- pp$gene_names
  cat(sprintf("    n=%d, p=%d\n", nrow(pp$Y), ncol(pp$Y)))
}

# --------------------------------------------------------------------------
# 4. K-CV for D3 and D4 (D1/D2 reuse K from existing globals.yml entries)
# --------------------------------------------------------------------------

cat("\n--- K-CV for DeSurv-aligned configs (D3, D4, D5) ---\n")
cfg <- yaml::read_yaml(YML_PATH); b <- cfg$benchmark  # re-read after preprocessing

for (dcfg in DESURV_CONFIGS[3:5]) {
  if (!is.null(b[[dcfg$k_key]])) {
    cat(sprintf("  %s: K=%d already set — skipping CV\n", dcfg$id, b[[dcfg$k_key]]))
    next
  }
  cat(sprintf("  Running K-CV for %s (%s) ...\n", dcfg$id, dcfg$label))
  ckey    <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, dcfg$norm_method, sep = "_")
  Y_train <- preproc_cache[[ckey]]
  cid_cv  <- if (dcfg$cohort_id) cohort_labels else NULL

  set.seed(42L)
  cv_res  <- select_K_cv(
    Y          = Y_train,
    time       = time_train,
    status     = status_train,
    K_grid     = K_GRID,
    n_folds    = N_CV_FOLDS,
    seed       = 42L,
    verbose    = TRUE,
    model      = dcfg$model,
    max_iter   = MAX_ITER,
    prior_beta = PRIOR_BETA,
    alpha      = ALPHA,
    cohort_id  = cid_cv,
    sign_correction = FALSE
  )
  K_opt   <- cv_res$K_opt
  K_final <- max(K_opt, K_MIN_BIOLOGICAL)
  cat(sprintf("  %s: K_opt=%d -> K_final=%d (biological floor K>=%d)\n",
              dcfg$id, K_opt, K_final, K_MIN_BIOLOGICAL))
  set_key(YML_PATH, dcfg$k_key, as.character(K_final))
}

cfg <- yaml::read_yaml(YML_PATH); b <- cfg$benchmark  # re-read after K-CV writes

# --------------------------------------------------------------------------
# 5. Fit all 4 configurations
# --------------------------------------------------------------------------

cat("\n=== Fitting 5 configurations ===\n\n")
fits <- list()

for (dcfg in DESURV_CONFIGS) {
  ckey      <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, dcfg$norm_method, sep = "_")
  Y_train   <- preproc_cache[[ckey]]
  K         <- b[[dcfg$k_key]]
  cohort_id <- if (dcfg$cohort_id) cohort_labels else NULL

  cat(sprintf("--- %s [%s] K=%d cohort_id=%s ---\n",
              dcfg$id, dcfg$label, K, !is.null(cohort_id)))
  set.seed(42L)

  fit <- suppressMessages(
    if (dcfg$model == "LB")
      fit_supervised_mf_modular(
        Y_train, time_train, status_train,
        K = K, max_iter = MAX_ITER, alpha = ALPHA,
        prior_beta = PRIOR_BETA, verbose = TRUE,
        cohort_id = cohort_id, sigma_F_cohort = SIGMA_COH)
    else
      fit_cox_on_yf(
        Y_train, time_train, status_train,
        K = K, max_iter = MAX_ITER, alpha = ALPHA,
        prior_beta = PRIOR_BETA, verbose = TRUE,
        cohort_id = cohort_id, sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fits[[dcfg$id]] <- fit
}

# --------------------------------------------------------------------------
# 6. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows     <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  # top_n=NULL: keep all external genes; intersection with train_genes controls
  # the final gene set. This prevents the external-cohort top-N filter from
  # discarding genes that happen to be in the training gene set.
  #
  # rank_transform=FALSE, per_platform_standardize=TRUE: matches the training
  # preprocessing exactly (Section 3 above calls preprocess_merged_cohorts with
  # the same two settings for every DESURV_CONFIGS entry, including D1/D2).
  # Before Phase 1c this was the reverse -- external cohorts were rank-
  # transformed and never per-platform z-standardized, while training was
  # per-platform z-standardized and never rank-transformed. See DECISIONS.md.
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE
  )

  for (dcfg in DESURV_CONFIGS) {
    ckey        <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, dcfg$norm_method, sep = "_")
    train_genes <- gene_set_cache[[ckey]]
    fit         <- fits[[dcfg$id]]
    K           <- b[[dcfg$k_key]]

    common    <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) {
      cat(sprintf("    Skipping %s x %s: only %d common genes\n",
                  ext_cohort, dcfg$id, length(common)))
      next
    }

    Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    train_idx <- match(common, train_genes)
    EF_sub    <- fit$EF[train_idx, , drop = FALSE]

    pred <- if (dcfg$model == "LB")
      predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)
    else
      predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)

    c_val <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      model          = dcfg$id,
      label          = dcfg$label,
      model_type     = dcfg$model,
      cohort         = ext_cohort,
      c_index        = round(c_val, 4),
      K              = K,
      k_eff          = sum(abs(fit$EBeta) > BETA_THRESH),
      beta_max       = round(max(abs(fit$EBeta)), 4),
      top_n          = dcfg$top_n,
      sel_method     = dcfg$sel_method,
      per_cohort     = dcfg$per_cohort,
      n_common_genes = length(common),
      stringsAsFactors = FALSE
    )
  }
}

# --------------------------------------------------------------------------
# 7. Save and report
# --------------------------------------------------------------------------

if (length(results_rows) == 0)
  stop("No external validation rows: all cohort x config pairs had < 100 common genes.")

results <- do.call(rbind, results_rows)
out_csv <- file.path(OUT_DIR, "desurv_comparison_results.csv")
write.csv(results, out_csv, row.names = FALSE)

# Save fit objects for post-hoc analysis (factor loadings, pathway enrichment).
saveRDS(fits, file.path(OUT_DIR, "desurv_comparison_fits.rds"))

cat(sprintf("\n=== Results saved: %s ===\n\n", out_csv))
cat("Mean C-index by configuration:\n")
agg <- aggregate(c_index ~ model + label + K, data = results, FUN = mean)
agg$k_eff <- sapply(agg$model, function(m) unique(results$k_eff[results$model == m]))
agg <- agg[order(agg$c_index, decreasing = TRUE), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  %s (%s): mean C=%.3f | K=%d | K_eff=%d\n",
              agg$model[i], agg$label[i], agg$c_index[i], agg$K[i], agg$k_eff[i]))
}
cat("\nPer-cohort C-index:\n")
for (m in unique(results$model)) {
  sub <- results[results$model == m, ]
  cat(sprintf("  %s: %s\n", m,
              paste(sprintf("%s=%.3f", sub$cohort, sub$c_index), collapse = "  ")))
}
