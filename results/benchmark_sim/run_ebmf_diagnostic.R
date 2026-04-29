# ============================================================
# Script: run_ebmf_diagnostic.R
# Purpose: Unsupervised factor analysis diagnostic on the merged TCGA_PAAD +
#          CPTAC v2-preprocessed training data.
#
#          This is a "data vs. model" diagnostic. SSBMF uses a joint genomics +
#          survival objective, so a zero-β result could mean either:
#            (a) No survival signal exists in the merged data  (data problem)
#            (b) Survival signal exists but the model can't find it  (model problem)
#
#          This script answers the question by running two unsupervised methods
#          (EBMF via flashier, PCA) on the same merged training data with NO
#          survival objective. For each factor, a univariate Cox model is fitted
#          to test for association with overall survival.
#
#          Decision logic (post-run):
#            EBMF factors associate with survival → data has signal; SSBMF has
#              a model/objective-weighting problem → investigate lambda / prior
#            EBMF factors split cleanly by cohort, no survival association →
#              batch dominates all variance → preprocessing fix or modality-
#              specific training is required
#            EBMF factors neither batch-split nor survival-associated →
#              the merged gene set may lack prognostic information entirely
#
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-29
# Dependencies: flashier (>= 1.0.0), survival, pheatmap
# Inputs:  results/benchmark_sim/outputs/real_data/merged/v2_point_normal/
#            tables/final_model.rds   (reuses saved training data metadata)
#          code/preprocess_desurv.R  (re-runs preprocess_merged_cohorts if needed)
# Outputs: results/benchmark_sim/outputs/ebmf_diagnostic/
#            tables/ebmf_cox_summary.csv
#            tables/pca_cox_summary.csv
#            figures/ebmf_loading_heatmap.{pdf,png}
#            figures/pca_loading_heatmap.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(flashier)
  library(survival)
  library(pheatmap)
})

# ── resolve repo root ──────────────────────────────────────────────────────────
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R")) {
  if (file.exists("../../code/fit_modular.R")) setwd("../..")
}

# ── source preprocessing helpers ───────────────────────────────────────────────
source("code/preprocess_desurv.R")
source("results/benchmark_sim/run_phase1_diagnostics.R")  # for plot_cohort_loading_heatmap

# ── helper: load PDAC raw data (mirrors run_ssbmf_benchmark.R) ────────────────
PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT",
                             unset = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/PDAC_data")

if (!exists("load_pdac_raw")) {
  # Source the benchmark runner in a minimal way just to get load_pdac_raw().
  # The if (sys.nframe() == 0) guard keeps the entry-point from running.
  source("results/benchmark_sim/run_ssbmf_benchmark.R")
}

PLATFORM_LOG_TRANSFORM <- c(TCGA_PAAD = TRUE, CPTAC = FALSE)

# ── output directories ─────────────────────────────────────────────────────────
OUT_ROOT  <- "results/benchmark_sim/outputs/ebmf_diagnostic"
TABLE_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR   <- file.path(OUT_ROOT, "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

# ── parameters ─────────────────────────────────────────────────────────────────
K_EBMF     <- 20    # request up to 20 factors; flashier may converge on fewer
TOP_N      <- 2000  # must match the v2 SSBMF run for a fair comparison
COX_ALPHA  <- 0.05  # significance threshold for flagging survival association

# ==============================================================================
# 1. Load and preprocess merged training data (v2 pipeline)
# ==============================================================================
cat("=== EBMF Diagnostic: loading and preprocessing merged training data ===\n")

if (!dir.exists(PDAC_DATA_ROOT)) {
  stop(sprintf(
    "PDAC data not found at '%s'. Set PDAC_DATA_ROOT env var to the data directory.",
    PDAC_DATA_ROOT))
}

train_cohorts <- c("TCGA_PAAD", "CPTAC")
train_raw <- lapply(setNames(train_cohorts, train_cohorts), function(ds) {
  cat(sprintf("  Loading %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})

cat("  Running v2 preprocessing (intersect → log2 → QN → top-2000 → rank) ...\n")
merged <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM,
  top_n               = TOP_N,
  rank_transform      = TRUE
)

Y             <- merged$Y             # n × p rank-transformed matrix
cohort_labels <- merged$dataset_labels
time_all      <- unlist(lapply(train_cohorts, function(ds) train_raw[[ds]]$time))
status_all    <- unlist(lapply(train_cohorts, function(ds) train_raw[[ds]]$status))

cat(sprintf("  Training data: n=%d, p=%d, event_rate=%.1f%%\n",
            nrow(Y), ncol(Y), 100 * mean(status_all)))
cat(sprintf("  Cohort counts: %s\n",
            paste(sprintf("%s=%d", names(table(cohort_labels)),
                          as.integer(table(cohort_labels))), collapse = ", ")))

# ==============================================================================
# 2. EBMF via flashier
# ==============================================================================
cat(sprintf("\n=== Fitting EBMF (flashier, K_max=%d, var_type=2) ===\n", K_EBMF))

# var_type = 2 (per-column / per-gene residual variance) is appropriate for
# genomics: each gene has its own noise level, so a single global σ² would
# under-shrink low-noise genes and over-shrink high-noise ones.
#
# greedy + backfit: greedy adds factors one at a time, backfit refines jointly.
# K_max controls the upper bound; flashier stops early if additional factors
# explain negligible variance.
flash_fit <- flash(
  Y,
  Kmax       = K_EBMF,
  var_type   = 2,
  greedy_Kmax = K_EBMF,
  backfit    = TRUE,
  verbose    = 1
)

K_ebmf <- flash_fit$n_factors
cat(sprintf("  EBMF converged with K=%d active factors.\n", K_ebmf))

# Extract L (n × K) loadings via LDF decomposition (type "2" rescales so
# columns of L have unit L2 norm; D absorbs the scale).  This puts all factors
# on a comparable scale for the Cox analysis.
ldf_res  <- ldf(flash_fit, type = "2")
L_ebmf   <- ldf_res$L   # n × K, unit-norm columns
D_ebmf   <- ldf_res$D   # length-K scale vector

saveRDS(list(flash_fit = flash_fit, L = L_ebmf, D = D_ebmf,
             cohort_labels = cohort_labels, time = time_all, status = status_all),
        file.path(TABLE_DIR, "ebmf_fit.rds"))

# ==============================================================================
# 3. Univariate Cox per EBMF factor
# ==============================================================================
cat("\n=== Univariate Cox per EBMF factor ===\n")

#' Fit a univariate Cox model for one factor column and extract key scalars.
#'
#' @param scores  numeric vector of length n — subject loadings for factor k
#' @param time    numeric vector of length n — follow-up times
#' @param status  integer vector of length n — event indicators (1 = event)
#' @param k       factor index (for labelling)
#' @return data.frame row: Factor, coef, HR, z_stat, p_value, c_index
fit_univariate_cox <- function(scores, time, status, k) {
  fit    <- tryCatch(
    coxph(Surv(time, status) ~ scores),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(Factor = paste0("EBMF", k), coef = NA_real_,
                      HR = NA_real_, z_stat = NA_real_,
                      p_value = NA_real_, c_index = NA_real_,
                      stringsAsFactors = FALSE))
  }
  s <- summary(fit)
  data.frame(
    Factor   = paste0("EBMF", k),
    coef     = round(coef(fit)[["scores"]], 4),
    HR       = round(s$conf.int[1, "exp(coef)"], 4),
    z_stat   = round(s$coefficients[1, "z"], 4),
    p_value  = round(s$coefficients[1, "Pr(>|z|)"], 6),
    c_index  = round(s$concordance[["C"]], 4),
    stringsAsFactors = FALSE
  )
}

ebmf_cox_rows <- lapply(seq_len(K_ebmf), function(k) {
  fit_univariate_cox(L_ebmf[, k], time_all, status_all, k)
})
ebmf_cox_df <- do.call(rbind, ebmf_cox_rows)
ebmf_cox_df$Significant <- ebmf_cox_df$p_value < COX_ALPHA

write.csv(ebmf_cox_df, file.path(TABLE_DIR, "ebmf_cox_summary.csv"),
          row.names = FALSE)

cat(sprintf(
  "  EBMF factors with p < %.2f: %d / %d\n",
  COX_ALPHA,
  sum(ebmf_cox_df$Significant, na.rm = TRUE),
  K_ebmf
))
print(ebmf_cox_df)

# ==============================================================================
# 4. PCA baseline
# ==============================================================================
cat(sprintf("\n=== PCA baseline (K=%d PCs) ===\n", min(K_EBMF, nrow(Y) - 1)))

K_pca    <- min(K_EBMF, nrow(Y) - 1)
pca_fit  <- prcomp(Y, rank. = K_pca, center = TRUE, scale. = FALSE)
L_pca    <- pca_fit$x  # n × K_pca score matrix (unit-variance columns)

pca_cox_rows <- lapply(seq_len(K_pca), function(k) {
  fit_univariate_cox(L_pca[, k], time_all, status_all, k)
})
pca_cox_df <- do.call(rbind, pca_cox_rows)
pca_cox_df$Factor      <- paste0("PC", seq_len(K_pca))
pca_cox_df$Significant <- pca_cox_df$p_value < COX_ALPHA

write.csv(pca_cox_df, file.path(TABLE_DIR, "pca_cox_summary.csv"),
          row.names = FALSE)

cat(sprintf(
  "  PCA factors with p < %.2f: %d / %d\n",
  COX_ALPHA,
  sum(pca_cox_df$Significant, na.rm = TRUE),
  K_pca
))
print(pca_cox_df[, c("Factor", "coef", "HR", "z_stat", "p_value", "c_index", "Significant")])

# ==============================================================================
# 5. Cohort-stratified loading heatmaps
# ==============================================================================
cat("\n=== Generating loading heatmaps ===\n")

# Reuse plot_cohort_loading_heatmap from run_phase1_diagnostics.R.
# EBeta is not meaningful here (unsupervised); pass a zero vector so the
# Beta_Status annotation column shows "Shrunk" for all factors (gray) — this
# cleanly communicates that no survival regression was run.

# ── EBMF heatmap ──
plot_cohort_loading_heatmap(
  EL            = L_ebmf,
  EBeta         = rep(0, K_ebmf),
  cohort_labels = cohort_labels,
  out_stub      = file.path(FIG_DIR, "ebmf_loading_heatmap"),
  beta_thresh   = COX_ALPHA   # repurposed: column annotation shows "Active" if
                               # the factor is Cox-significant at COX_ALPHA
)

# Override the column annotation to show Cox significance instead of β status.
# pheatmap is already written above; produce a second version with p-value
# annotation via a named vector trick: we pass |−log10(p)| as the EBeta
# substitute so the Abs_Beta gradient strip reflects Cox significance.
neg_log_p_ebmf <- -log10(pmax(ebmf_cox_df$p_value, 1e-10))
plot_cohort_loading_heatmap(
  EL            = L_ebmf,
  EBeta         = ifelse(ebmf_cox_df$Significant, neg_log_p_ebmf, 0),
  cohort_labels = cohort_labels,
  out_stub      = file.path(FIG_DIR, "ebmf_loading_heatmap_coxannot"),
  beta_thresh   = -log10(COX_ALPHA)   # threshold = -log10(0.05) ≈ 1.3
)

# ── PCA heatmap (top K_pca PCs, same format) ──
plot_cohort_loading_heatmap(
  EL            = L_pca,
  EBeta         = ifelse(pca_cox_df$Significant,
                         -log10(pmax(pca_cox_df$p_value, 1e-10)), 0),
  cohort_labels = cohort_labels,
  out_stub      = file.path(FIG_DIR, "pca_loading_heatmap_coxannot"),
  beta_thresh   = -log10(COX_ALPHA)
)

cat(sprintf("  EBMF heatmap    : %s/ebmf_loading_heatmap.{pdf,png}\n",   FIG_DIR))
cat(sprintf("  EBMF Cox heatmap: %s/ebmf_loading_heatmap_coxannot.{pdf,png}\n", FIG_DIR))
cat(sprintf("  PCA Cox heatmap : %s/pca_loading_heatmap_coxannot.{pdf,png}\n",  FIG_DIR))

# ==============================================================================
# 6. Summary print
# ==============================================================================
cat("\n=== EBMF Diagnostic Summary ===\n")
cat(sprintf("  EBMF K         : %d factors\n",  K_ebmf))
cat(sprintf("  PCA K          : %d components\n", K_pca))
cat(sprintf("  EBMF Cox-sig   : %d / %d (p < %.2f)\n",
            sum(ebmf_cox_df$Significant, na.rm = TRUE), K_ebmf, COX_ALPHA))
cat(sprintf("  PCA Cox-sig    : %d / %d (p < %.2f)\n",
            sum(pca_cox_df$Significant,  na.rm = TRUE), K_pca,  COX_ALPHA))

cat("\n  Interpretation guide:\n")
cat("  ┌─ EBMF Cox-sig > 0  AND  factor loads uniformly on both cohorts\n")
cat("  │  → Merged data has survival signal; SSBMF is an objective-weighting issue\n")
cat("  │  → Next step: increase lambda, debug β prior\n")
cat("  ├─ EBMF Cox-sig > 0  BUT  factor splits by cohort\n")
cat("  │  → Batch signal correlates with survival in training data (confounding)\n")
cat("  │  → TCGA-only training may give cleaner β estimates\n")
cat("  └─ EBMF Cox-sig = 0  entirely\n")
cat("     → No detectable survival signal regardless of method\n")
cat("     → Consider single-cohort training or feature engineering\n")

cat("\nDone.\n")
