# ============================================================
# Script:       run_YFB_benchmark.R
# Purpose:      Cluster B (YFB) benchmark — side-by-side comparison of
#               prior_beta="point_normal" vs "normal" on synthetic validation,
#               PDAC training (merged TCGA_PAAD + CPTAC), and 5 external cohorts.
#               Replaces run_cox_on_yf_benchmark.R.
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-05-04
# Dependencies: code/fit_cox_on_yf.R, code/predict_cox_on_yf.R,
#               code/train_test_split.R, code/preprocess_desurv.R,
#               results/benchmark_sim/run_ssbmf_benchmark.R
# Usage:        Rscript results/benchmark_sim/run_YFB_benchmark.R [--quick]
# ============================================================

# --------------------------------------------------------------------------
# Setup — working directory and config
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/fit_cox_on_yf.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_cox_on_yf.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
  library(ebnm)
})

cfg <- yaml::read_yaml("config/globals.yml")

# Source data-loading hub (load_pdac_raw, preprocess_merged_cohorts, EXTERNAL_COHORTS,
# PLATFORM_LOG_TRANSFORM, PDAC_DATA_ROOT, generate_synthetic_benchmark_data).
# The hub also sources fit_modular.R via tryCatch — we don't need it here, but it's harmless.
suppressMessages(source("results/benchmark_sim/run_ssbmf_benchmark.R"))

# Source Cluster B fit and predict functions via tryCatch so their runner blocks
# (which require real_Y to be set) do not execute during source().
tryCatch(source("code/fit_cox_on_yf.R"),   error = function(e) invisible(NULL))
tryCatch(source("code/predict_cox_on_yf.R"), error = function(e) invisible(NULL))

source("code/train_test_split.R")
source("code/preprocess_desurv.R")

# Benchmark defaults from globals.yml
K            <- if (QUICK_MODE) 5 else cfg$benchmark$k_pdac
K_SYN        <- if (QUICK_MODE) 5 else cfg$benchmark$k_pdac_synthetic
ALPHA        <- cfg$benchmark$alpha
LAMBDA       <- cfg$benchmark$lambda
N_BURNIN     <- cfg$benchmark$n_burnin
NORM_AB      <- cfg$benchmark$normalize_ab
COX_WARMSTART <- cfg$benchmark$cox_warmstart
PRIORS       <- cfg$benchmark$prior_beta_compare
BETA_THRESH  <- cfg$k_selection$beta_threshold
MAX_ITER     <- if (QUICK_MODE) 30 else cfg$cavi$max_iter
ALPHA_F      <- 0   # Cluster B baseline: F update is pure-genomics (see DECISIONS.md 2026-04-30)

cat("=== YFB Benchmark (Cluster B — eta = YF·beta) ===\n")
cat(sprintf("    K=%d (synthetic K=%d) | alpha=%.2f | alpha_F=%.2f | lambda=%.2f | N_burnin=%d\n",
            K, K_SYN, ALPHA, ALPHA_F, LAMBDA, N_BURNIN))
cat(sprintf("    cox_warmstart=%s | normalize_AB=%s\n", COX_WARMSTART, NORM_AB))
cat(sprintf("    Priors: %s\n", paste(PRIORS, collapse = " vs ")))
cat(sprintf("    Quick mode: %s\n\n", QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/YFB_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results_rows <- list()

# ============================================================
# Section 1 — Synthetic Validation
# ============================================================
cat("--- Section 1: Synthetic Validation ---\n")

syn <- generate_synthetic_benchmark_data(
  n = cfg$synthetic$n, p = cfg$synthetic$p,
  K_true = cfg$synthetic$k_true, seed = cfg$synthetic$seed
)
cat(sprintf("    DGP: n=%d, p=%d, K_true=%d, censoring=%.1f%%\n",
            syn$n, syn$p, syn$K_true, 100 * syn$censoring_rate))

set.seed(cfg$evaluation$holdout_seed)
split_idx <- stratified_split(syn$status,
                              test_frac = cfg$evaluation$holdout_frac,
                              seed      = cfg$evaluation$holdout_seed)
Y_tr <- syn$Y[split_idx$train, ]
Y_te <- syn$Y[split_idx$test,  ]
t_tr <- syn$time[split_idx$train];   s_tr <- syn$status[split_idx$train]
t_te <- syn$time[split_idx$test];    s_te <- syn$status[split_idx$test]

for (pr in PRIORS) {
  cat(sprintf("    Fitting YFB synthetic — prior_beta='%s' ...\n", pr))
  fit <- fit_cox_on_yf(
    Y_tr, t_tr, s_tr,
    K = K_SYN, max_iter = MAX_ITER, tol = cfg$cavi$tol,
    prior_LF = "point_exponential", prior_beta = pr,
    alpha = ALPHA, lambda = LAMBDA, N_burnin = N_BURNIN,
    cox_warmstart = COX_WARMSTART, normalize_AB = NORM_AB,
    verbose = FALSE
  )
  pred  <- predict_cox_on_yf(Y_te, fit$EF, fit$EBeta)
  c_idx <- as.numeric(concordance(Surv(t_te, s_te) ~ pred$risk_scores)$concordance)
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("      C-index=%.4f | K_eff=%d | EBeta: %s\n",
              c_idx, k_eff, paste(sprintf("%.4f", fit$EBeta), collapse = " ")))
  results_rows[[length(results_rows) + 1]] <- data.frame(
    section    = "1_synthetic",
    model      = "YFB",
    prior_beta = pr,
    cohort     = "synthetic_holdout",
    c_index    = round(c_idx, 4),
    k_eff      = k_eff,
    beta_max   = round(max(abs(fit$EBeta)), 4),
    n_iters    = fit$history$n_iter,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Section 2 — PDAC Training (merged TCGA_PAAD + CPTAC)
# ============================================================
cat("\n--- Section 2: PDAC Training (merged TCGA_PAAD + CPTAC) ---\n")

pdac_available <- dir.exists(PDAC_DATA_ROOT)
if (!pdac_available) {
  cat("    [SKIP] PDAC_DATA_ROOT not found — set PDAC_DATA_ROOT env var to run.\n")
  fit_yfb <- list()
} else {
  cat("    Loading training cohorts ...\n")
  train_cohorts <- c("TCGA_PAAD", "CPTAC")
  train_raw <- lapply(setNames(train_cohorts, train_cohorts), function(ds) {
    cat(sprintf("      Loading %s ...\n", ds))
    load_pdac_raw(ds, PDAC_DATA_ROOT)
  })

  cat("    Preprocessing (v2: intersect-first + QN) ...\n")
  log_flags <- PLATFORM_LOG_TRANSFORM[train_cohorts]
  merged    <- preprocess_merged_cohorts(
    cohort_raw_list     = train_raw,
    log_transform_flags = log_flags,
    top_n               = cfg$preprocessing$top_n_genes,
    rank_transform      = TRUE
  )
  Y_train      <- merged$Y
  time_train   <- unlist(lapply(train_cohorts, function(ds) train_raw[[ds]]$time))
  status_train <- unlist(lapply(train_cohorts, function(ds) train_raw[[ds]]$status))
  train_genes  <- merged$gene_names

  cat(sprintf("    Training matrix: n=%d, p=%d\n", nrow(Y_train), ncol(Y_train)))

  fit_yfb <- list()
  for (pr in PRIORS) {
    cat(sprintf("    Fitting YFB PDAC — prior_beta='%s' (K=%d, max_iter=%d) ...\n",
                pr, K, MAX_ITER))
    fit <- fit_cox_on_yf(
      Y_train, time_train, status_train,
      K = K, max_iter = MAX_ITER, tol = cfg$cavi$tol,
      prior_LF = "point_exponential", prior_beta = pr,
      alpha = ALPHA, lambda = LAMBDA, N_burnin = N_BURNIN,
      cox_warmstart = COX_WARMSTART, normalize_AB = NORM_AB,
      verbose = TRUE
    )
    fit_yfb[[pr]] <- list(fit = fit, train_genes = train_genes)

    k_eff    <- sum(abs(fit$EBeta) > BETA_THRESH)
    cat(sprintf("      Converged iter=%d | K_eff=%d | EBeta (all): %s\n",
                fit$history$n_iter, k_eff,
                paste(sprintf("%.4f", fit$EBeta), collapse = " ")))
    results_rows[[length(results_rows) + 1]] <- data.frame(
      section    = "2_pdac_train",
      model      = "YFB",
      prior_beta = pr,
      cohort     = "merged_TCGA_CPTAC",
      c_index    = NA_real_,
      k_eff      = k_eff,
      beta_max   = round(max(abs(fit$EBeta)), 4),
      n_iters    = fit$history$n_iter,
      stringsAsFactors = FALSE
    )
  }
}

# ============================================================
# Section 3 — External Validation (5 held-out PDAC cohorts)
# ============================================================
cat("\n--- Section 3: External Validation ---\n")

if (!pdac_available || length(fit_yfb) == 0) {
  cat("    [SKIP] PDAC data unavailable — skipping external validation.\n")
} else {
  for (pr in PRIORS) {
    cat(sprintf("  [prior_beta='%s']\n", pr))
    fl       <- fit_yfb[[pr]]
    fit_obj  <- fl$fit
    tr_genes <- fl$train_genes

    for (cohort in EXTERNAL_COHORTS) {
      tryCatch({
        raw     <- load_pdac_raw(cohort, PDAC_DATA_ROOT)
        ext_pre <- preprocess_desurv_cohort(
          Y             = raw$Y,
          gene_names    = raw$gene_names,
          top_n         = NULL,
          log_transform = PLATFORM_LOG_TRANSFORM[[cohort]],
          cohort_name   = cohort
        )
        common <- intersect(tr_genes, ext_pre$gene_names)
        if (length(common) < 10) {
          cat(sprintf("    [SKIP] %s: only %d common genes\n", cohort, length(common)))
          next
        }
        tr_idx  <- match(common, tr_genes)
        ext_idx <- match(common, ext_pre$gene_names)
        Y_ext   <- ext_pre$Y[, ext_idx, drop = FALSE]
        EF_sub  <- fit_obj$EF[tr_idx, , drop = FALSE]

        # Cluster B prediction: Y_ext · EF_sub · beta_tilde
        pred  <- predict_cox_on_yf(Y_ext, EF_sub, fit_obj$EBeta)
        c_idx <- as.numeric(concordance(Surv(raw$time, raw$status) ~ pred$risk_scores)$concordance)
        cat(sprintf("    %s (n=%d, genes=%d): C-index=%.4f\n",
                    cohort, raw$n, length(common), c_idx))
        results_rows[[length(results_rows) + 1]] <- data.frame(
          section    = "3_external",
          model      = "YFB",
          prior_beta = pr,
          cohort     = cohort,
          c_index    = round(c_idx, 4),
          k_eff      = sum(abs(fit_obj$EBeta) > BETA_THRESH),
          beta_max   = round(max(abs(fit_obj$EBeta)), 4),
          n_iters    = fit_obj$history$n_iter,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        cat(sprintf("    [ERROR] %s: %s\n", cohort, conditionMessage(e)))
      })
    }
  }
}

# ============================================================
# Section 4 — Summary Table
# ============================================================
cat("\n--- Section 4: Summary ---\n")

if (length(results_rows) > 0) {
  summary_df <- do.call(rbind, results_rows)
  csv_path   <- file.path(OUT_DIR, "YFB_benchmark_results.csv")
  write.csv(summary_df, csv_path, row.names = FALSE)
  cat(sprintf("  Results saved to: %s\n\n", csv_path))

  ext_rows <- summary_df[summary_df$section == "3_external", ]
  if (nrow(ext_rows) > 0) {
    cat("  External C-index comparison (point_normal vs normal):\n")
    ext_wide <- reshape(
      ext_rows[, c("cohort", "prior_beta", "c_index")],
      idvar = "cohort", timevar = "prior_beta", direction = "wide"
    )
    print(ext_wide, row.names = FALSE)
  }
} else {
  cat("  No results collected.\n")
}

cat("\nYFB benchmark complete.\n")
