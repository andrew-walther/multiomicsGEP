# ============================================================
# Script:       run_LB_benchmark.R
# Purpose:      Cluster A (LB) benchmark — side-by-side comparison of
#               prior_beta="point_normal" vs "normal" on synthetic validation,
#               PDAC training (merged TCGA_PAAD + CPTAC), and 5 external cohorts.
#               Replaces run_cluster_a_smoke.R + run_cluster_a_external.R.
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-05-04
# Dependencies: code/fit_modular.R, code/predict.R, code/train_test_split.R,
#               code/preprocess_desurv.R, results/benchmark_sim/benchmark_helpers.R
# Usage:        Rscript results/benchmark_sim/run_LB_benchmark.R [--quick] [--train-mode merged|tcga_only|cptac_only]
# ============================================================

# --------------------------------------------------------------------------
# Setup — working directory and config
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args   # smaller K + fewer iters for smoke testing

TRAIN_MODE <- "merged"
if ("--train-mode" %in% args) {
  TRAIN_MODE <- args[which(args == "--train-mode") + 1]
}
stopifnot(TRAIN_MODE %in% c("merged", "tcga_only", "cptac_only"))
RUN_SYNTHETIC <- (TRAIN_MODE == "merged")

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")  # PDAC constants + load_pdac_raw + generate_synthetic_benchmark_data
tryCatch(source("code/fit_modular.R"), error = function(e) invisible(NULL))  # fit_supervised_mf_modular
source("code/predict.R")
source("code/train_test_split.R")
source("code/preprocess_desurv.R")
source("code/select_alpha_cv.R")  # select_alpha_cv() — LB-only alpha CV via 5-fold stratified CV

# Benchmark defaults from globals.yml
K          <- if (QUICK_MODE) 5 else if (TRAIN_MODE == "merged") cfg$benchmark$k_pdac else cfg$benchmark$k_pdac_single
K_SYN      <- if (QUICK_MODE) 5 else cfg$benchmark$k_pdac_synthetic
ALPHA      <- cfg$benchmark$alpha
LAMBDA     <- cfg$benchmark$lambda
N_BURNIN   <- cfg$benchmark$n_burnin
NORM_AB    <- cfg$benchmark$normalize_ab
PRIORS     <- cfg$benchmark$prior_beta_compare
BETA_THRESH <- cfg$k_selection$beta_threshold
MAX_ITER   <- if (QUICK_MODE) 30 else cfg$cavi$max_iter

cat("=== LB Benchmark (Cluster A — eta = L·beta) ===\n")
cat(sprintf("    K=%d (synthetic K=%d) | alpha=%s | lambda=%.2f | N_burnin=%d\n",
            K, K_SYN,
            if (QUICK_MODE) sprintf("%.2f (fixed)", ALPHA) else sprintf("%.2f (CV-selected per mode)", ALPHA),
            LAMBDA, N_BURNIN))
cat(sprintf("    Priors: %s\n", paste(PRIORS, collapse = " vs ")))
cat(sprintf("    Train mode: %s | Quick mode: %s\n\n", TRAIN_MODE, QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/LB_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Accumulate one-row results for each (section, prior) combination
results_rows <- list()
ebeta_rows   <- list()

# ============================================================
# Section 1 — Synthetic Validation
# ============================================================
if (RUN_SYNTHETIC) {
  cat("--- Section 1: Synthetic Validation ---\n")

  syn <- generate_synthetic_benchmark_data(
    n = cfg$synthetic$n, p = cfg$synthetic$p,
    K_true = cfg$synthetic$k_true, seed = cfg$synthetic$seed
  )
  cat(sprintf("    DGP: n=%d, p=%d, K_true=%d, censoring=%.1f%%\n",
              syn$n, syn$p, syn$K_true, 100 * syn$censoring_rate))

  # 80/20 holdout split
  set.seed(cfg$evaluation$holdout_seed)
  split_idx <- stratified_split(syn$status,
                                test_frac = cfg$evaluation$holdout_frac,
                                seed      = cfg$evaluation$holdout_seed)
  Y_tr <- syn$Y[split_idx$train, ]
  Y_te <- syn$Y[split_idx$test,  ]
  t_tr <- syn$time[split_idx$train];   s_tr <- syn$status[split_idx$train]
  t_te <- syn$time[split_idx$test];    s_te <- syn$status[split_idx$test]

  for (pr in PRIORS) {
    cat(sprintf("    Fitting LB synthetic — prior_beta='%s' ...\n", pr))
    fit <- fit_supervised_mf_modular(
      Y_tr, t_tr, s_tr,
      K = K_SYN, max_iter = MAX_ITER, tol = cfg$cavi$tol,
      prior_LF = "point_exponential", prior_beta = pr,
      alpha = ALPHA, lambda = LAMBDA, N_burnin = N_BURNIN,
      normalize_AB = NORM_AB, verbose = FALSE
    )
    pred     <- predict_supervised_mf(Y_te, fit$EF, fit$EBeta)
    c_idx    <- as.numeric(concordance(Surv(t_te, s_te) ~ pred$risk_scores)$concordance)
    k_eff    <- sum(abs(fit$EBeta) > BETA_THRESH)
    beta_rng <- range(fit$EBeta)
    cat(sprintf("      C-index=%.4f | K_eff=%d | EBeta range [%.3e, %.3e]\n",
                c_idx, k_eff, beta_rng[1], beta_rng[2]))
    results_rows[[length(results_rows) + 1]] <- data.frame(
      train_mode = TRAIN_MODE,
      section    = "1_synthetic",
      model      = "LB",
      prior_beta = pr,
      cohort     = "synthetic_holdout",
      c_index    = round(c_idx, 4),
      k_eff      = k_eff,
      beta_max   = round(max(abs(fit$EBeta)), 4),
      n_iters    = fit$history$n_iter,
      stringsAsFactors = FALSE
    )
  }
} else {
  cat("--- Section 1: Skipped (synthetic only runs for --train-mode merged) ---\n")
}

# ============================================================
# Section 2 — PDAC Training
# ============================================================
cat(sprintf("\n--- Section 2: PDAC Training (%s) ---\n", TRAIN_MODE))

pdac_available <- dir.exists(PDAC_DATA_ROOT)
if (!pdac_available) {
  cat("    [SKIP] PDAC_DATA_ROOT not found — set PDAC_DATA_ROOT env var to run.\n")
  fit_lb <- list()  # placeholder so Sections 3/4 can be skipped gracefully
} else {
  if (TRAIN_MODE == "merged") {
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
    train_cohort_label <- "merged_TCGA_CPTAC"
  } else {
    ds <- if (TRAIN_MODE == "tcga_only") "TCGA_PAAD" else "CPTAC"
    cat(sprintf("    Loading single cohort: %s ...\n", ds))
    raw <- load_pdac_raw(ds, PDAC_DATA_ROOT)
    pre <- preprocess_desurv_cohort(
      Y             = raw$Y,
      gene_names    = raw$gene_names,
      top_n         = cfg$preprocessing$top_n_genes,
      log_transform = PLATFORM_LOG_TRANSFORM[[ds]],
      cohort_name   = ds
    )
    Y_train      <- pre$Y
    time_train   <- raw$time
    status_train <- raw$status
    train_genes  <- pre$gene_names
    train_cohort_label <- ds
  }

  cat(sprintf("    Training matrix: n=%d, p=%d\n", nrow(Y_train), ncol(Y_train)))

  # Alpha CV — skip in QUICK_MODE to keep smoke runs fast
  if (!QUICK_MODE) {
    cat("    Running alpha CV (5-fold) ...\n")
    cv_res <- select_alpha_cv(
      Y_train, time_train, status_train,
      alpha_grid   = cfg$cavi$alpha_grid,
      n_folds      = cfg$cavi$n_cv_folds,
      K_max        = K,
      max_iter     = MAX_ITER,
      tol          = cfg$cavi$tol,
      prior_LF     = "point_exponential",
      lambda       = LAMBDA,
      N_burnin     = N_BURNIN,
      normalize_AB = NORM_AB,
      verbose      = FALSE
    )
    ALPHA <- cv_res$alpha_opt
    cat(sprintf("    Alpha CV selected: %.2f (rule: %s)\n", ALPHA, cv_res$selection_rule))
    cat("    CV table:\n")
    print(cv_res$cv_table)
  } else {
    cat(sprintf("    Alpha CV skipped (QUICK_MODE) — using fixed alpha=%.2f\n", ALPHA))
  }

  fit_lb <- list()   # store one fit per prior for Section 3
  for (pr in PRIORS) {
    cat(sprintf("    Fitting LB PDAC — prior_beta='%s' (K=%d, max_iter=%d) ...\n",
                pr, K, MAX_ITER))
    fit <- fit_supervised_mf_modular(
      Y_train, time_train, status_train,
      K = K, max_iter = MAX_ITER, tol = cfg$cavi$tol,
      prior_LF = "point_exponential", prior_beta = pr,
      alpha = ALPHA, lambda = LAMBDA, N_burnin = N_BURNIN,
      normalize_AB = NORM_AB, verbose = TRUE
    )
    fit_lb[[pr]] <- list(fit = fit, train_genes = train_genes)

    k_eff    <- sum(abs(fit$EBeta) > BETA_THRESH)
    beta_rng <- range(fit$EBeta)
    cat(sprintf("      Converged iter=%d | K_eff=%d | EBeta: [%.3e, %.3e]\n",
                fit$history$n_iter, k_eff, beta_rng[1], beta_rng[2]))
    cat(sprintf("      EBeta (all): %s\n",
                paste(sprintf("%.4f", fit$EBeta), collapse = " ")))
    results_rows[[length(results_rows) + 1]] <- data.frame(
      train_mode = TRAIN_MODE,
      section    = "2_pdac_train",
      model      = "LB",
      prior_beta = pr,
      cohort     = train_cohort_label,
      c_index    = NA_real_,
      k_eff      = k_eff,
      beta_max   = round(max(abs(fit$EBeta)), 4),
      n_iters    = fit$history$n_iter,
      stringsAsFactors = FALSE
    )
    ebeta_rows[[length(ebeta_rows) + 1]] <- data.frame(
      train_mode = TRAIN_MODE,
      section    = "2_pdac_train",
      prior_beta = pr,
      cohort     = NA_character_,
      factor_k   = seq_along(fit$EBeta),
      ebeta      = fit$EBeta,
      stringsAsFactors = FALSE
    )
  }
}

# ============================================================
# Section 3 — External Validation (5 held-out PDAC cohorts)
# ============================================================
cat("\n--- Section 3: External Validation ---\n")

if (!pdac_available || length(fit_lb) == 0) {
  cat("    [SKIP] PDAC data unavailable — skipping external validation.\n")
} else {
  for (pr in PRIORS) {
    cat(sprintf("  [prior_beta='%s']\n", pr))
    fl       <- fit_lb[[pr]]
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
        # Align to training genes
        common <- intersect(tr_genes, ext_pre$gene_names)
        if (length(common) < 10) {
          cat(sprintf("    [SKIP] %s: only %d common genes\n", cohort, length(common)))
          next
        }
        tr_idx  <- match(common, tr_genes)
        ext_idx <- match(common, ext_pre$gene_names)
        Y_ext   <- ext_pre$Y[, ext_idx, drop = FALSE]
        EF_sub  <- fit_obj$EF[tr_idx, , drop = FALSE]

        pred  <- predict_supervised_mf(Y_ext, EF_sub, fit_obj$EBeta)
        c_idx <- as.numeric(concordance(Surv(raw$time, raw$status) ~ pred$risk_scores)$concordance)
        cat(sprintf("    %s (n=%d, genes=%d): C-index=%.4f\n",
                    cohort, raw$n, length(common), c_idx))
        results_rows[[length(results_rows) + 1]] <- data.frame(
          train_mode = TRAIN_MODE,
          section    = "3_external",
          model      = "LB",
          prior_beta = pr,
          cohort     = cohort,
          c_index    = round(c_idx, 4),
          k_eff      = sum(abs(fit_obj$EBeta) > BETA_THRESH),
          beta_max   = round(max(abs(fit_obj$EBeta)), 4),
          n_iters    = fit_obj$history$n_iter,
          stringsAsFactors = FALSE
        )
        ebeta_rows[[length(ebeta_rows) + 1]] <- data.frame(
          train_mode = TRAIN_MODE,
          section    = "3_external",
          prior_beta = pr,
          cohort     = cohort,
          factor_k   = seq_along(fit_obj$EBeta),
          ebeta      = fit_obj$EBeta,
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
  summary_df     <- do.call(rbind, results_rows)
  csv_path       <- file.path(OUT_DIR, sprintf("LB_benchmark_results_%s.csv", TRAIN_MODE))
  ebeta_csv_path <- file.path(OUT_DIR, sprintf("LB_ebeta_detail_%s.csv",      TRAIN_MODE))

  write.csv(summary_df, csv_path, row.names = FALSE)
  cat(sprintf("  Results saved to: %s\n", csv_path))

  if (length(ebeta_rows) > 0) {
    ebeta_df <- do.call(rbind, ebeta_rows)
    write.csv(ebeta_df, ebeta_csv_path, row.names = FALSE)
    cat(sprintf("  EBeta detail saved to: %s\n", ebeta_csv_path))
  }
  cat("\n")

  # Side-by-side external C-index table
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

cat("\nLB benchmark complete.\n")
