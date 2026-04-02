# ==============================================================================
# Script:       run_prior_k_comparison.R
# Purpose:      Compare point_normal vs point_laplace priors across all 7 PDAC
#               datasets, and run auto_prune_K to determine K_effective per
#               dataset.  Results feed the Prior Comparison and K Selection
#               sections of factor_modular_sim_report_PDAC.qmd.
#
#               Outputs (saved to results/tables/):
#                 PDAC_cross_dataset/prior_comparison.csv    — all 7 × 3 rows
#                 PDAC_cross_dataset/k_selection_summary.csv — all 7 × K_max cols
#                 {dataset}_pl/factor_summary_table.csv      — per-dataset PL fit
#                 {dataset}_pl/holdout_cindex.csv            — per-dataset PL holdout
#                 {dataset}_Keff/factor_summary_table.csv    — per-dataset Keff fit
#                 {dataset}_Keff/holdout_cindex.csv          — per-dataset Keff holdout
#
# Usage:
#   Rscript results/modular_sim_factor/run_prior_k_comparison.R
#
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-03-31
# Dependencies: code/fit_modular.R (sources update_*.R), code/predict.R,
#               code/train_test_split.R, code/select_K.R, survival, ebnm, limma
# ==============================================================================

# ==============================================================================
# Setup
# ==============================================================================

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/update_L.R")) {
  # Already at repo root
} else if (file.exists("../code/update_L.R")) {
  setwd("..")
} else {
  stop("Cannot find repo root. Run from project root or set REPO_ROOT env var.")
}

suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
source("code/predict.R")
source("code/train_test_split.R")
source("code/select_K.R")

# PDAC data root
pdac_data_root <- Sys.getenv("PDAC_DATA_ROOT", unset = path.expand(
  paste0("~/Library/CloudStorage/",
         "OneDrive-UniversityofNorthCarolinaatChapelHill/",
         "UNC Dissertation (Liu)/PDAC_data")
))

ALL_DATASETS <- c("TCGA_PAAD", "CPTAC", "Dijk", "Moffitt_GEO_array",
                   "PACA_AU_array", "PACA_AU_seq", "Puleo_array")

PLATFORM_MAP <- c(
  TCGA_PAAD         = "RNA-seq",
  CPTAC             = "Proteomics",
  Dijk              = "RNA-seq",
  Moffitt_GEO_array = "Microarray",
  PACA_AU_array     = "Microarray",
  PACA_AU_seq       = "RNA-seq",
  Puleo_array       = "Microarray"
)

K_MAX   <- 10     # max K for auto_prune
K_FIXED <- 5      # fixed K for prior comparison

# ==============================================================================
# Data Loading Helpers
# (extracted from run_factor_modular_simulation.R to keep this script
#  self-contained; must stay in sync with any upstream changes)
# ==============================================================================

filter_top_genes <- function(Y, gene_names, top_n) {
  if (is.null(top_n) || top_n >= ncol(Y))
    return(list(Y = Y, gene_names = gene_names))
  gene_var <- apply(Y, 2, var)
  keep_idx <- order(gene_var, decreasing = TRUE)[seq_len(top_n)]
  list(Y = Y[, keep_idx, drop = FALSE], gene_names = gene_names[keep_idx])
}

load_real_data <- function(dataset_name, pdac_root, top_n = 5000) {
  if (!dir.exists(pdac_root))
    stop(sprintf("PDAC data root not found: %s", pdac_root))

  tmp_wd    <- tempfile("pdac_wd_")
  dir.create(tmp_wd, showWarnings = FALSE)
  data_link <- file.path(tmp_wd, "data")
  file.symlink(pdac_root, data_link)

  old_wd <- getwd()
  on.exit({ setwd(old_wd); unlink(data_link) }, add = TRUE)

  setwd(tmp_wd)
  source(file.path(pdac_root, "load_data_internal.R"), local = TRUE)
  result <- load_data_internal(dataset_name)
  setwd(old_wd)

  keeps <- which(result$sampInfo$keep == 1)
  Y <- t(result$ex[, keeps])
  rownames(Y) <- NULL; colnames(Y) <- NULL

  time   <- result$sampInfo$time[keeps]
  status <- as.integer(result$sampInfo$event[keeps])

  fi <- result$featInfo
  if (is.data.frame(fi) && "SYMBOL" %in% names(fi)) {
    gene_names <- fi$SYMBOL
  } else if (is.character(fi)) {
    gene_names <- fi
  } else {
    gene_names <- rownames(result$ex)
  }
  if (length(gene_names) != ncol(Y)) {
    gene_names <- seq_len(ncol(Y))
    warning("gene_names length mismatch; falling back to integer indices.")
  }

  filtered   <- filter_top_genes(Y, gene_names, top_n)
  Y          <- filtered$Y
  gene_names <- filtered$gene_names

  # Column-centre (critical for CAVI initialisation)
  Y <- sweep(Y, 2, colMeans(Y), "-")

  stopifnot(is.numeric(Y), !anyNA(Y))
  stopifnot(length(time)   == nrow(Y), !anyNA(time),   all(time > 0))
  stopifnot(length(status) == nrow(Y), !anyNA(status), all(status %in% c(0L, 1L)))

  list(Y = Y, time = time, status = status,
       gene_names = gene_names, n = nrow(Y), p = ncol(Y),
       dataset_name = dataset_name)
}

# ==============================================================================
# Fit + Evaluate Helper
#
# Fits the model with a given prior and K, optionally with holdout evaluation.
# Returns a named list of summary metrics.
# ==============================================================================

#' Fit SBMF and compute C-index metrics for one (dataset, prior, K) combination.
#'
#' @param d           result of load_real_data()
#' @param prior       character: "point_normal" or "point_laplace"
#' @param K_val       integer: number of latent factors
#' @param holdout     logical: compute 80/20 stratified holdout C-index?
#' @param seed        integer: random seed for holdout split
#' @param table_dir   character or NULL: directory to write per-dataset tables
#'
#' @return named list:
#'   Prior, K_fit, C_InSample, C_Train, C_Test, ELBO_Final, N_Iters,
#'   Mean_Sparsity_Pct, Converged
fit_and_evaluate <- function(d, prior, K_val, holdout = TRUE,
                              seed = 42, table_dir = NULL) {
  Y      <- d$Y
  time   <- d$time
  status <- d$status

  cat(sprintf("  [%s | K=%d | prior=%s]\n", d$dataset_name, K_val, prior))

  # --------------------------------------------------------------------------
  # Full-data in-sample fit (for C_InSample)
  # --------------------------------------------------------------------------
  res_full <- fit_supervised_mf_modular(
    Y, time, status, K = K_val, max_iter = 300, tol = 1e-3,
    prior_family = prior, init_method = "svd", verbose = FALSE
  )
  EL_full    <- res_full$EL
  EF_full    <- res_full$EF
  EBeta_full <- res_full$EBeta
  history    <- res_full$history

  fit_cox_full   <- coxph(Surv(time, status) ~ EL_full)
  c_in_sample    <- round(summary(fit_cox_full)$concordance[1], 4)

  elbo_final <- round(tail(history$elbo_full, 1), 2)   # full ELBO (all 5 terms)
  n_iters    <- history$n_iter
  converged  <- history$converged

  # Sparsity: fraction of EF weights with |value| < 1e-6
  mean_sparsity <- round(100 * mean(abs(EF_full) < 1e-6), 2)

  # Factor summary table
  K_eff_cols <- ncol(EL_full)
  if (!is.null(table_dir)) {
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    total_var  <- sum(Y^2)
    pve_pct    <- sapply(seq_len(K_eff_cols), function(k)
      round(sum(outer(EL_full[, k], EF_full[, k])^2) / total_var * 100, 3))
    nonzero_pct <- round(colMeans(abs(EF_full) > 1e-6) * 100, 2)
    summary_tab <- data.frame(
      Factor      = seq_len(K_eff_cols),
      Beta        = round(EBeta_full, 4),
      NonZero_Pct = nonzero_pct,
      PVE_Pct     = pve_pct,
      stringsAsFactors = FALSE
    )
    write.csv(summary_tab,
              file.path(table_dir, "factor_summary_table.csv"), row.names = FALSE)
    conv_df <- data.frame(
      Iteration  = seq_along(history$rmse),
      RMSE       = history$rmse,
      ELBO_Proxy = history$elbo_proxy,
      ELBO_Full  = history$elbo_full
    )
    write.csv(conv_df,
              file.path(table_dir, "convergence_history.csv"), row.names = FALSE)
  }

  # --------------------------------------------------------------------------
  # Hold-out evaluation (80/20 stratified split)
  # --------------------------------------------------------------------------
  c_train <- NA_real_
  c_test  <- NA_real_
  n_test  <- NA_integer_

  if (holdout) {
    sp <- tryCatch(
      stratified_split(status, test_frac = 0.2, seed = seed),
      error = function(e) {
        cat(sprintf("    [Hold-out] Split failed: %s\n", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(sp)) {
      Y_train      <- Y[sp$train_idx, , drop = FALSE]
      Y_test       <- Y[sp$test_idx,  , drop = FALSE]
      time_train   <- time[sp$train_idx]
      time_test    <- time[sp$test_idx]
      status_train <- status[sp$train_idx]
      status_test  <- status[sp$test_idx]
      n_test       <- sp$n_test

      # Fit on training data
      res_train <- fit_supervised_mf_modular(
        Y_train, time_train, status_train, K = K_val, max_iter = 300, tol = 1e-3,
        prior_family = prior, init_method = "svd", verbose = FALSE
      )
      EF_train    <- res_train$EF
      EBeta_train <- res_train$EBeta
      EL_train    <- res_train$EL

      # In-sample (train) C-index
      fit_cox_tr <- coxph(Surv(time_train, status_train) ~ EL_train)
      c_train    <- round(summary(fit_cox_tr)$concordance[1], 4)

      # Project test patients and evaluate
      test_pred  <- tryCatch(
        predict_supervised_mf(Y_test, EF_train, EBeta_train),
        error = function(e) {
          cat(sprintf("    [Hold-out] Predict failed: %s\n", conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(test_pred) && sum(status_test) >= 2) {
        fit_cox_te <- tryCatch(
          coxph(Surv(time_test, status_test) ~ test_pred$risk_scores),
          error = function(e) NULL
        )
        if (!is.null(fit_cox_te))
          c_test <- round(summary(fit_cox_te)$concordance[1], 4)
      } else {
        cat(sprintf("    [Hold-out] n_test=%d with %d events — skipping test C-index\n",
                    n_test, sum(status_test)))
      }

      if (!is.null(table_dir)) {
        holdout_df <- data.frame(
          Method     = c("Null (0.5)", "Supervised (train)", "Supervised (test)"),
          Prior      = prior,
          K_fit      = K_val,
          C_Index    = c(0.5, c_train, c_test),
          Split      = c("—", "train", "test"),
          n_patients = c(NA, sp$n_train, sp$n_test),
          n_events   = c(NA, sum(status_train), sum(status_test)),
          stringsAsFactors = FALSE
        )
        write.csv(holdout_df,
                  file.path(table_dir, "holdout_cindex.csv"), row.names = FALSE)
      }
    }
  }

  list(
    Prior            = prior,
    K_fit            = K_val,
    C_InSample       = c_in_sample,
    C_Train          = c_train,
    C_Test           = c_test,
    ELBO_Final       = elbo_final,
    N_Iters          = n_iters,
    Mean_Sparsity_Pct = mean_sparsity,
    Converged        = converged,
    n_test           = n_test
  )
}

# ==============================================================================
# Main Loop
# ==============================================================================

prior_rows  <- list()   # accumulate prior comparison rows
k_sel_rows  <- list()   # accumulate K selection rows

for (dsname in ALL_DATASETS) {

  cat(sprintf("\n========== %s ==========\n", dsname))

  # --- Load data ---
  d <- tryCatch(
    load_real_data(dsname, pdac_data_root, top_n = 5000),
    error = function(e) {
      cat(sprintf("  ERROR loading %s: %s\n", dsname, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(d)) next

  cat(sprintf("  Loaded: n=%d, p=%d, censoring=%.1f%%\n",
              d$n, d$p, 100 * mean(d$status == 0)))

  # --------------------------------------------------------------------------
  # 1. Read existing point_normal K=5 results
  #    (already computed in the main runner with holdout_eval=TRUE)
  # --------------------------------------------------------------------------
  base_dir <- sprintf("results/tables/%s", dsname)

  pn_row <- tryCatch({
    ho   <- read.csv(file.path(base_dir, "holdout_cindex.csv"),
                     stringsAsFactors = FALSE)
    conv <- read.csv(file.path(base_dir, "convergence_history.csv"),
                     stringsAsFactors = FALSE)
    fsm  <- read.csv(file.path(base_dir, "factor_summary_table.csv"),
                     stringsAsFactors = FALSE)

    c_tr  <- ho$C_Index[grepl("train", ho$Split, ignore.case = TRUE)]
    c_te  <- ho$C_Index[grepl("test",  ho$Split, ignore.case = TRUE) &
                        !grepl("PCA",  ho$Method, ignore.case = TRUE)]
    # Some tables have PCA test row; pick Supervised test
    if (length(c_te) == 0)
      c_te <- NA_real_
    else
      c_te <- c_te[length(c_te)]

    elbo_val     <- tail(conv$ELBO_Proxy, 1)
    n_iters_val  <- nrow(conv)
    sparsity_val <- round(mean(100 - fsm$NonZero_Pct), 2)  # % sparse = 100 - nonzero

    list(
      Prior            = "point_normal",
      K_fit            = nrow(fsm),
      C_InSample       = NA_real_,   # in-sample not in holdout CSV; set to NA
      C_Train          = round(c_tr, 4),
      C_Test           = round(as.numeric(c_te), 4),
      ELBO_Final       = round(elbo_val, 2),
      N_Iters          = n_iters_val,
      Mean_Sparsity_Pct = sparsity_val,
      Converged        = TRUE,
      n_test           = NA_integer_
    )
  }, error = function(e) {
    cat(sprintf("  Could not read existing PN results for %s: %s\n",
                dsname, conditionMessage(e)))
    NULL
  })

  # If read failed, re-run point_normal
  if (is.null(pn_row)) {
    cat(sprintf("  Re-running point_normal K=%d...\n", K_FIXED))
    pn_row <- tryCatch(
      fit_and_evaluate(d, prior = "point_normal", K_val = K_FIXED, holdout = TRUE,
                       table_dir = base_dir),
      error = function(e) { cat(sprintf("  FAILED: %s\n", conditionMessage(e))); NULL }
    )
  }

  if (!is.null(pn_row)) {
    pn_row$Dataset  <- dsname
    pn_row$Platform <- PLATFORM_MAP[dsname]
    pn_row$n        <- d$n
    prior_rows[[length(prior_rows) + 1]] <- pn_row
  }

  # --------------------------------------------------------------------------
  # 2. Fit point_laplace K=5 with holdout
  # --------------------------------------------------------------------------
  pl_dir <- sprintf("results/tables/%s_pl", dsname)
  cat(sprintf("  Fitting point_laplace K=%d...\n", K_FIXED))

  pl_row <- tryCatch(
    fit_and_evaluate(d, prior = "point_laplace", K_val = K_FIXED, holdout = TRUE,
                     table_dir = pl_dir),
    error = function(e) {
      cat(sprintf("  point_laplace fit FAILED: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(pl_row)) {
    pl_row$Dataset  <- dsname
    pl_row$Platform <- PLATFORM_MAP[dsname]
    pl_row$n        <- d$n
    prior_rows[[length(prior_rows) + 1]] <- pl_row
  }

  # --------------------------------------------------------------------------
  # 3. Auto-prune K (K_max=10, point_normal) + K_effective refit with holdout
  # --------------------------------------------------------------------------
  cat(sprintf("  Running auto_prune_K (K_max=%d)...\n", K_MAX))

  prune_res <- tryCatch(
    auto_prune_K(d$Y, d$time, d$status, K_max = K_MAX,
                 prior_family = "point_normal", init_method = "svd",
                 max_iter = 200, tol = 1e-3, verbose = FALSE),
    error = function(e) {
      cat(sprintf("  auto_prune_K FAILED: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(prune_res)) {
    K_eff <- max(1L, prune_res$K_effective)
    cat(sprintf("  K_effective = %d / %d\n", K_eff, K_MAX))

    # Save K selection diagnostics to base dir
    k_sel_df <- data.frame(
      Factor   = seq_len(K_MAX),
      Abs_Beta = round(prune_res$beta, 4),
      PVE_Pct  = round(prune_res$pve * 100, 3),
      Active   = prune_res$active,
      stringsAsFactors = FALSE
    )
    k_sel_path <- file.path(base_dir, "k_selection_pve.csv")
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(k_sel_df, k_sel_path, row.names = FALSE)
    cat(sprintf("  Saved K selection table: %s\n", k_sel_path))

    # Accumulate K selection summary row
    k_row <- c(
      list(Dataset = dsname, K_max = K_MAX, K_effective = K_eff),
      setNames(as.list(round(prune_res$pve * 100, 2)),
               paste0("PVE_F", seq_len(K_MAX))),
      setNames(as.list(round(prune_res$beta, 4)),
               paste0("AbsBeta_F", seq_len(K_MAX)))
    )
    k_sel_rows[[length(k_sel_rows) + 1]] <- k_row

    # Save PVE scree figure
    fig_dir <- sprintf("results/figures/%s", dsname)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    for (ext in c("pdf", "png")) {
      fpath <- file.path(fig_dir, paste0("figK_pve_scree.", ext))
      if (ext == "pdf") pdf(fpath, width = 7, height = 5)
      else              png(fpath, width = 700, height = 500, res = 120)
      par(mar = c(5, 5, 4, 2))
      bar_cols <- ifelse(prune_res$active, "#1f77b4", "#aec7e8")
      barplot(prune_res$pve * 100,
              names.arg = paste0("F", seq_len(K_MAX)),
              col = bar_cols, border = "white",
              main = sprintf("K Selection: PVE per Factor (K_max=%d)\n%s", K_MAX, dsname),
              xlab = "Factor", ylab = "PVE (%)",
              cex.lab = 1.1, cex.main = 1.2, cex.names = 0.9)
      abline(h = 1.0, col = "#d62728", lty = 2, lwd = 1.5)
      legend("topright",
             legend = c(sprintf("Active (K_eff=%d)", K_eff), "Pruned", "1% threshold"),
             fill = c("#1f77b4", "#aec7e8", NA), border = NA,
             lty = c(NA, NA, 2), lwd = c(NA, NA, 1.5),
             col = c(NA, NA, "#d62728"), bty = "n")
      dev.off()
    }
    cat(sprintf("  PVE scree plots saved to %s\n", fig_dir))

    # Refit with K_effective and holdout
    if (K_eff != K_FIXED) {
      keff_dir <- sprintf("results/tables/%s_Keff", dsname)
      cat(sprintf("  Refitting with K_effective=%d...\n", K_eff))
      keff_row <- tryCatch(
        fit_and_evaluate(d, prior = "point_normal", K_val = K_eff, holdout = TRUE,
                         table_dir = keff_dir),
        error = function(e) {
          cat(sprintf("  K_effective refit FAILED: %s\n", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(keff_row)) {
        keff_row$Dataset  <- dsname
        keff_row$Platform <- PLATFORM_MAP[dsname]
        keff_row$n        <- d$n
        keff_row$Prior    <- "point_normal (K_eff)"
        prior_rows[[length(prior_rows) + 1]] <- keff_row
      }
    } else {
      cat(sprintf("  K_effective == K_fixed=%d — no separate refit needed\n", K_eff))
    }

    # --------------------------------------------------------------------------
    # 4. Fit point_laplace at K_effective (the missing experiment:
    #    best prior family at the data-driven number of factors)
    # --------------------------------------------------------------------------
    pl_keff_dir <- sprintf("results/tables/%s_pl_Keff", dsname)
    cat(sprintf("  Fitting point_laplace K=%d (K_eff)...\n", K_eff))
    pl_keff_row <- tryCatch(
      fit_and_evaluate(d, prior = "point_laplace", K_val = K_eff, holdout = TRUE,
                       table_dir = pl_keff_dir),
      error = function(e) {
        cat(sprintf("  point_laplace K_eff fit FAILED: %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(pl_keff_row)) {
      pl_keff_row$Dataset  <- dsname
      pl_keff_row$Platform <- PLATFORM_MAP[dsname]
      pl_keff_row$n        <- d$n
      pl_keff_row$Prior    <- "point_laplace (K_eff)"
      prior_rows[[length(prior_rows) + 1]] <- pl_keff_row
    }
  }
}

# ==============================================================================
# Save Combined Summary Tables
# ==============================================================================

cross_dir <- "results/tables/PDAC_cross_dataset"
dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)

# --- Prior comparison table ---
if (length(prior_rows) > 0) {
  col_order <- c("Dataset", "Platform", "n", "Prior", "K_fit",
                 "C_Train", "C_Test", "ELBO_Final", "N_Iters",
                 "Mean_Sparsity_Pct", "Converged")

  prior_df <- do.call(rbind, lapply(prior_rows, function(r) {
    row <- as.data.frame(r[intersect(names(r), col_order)],
                         stringsAsFactors = FALSE)
    # Fill in any missing columns with NA
    for (col in col_order) {
      if (!col %in% names(row)) row[[col]] <- NA
    }
    row[, col_order]
  }))

  pc_path <- file.path(cross_dir, "prior_comparison.csv")
  write.csv(prior_df, pc_path, row.names = FALSE)
  cat(sprintf("\nPrior comparison table saved: %s\n", pc_path))
  print(prior_df)
}

# --- K selection summary table ---
if (length(k_sel_rows) > 0) {
  k_sel_df <- do.call(rbind, lapply(k_sel_rows, function(r) {
    as.data.frame(r, stringsAsFactors = FALSE)
  }))
  ks_path <- file.path(cross_dir, "k_selection_summary.csv")
  write.csv(k_sel_df, ks_path, row.names = FALSE)
  cat(sprintf("K selection summary saved: %s\n", ks_path))
  print(k_sel_df[, c("Dataset", "K_max", "K_effective",
                      paste0("PVE_F", 1:K_MAX))])
}

cat("\n=== run_prior_k_comparison.R complete ===\n")
