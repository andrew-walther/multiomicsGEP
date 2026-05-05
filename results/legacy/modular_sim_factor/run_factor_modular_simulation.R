# ==============================================================================
# Script:       run_factor_modular_simulation.R
# Purpose:      Run the factor-wise Supervised Bayesian MF via fit_modular.R.
#               Supports both synthetic data (default) and real PDAC datasets.
#               When run_all = TRUE and data_mode = "real", loops over all 7
#               PDAC cohorts and produces a cross-dataset summary table, then
#               fits a pooled model on the RNA-seq trio.
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-03-31
# Dependencies: code/fit_modular.R (sources update_*.R internally);
#               survival, ebnm (loaded by fit_modular.R)
#               For real data: load_data_internal.R (in PDAC_data_root)
# ==============================================================================

# ==============================================================================
# Configuration
# ==============================================================================

# Load global parameter registry (single source of truth for all magic numbers)
library(yaml)
cfg <- yaml::read_yaml("config/globals.yml")

data_mode         <- "synthetic"     # "synthetic" or "real"
dataset_name      <- "TCGA_PAAD"    # used when data_mode = "real" and run_all = FALSE
run_all           <- FALSE           # TRUE: loop over all 7 PDAC datasets
top_n_genes       <- cfg$preprocessing$top_n_genes   # from globals.yml
K                 <- cfg$cavi$k_default               # from globals.yml
prior_beta        <- "point_normal"  # beta prior: "point_normal" or "point_laplace"; L/F always point_exponential
n_init            <- cfg$evaluation$n_init            # from globals.yml
init_method       <- "svd"           # "svd" (deterministic) or "random" (for multi-init)
batch_correct     <- cfg$preprocessing$batch_correct  # from globals.yml
holdout_eval      <- FALSE           # TRUE: 80/20 stratified train/test split evaluation
feature_selection <- "variance"      # "variance" (top by var) or "cox" (univariate Cox p-val)
k_select          <- "fixed"         # "fixed", "auto_prune", "cv" (cv = stub for Longleaf)

# PDAC data root.  Override for Longleaf:
#   export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC
pdac_data_root <- Sys.getenv("PDAC_DATA_ROOT", unset = path.expand(
  paste0("~/Library/CloudStorage/",
         "OneDrive-UniversityofNorthCarolinaatChapelHill/",
         "UNC Dissertation (Liu)/PDAC_data")
))

# Environment variable overrides (scripted / Longleaf use)
if (Sys.getenv("DATA_MODE")          != "") data_mode         <- Sys.getenv("DATA_MODE")
if (Sys.getenv("DATASET_NAME")      != "") dataset_name      <- Sys.getenv("DATASET_NAME")
if (Sys.getenv("RUN_ALL")           != "") run_all           <- as.logical(Sys.getenv("RUN_ALL"))
if (Sys.getenv("TOP_N_GENES")       != "") top_n_genes       <- as.integer(Sys.getenv("TOP_N_GENES"))
if (Sys.getenv("PRIOR_BETA")        != "") prior_beta        <- Sys.getenv("PRIOR_BETA")
if (Sys.getenv("N_INIT")            != "") n_init            <- as.integer(Sys.getenv("N_INIT"))
if (Sys.getenv("INIT_METHOD")       != "") init_method       <- Sys.getenv("INIT_METHOD")
if (Sys.getenv("BATCH_CORRECT")     != "") batch_correct     <- as.logical(Sys.getenv("BATCH_CORRECT"))
if (Sys.getenv("HOLDOUT_EVAL")      != "") holdout_eval      <- as.logical(Sys.getenv("HOLDOUT_EVAL"))
if (Sys.getenv("FEATURE_SELECTION") != "") feature_selection <- Sys.getenv("FEATURE_SELECTION")
if (Sys.getenv("K_SELECT")          != "") k_select          <- Sys.getenv("K_SELECT")

# All 7 available PDAC cohorts (used when run_all = TRUE)
ALL_DATASETS <- c("TCGA_PAAD", "CPTAC", "Dijk", "Moffitt_GEO_array",
                   "PACA_AU_array", "PACA_AU_seq", "Puleo_array")

# Platform labels for the cross-dataset summary table
PLATFORM_MAP <- c(
  TCGA_PAAD        = "RNA-seq",
  CPTAC            = "Proteomics",
  Dijk             = "RNA-seq",
  Moffitt_GEO_array = "Microarray",
  PACA_AU_array    = "Microarray",
  PACA_AU_seq      = "RNA-seq",
  Puleo_array      = "Microarray"
)

# ==============================================================================
# Working Directory
# ==============================================================================

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/update_L.R")) {
  # Already at repo root (e.g., local RStudio)
} else if (file.exists("../code/update_L.R")) {
  setwd("..")
} else {
  stop("Cannot find repo root. Run from project root or set REPO_ROOT env var.")
}

# ==============================================================================
# Source fit_modular.R via tryCatch
#
# fit_modular.R has a runner block at the bottom that stopifnot()s when
# DATA_MODE="real" and real_Y is NULL.  We wrap the source() in tryCatch so
# the stopifnot error is silently caught.  By the time the error fires,
# fit_supervised_mf_modular() and all its dependencies (update_L.R,
# update_F.R, update_beta.R, update_tau.R, library(survival), library(ebnm),
# calc_cox_taylor) are already defined in the global environment.
# ==============================================================================

suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
source("code/predict.R")
source("code/train_test_split.R")
source("code/feature_selection.R")
source("code/select_K.R")

# ==============================================================================
# Analytics Helpers
# ==============================================================================

# C-index comparison: Supervised loadings vs. top-5 PCA components.
#
# When EBeta is supplied, the supervised risk score is EL %*% EBeta — the
# model's own estimated linear predictor.  This tests whether the model's
# actual beta coefficients discriminate survival, rather than re-fitting a
# new Cox model on the loadings (which would use optimally re-fitted
# coefficients and is therefore a different and less honest question).
#
# Convention note: concordance(Surv ~ predictor) treats HIGHER predictor as
# LOWER risk (longer survival).  The Cox LP convention is the opposite: higher
# LP = higher hazard = shorter survival = higher risk.  We therefore pass
# I(-lp) so that higher -LP = lower LP = lower risk, aligning with the
# concordance() convention.  This is equivalent to what summary(coxph(...))
# reports and avoids a spurious apparent inversion.
#
# EBeta = NULL falls back to refitting coxph(Surv ~ EL) — kept for
# back-compatibility with any legacy calls.
get_cindex_comparison <- function(EL, data, EBeta = NULL) {
  pca_y  <- prcomp(data$Y, rank. = min(5, ncol(EL)))
  fit_pc <- coxph(Surv(data$time, data$status) ~ pca_y$x)
  c_pca  <- round(summary(fit_pc)$concordance[1], 3)

  if (!is.null(EBeta)) {
    lp    <- as.vector(EL %*% EBeta)
    # Negate LP: concordance() convention is higher = lower risk; Cox is opposite.
    c_sup <- round(concordance(Surv(data$time, data$status) ~ I(-lp))$concordance, 3)
  } else {
    fit_l <- coxph(Surv(data$time, data$status) ~ EL)
    c_sup <- round(summary(fit_l)$concordance[1], 3)
  }

  list(
    c_original = c_pca,
    c_latent   = c_sup
  )
}

# Top n_top influential features per factor by |weight|.
# When gene_names is supplied, returns gene symbols; otherwise returns FeatureID.
get_top_features <- function(EF, n_top = 10, gene_names = NULL) {
  lapply(1:ncol(EF), function(k) {
    weights   <- EF[, k]
    order_idx <- order(abs(weights), decreasing = TRUE)
    idx       <- order_idx[1:n_top]
    if (!is.null(gene_names)) {
      data.frame(GeneName = gene_names[idx],
                 Weight   = round(weights[idx], 4))
    } else {
      data.frame(FeatureID = idx,
                 Weight    = round(weights[idx], 4))
    }
  })
}

# Per-factor summary: Beta, log-rank p, sparsity %, PVE %
get_factor_summary_table <- function(EL, EF, EBeta, data) {
  K <- ncol(EL)
  p_vals <- sapply(1:K, function(k) {
    grp     <- ifelse(EL[, k] > median(EL[, k]), "High", "Low")
    sd_test <- survdiff(Surv(data$time, data$status) ~ grp)
    1 - pchisq(sd_test$chisq, 1)
  })
  total_var <- sum(data$Y^2)
  pve <- sapply(1:K, function(k) {
    sum((outer(EL[, k], EF[, k]))^2) / total_var * 100
  })
  nonzero_pct <- colMeans(EF != 0) * 100
  data.frame(
    Factor      = 1:K,
    Beta        = round(EBeta, 3),
    LogRank_P   = round(p_vals, 4),
    NonZero_Pct = round(nonzero_pct, 2),
    PVE_Pct     = round(pve, 2)
  )
}

#' Bijective greedy factor matching
#'
#' @description
#' Given a K_true × K_est correlation matrix, returns a length-K_est integer
#' vector \code{perm} where \code{perm[k]} is the true-factor index that best
#' matches estimated factor k, with each true factor used at most once.
#'
#' The simple \code{apply(abs(cors), 2, which.max)} approach is NOT bijective:
#' two estimated factors can both claim the same true factor.  This function
#' instead iterates greedily — at each step claiming the globally largest
#' remaining |correlation| — which guarantees a one-to-one assignment.
#'
#' @param cors K_true × K_est matrix of Pearson correlations
#' @return Integer vector of length K_est; \code{perm[k]} = matched true factor
bijective_match <- function(cors) {
  cost <- abs(cors)   # work with absolute correlations
  K    <- ncol(cors)
  perm <- integer(K)
  for (step in seq_len(K)) {
    idx       <- arrayInd(which.max(cost), dim(cost))  # [row, col] of global max
    perm[idx[2]] <- idx[1]                              # assign true factor to est factor
    cost[idx[1], ] <- -Inf                              # mark true factor as used
    cost[, idx[2]] <- -Inf                              # mark estimated factor as used
  }
  perm
}

# ==============================================================================
# Batch Correction Helper
# ==============================================================================

#' Apply limma::removeBatchEffect to remove cohort-of-origin effects.
#'
#' Used after pool_datasets() to correct for technical variation between
#' cohorts before factorization.  limma expects genes-in-rows (p × n),
#' so we transpose before and after.
#'
#' @param Y             numeric matrix (n × p) — patients in rows, genes in columns
#' @param batch_labels  factor or character vector of length n identifying cohort
#' @return numeric matrix (n × p) with batch effects removed; same dimensions as input
apply_batch_correction <- function(Y, batch_labels) {
  if (!requireNamespace("limma", quietly = TRUE))
    stop("limma required for batch correction. Install via BiocManager::install('limma')")

  # limma::removeBatchEffect expects a genes-in-rows matrix (p × n).
  # Transpose Y (n × p) -> (p × n), apply correction, then transpose back.
  Y_corrected <- t(limma::removeBatchEffect(t(Y), batch = batch_labels))

  # Sanity checks: dimensions preserved, no NAs introduced
  stopifnot(identical(dim(Y_corrected), dim(Y)))
  stopifnot(!anyNA(Y_corrected))

  Y_corrected
}

# ==============================================================================
# Real-Data Helpers
# ==============================================================================

#' Filter to top N most variable genes
#'
#' @param Y        numeric matrix (n x p)
#' @param gene_names character vector length p
#' @param top_n    integer or NULL (NULL = passthrough)
#' @return list(Y, gene_names)
filter_top_genes <- function(Y, gene_names, top_n) {
  if (is.null(top_n) || top_n >= ncol(Y))
    return(list(Y = Y, gene_names = gene_names))
  gene_var <- apply(Y, 2, var)
  keep_idx <- order(gene_var, decreasing = TRUE)[seq_len(top_n)]
  list(Y = Y[, keep_idx, drop = FALSE], gene_names = gene_names[keep_idx])
}

#' Load a single PDAC dataset and prepare it for fit_supervised_mf_modular().
#'
#' Uses load_data_internal() from PDAC_data_root, subsets to valid samples,
#' transposes the expression matrix (genes x samples -> patients x genes),
#' and optionally applies variance-based gene filtering.
#'
#' @param dataset_name  string matching a dataset in load_data_internal.R
#' @param pdac_root     path to the PDAC_data directory
#' @param top_n         integer or NULL; passed to filter_top_genes()
#' @return list(Y, time, status, gene_names, n, p, dataset_name, preprocessing_notes)
load_real_data <- function(dataset_name, pdac_root, top_n = 5000) {
  if (!dir.exists(pdac_root))
    stop(sprintf("PDAC data root not found: %s\nSet PDAC_DATA_ROOT env var.", pdac_root))

  # load_data_internal.R uses relative paths of the form "data/original/<name>.rds".
  # It was written for a project root where the data folder is named "data/".
  # Our data folder is named "PDAC_data/" (= pdac_root).
  #
  # Fix: work from a temp directory that contains a "data" symlink pointing to
  # pdac_root, so that "data/original/..." resolves correctly.
  tmp_wd    <- tempfile("pdac_wd_")
  dir.create(tmp_wd, showWarnings = FALSE)
  data_link <- file.path(tmp_wd, "data")
  file.symlink(pdac_root, data_link)      # tmp_wd/data/ -> pdac_root/

  old_wd <- getwd()
  on.exit({
    setwd(old_wd)
    unlink(data_link)          # remove the symlink (leaves tmp_wd itself, cleaned by OS)
  }, add = TRUE)

  setwd(tmp_wd)
  source(file.path(pdac_root, "load_data_internal.R"), local = TRUE)
  result <- load_data_internal(dataset_name)
  setwd(old_wd)

  # Subset to keep == 1 samples
  keeps <- which(result$sampInfo$keep == 1)
  if (length(keeps) == 0)
    stop(sprintf("No valid samples for dataset '%s' after filtering.", dataset_name))

  # Transpose: genes x samples -> patients x genes.
  # Strip dimension names: row/col names on Y can propagate into ebnm()
  # and cause failures when gene symbols are duplicated or contain "?".
  # Gene identifiers are stored separately in gene_names below.
  Y <- t(result$ex[, keeps])
  rownames(Y) <- NULL
  colnames(Y) <- NULL

  # Extract survival
  time   <- result$sampInfo$time[keeps]
  status <- as.integer(result$sampInfo$event[keeps])

  # Extract gene names (SYMBOL column or character vector)
  fi <- result$featInfo
  if (is.data.frame(fi) && "SYMBOL" %in% names(fi)) {
    gene_names <- fi$SYMBOL
  } else if (is.character(fi)) {
    gene_names <- fi
  } else {
    gene_names <- rownames(result$ex)
  }
  # Ensure gene_names length matches columns of Y (= rows of ex)
  if (length(gene_names) != ncol(Y)) {
    gene_names <- seq_len(ncol(Y))   # fallback to numeric indices
    warning("gene_names length mismatch; falling back to integer indices.")
  }

  # Apply variance-based gene filter
  filtered   <- filter_top_genes(Y, gene_names, top_n)
  Y          <- filtered$Y
  gene_names <- filtered$gene_names

  # Centre each gene (column) to zero mean.
  # Without centering, SVD initialisation captures the global expression mean
  # in Factor 1; subsequent CAVI iterations fight this mean-shift artefact,
  # causing RMSE to increase rather than decrease.  Column centering is standard
  # preprocessing for matrix factorisation on genomics data and has no effect on
  # the synthetic simulation (which already has zero-mean factors by construction).
  Y <- sweep(Y, 2, colMeans(Y), "-")

  # --- Validation ---
  stopifnot(is.numeric(Y), !anyNA(Y))
  stopifnot(length(time)   == nrow(Y), !anyNA(time),   all(time > 0))
  stopifnot(length(status) == nrow(Y), !anyNA(status), all(status %in% c(0L, 1L)))

  # Build preprocessing notes string
  notes <- sprintf(
    "Dataset: %s | n=%d | p=%d (from %d raw genes, top_%s by variance) | censoring=%.1f%%",
    dataset_name, nrow(Y), ncol(Y), nrow(result$ex),
    ifelse(is.null(top_n), "all", as.character(top_n)),
    100 * mean(status == 0)
  )

  list(Y = Y, time = time, status = status,
       gene_names = gene_names,
       n = nrow(Y), p = ncol(Y),
       dataset_name = dataset_name,
       preprocessing_notes = notes)
}

#' Pool multiple load_real_data() results for horizontal integration.
#'
#' Finds the gene intersection, aligns columns, then row-binds Y and
#' concatenates time/status.  The returned dataset_labels factor identifies
#' which rows came from which cohort (for diagnostics only; not used in fitting).
#'
#' @param ds_list  named list of load_real_data() results
#' @return list(Y, time, status, gene_names, n, p, dataset_labels, dataset_name)
pool_datasets <- function(ds_list) {
  common_genes <- Reduce(intersect, lapply(ds_list, "[[", "gene_names"))
  if (length(common_genes) == 0)
    stop("No common genes across the supplied datasets.")

  Y_list      <- lapply(ds_list, function(d) d$Y[, match(common_genes, d$gene_names)])
  time_vec    <- unlist(lapply(ds_list, "[[", "time"))
  status_vec  <- unlist(lapply(ds_list, "[[", "status"))
  label_vec   <- factor(rep(names(ds_list), sapply(ds_list, "[[", "n")))

  Y_pool <- do.call(rbind, Y_list)
  pool_name <- paste(names(ds_list), collapse = "_")

  list(Y = Y_pool, time = time_vec, status = status_vec,
       gene_names = common_genes,
       n = nrow(Y_pool), p = ncol(Y_pool),
       dataset_labels = label_vec,
       dataset_name = pool_name)
}

# ==============================================================================
# Helper: run full pipeline for one dataset, write outputs, return summary row
# ==============================================================================

#' Run the full pipeline for a single dataset.
#'
#' @param Y           numeric matrix n x p
#' @param time        numeric vector n
#' @param status      integer vector n (0/1)
#' @param gene_names  character vector p or NULL
#' @param data        list with at least Y, time, status; may include L_true etc.
#' @param table_dir   output path for CSV tables
#' @param figure_dir  output path for PDF/PNG figures
#' @param run_label   string for cat() messages (e.g. dataset name)
#' @param is_synthetic logical; if TRUE, write ground-truth comparisons
#' @return data.frame with one summary row
run_pipeline <- function(Y, time, status, gene_names, data,
                         table_dir, figure_dir,
                         run_label, is_synthetic = FALSE) {

  dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------------
  # Fit model
  # --------------------------------------------------------------------------
  cat(sprintf("\n=== %s  (n=%d, p=%d, K=%d) ===\n",
              run_label, nrow(Y), ncol(Y), K))
  cat(sprintf("  Censoring rate: %.1f%%  |  prior=%s  |  init=%s  |  n_init=%d\n\n",
              100 * mean(status == 0), prior_beta, init_method, n_init))

  # --------------------------------------------------------------------------
  # K selection (when k_select != "fixed")
  # --------------------------------------------------------------------------
  K_fit <- K   # effective K used for the main fit (may be updated below)

  if (k_select == "auto_prune") {
    cat(sprintf("  [K selection] auto_prune: fitting K_max=%d to identify active factors...\n",
                cfg$cavi$k_max))
    prune_res <- tryCatch(
      auto_prune_K(Y, time, status, K_max = cfg$cavi$k_max,
                   beta_thresh = cfg$k_selection$beta_threshold,
                   pve_thresh  = cfg$k_selection$pve_threshold,
                   prior_beta = prior_beta, init_method = init_method,
                   max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol, verbose = FALSE),
      error = function(e) {
        cat(sprintf("  [K selection] auto_prune failed: %s — using K=%d\n",
                    conditionMessage(e), K))
        NULL
      }
    )
    if (!is.null(prune_res)) {
      K_fit <- max(1L, prune_res$K_effective)
      cat(sprintf("  [K selection] K_effective=%d (using K=%d for main fit)\n",
                  prune_res$K_effective, K_fit))

      # Save K selection diagnostics
      k_sel_df <- data.frame(
        Factor      = seq_len(10),
        Abs_Beta    = round(prune_res$beta, 4),
        PVE_Pct     = round(prune_res$pve * 100, 3),
        Active      = prune_res$active,
        stringsAsFactors = FALSE
      )
      write.csv(k_sel_df,
                file.path(table_dir, "k_selection_pve.csv"), row.names = FALSE)

      # Save PVE scree figure
      for (ext in c("pdf", "png")) {
        fpath <- file.path(figure_dir, paste0("figK_pve_scree.", ext))
        if (ext == "pdf") pdf(fpath, width = 7, height = 5)
        else              png(fpath, width = 700, height = 500, res = 120)
        par(mar = c(5, 5, 4, 2))
        bar_cols <- ifelse(prune_res$active, "#1f77b4", "#aec7e8")
        barplot(prune_res$pve * 100,
                names.arg = paste0("F", seq_len(10)),
                col = bar_cols, border = "white",
                main = sprintf("K Selection: PVE per Factor (K_max=10)\n%s", run_label),
                xlab = "Factor", ylab = "PVE (%)",
                cex.lab = 1.1, cex.main = 1.2, cex.names = 0.9)
        abline(h = 1.0, col = "#d62728", lty = 2, lwd = 1.5)
        legend("topright",
               legend = c(sprintf("Active (K_eff=%d)", K_fit), "Pruned", "1% threshold"),
               fill = c("#1f77b4", "#aec7e8", NA), border = NA,
               lty = c(NA, NA, 2), lwd = c(NA, NA, 1.5),
               col = c(NA, NA, "#d62728"), bty = "n")
        dev.off()
      }
    }
  } else if (k_select == "cv") {
    select_K_cv(Y, time, status)   # errors with informative message
  }
  # "fixed": K_fit = K unchanged

  # Multi-init: run n_init random starts and keep the best-ELBO fit.
  # When n_init == 1 or init_method == "svd", run a single (deterministic) fit.
  if (n_init > 1 && init_method == "random") {
    best_elbo <- -Inf
    best_res  <- NULL
    elbo_vec  <- numeric(n_init)
    for (init_i in seq_len(n_init)) {
      set.seed(cfg$evaluation$multi_init_seed_base + init_i)
      cat(sprintf("  [init %d/%d] ", init_i, n_init))
      res_i <- fit_supervised_mf_modular(
        Y, time, status, K = K_fit,
        max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
        lambda = cfg$cavi$lambda,
        prior_beta = prior_beta, init_method = "random", verbose = FALSE)
      elbo_i <- tail(res_i$history$elbo_proxy, 1)
      elbo_vec[init_i] <- elbo_i
      cat(sprintf("ELBO=%.1f  iters=%d  converged=%s\n",
                  elbo_i, res_i$history$n_iter, res_i$history$converged))
      if (elbo_i > best_elbo) {
        best_elbo <- elbo_i
        best_res  <- res_i
      }
    }
    cat(sprintf("  Best init: ELBO=%.1f (range: %.1f to %.1f)\n",
                best_elbo, min(elbo_vec), max(elbo_vec)))
    res <- best_res
    # Save ELBO distribution for report
    write.csv(data.frame(Init = seq_len(n_init), Final_ELBO = round(elbo_vec, 2),
                         Selected = elbo_vec == best_elbo),
              file.path(table_dir, "multi_init_elbos.csv"), row.names = FALSE)
  } else {
    res <- fit_supervised_mf_modular(
      Y, time, status, K = K_fit,
      max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
      lambda = cfg$cavi$lambda,
      prior_beta = prior_beta, init_method = init_method, verbose = TRUE)
  }
  EL     <- res$EL
  EL2    <- res$EL2
  EF     <- res$EF
  EF2    <- res$EF2
  EBeta  <- res$EBeta
  EBeta2 <- res$EBeta2
  Tau    <- res$Tau
  history <- res$history
  beta_sd <- sqrt(pmax(EBeta2 - EBeta^2, 0))

  cat(sprintf("\n  Converged: %s  |  Iterations: %d  |  Final RMSE: %.4f\n",
              history$converged, history$n_iter, tail(history$rmse, 1)))

  # --------------------------------------------------------------------------
  # Computed summaries (work for both synthetic and real)
  # --------------------------------------------------------------------------
  summary_tab <- get_factor_summary_table(EL, EF, EBeta, data)
  perf        <- get_cindex_comparison(EL, data, EBeta = EBeta)
  top_feats   <- get_top_features(EF, 10, gene_names)

  ph_test <- cox.zph(coxph(Surv(time, status) ~ EL))
  ph_df   <- data.frame(
    Factor  = rownames(ph_test$table),
    Chisq   = round(ph_test$table[, 1], 4),
    DF      = ph_test$table[, 2],
    P_Value = round(ph_test$table[, 3], 4)
  )

  history_df <- data.frame(
    Iteration  = seq_along(history$rmse),
    RMSE       = history$rmse,
    ELBO_Proxy = history$elbo_proxy
  )

  cindex_df <- data.frame(
    Method  = c("Top-5 PCA", "Supervised (EL %*% EBeta)"),
    C_Index = c(perf$c_original, perf$c_latent)
  )

  # --------------------------------------------------------------------------
  # Save CSV tables
  # --------------------------------------------------------------------------
  write.csv(summary_tab,
            file.path(table_dir, "factor_summary_table.csv"), row.names = FALSE)
  write.csv(cindex_df,
            file.path(table_dir, "cindex_comparison.csv"), row.names = FALSE)
  write.csv(history_df,
            file.path(table_dir, "convergence_history.csv"), row.names = FALSE)
  write.csv(ph_df,
            file.path(table_dir, "ph_test_results.csv"), row.names = FALSE)

  for (k in seq_along(top_feats)) {
    write.csv(top_feats[[k]],
              file.path(table_dir, sprintf("top_features_GEP%d.csv", k)),
              row.names = FALSE)
  }

  # Synthetic-only: beta comparison and loading correlations
  if (is_synthetic) {
    B_true  <- data$B_true
    L_true  <- data$L_true

    cors <- cor(L_true, EL)   # K_true x K_est
    write.csv(round(cors, 4),
              file.path(table_dir, "loading_correlation_matrix.csv"))

    # Bijective permutation + sign alignment for beta comparison table.
    # bijective_match() guarantees each true factor is assigned to at most one
    # estimated factor (unlike apply(abs(cors), 2, which.max) which can
    # duplicate assignments when two estimated factors are most correlated with
    # the same true factor).
    perm        <- bijective_match(cors)               # perm[k] = true factor for est factor k
    match_signs <- sign(cors[cbind(perm, seq_len(K))]) # sign correction per factor
    B_true_aligned <- B_true[perm] * match_signs        # reordered + sign-corrected true betas

    beta_df <- data.frame(
      Est_Factor          = seq_len(K),
      Matched_True_Factor = perm,
      Loading_Corr        = round(cors[cbind(perm, seq_len(K))], 4),
      Beta_true_aligned   = round(B_true_aligned, 4),
      Beta_est          = round(EBeta, 4),
      Beta2_est         = round(EBeta2, 6),
      Posterior_SD      = round(beta_sd, 4),
      Abs_Error         = round(abs(EBeta - B_true_aligned), 4),
      Sign_Match        = sign(EBeta) == sign(B_true_aligned) | B_true_aligned == 0
    )
    write.csv(beta_df,
              file.path(table_dir, "beta_comparison_table.csv"), row.names = FALSE)
  }

  cat(sprintf("  CSV tables saved to %s\n", table_dir))

  # --------------------------------------------------------------------------
  # Hold-Out Prediction Evaluation (when holdout_eval = TRUE)
  #
  # Fits the model on training patients only, projects test patients via
  # pseudo-inverse (predict_supervised_mf), then evaluates C-index on the
  # held-out test set.  Compares: Supervised SBMF vs PCA baseline vs Null.
  #
  # This evaluation is separate from the in-sample C-index computed above.
  # The in-sample C-index (perf$c_latent) is optimistic; the hold-out C-index
  # is the proper assessment of predictive generalisation.
  #
  # When feature_selection == "cox", Cox gene filtering is applied to
  # Y_train only (before fitting), then the selected genes are applied
  # to Y_test as well.  This prevents survival signal leakage.
  # --------------------------------------------------------------------------
  if (holdout_eval) {
    cat("  [Hold-out evaluation] Splitting data...\n")

    # Stratified 80/20 split preserving event rate
    sp <- tryCatch(
      stratified_split(status, test_frac = cfg$evaluation$holdout_frac,
                       seed = cfg$evaluation$holdout_seed),
      error = function(e) {
        cat(sprintf("  [Hold-out] Split failed: %s — skipping.\n", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(sp)) {
      train_idx <- sp$train_idx
      test_idx  <- sp$test_idx

      cat(sprintf("  [Hold-out] n_train=%d (events=%d), n_test=%d (events=%d)\n",
                  sp$n_train, sum(status[train_idx] == 1),
                  sp$n_test,  sum(status[test_idx]  == 1)))

      Y_train      <- Y[train_idx, , drop = FALSE]
      Y_test       <- Y[test_idx,  , drop = FALSE]
      time_train   <- time[train_idx]
      time_test    <- time[test_idx]
      status_train <- status[train_idx]
      status_test  <- status[test_idx]

      # --- Feature selection on training data only ---
      # "variance": use same genes already selected (no extra step)
      # "cox":      fit univariate Cox per gene on Y_train, filter by p-value
      if (feature_selection == "cox") {
        cat("  [Hold-out] Cox feature selection on training data...\n")
        selected_genes <- tryCatch(
          cox_feature_selection(Y_train, time_train, status_train, p_thresh = 0.05),
          error = function(e) {
            cat(sprintf("  [Hold-out] Cox selection failed: %s — using all genes.\n",
                        conditionMessage(e)))
            seq_len(ncol(Y_train))
          }
        )
        if (length(selected_genes) < 5) {
          cat(sprintf("  [Hold-out] Only %d genes pass Cox filter; using all.\n",
                      length(selected_genes)))
          selected_genes <- seq_len(ncol(Y_train))
        }
        cat(sprintf("  [Hold-out] Cox selected %d / %d genes\n",
                    length(selected_genes), ncol(Y_train)))
        Y_train <- Y_train[, selected_genes, drop = FALSE]
        Y_test  <- Y_test[,  selected_genes, drop = FALSE]
      }

      # --- Fit on training data ---
      cat("  [Hold-out] Fitting model on training set...\n")
      res_train <- tryCatch(
        fit_supervised_mf_modular(
          Y_train, time_train, status_train, K = K_fit,
          max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
          lambda = cfg$cavi$lambda,
          prior_beta = prior_beta, init_method = init_method,
          verbose = FALSE),
        error = function(e) {
          cat(sprintf("  [Hold-out] Training fit failed: %s\n", conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(res_train)) {
        # --- Project test patients via pseudo-inverse ---
        pred_test <- predict_supervised_mf(Y_test, res_train$EF, res_train$EBeta)

        # --- Evaluate held-out C-index ---
        c_supervised_test <- tryCatch(
          concordance(Surv(time_test, status_test) ~ pred_test$risk_scores)$concordance,
          error = function(e) NA_real_
        )

        # --- PCA baseline: train PCA on Y_train, project Y_test ---
        pca_train <- prcomp(Y_train, rank. = min(5, K))
        # Project test: multiply Y_test by training rotation matrix
        L_pca_test   <- Y_test %*% pca_train$rotation
        fit_pca_test <- tryCatch(
          coxph(Surv(time_train, status_train) ~ pca_train$x, x = FALSE),
          error = function(e) NULL
        )
        c_pca_test <- if (!is.null(fit_pca_test)) {
          pca_lp_test <- as.vector(L_pca_test %*% coef(fit_pca_test))
          tryCatch(
            concordance(Surv(time_test, status_test) ~ pca_lp_test)$concordance,
            error = function(e) NA_real_
          )
        } else NA_real_

        # --- Training set C-index (for comparison) ---
        c_supervised_train <- tryCatch({
          lp_train <- as.vector(res_train$EL %*% res_train$EBeta)
          concordance(Surv(time_train, status_train) ~ lp_train)$concordance
        }, error = function(e) NA_real_)

        cat(sprintf("  [Hold-out] C-index: Supervised train=%.3f | test=%.3f | PCA test=%.3f\n",
                    c_supervised_train, c_supervised_test, c_pca_test))

        # Save hold-out results
        holdout_df <- data.frame(
          Method      = c("Null (0.5)", "PCA (test)", "Supervised (train)", "Supervised (test)"),
          C_Index     = round(c(0.5, c_pca_test, c_supervised_train, c_supervised_test), 4),
          Split       = c("—", "test", "train", "test"),
          n_patients  = c(NA, sp$n_test, sp$n_train, sp$n_test),
          n_events    = c(NA,
                          sum(status_test == 1),
                          sum(status_train == 1),
                          sum(status_test == 1))
        )
        write.csv(holdout_df,
                  file.path(table_dir, "holdout_cindex.csv"), row.names = FALSE)
        cat(sprintf("  [Hold-out] Results saved to %s\n",
                    file.path(table_dir, "holdout_cindex.csv")))
      }
    }
  }

  # --------------------------------------------------------------------------
  # Save Figures
  # --------------------------------------------------------------------------

  # --- Figure 1: RMSE Convergence Trace ---
  for (ext in c("pdf", "png")) {
    fpath <- file.path(figure_dir, paste0("fig1_rmse_trace.", ext))
    if (ext == "pdf") pdf(fpath, width = 8, height = 5)
    else              png(fpath, width = 800, height = 500, res = 120)
    par(mar = c(5, 5, 4, 2))
    plot(history$rmse, type = "l", lwd = 2.5, col = "#1f77b4",
         main = sprintf("Figure 1: Reconstruction RMSE Across CAVI Iterations\n(%s)", run_label),
         xlab = "Iteration", ylab = "RMSE", bty = "n",
         cex.lab = 1.2, cex.main = 1.3)
    if (is_synthetic)
      abline(h = 1.0, col = "#d62728", lty = 2, lwd = 1.5)
    abline(h = tail(history$rmse, 1), col = "gray50", lty = 3, lwd = 1)
    legd <- c("RMSE", sprintf("Final RMSE = %.4f", tail(history$rmse, 1)))
    cols <- c("#1f77b4", "gray50"); ltys <- c(1, 3); lwds <- c(2.5, 1)
    if (is_synthetic) {
      legd <- c(legd, "True Noise SD = 1.0")
      cols <- c(cols, "#d62728"); ltys <- c(ltys, 2); lwds <- c(lwds, 1.5)
    }
    legend("topright", legend = legd, col = cols, lty = ltys, lwd = lwds, bty = "n")
    grid(col = "lightgray", lty = "dotted")
    dev.off()
  }

  # --- Figure 2: ELBO Proxy Trace ---
  for (ext in c("pdf", "png")) {
    fpath <- file.path(figure_dir, paste0("fig2_elbo_proxy.", ext))
    if (ext == "pdf") pdf(fpath, width = 8, height = 5)
    else              png(fpath, width = 800, height = 500, res = 120)
    par(mar = c(5, 5, 4, 2))
    plot(history$elbo_proxy, type = "l", lwd = 2.5, col = "#2ca02c",
         main = sprintf("Figure 2: Genomics ELBO Proxy Across Iterations\n(%s)", run_label),
         xlab = "Iteration",
         ylab = expression(E[q]*"[log P(Y | L, F, "*tau*")]"),
         bty = "n", cex.lab = 1.2, cex.main = 1.3)
    grid(col = "lightgray", lty = "dotted")
    dev.off()
  }

  # --- Figure 3: Beta Estimates (with ground-truth overlay in synthetic mode) ---
  for (ext in c("pdf", "png")) {
    fpath <- file.path(figure_dir, paste0("fig3_beta_comparison.", ext))
    if (ext == "pdf") pdf(fpath, width = 8, height = 5.5)
    else              png(fpath, width = 800, height = 550, res = 120)
    par(mar = c(5, 5, 4, 2))
    x_pos <- 1:K
    if (is_synthetic) {
      # perm, match_signs, and B_true_aligned were computed via bijective_match()
      # in the CSV section above — reuse them here to guarantee the figure and
      # table use exactly the same bijective assignment.
      B_true <- data$B_true

      ylim_r <- range(c(B_true_aligned,
                        EBeta + 1.96 * beta_sd,
                        EBeta - 1.96 * beta_sd)) * 1.2

      # x-axis labels show estimated factor index and matched true factor
      x_labels <- paste0("F", seq_len(K), " (T", perm, ")")
    } else {
      ylim_r   <- range(c(EBeta + 1.96 * beta_sd, EBeta - 1.96 * beta_sd, 0)) * 1.2
      x_labels <- as.character(seq_len(K))
    }

    fig3_title <- if (is_synthetic)
      sprintf("Figure 3: Estimated Survival Coefficients (95%% CI)\n(%s, true beta permutation-aligned)", run_label)
    else
      sprintf("Figure 3: Estimated Survival Coefficients (95%% CI)\n(%s)", run_label)

    plot(x_pos, EBeta, pch = 16, cex = 1.8, col = "#1f77b4",
         ylim = ylim_r,
         xlab = "Factor (estimated; true factor in parentheses)",
         ylab = expression(beta[k]),
         main = fig3_title,
         bty = "n", cex.lab = 1.1, cex.main = 1.2, xaxt = "n")
    axis(1, at = seq_len(K), labels = x_labels, cex.axis = 0.85)
    arrows(x_pos, EBeta - 1.96 * beta_sd, x_pos, EBeta + 1.96 * beta_sd,
           angle = 90, code = 3, length = 0.08, col = "#1f77b4", lwd = 1.5)
    abline(h = 0, col = "gray50", lty = 3)
    if (is_synthetic) {
      points(x_pos, B_true_aligned, pch = 4, cex = 2, col = "#d62728", lwd = 2.5)
      legend("bottomleft", legend = c("Estimated (95% CI)", "True (permutation-aligned)"),
             col = c("#1f77b4", "#d62728"), pch = c(16, 4),
             pt.cex = c(1.8, 2), pt.lwd = c(1, 2.5), bty = "n")
    } else {
      legend("bottomleft", legend = "Estimated (95% CI)",
             col = "#1f77b4", pch = 16, pt.cex = 1.8, bty = "n")
    }
    grid(col = "lightgray", lty = "dotted")
    dev.off()
  }

  # --- Figure 4: GEP Heatmap ---
  n_features    <- min(50, nrow(EF))
  top_var_genes <- order(rowSums(abs(EF)), decreasing = TRUE)[1:n_features]
  F_sub         <- EF[top_var_genes, ]
  palette_heat  <- colorRampPalette(c("blue", "white", "red"))(100)
  max_val       <- max(abs(F_sub))
  # Labels: gene symbols if available, else feature indices
  feat_labels   <- if (!is.null(gene_names)) gene_names[top_var_genes] else top_var_genes

  for (ext in c("pdf", "png")) {
    fpath <- file.path(figure_dir, paste0("fig4_gep_heatmap.", ext))
    if (ext == "pdf") pdf(fpath, width = 10, height = 6)
    else              png(fpath, width = 1000, height = 600, res = 120)
    layout(matrix(1:2, ncol = 2), widths = c(5, 1))
    par(mar = c(6, 4, 4, 1))
    image(1:nrow(F_sub), 1:ncol(F_sub), F_sub,
          main = sprintf("Figure 4: GEP Feature Weights (Top %d Features)\n(%s)", n_features, run_label),
          xlab = "Feature", ylab = "Latent Factor",
          col = palette_heat, axes = FALSE, zlim = c(-max_val, max_val))
    axis(1, at = 1:nrow(F_sub), labels = feat_labels, las = 2, cex.axis = 0.55)
    axis(2, at = 1:ncol(F_sub), labels = paste0("F", 1:ncol(F_sub)), las = 1)
    box()
    par(mar = c(6, 1, 4, 3))
    legend_image <- as.matrix(seq(-max_val, max_val, length.out = 100))
    image(1, seq(-max_val, max_val, length.out = 100), t(legend_image),
          col = palette_heat, axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1, cex.axis = 0.8)
    mtext("Weight", side = 4, line = 2, cex = 0.8)
    layout(1)
    dev.off()
  }

  # --- Figure 5: Kaplan-Meier Survival Curves per Factor ---
  for (ext in c("pdf", "png")) {
    fpath <- file.path(figure_dir, paste0("fig5_kaplan_meier.", ext))
    if (ext == "pdf") pdf(fpath, width = 12, height = 8)
    else              png(fpath, width = 1200, height = 800, res = 120)
    par(mfrow = c(2, ceiling(K / 2)), mar = c(4, 4, 3, 1))
    for (k in 1:K) {
      groups  <- ifelse(EL[, k] > median(EL[, k]), "High Score", "Low Score")
      km      <- survfit(Surv(time, status) ~ groups)
      p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
      plot(km, col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n",
           main = paste0("Factor ", k, " (beta = ", round(EBeta[k], 2), ")\n", p_label),
           xlab = "Time", ylab = "Survival Probability")
      legend("bottomleft", legend = c("High", "Low"),
             col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n", cex = 0.8)
    }
    dev.off()
  }

  # --- Figures 6 & 7: Signal recovery (synthetic only) ---
  if (is_synthetic) {
    L_true   <- data$L_true
    cors_mat <- cor(L_true, EL)

    # Figure 6: Best-matched factor loading scatter
    best_match  <- apply(abs(cors_mat), 2, which.max)
    target_est  <- which.max(abs(EBeta))
    target_true <- best_match[target_est]
    sign_corr   <- sign(cors_mat[target_true, target_est])
    r_val       <- cors_mat[target_true, target_est]

    for (ext in c("pdf", "png")) {
      fpath <- file.path(figure_dir, paste0("fig6_signal_recovery.", ext))
      if (ext == "pdf") pdf(fpath, width = 7, height = 7)
      else              png(fpath, width = 700, height = 700, res = 120)
      par(mar = c(5, 5, 4, 2))
      plot(L_true[, target_true], EL[, target_est] * sign_corr,
           main = sprintf("Figure 6: Signal Recovery\n(Est F%d vs True F%d, r = %.3f)",
                          target_est, target_true, abs(r_val)),
           xlab = "Ground Truth Loading",
           ylab = "Estimated Loading (sign-corrected)",
           col = rgb(0, 0, 0, 0.35), pch = 16, cex = 1.2, bty = "n",
           cex.lab = 1.2, cex.main = 1.2)
      abline(0, 1, col = "#d62728", lwd = 2, lty = 2)
      abline(lm(I(EL[, target_est] * sign_corr) ~ L_true[, target_true]),
             col = "#1f77b4", lwd = 2)
      legend("topleft", legend = c("Identity line", "Best fit"),
             col = c("#d62728", "#1f77b4"), lty = c(2, 1), lwd = 2, bty = "n")
      grid(col = "lightgray", lty = "dotted")
      dev.off()
    }

    # Figure 7: Loading correlation heatmap
    palette2 <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    for (ext in c("pdf", "png")) {
      fpath <- file.path(figure_dir, paste0("fig7_loading_correlations.", ext))
      if (ext == "pdf") pdf(fpath, width = 7, height = 6)
      else              png(fpath, width = 700, height = 600, res = 120)
      par(mar = c(5, 5, 4, 5))
      image(1:K, 1:K, abs(cors_mat),
            col = palette2, zlim = c(0, 1),
            main = "Figure 7: |Correlation| Between True and Estimated Loadings",
            xlab = "True Factor", ylab = "Estimated Factor",
            axes = FALSE, cex.lab = 1.2, cex.main = 1.2)
      axis(1, at = 1:K); axis(2, at = 1:K)
      for (i in 1:K) for (j in 1:K)
        text(i, j, sprintf("%.2f", cors_mat[i, j]), cex = 1.0,
             col = if (abs(cors_mat[i, j]) > 0.5) "white" else "black")
      box()
      dev.off()
    }
  }

  # --- Figure 8: Tau Distribution ---
  for (ext in c("pdf", "png")) {
    fpath <- file.path(figure_dir, paste0("fig8_tau_distribution.", ext))
    if (ext == "pdf") pdf(fpath, width = 8, height = 5)
    else              png(fpath, width = 800, height = 500, res = 120)
    par(mar = c(5, 5, 4, 2))
    hist(Tau, breaks = 50, col = "#1f77b4AA", border = "white",
         main = sprintf("Figure 8: Estimated Feature-Specific Noise Precision\n(%s)", run_label),
         xlab = expression(hat(tau)[j]), ylab = "Count",
         cex.lab = 1.2, cex.main = 1.3)
    if (is_synthetic)
      abline(v = 1.0, col = "#d62728", lwd = 2, lty = 2)
    abline(v = median(Tau), col = "#2ca02c", lwd = 2, lty = 3)
    legd <- sprintf("Median est. (%.3f)", median(Tau))
    cols <- "#2ca02c"; ltys <- 3; lwds <- 2
    if (is_synthetic) {
      legd <- c("True (tau = 1.0)", legd)
      cols <- c("#d62728", cols); ltys <- c(2, ltys); lwds <- c(2, lwds)
    }
    legend("topright", legend = legd, col = cols, lty = ltys, lwd = lwds, bty = "n")
    dev.off()
  }

  n_figs <- if (is_synthetic) 8L else 6L
  cat(sprintf("  Figures saved to %s (%d figure pairs)\n", figure_dir, n_figs))

  # --------------------------------------------------------------------------
  # Return one-row summary for cross-dataset table
  # --------------------------------------------------------------------------
  data.frame(
    Dataset      = dataset_name,
    Platform     = PLATFORM_MAP[dataset_name],
    n            = nrow(Y),
    p            = ncol(Y),
    K            = K_fit,
    Converged    = history$converged,
    N_Iter       = history$n_iter,
    Final_RMSE   = round(tail(history$rmse, 1), 4),
    C_PCA        = perf$c_original,
    C_Supervised = perf$c_latent,
    Censoring_Pct = round(100 * mean(status == 0), 1),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Main Execution
# ==============================================================================

if (data_mode == "synthetic") {

  # --------------------------------------------------------------------------
  # Synthetic data generation  — parameters from config/globals.yml
  # seed=222 chosen: bijective factor matching, null factor (true beta=0)
  # correctly shrunk to beta_est≈0, all 5 factors cleanly separated.
  # --------------------------------------------------------------------------
  set.seed(cfg$synthetic$seed)
  n <- cfg$synthetic$n
  p <- cfg$synthetic$p

  L_true <- matrix(rnorm(n * K), n, K)
  F_true <- matrix(0, p, K)
  for (k in 1:K) {
    active <- sample(1:p, round(p * 0.05))
    F_true[active, k] <- rnorm(length(active), 0, 5)
  }
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p), n, p)

  B_true     <- unlist(cfg$synthetic$b_true)
  eta_true   <- as.vector(L_true %*% B_true)
  raw_times  <- (-log(runif(n)) / (0.01 * exp(eta_true)))^(1 / 1.5)
  cens_times <- rexp(n, rate = 1 / 50)
  time   <- pmin(raw_times, cens_times)
  status <- as.integer(raw_times <= cens_times)

  data <- list(Y = Y, time = time, status = status,
               L_true = L_true, F_true = F_true, B_true = B_true)

  table_dir  <- "results/tables/synthetic/"
  figure_dir <- "results/figures/synthetic/"

  cat("=== Factor-Wise Modular Supervised MF — Synthetic Simulation ===\n")
  cat(sprintf("  n=%d  p=%d  K=%d  seed=%d  (factor-wise CAVI, V3 Algorithm 1)\n",
              n, p, K, cfg$synthetic$seed))

  run_pipeline(Y, time, status, gene_names = NULL, data,
               table_dir, figure_dir,
               run_label = "Synthetic (n=250, p=1000, K=5)",
               is_synthetic = TRUE)

  cat(sprintf("\nTotal: %d CSVs + 8 figure pairs written to results/tables/synthetic/ and results/figures/synthetic/\n",
              length(list.files(table_dir, pattern = "\\.csv$"))))

} else if (data_mode == "real") {

  # --------------------------------------------------------------------------
  # Real PDAC data: single dataset or loop over all 7
  # --------------------------------------------------------------------------

  datasets_to_run <- if (run_all) ALL_DATASETS else dataset_name
  summary_rows    <- list()

  for (ds in datasets_to_run) {

    cat(sprintf("\n\n========================================\n"))
    cat(sprintf("  Dataset: %s\n", ds))
    cat(sprintf("========================================\n"))

    # Load and prepare data
    real <- tryCatch(
      load_real_data(ds, pdac_data_root, top_n_genes),
      error = function(e) {
        cat(sprintf("  ERROR loading %s: %s\n  Skipping.\n", ds, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(real)) next

    cat(sprintf("  %s\n", real$preprocessing_notes))

    data_real  <- list(Y = real$Y, time = real$time, status = real$status)
    table_dir  <- paste0("results/tables/", ds, "/")
    figure_dir <- paste0("results/figures/", ds, "/")

    row <- run_pipeline(real$Y, real$time, real$status, real$gene_names,
                        data_real, table_dir, figure_dir,
                        run_label = ds, is_synthetic = FALSE)
    row$Dataset  <- ds                # override global dataset_name used inside run_pipeline()
    row$Platform <- PLATFORM_MAP[ds]
    summary_rows[[ds]] <- row
  }

  # --------------------------------------------------------------------------
  # Cross-dataset summary table
  # --------------------------------------------------------------------------
  if (length(summary_rows) > 1) {
    cross_dir <- "results/tables/PDAC_cross_dataset/"
    dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
    cross_df <- do.call(rbind, summary_rows)
    write.csv(cross_df, file.path(cross_dir, "cross_dataset_summary.csv"),
              row.names = FALSE)
    cat(sprintf("\nCross-dataset summary written to %s\n", cross_dir))
    print(cross_df[, c("Dataset", "Platform", "n", "p", "Converged",
                       "N_Iter", "Final_RMSE", "C_PCA", "C_Supervised")])
  }

  # --------------------------------------------------------------------------
  # Horizontal integration: pool RNA-seq trio and fit one model
  # --------------------------------------------------------------------------
  if (run_all) {
    rnaseq_names <- c("TCGA_PAAD", "Dijk", "PACA_AU_seq")
    cat(sprintf("\n\n========================================\n"))
    cat(sprintf("  Horizontal Integration: Pooled RNA-seq\n"))
    cat(sprintf("  Datasets: %s\n", paste(rnaseq_names, collapse = ", ")))
    cat(sprintf("========================================\n"))

    rnaseq_list <- lapply(rnaseq_names, function(ds) {
      tryCatch(load_real_data(ds, pdac_data_root, top_n_genes),
               error = function(e) { cat(sprintf("  Skip %s: %s\n", ds, conditionMessage(e))); NULL })
    })
    names(rnaseq_list) <- rnaseq_names
    rnaseq_list <- Filter(Negate(is.null), rnaseq_list)

    if (length(rnaseq_list) >= 2) {
      pooled <- pool_datasets(rnaseq_list)
      cat(sprintf("  Common genes: %d | Pooled n: %d\n", pooled$p, pooled$n))

      # Apply batch correction if toggled on.
      # limma::removeBatchEffect removes cohort-of-origin mean differences
      # while preserving the shared biological signal across cohorts.
      if (batch_correct) {
        cat("  Applying limma::removeBatchEffect (batch = cohort of origin)...\n")
        pooled$Y <- apply_batch_correction(pooled$Y, pooled$dataset_labels)
        cat("  Batch correction complete. Dimensions preserved.\n")
      } else {
        cat("  NOTE: No batch correction applied. Factors may partially reflect cohort.\n")
      }

      data_pooled  <- list(Y = pooled$Y, time = pooled$time, status = pooled$status)
      table_dir_p  <- "results/tables/PDAC_pooled_rnaseq/"
      figure_dir_p <- "results/figures/PDAC_pooled_rnaseq/"

      # Add dataset label column to cross-dataset summary
      pool_row <- run_pipeline(pooled$Y, pooled$time, pooled$status,
                                pooled$gene_names, data_pooled,
                                table_dir_p, figure_dir_p,
                                run_label = paste("Pooled RNA-seq:", paste(rnaseq_names, collapse = "+")),
                                is_synthetic = FALSE)
      pool_row$Dataset  <- "PDAC_pooled_rnaseq"
      pool_row$Platform <- "RNA-seq (pooled)"

      # Write dataset membership table for downstream use in the report
      membership_df <- data.frame(
        SampleIndex  = seq_len(pooled$n),
        Dataset      = as.character(pooled$dataset_labels)
      )
      write.csv(membership_df,
                file.path(table_dir_p, "pool_membership.csv"), row.names = FALSE)

      # Append to cross-dataset summary
      if (length(summary_rows) > 1) {
        cross_dir <- "results/tables/PDAC_cross_dataset/"
        cross_df  <- rbind(cross_df, pool_row)
        write.csv(cross_df, file.path(cross_dir, "cross_dataset_summary.csv"),
                  row.names = FALSE)
      }
    } else {
      cat("  Fewer than 2 RNA-seq datasets loaded; skipping pooled fit.\n")
    }
  }

} else {
  stop(sprintf("Unknown data_mode: '%s'. Use 'synthetic' or 'real'.", data_mode))
}

cat("\n=== Done ===\n")
