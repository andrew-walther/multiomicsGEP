# ==============================================================================
# Script:       run_cox_on_yf_benchmark.R
# Purpose:      Cluster B (Cox-on-YF) benchmark — synthetic DGP + PDAC real-data
#               C-index evaluation parallel to run_ssbmf_benchmark.R (Cluster A).
#
# MODEL:        eta_i = (y_i · F) · beta_tilde   [Cluster B]
#               vs.  eta_i = l_i · beta           [Cluster A, fit_modular.R]
#
# ALPHA:        Fixed alpha=0.5 for beta update (F update always uses alpha_F=0).
#               Alpha CV for Cluster B is deferred — select_alpha_cv.R hard-codes
#               fit_supervised_mf_modular; a Cluster B CV wrapper is a future task.
#               Fixed alpha matches the alpha used in Cluster A runs for comparison.
#
# USAGE:
#   Rscript results/benchmark_sim/run_cox_on_yf_benchmark.R
#   Rscript results/benchmark_sim/run_cox_on_yf_benchmark.R --quick
#
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-05-04
# Dependencies: code/fit_cox_on_yf.R, code/predict_cox_on_yf.R,
#               code/preprocess_desurv.R, code/train_test_split.R
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/fit_cox_on_yf.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_cox_on_yf.R")) {
  setwd("../..")
} else if (file.exists("../code/fit_cox_on_yf.R")) {
  setwd("..")
}

cfg <- yaml::read_yaml("config/globals.yml")

source("code/train_test_split.R")
source("code/predict_cox_on_yf.R")
source("code/preprocess_desurv.R")
suppressMessages(tryCatch(
  source("code/fit_cox_on_yf.R"),
  error = function(e) invisible(NULL)
))

# Source run_ssbmf_benchmark.R to get load_pdac_raw, preprocess_merged_cohorts,
# and other shared PDAC utilities. Suppress the runner block (DATA_MODE is not set there).
suppressMessages(tryCatch(
  source("results/benchmark_sim/run_ssbmf_benchmark.R"),
  error = function(e) invisible(NULL)
))

# Platform-specific log2+1 transform flag (same as run_ssbmf_benchmark.R)
PLATFORM_LOG_TRANSFORM <- c(
  TCGA_PAAD         = TRUE,
  CPTAC             = FALSE,
  Dijk              = TRUE,
  Moffitt_GEO_array = FALSE,
  PACA_AU_array     = FALSE,
  PACA_AU_seq       = TRUE,
  Puleo_array       = FALSE
)

PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT", unset = path.expand(
  paste0("~/Library/CloudStorage/",
         "OneDrive-UniversityofNorthCarolinaatChapelHill/",
         "UNC Dissertation (Liu)/PDAC_data")
))

TRAINING_COHORTS <- c("TCGA_PAAD", "CPTAC")
EXTERNAL_COHORTS <- c("Dijk", "Moffitt_GEO_array", "PACA_AU_array", "PACA_AU_seq", "Puleo_array")

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Synthetic DGP (identical to run_ssbmf_benchmark.R)
# ==============================================================================

calibrate_censor_scale <- function(event_times, base_censor, target = 0.30, n_iter = 40) {
  lo <- 1e-3; hi <- max(event_times) * 10
  for (i in seq_len(n_iter)) {
    mid <- sqrt(lo * hi)
    if (mean(base_censor * mid < event_times) > target) lo <- mid else hi <- mid
  }
  hi
}

generate_synthetic_benchmark_data <- function(n = 300, p = 1000, K_true = 5,
                                               seed = 222, target_censoring = 0.30) {
  set.seed(seed)
  beta_true <- as.numeric(cfg$synthetic$b_true)
  if (length(beta_true) != K_true)
    stop(sprintf("cfg$synthetic$b_true has length %d but K_true=%d.", length(beta_true), K_true))

  L_true <- matrix(rexp(n * K_true, rate = 1), n, K_true)
  signal_scale <- 0.25
  F_true <- matrix(0, p, K_true)
  for (k in seq_len(K_true)) {
    active <- sample.int(p, size = max(1L, round(0.05 * p)))
    F_true[active, k] <- signal_scale
  }
  tau_true <- rgamma(p, shape = 2, rate = 2)
  E <- sweep(matrix(rnorm(n * p), n, p), 2, sqrt(tau_true), "/") * signal_scale
  Y <- L_true %*% t(F_true) + E

  eta_true <- as.vector(L_true %*% beta_true)
  event_times <- (-log(runif(n)) / (0.01 * exp(eta_true)))^(1 / 1.5)
  censor_base <- rexp(n, rate = 1)
  censor_times <- censor_base * calibrate_censor_scale(event_times, censor_base,
                                                        target = target_censoring)
  time   <- pmin(event_times, censor_times)
  status <- as.integer(event_times <= censor_times)

  list(Y = Y, time = time, status = status,
       L_true = L_true, F_true = F_true, beta_true = beta_true,
       censoring_rate = mean(status == 0))
}

# ==============================================================================
# Main benchmark function
# ==============================================================================

run_cox_on_yf_benchmark <- function(
    output_root      = "results/benchmark_sim/outputs_cox_on_yf",
    synthetic_n      = 300,
    synthetic_p      = 1000,
    K_true           = 5,
    K_max            = 15,
    alpha            = 0.5,   # fixed; no CV for Cluster B yet
    N_burnin         = 5,
    max_iter         = 80,
    tol              = 1e-4,
    seed             = cfg$synthetic$seed,
    holdout_frac     = 0.2,
    prior_beta       = "point_normal",
    run_real_data    = TRUE,
    quick            = FALSE) {

  if (quick) {
    synthetic_n <- 120; synthetic_p <- 300
    K_true <- 5; K_max <- 8
    max_iter <- 35; tol <- 1e-4
  }

  table_dir  <- file.path(output_root, "tables")
  figure_dir <- file.path(output_root, "figures")
  ensure_dir(table_dir); ensure_dir(figure_dir)

  cat("=== Cox-on-YF (Cluster B) Benchmark ===\n")
  cat(sprintf("  n=%d p=%d K_true=%d K_max=%d alpha=%.2f alpha_F=0\n",
              synthetic_n, synthetic_p, K_true, K_max, alpha))

  # ---- Synthetic benchmark ----
  data <- generate_synthetic_benchmark_data(
    n = synthetic_n, p = synthetic_p, K_true = K_true, seed = seed
  )
  cat(sprintf("  Censoring rate: %.1f%%\n\n", 100 * data$censoring_rate))

  split     <- stratified_split(data$status, test_frac = holdout_frac, seed = seed)
  train_idx <- split$train_idx
  test_idx  <- split$test_idx

  Y_train <- data$Y[train_idx, , drop = FALSE]
  Y_test  <- data$Y[test_idx,  , drop = FALSE]

  cat("Fitting Cluster B (Cox-on-YF) on synthetic training data...\n")
  fit_synth <- fit_cox_on_yf(
    Y_train, data$time[train_idx], data$status[train_idx],
    K        = K_max,
    alpha    = alpha,
    N_burnin = N_burnin,
    max_iter = max_iter,
    tol      = tol,
    prior_beta = prior_beta,
    verbose  = FALSE
  )
  pred_synth <- predict_cox_on_yf(Y_test, fit_synth$EF, fit_synth$EBeta)
  cindex_synth <- as.numeric(concordance(
    Surv(data$time[test_idx], data$status[test_idx]) ~ I(-pred_synth$risk_scores)
  )$concordance)

  # PCA baseline
  pca_train  <- prcomp(Y_train, rank. = min(5, K_max))
  pca_scores <- predict(pca_train, newdata = Y_train)
  pca_cox    <- coxph(Surv(data$time[train_idx], data$status[train_idx]) ~ pca_scores)
  pca_test   <- predict(pca_train, newdata = Y_test)
  pca_lp     <- as.vector(pca_test %*% coef(pca_cox))
  cindex_pca <- as.numeric(concordance(
    Surv(data$time[test_idx], data$status[test_idx]) ~ I(-pca_lp)
  )$concordance)

  cat(sprintf("  Synthetic holdout C-index: Cox-on-YF=%.3f  PCA=%.3f\n",
              cindex_synth, cindex_pca))
  cat(sprintf("  Converged: %s  Iterations: %d  K_eff(|beta|>0.001): %d\n",
              fit_synth$history$converged, fit_synth$history$n_iter,
              sum(abs(fit_synth$EBeta) > 0.001)))
  cat(sprintf("  EBeta range: [%.4e, %.4e]\n\n",
              min(fit_synth$EBeta), max(fit_synth$EBeta)))

  synth_summary <- data.frame(
    Method        = c("Cox-on-YF (Cluster B)", "PCA"),
    Holdout_CIndex = round(c(cindex_synth, cindex_pca), 4),
    Alpha         = c(alpha, NA_real_),
    K_max         = K_max,
    N_iter        = c(fit_synth$history$n_iter, NA_integer_),
    Converged     = c(fit_synth$history$converged, NA)
  )
  write.csv(synth_summary, file.path(table_dir, "synthetic_holdout_cindex.csv"),
            row.names = FALSE)

  # ---- ELBO trace ----
  elbo_df <- data.frame(
    Iteration  = seq_along(fit_synth$history$elbo_full),
    ELBO_Full  = fit_synth$history$elbo_full,
    ELBO_Proxy = fit_synth$history$elbo_proxy,
    RMSE       = fit_synth$history$rmse
  )
  write.csv(elbo_df, file.path(table_dir, "synthetic_elbo_trace.csv"), row.names = FALSE)

  # ---- ELBO plot ----
  pdf(file.path(figure_dir, "synthetic_elbo_trace.pdf"), width = 8, height = 5)
  par(mar = c(5, 5, 4, 2))
  plot(elbo_df$ELBO_Full, type = "l", lwd = 2, col = "#1F77B4",
       xlab = "Iteration", ylab = "ELBO",
       main = sprintf("Cox-on-YF ELBO Trace (synthetic n=%d, p=%d, alpha=%.2f)",
                      synthetic_n, synthetic_p, alpha), bty = "n")
  lines(elbo_df$ELBO_Proxy, col = "#2CA02C", lwd = 2, lty = 2)
  legend("bottomright", legend = c("Full ELBO", "Genomics proxy"),
         col = c("#1F77B4", "#2CA02C"), lwd = 2, lty = c(1, 2), bty = "n")
  dev.off()

  # ---- C-index bar chart ----
  pdf(file.path(figure_dir, "synthetic_holdout_cindex.pdf"), width = 6, height = 4.5)
  par(mar = c(4, 5, 4, 1))
  bp <- barplot(synth_summary$Holdout_CIndex,
                names.arg = synth_summary$Method,
                col = c("#FF7F0E", "#1F77B4"),
                ylim = c(0, 1), main = "Held-Out C-index: Cox-on-YF vs PCA",
                ylab = "C-index", las = 2, cex.names = 0.85, bty = "n")
  abline(h = 0.5, col = "#666666", lty = 2)
  text(bp, synth_summary$Holdout_CIndex + 0.03,
       sprintf("%.3f", synth_summary$Holdout_CIndex), cex = 0.9, font = 2)
  dev.off()

  results <- list(fit_synth = fit_synth, synth_summary = synth_summary)

  # ==============================================================================
  # Real PDAC Data Validation
  # ==============================================================================
  if (!run_real_data) return(invisible(results))

  if (!dir.exists(PDAC_DATA_ROOT)) {
    cat(sprintf("PDAC data root not found: %s\nSkipping real-data validation.\n", PDAC_DATA_ROOT))
    return(invisible(results))
  }

  cat("=== Real PDAC Data Validation ===\n")

  # ---- Load raw training cohorts (same pattern as run_cluster_a_smoke.R) ----
  # load_pdac_raw() sources load_data_internal.R from PDAC_DATA_ROOT.
  train_raw <- lapply(setNames(TRAINING_COHORTS, TRAINING_COHORTS), function(ds) {
    tryCatch(load_pdac_raw(ds, PDAC_DATA_ROOT),
             error = function(e) { cat(sprintf("  [skip] %s\n", ds)); NULL })
  })
  train_raw <- Filter(Negate(is.null), train_raw)

  if (length(train_raw) == 0) {
    cat("No training cohorts loaded. Skipping PDAC benchmark.\n")
    return(invisible(results))
  }

  cat(sprintf("  Loaded training cohorts: %s\n", paste(names(train_raw), collapse = " + ")))

  # Merge with v2 preprocessing (intersect → log2 → QN → top-2000 → rank)
  merged <- preprocess_merged_cohorts(
    cohort_raw_list    = train_raw,
    log_transform_flags = PLATFORM_LOG_TRANSFORM[names(train_raw)],
    top_n = 2000
  )

  # Survival labels for merged training set
  time_train   <- unlist(lapply(names(train_raw), function(ds) train_raw[[ds]]$time))
  status_train <- unlist(lapply(names(train_raw), function(ds) train_raw[[ds]]$status))

  cat(sprintf("  Training: n=%d p=%d (merged %s)\n",
              nrow(merged$Y), ncol(merged$Y),
              paste(names(train_raw), collapse = " + ")))

  # Fit Cluster B on merged training data
  cat("  Fitting Cox-on-YF on merged training cohorts...\n")
  fit_pdac <- fit_cox_on_yf(
    merged$Y, time_train, status_train,
    K        = K_max,
    alpha    = alpha,
    N_burnin = N_burnin,
    max_iter = max_iter,
    tol      = tol,
    prior_beta = prior_beta,
    verbose  = FALSE
  )
  cat(sprintf("  Converged: %s  Iterations: %d\n",
              fit_pdac$history$converged, fit_pdac$history$n_iter))
  cat(sprintf("  EBeta range: [%.4e, %.4e]  K_eff(|b|>0.001): %d\n\n",
              min(fit_pdac$EBeta), max(fit_pdac$EBeta),
              sum(abs(fit_pdac$EBeta) > 0.001)))

  # External validation on each cohort
  # Uses preprocess_desurv_cohort (single-cohort QN) + alignment to training gene set.
  # Same pattern as run_cluster_a_external.R eval_external().
  training_gene_names <- colnames(merged$Y)
  TOP_N <- 2000

  ext_results <- list()
  for (cohort in EXTERNAL_COHORTS) {
    res <- tryCatch({
      raw      <- load_pdac_raw(cohort, PDAC_DATA_ROOT)
      log_flag <- PLATFORM_LOG_TRANSFORM[[cohort]]
      if (is.null(log_flag)) log_flag <- FALSE
      preproc  <- preprocess_desurv_cohort(
        Y = raw$Y, gene_names = raw$gene_names,
        top_n = TOP_N, log_transform = log_flag,
        cohort_name = cohort
      )

      # Align to training gene order; fill missing genes with 0
      common_idx <- match(training_gene_names, preproc$gene_names)
      Y_ext <- matrix(0.0, nrow = raw$n, ncol = length(training_gene_names))
      present <- !is.na(common_idx)
      Y_ext[, present] <- preproc$Y[, common_idx[present], drop = FALSE]
      EF_aligned <- fit_pdac$EF  # already ordered by training_gene_names

      pred_ext <- predict_cox_on_yf(Y_ext, EF_aligned, fit_pdac$EBeta)

      cindex_ext <- if (sd(pred_ext$risk_scores) < 1e-10) NA_real_ else
        as.numeric(concordance(
          Surv(raw$time, raw$status) ~ I(-pred_ext$risk_scores)
        )$concordance)

      cat(sprintf("  %s: n=%d  C-index=%s  p_intersect=%d\n",
                  cohort, raw$n, format(round(cindex_ext, 3), nsmall=3), sum(present)))
      list(Cohort = cohort, n = raw$n, CIndex = round(cindex_ext, 4))
    }, error = function(e) {
      cat(sprintf("  [skip] %s: %s\n", cohort, conditionMessage(e)))
      NULL
    })
    if (!is.null(res)) ext_results[[cohort]] <- as.data.frame(res)
  }

  if (length(ext_results) > 0) {
    ext_df <- do.call(rbind, ext_results)
    rownames(ext_df) <- NULL
    write.csv(ext_df, file.path(table_dir, "pdac_external_cindex.csv"), row.names = FALSE)

    # Bar chart of external C-indices
    pdf(file.path(figure_dir, "pdac_external_cindex.pdf"), width = 8, height = 5)
    par(mar = c(7, 5, 4, 1))
    bp <- barplot(ext_df$CIndex, names.arg = ext_df$Cohort,
                  col = "#FF7F0E", ylim = c(0, 1),
                  main = "Cox-on-YF External C-index (PDAC cohorts)",
                  ylab = "C-index", las = 2, cex.names = 0.85, bty = "n")
    abline(h = 0.5, col = "#666666", lty = 2)
    text(bp, ext_df$CIndex + 0.03,
         sprintf("%.3f", ext_df$CIndex), cex = 0.85, font = 2)
    dev.off()

    cat("\nExternal C-index summary:\n")
    print(ext_df)
    results$ext_df <- ext_df
  }

  # ---- Training ELBO ----
  pdac_elbo_df <- data.frame(
    Iteration  = seq_along(fit_pdac$history$elbo_full),
    ELBO_Full  = fit_pdac$history$elbo_full,
    ELBO_Proxy = fit_pdac$history$elbo_proxy,
    RMSE       = fit_pdac$history$rmse
  )
  write.csv(pdac_elbo_df, file.path(table_dir, "pdac_elbo_trace.csv"), row.names = FALSE)

  cat("\n=== Cox-on-YF Benchmark Complete ===\n")
  cat(sprintf("  Output: %s\n", output_root))

  results$fit_pdac <- fit_pdac
  invisible(results)
}

# ==============================================================================
# Entry point: parse --quick flag and run
# ==============================================================================

args  <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

run_cox_on_yf_benchmark(quick = quick)
