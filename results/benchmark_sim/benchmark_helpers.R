# ============================================================
# Script:       benchmark_helpers.R
# Purpose:      Shared constants and data-loading functions for
#               run_LB_benchmark.R and run_YFB_benchmark.R.
#               Must be sourced after cfg <- yaml::read_yaml("config/globals.yml").
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-05-05
# Dependencies: config/globals.yml (cfg must exist in calling env)
# ============================================================

# --------------------------------------------------------------------------
# Constants derived from config
# --------------------------------------------------------------------------

PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT",
                              unset = path.expand(cfg$pdac$data_root_default))

# Named logical vector: TRUE = RNA-seq (apply log2(x+1)), FALSE = already normalised
PLATFORM_LOG_TRANSFORM <- unlist(cfg$pdac$platform_log_transform)

EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

# --------------------------------------------------------------------------
# load_pdac_raw
# --------------------------------------------------------------------------
#' Load a single PDAC cohort from the local data directory.
#'
#' Calls the cohort's internal load_data_internal() via a temporary working
#' directory symlink, applies the keep-sample filter, transposes to n×p,
#' deduplicates gene names, and returns a named list.
#'
#' @param dataset_name  Character. One of the cohort keys in globals.yml pdac section.
#' @param pdac_root     Character. Path to the PDAC_data directory.
#' @return Named list: Y (n×p), gene_names, time, status, n, p, dataset_name.
load_pdac_raw <- function(dataset_name, pdac_root) {
  if (!dir.exists(pdac_root))
    stop(sprintf("PDAC data root not found: %s\nSet PDAC_DATA_ROOT env var.", pdac_root))

  tmp_wd    <- tempfile("pdac_wd_")
  dir.create(tmp_wd, showWarnings = FALSE)
  data_link <- file.path(tmp_wd, "data")
  file.symlink(pdac_root, data_link)

  old_wd <- getwd()
  on.exit({
    setwd(old_wd)
    unlink(data_link)
  }, add = TRUE)

  setwd(tmp_wd)
  source(file.path(pdac_root, "load_data_internal.R"), local = TRUE)
  result <- load_data_internal(dataset_name)
  setwd(old_wd)

  keeps <- which(result$sampInfo$keep == 1)
  if (length(keeps) == 0)
    stop(sprintf("No valid samples for dataset '%s' after keep filter.", dataset_name))

  # Transpose: genes × samples → patients × genes (n × p)
  Y <- t(result$ex[, keeps])

  fi <- result$featInfo
  if (is.data.frame(fi) && "SYMBOL" %in% names(fi)) {
    gene_names <- fi$SYMBOL
  } else if (is.character(fi)) {
    gene_names <- fi
  } else {
    gene_names <- rownames(result$ex)
  }
  if (length(gene_names) != ncol(Y))
    gene_names <- paste0("Gene", seq_len(ncol(Y)))

  # Keep first occurrence of duplicated gene names
  dup_mask <- duplicated(gene_names)
  if (any(dup_mask)) {
    keep_cols  <- !dup_mask
    Y          <- Y[, keep_cols, drop = FALSE]
    gene_names <- gene_names[keep_cols]
  }

  colnames(Y) <- gene_names
  rownames(Y) <- NULL

  time   <- result$sampInfo$time[keeps]
  status <- as.integer(result$sampInfo$event[keeps])

  list(
    Y = Y, gene_names = gene_names,
    time = time, status = status,
    n = nrow(Y), p = ncol(Y),
    dataset_name = dataset_name
  )
}

# --------------------------------------------------------------------------
# generate_synthetic_benchmark_data
# --------------------------------------------------------------------------
#' Generate a synthetic benchmark dataset for validating CAVI fits.
#'
#' DGP: Y = L F' + E (point-exponential L/F, Gaussian noise scaled by τ).
#' Survival times drawn from a Weibull PH model with η = Lβ_true.
#' Censoring calibrated to target_censoring via bisection.
#'
#' Reads beta_true from cfg$synthetic$b_true; cfg must exist in calling env.
#'
#' @param n                 Number of subjects.
#' @param p                 Number of genes.
#' @param K_true            Number of true latent factors.
#' @param seed              RNG seed.
#' @param target_censoring  Target censoring rate (fraction of censored obs).
#' @return Named list: Y, time, status, L_true, F_true, beta_true, tau_true,
#'   censoring_rate, shape, scale0, K_true, n, p.
generate_synthetic_benchmark_data <- function(n = 300, p = 1000, K_true = 5,
                                              seed = 222, target_censoring = 0.30) {
  set.seed(seed)

  beta_true <- as.numeric(cfg$synthetic$b_true)
  if (length(beta_true) != K_true) {
    stop(sprintf(
      "cfg$synthetic$b_true has length %d but K_true=%d. Update globals.yml.",
      length(beta_true), K_true
    ))
  }

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
  shape  <- 1.5
  scale0 <- 0.01
  event_times  <- (-log(runif(n)) / (scale0 * exp(eta_true)))^(1 / shape)
  censor_base  <- rexp(n, rate = 1)
  censor_scale <- .calibrate_censor_scale(event_times, censor_base,
                                          target = target_censoring)
  censor_times <- censor_base * censor_scale
  time   <- pmin(event_times, censor_times)
  status <- as.integer(event_times <= censor_times)

  list(
    Y = Y, time = time, status = status,
    L_true = L_true, F_true = F_true,
    beta_true = beta_true, tau_true = tau_true,
    censoring_rate = mean(status == 0),
    shape = shape, scale0 = scale0,
    K_true = K_true, n = n, p = p
  )
}

# Internal helper — not exported
.calibrate_censor_scale <- function(event_times, base_censor, target = 0.30, n_iter = 40) {
  lo <- 1e-3
  hi <- max(event_times) * 10
  for (i in seq_len(n_iter)) {
    mid <- sqrt(lo * hi)
    censor_rate <- mean(base_censor * mid < event_times)
    if (censor_rate > target) lo <- mid else hi <- mid
  }
  hi
}
