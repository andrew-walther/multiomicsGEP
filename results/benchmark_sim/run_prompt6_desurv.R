# ============================================================
# Script: run_prompt6_desurv.R
# Purpose: Fit the author DeSurv model on the fixed S-11 training patients.
# Created: 2026-11-06
# Dependencies: DeSurv, dplyr, parallel, yaml
# ============================================================
# One fit per fork worker, with one BLAS thread each. Checkpoints containing
# training fits stay under ignored data/. Workers share the parent matrix until modified.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
if (dir.exists("data/round2_R_library"))
  .libPaths(c(normalizePath("data/round2_R_library"), .libPaths()))
stopifnot(identical(packageDescription("DeSurv")$RemoteSha,
                    "cc9d4c2486b3125704ab71eea61f35e040a40080"))
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("results/benchmark_sim/desurv_author_adapter.R")
source("results/benchmark_sim/run_desurv_h2h.R")

#' Load an author cohort with its aligned eligible sample metadata.
#' @param cohort Key from the benchmark configuration.
#' @return Original loader list containing ex, featInfo and sampInfo.
#' @examples
#' # prompt6_author_data("CPTAC")
prompt6_author_data <- function(cohort) {
  root <- normalizePath(PDAC_DATA_ROOT, mustWork = TRUE)
  loader <- file.path(root, "load_data_internal.R")
  stopifnot(unname(tools::md5sum(loader)) == "3792de03e1649c4fbc1987770905aab9")
  tmp <- tempfile()
  dir.create(tmp)
  old <- getwd()
  on.exit({setwd(old); unlink(tmp, recursive = TRUE)})
  stopifnot(file.symlink(root, file.path(tmp, "data")))
  setwd(tmp)
  source(loader, local = TRUE)
  load_data_internal(cohort)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  workers <- if (length(args)) as.integer(args[1]) else 2L
  training <- desurv_author_prepare(lapply(setNames(cfg$pdac$training_cohorts,
                                                    cfg$pdac$training_cohorts),
                                          prompt6_author_data))
  s11 <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  stopifnot(identical(as.numeric(training$sampInfo$time), as.numeric(s11$time)),
            identical(as.numeric(training$sampInfo$event), as.numeric(s11$status)))
  work <- "data/prompt6_desurv"
  dir.create(work, showWarnings = FALSE)
  seeds <- 1:100
  stems <- file.path(work, sprintf("seed_%03d", seeds))
  pending <- seeds[!file.exists(paste0(stems, "_status.csv"))]
  if (length(pending)) {
    parallel::mclapply(pending, desurv_seed_fit, training = training, out_dir = work,
                       mc.cores = workers, mc.preschedule = FALSE)
  }
  status <- do.call(rbind, lapply(paste0(stems, "_status.csv"), read.csv))
  stopifnot(identical(status$seed, seeds))
  keep <- is.na(status$error)
  stopifnot(any(keep))
  fits <- lapply(paste0(stems[keep], ".rds"), readRDS)
  init <- DeSurv::desurv_consensus_seed(fits, training$ex, ntop = 270L, k = 3L,
                                      min_frequency = ceiling(.3 * sum(keep)))
  p <- desurv_reported_controls(NULL)
  final <- DeSurv::desurv_fit(training$ex, s11$time, s11$status, k = p$k,
    alpha = p$alpha, lambda = p$lambda, nu = p$nu, lambdaW = 0, lambdaH = 0,
    W0 = init$W0, H0 = init$H0, beta0 = init$beta0, seed = NULL,
    tol = 1e-7, maxit = 4000L, verbose = FALSE)
  stopifnot(!isTRUE(final$nan_flag), all(is.finite(unlist(final$loss))))
  saveRDS(final, file.path(work, "final.rds"))
}
