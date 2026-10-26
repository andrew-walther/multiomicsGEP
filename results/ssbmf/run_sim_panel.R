# ============================================================
# Script: run_sim_panel.R
# Purpose: Run the registered w = 1 simulation panel of results/ssbmf/sim_design.md
#          (reduced local version: 5 datasets per cell) and write one row per
#          fit with the metrics for O1-O4 and O7.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: ssbmf (loaded from ./ssbmf with pkgload), survival, parallel
# Usage: Rscript results/ssbmf/run_sim_panel.R [cores] [cells]   (from the repository root)
#        e.g. Rscript results/ssbmf/run_sim_panel.R 8 A,A2,B,C,C2,C3
# Memory: the two-block cells hold 200 x 12,000 matrices (about 20 MB each,
#         several copies per fit); 8 parallel fits use well under 8 GB.
# ============================================================

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")   # one BLAS thread per worker
load_ssbmf <- function() {
  lib <- Sys.getenv("SSBMF_LIB")
  if (nzchar(lib)) library(ssbmf, lib.loc = lib) else suppressMessages(pkgload::load_all("ssbmf", quiet = TRUE))
}
load_ssbmf()
library(survival)

args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 8L
cells_run <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("A", "A2", "B", "C", "C2", "C3")
out_dir <- "results/ssbmf/outputs/sim_panel"
dir.create(file.path(out_dir, "jobs"), recursive = TRUE, showWarnings = FALSE)

# Cells (sim_design.md) ----
cells <- list(
  A  = list(p = c(expression = 2000), noise = "lower", lowvar = FALSE),
  A2 = list(p = c(expression = 2000), noise = "higher", lowvar = FALSE),
  B  = list(p = c(expression = 2000, methylation = 10000), noise = "lower", lowvar = FALSE),
  C  = list(p = c(expression = 2000), noise = "lower", lowvar = TRUE),
  C2 = list(p = c(expression = 2000, methylation = 10000), noise = "lower", lowvar = TRUE),
  C3 = list(p = c(expression = 2000), noise = "lower", lowvar = TRUE, lowvar_pve = 0.0025)
)
K_true <- 4L
K_grid <- c(4L, 6L, 8L)
seeds <- 1:5
arms <- c("ssbmf_point_laplace", "ssbmf_point_exponential", "two_step")

# Metrics ----

#' Evaluate one fitted arm against the truth
#'
#' @param Fhat concatenated fitted loadings (features x K)
#' @param bhat,lower,upper coefficients and 90% interval per SD of fitted score
#' @param active logical, program survival-active
#' @param val_risk validation risk score
#' @param sim simulated data (truth)
#' @return one-row data.frame plus per-true-program columns
evaluate <- function(Fhat, bhat, lower, upper, active, val_risk, sim) {
  Ft <- do.call(rbind, sim$truth$F)
  r <- suppressWarnings(cor(Fhat, Ft)); r[!is.finite(r)] <- 0
  pve_t <- apply(sim$truth$pve, 1, max)
  out <- list()
  for (t in seq_len(ncol(Ft))) {
    k <- which.max(abs(r[, t])); s <- sign(r[k, t])
    b_al <- s * bhat[k]; lo <- if (s > 0) lower[k] else -upper[k]; hi <- if (s > 0) upper[k] else -lower[k]
    out[[t]] <- data.frame(true_program = t, pve_true = pve_t[t], beta_true = sim$truth$beta[t],
                           best_r = abs(r[k, t]), matched_k = k, beta_aligned = b_al,
                           lower_aligned = lo, upper_aligned = hi,
                           covered = lo <= sim$truth$beta[t] & sim$truth$beta[t] <= hi,
                           active = active[k], sign_ok = sign(b_al) == sign(sim$truth$beta[t]))
  }
  prog <- do.call(rbind, out)
  prog$val_C <- unname(survival::concordance(Surv(sim$validation$time, sim$validation$status) ~ val_risk,
                                             reverse = TRUE)$concordance)
  prog
}

fit_arm <- function(arm, dd, sim, K) {
  t0 <- proc.time()[["elapsed"]]
  if (arm == "two_step") {
    fit <- ssbmf_fit(dd, K, ssbmf_control(w_surv = 0, prior_F = "point_laplace"))
    keep <- fit$roles$role != "inactive"
    Z <- predict(fit, type = "scores")
    Zs <- sweep(Z, 2, ssbmf:::.canon_scale(fit), "/")      # per SD of score
    cx <- coxph(Surv(dd$time, dd$status) ~ Zs[, keep, drop = FALSE])
    b <- se <- numeric(K); b[keep] <- coef(cx); se[keep] <- sqrt(diag(vcov(cx)))
    b[!is.finite(b)] <- 0; se[!is.finite(se)] <- Inf
    pval <- 2 * pnorm(-abs(b / se))
    Zv <- sweep(predict(fit, newdata = sim$validation$X, type = "scores"), 2, ssbmf:::.canon_scale(fit), "/")
    val_risk <- drop(Zv %*% b)
    ev <- evaluate(do.call(rbind, lapply(fit$F, `[[`, "mean")), b, b - 1.645 * se, b + 1.645 * se,
                   keep & pval <= 0.10, val_risk, sim)
    share <- NA_real_
  } else {
    fam <- sub("ssbmf_", "", arm)
    fit <- ssbmf_fit(dd, K, ssbmf_control(prior_F = fam))
    cf <- coef(fit, level = 0.90)
    ev <- evaluate(do.call(rbind, lapply(fit$F, `[[`, "mean")), cf$mean, cf$lower, cf$upper,
                   fit$roles$role == "survival_active", predict(fit, newdata = sim$validation$X), sim)
    share <- fit$surv_share
  }
  ev$arm <- arm; ev$K_init <- K
  ev$converged <- fit$converged; ev$sweeps <- fit$sweeps
  ev$step_decreases <- fit$step_decreases; ev$worst_relative_step <- fit$worst_relative_step
  ev$rejections <- sum(fit$rejections); ev$halvings <- fit$halvings
  ev$n_survival_active <- sum(fit$roles$role == "survival_active")
  ev$n_reconstruction_only <- sum(fit$roles$role == "reconstruction_only")
  ev$surv_share_matched <- if (is.matrix(share)) vapply(ev$matched_k, function(k) max(share[k, ]), 0) else NA_real_
  ev$runtime_s <- proc.time()[["elapsed"]] - t0
  ev
}

# Jobs ----
jobs <- expand.grid(cell = cells_run, seed = seeds, K = K_grid, stringsAsFactors = FALSE)
jobs <- jobs[order(-nchar(jobs$cell), jobs$cell, jobs$seed, jobs$K), ]  # two-block cells first

run_job <- function(i) {
  jb <- jobs[i, ]
  f_out <- file.path(out_dir, "jobs", sprintf("%s_seed%d_K%d.rds", jb$cell, jb$seed, jb$K))
  if (file.exists(f_out)) return(readRDS(f_out))
  cl <- cells[[jb$cell]]
  sim <- ssbmf_simulate(n = 200, p = cl$p, K = K_true, low_variance_prognostic = cl$lowvar,
                        low_variance_pve = if (is.null(cl$lowvar_pve)) 0.01 else cl$lowvar_pve,
                        feature_correlation = "block", noise = cl$noise, seed = jb$seed)
  dd <- ssbmf_data(sim$X, sim$time, sim$status)
  res <- do.call(rbind, lapply(arms, function(a) fit_arm(a, dd, sim, jb$K)))   # errors fail the job
  res$cell <- jb$cell; res$seed <- jb$seed
  saveRDS(res, f_out)
  res
}

t_start <- Sys.time()
# Socket workers (forked workers crash in the macOS Accelerate BLAS)
cl <- parallel::makePSOCKcluster(cores)
parallel::clusterCall(cl, function(wd) { setwd(wd); NULL }, getwd())
parallel::clusterExport(cl, c("load_ssbmf", "cells", "K_true", "arms", "jobs", "out_dir", "evaluate", "fit_arm", "run_job"))
parallel::clusterEvalQ(cl, { load_ssbmf(); library(survival); NULL })
rows <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(i) try(run_job(i), silent = TRUE))
parallel::stopCluster(cl)
failed <- vapply(rows, function(x) inherits(x, "try-error"), TRUE)
if (any(failed)) message(paste(unlist(rows[failed]), collapse = "\n"))
if (any(failed)) message("Jobs that failed: ", paste(which(failed), collapse = ", "))
all_rows <- do.call(rbind, rows[!failed])
write.csv(all_rows, file.path(out_dir, "sim_panel_rows.csv"), row.names = FALSE)
message(sprintf("Panel done: %d jobs, %d rows, %.1f min", nrow(jobs), nrow(all_rows),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
