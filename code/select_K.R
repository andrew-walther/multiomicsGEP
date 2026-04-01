# ============================================================
# Script: select_K.R
# Purpose: Data-driven selection of the number of latent factors K
#          for the SBMF model.
#          Option A: fit large K and count active (non-pruned) factors.
#          Option B: cross-validated C-index (stub for Longleaf HPC).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-01
# Dependencies: code/fit_modular.R (must be sourced first)
# ============================================================

# ============================================================
# compute_pve() ----
# ============================================================

#' Compute per-factor proportion of variance explained (PVE).
#'
#' PVE_k = ||EL[,k] %*% t(EF[,k])||_F^2 / ||Y||_F^2
#'
#' Note: this uses the posterior means EL and EF, not the full
#' posterior (which would require second moments and a larger calculation).
#' It is a useful diagnostic but not identical to the model ELBO.
#'
#' @param res     list with $EL (n × K), $EF (p × K) — from fit_supervised_mf_modular()
#' @param Y       numeric matrix (n × p) — the original data
#'
#' @return numeric vector of length K: PVE per factor (values in [0, 1])
compute_pve <- function(res, Y) {
  K         <- ncol(res$EL)
  total_var <- sum(Y^2)
  if (total_var == 0) return(rep(0, K))
  sapply(seq_len(K), function(k) {
    sum(outer(res$EL[, k], res$EF[, k])^2) / total_var
  })
}

# ============================================================
# auto_prune_K() ----
# ============================================================

#' Select K by fitting a large model and counting active factors.
#'
#' Fits the SBMF model with K = K_max latent factors.  The point-normal
#' prior drives uninformative factors toward exactly zero, producing
#' sparse β and F columns.  An "active" factor is one where either:
#'   - |β_k| > beta_thresh (the factor has prognostic signal), OR
#'   - PVE_k > pve_thresh (the factor explains >1% of genomic variance)
#'
#' The count of active factors is the recommended effective K.  This is
#' analogous to the EBMF/flashr approach to rank selection.
#'
#' **Computational cost:** One fit at K_max instead of one fit per K
#' candidate — much cheaper than cross-validation.
#'
#' @param Y          numeric matrix (n × p)
#' @param time       numeric vector (n)
#' @param status     integer vector (n)
#' @param K_max      integer: maximum K to fit (default 10)
#' @param beta_thresh  numeric: |β_k| threshold for "active" (default 0.05)
#' @param pve_thresh   numeric: PVE_k threshold for "active" (default 0.01 = 1%)
#' @param ...        additional arguments passed to fit_supervised_mf_modular()
#'
#' @return Named list:
#'   $K_effective   integer: number of active factors
#'   $pve           numeric vector (K_max): PVE per factor
#'   $beta          numeric vector (K_max): |EBeta| per factor
#'   $active        logical vector (K_max): TRUE if factor is active
#'   $fit           full fit object from fit_supervised_mf_modular()
auto_prune_K <- function(Y, time, status, K_max = 10,
                          beta_thresh = 0.05, pve_thresh = 0.01, ...) {

  cat(sprintf("  [auto_prune_K] Fitting K=%d to identify active factors...\n", K_max))

  fit <- fit_supervised_mf_modular(Y, time, status, K = K_max, ...)

  pve     <- compute_pve(fit, Y)
  ab_beta <- abs(fit$EBeta)

  # A factor is "active" if it either has prognostic signal (|β| > thresh)
  # or explains non-trivial genomic variance (PVE > 1%)
  active <- (ab_beta > beta_thresh) | (pve > pve_thresh)

  K_eff <- sum(active)
  cat(sprintf("  [auto_prune_K] K_effective = %d / %d (beta_thresh=%.2f, pve_thresh=%.3f)\n",
              K_eff, K_max, beta_thresh, pve_thresh))
  cat(sprintf("  [auto_prune_K] |beta|: [%s]\n",
              paste(sprintf("%.3f", ab_beta), collapse = ", ")))
  cat(sprintf("  [auto_prune_K] PVE%%:   [%s]\n",
              paste(sprintf("%.1f%%", pve * 100), collapse = ", ")))

  list(
    K_effective = K_eff,
    pve         = pve,
    beta        = ab_beta,
    active      = active,
    fit         = fit
  )
}

# ============================================================
# select_K_cv() stub ----
# ============================================================

#' K selection via cross-validated held-out C-index (STUB — Longleaf HPC only).
#'
#' For each candidate K in K_grid:
#'   1. Do n_folds-fold CV: fit on (n_folds-1)/n_folds of data
#'   2. Evaluate C-index on the held-out fold via predict_supervised_mf()
#'   3. Select K that maximises mean held-out C-index
#'
#' This is the most principled K selection approach but is compute-intensive:
#'   n_folds × |K_grid| model fits, each up to max_iter=300 iterations.
#' Intended for Longleaf HPC with SLURM parallelisation.
#'
#' **Status:** Interface stub only.  Call auto_prune_K() for immediate use.
#'
#' @param Y       numeric matrix (n × p)
#' @param time    numeric vector (n)
#' @param status  integer vector (n)
#' @param K_grid  integer vector: candidate K values (default c(2, 3, 5, 8, 10))
#' @param n_folds integer: number of cross-validation folds (default 5)
#' @param ...     additional arguments passed to fit_supervised_mf_modular()
#'
#' @return (not implemented — errors with informative message)
select_K_cv <- function(Y, time, status,
                         K_grid  = c(2, 3, 5, 8, 10),
                         n_folds = 5, ...) {
  stop(paste0(
    "select_K_cv() is not yet implemented.\n",
    "This function is intended for Longleaf HPC with SLURM parallelisation.\n",
    "Use k_select = 'auto_prune' for immediate local use.\n",
    "See project_future_direction.md for K selection roadmap."
  ))
}
