# ============================================================
# Script: update_F_cohort.R
# Purpose: Closed-form Normal conjugate update for cohort indicator
#          F rows in the augmented SSBMF model (cohort-cols-L branch).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-22
# Dependencies: none (pure R arithmetic)
# ============================================================

# ==============================================================================
# Mathematical background
# ==============================================================================
#
# In the augmented model, L is extended with C-1 fixed binary cohort indicator
# columns (corner-point encoding).  The corresponding F rows f_c are estimated
# from the genomics likelihood alone — cohort columns never enter the Cox term,
# so beta_cohort = 0 by construction.
#
# For a single binary cohort column c (n-vector, {0,1}):
#
#   Likelihood (patients in cohort B, i.e. c_i = 1):
#     R_ij | f_cj ~ N(f_cj, tau_j^{-1})     for each gene j = 1..p
#
#   where R_ij = Y_ij - sum_{k=1}^{K} E[l_ik] * E[f_kj]
#   is the global-factor residual (Y minus the K biological factors).
#
#   Prior:  f_cj ~ N(0, sigma_F_cohort^2)
#
#   Posterior precision:  A_cj = n_c * tau_j + 1 / sigma_F_cohort^2   [p-vector]
#   Posterior mean:       E[f_cj] = B_cj / A_cj                       [p-vector]
#   where               B_cj = tau_j * sum_{i: c_i=1} R_ij            [p-vector]
#
#   Posterior second moment:  E[f_cj^2] = 1/A_cj + (B_cj/A_cj)^2    [p-vector]
#
# NOTE: c_i is OBSERVED, not a latent variable — only f_c is estimated.
#       EL2 for the cohort column equals L_cohort (binary: 0^2=0, 1^2=1).
#       EF2_cohort = 1/A + EF_c^2 (NEVER just EF_c^2 — that zeros posterior variance).
# ==============================================================================

# ==============================================================================
# Single-column update
# ==============================================================================

#' Update q(f_c) for a single cohort indicator column
#'
#' Closed-form Normal conjugate update for the platform loading vector f_c.
#' Unlike the biological factor updates (which use EBNM), this is exact:
#' no numerical optimisation is required.
#'
#' @param c_vec          n-vector, binary {0,1}: patients in this cohort (c_i = 1)
#' @param R_cohort       n x p matrix: global-factor residual Y - EL[,1:K] %*% t(EF[,1:K])
#' @param Tau            p-vector: noise precisions tau_j (> 0)
#' @param sigma_F_cohort scalar: prior SD (default 1.0).  sigma^2 = 1.0 is calibrated
#'                       for column-centred unit-variance data; run the pre-fitting
#'                       diagnostic to confirm median per-gene cohort offset < 2*sigma.
#'
#' @return Named list:
#'   $EF_c  -- p-vector: posterior means  E[f_cj] = B_cj / A_cj
#'   $EF2_c -- p-vector: posterior second moments  1/A_cj + (B_cj/A_cj)^2
#'   $A_c   -- p-vector: posterior precisions A_cj (needed for compute_var_term and ELBO)
#'
#' @seealso update_F_cohort_all, compute_normal_kl, code/update_F.R
update_F_cohort_col <- function(c_vec, R_cohort, Tau, sigma_F_cohort = 1.0) {
  n_c   <- sum(c_vec)                                  # scalar: patients in cohort B
  B_c   <- colSums(R_cohort * c_vec) * Tau             # p-vector: weighted residual sum
  A_c   <- n_c * Tau + 1.0 / sigma_F_cohort^2         # p-vector: posterior precision
  EF_c  <- B_c / A_c                                   # p-vector: posterior mean
  EF2_c <- 1.0 / A_c + EF_c^2                         # p-vector: posterior second moment
  list(EF_c = EF_c, EF2_c = EF2_c, A_c = A_c)
}

# ==============================================================================
# Multi-column wrapper
# ==============================================================================

#' Update q(f_c) for all C-1 cohort indicator columns
#'
#' Loops over columns of L_cohort, calling update_F_cohort_col for each, and
#' assembles results into p x (C-1) matrices.
#'
#' For C = 2 (one cohort column, e.g. TCGA + CPTAC) this calls
#' update_F_cohort_col exactly once.  For C > 2, this is a Gauss-Seidel
#' approximation: each column's update ignores current-iteration estimates
#' of the other cohort columns.  This converges to the correct fixed point but
#' more slowly.  The approximation is exact for C = 2.
#'
#' @param L_cohort       n x (C-1) matrix: fixed binary indicator columns from
#'                       model.matrix(~ factor(cohort_id))[, -1, drop = FALSE]
#' @param R_cohort       n x p matrix: global-factor residual passed to each column update
#' @param Tau            p-vector: noise precisions tau_j (> 0)
#' @param sigma_F_cohort scalar: prior SD (default 1.0)
#'
#' @return Named list:
#'   $EF_cohort  -- p x (C-1) matrix: posterior means
#'   $EF2_cohort -- p x (C-1) matrix: posterior second moments
#'   $A_cohort   -- p x (C-1) matrix: posterior precisions
#'
#' @seealso update_F_cohort_col, compute_normal_kl
update_F_cohort_all <- function(L_cohort, R_cohort, Tau, sigma_F_cohort = 1.0) {
  C_cols      <- ncol(L_cohort)
  p           <- ncol(R_cohort)
  EF_cohort   <- matrix(0.0, p, C_cols)
  EF2_cohort  <- matrix(0.0, p, C_cols)
  A_cohort    <- matrix(0.0, p, C_cols)

  for (c in seq_len(C_cols)) {
    res <- update_F_cohort_col(L_cohort[, c], R_cohort, Tau, sigma_F_cohort)
    EF_cohort[,  c] <- res$EF_c
    EF2_cohort[, c] <- res$EF2_c
    A_cohort[,   c] <- res$A_c
  }

  list(EF_cohort = EF_cohort, EF2_cohort = EF2_cohort, A_cohort = A_cohort)
}
