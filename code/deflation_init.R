# ============================================================
# Script:  deflation_init.R
# Purpose: Deflation-style CAVI initialization -- seeds factor 1 from a
#          rank-1 SVD of Y, subtracts its contribution, takes a rank-1 SVD
#          of the residual for factor 2, and so on through K factors.
#          Alternative to a single batch K-factor SVD (the existing "svd"
#          init_method in fit_modular.R / fit_cox_on_yf.R), intended to fix
#          the CAVI factor-collapse failure mode documented in DECISIONS.md
#          2026-07-12/13: batch SVD-init columns from near-tied singular
#          values can start out nearly symmetric in amplitude, which the
#          joint CAVI updates can then drive to a shared degenerate
#          (near-zero) fixed point. Deflation avoids this by construction --
#          each factor is fit to the residual after removing prior factors,
#          so successive factors always target strictly less remaining
#          signal, an approach analogous to EBMF's own greedy,
#          factor-at-a-time fitting (ROADMAP.md 2026-07-12 entry on why
#          EBMF avoids this collapse mode).
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Dependencies: none
# ============================================================

#' Deflation-Style Rank-K SVD Initialization
#'
#' Sequentially extracts K rank-1 SVD factors from Y, deflating (subtracting)
#' each factor's contribution before extracting the next.
#'
#' Returns the SIGNED loadings/factors -- callers apply their own
#' non-negativity convention (`abs()` in fit_cox_on_yf.R, `pmax(.,0)` in
#' fit_modular.R) to match their model's point_exponential-prior init,
#' exactly as each caller already does for its "svd" init_method. This
#' function does not itself choose a convention, so it stays usable by both.
#'
#' @param Y n x p numeric matrix (already residualized against any cohort
#'          effects -- i.e. the same `Y_for_svd` used by "svd" init).
#' @param K integer: number of factors to extract (1 <= K <= min(nrow(Y), ncol(Y))).
#'
#' @return list(EL = n x K matrix, EF = p x K matrix) of SIGNED loadings/
#'   factors. Column k is the k-th rank-1 SVD of the residual after removing
#'   columns 1..(k-1)'s contribution (EL %*% t(EF) reconstructs the greedy
#'   K-term deflation approximation of Y, matching a clean rank-K Y exactly).
deflation_svd_init <- function(Y, K) {
  n <- nrow(Y); p <- ncol(Y)
  if (!is.numeric(K) || length(K) != 1 || K != round(K) ||
      K < 1 || K > min(n, p)) {
    stop(sprintf(
      "K must be a positive integer <= min(nrow(Y), ncol(Y)) = %d; got %s.",
      min(n, p), deparse(K)
    ))
  }

  EL <- matrix(0.0, n, K)
  EF <- matrix(0.0, p, K)
  Y_resid <- Y

  for (k in seq_len(K)) {
    sv  <- svd(Y_resid, nu = 1, nv = 1)
    d1  <- sqrt(max(sv$d[1], 0))
    l_k <- sv$u[, 1] * d1
    f_k <- sv$v[, 1] * d1
    EL[, k] <- l_k
    EF[, k] <- f_k
    Y_resid <- Y_resid - outer(l_k, f_k)
  }

  list(EL = EL, EF = EF)
}
