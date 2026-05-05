# =============================================================================
# code/update_L_surv_YFB.R
#
# q(L) update for the Cox-on-YF model (Cluster B: eta = YFB)
#
# MODEL DIFFERENCE: Under eta = YFB, L appears ONLY in the genomics likelihood
# (Y ≈ LF'). L does NOT appear in the Cox likelihood — that role is now played
# by F via ZF = Y·F. Consequently, this is a PURE-GENOMICS update (no survival
# terms). Compare to update_L.R (Cluster A), which is DUAL-SOURCE (L appears
# in both genomics and Cox likelihoods).
#
# DERIVATION REFERENCE:
#   derivations/cox_on_YF/qL_unsupervised_derivation.tex
#
# MATHEMATICAL SUMMARY:
#   A_L = sum_j(tau_j * E[f^2_{jk}])   [scalar, same for all i]
#   B_L[i] = (R_k %*% (Tau * EF_k))[i] [n-vector]
#   EBNM(x = B_L/A_L, s = 1/sqrt(A_L)) <- VECTOR call (n observations)
#
# DESIGN PATTERN:
#   Part of the Cluster B modular update series:
#     update_beta.R  ->  update_L_surv_YFB.R  ->  update_F_surv_YFB.R  ->  update_tau.R
#   Each exposes a _k() (single factor) and _all() (full loop) function.
#
# NOTE: compute_R_k() is defined in update_L.R (Cluster A) and re-exported from
#       fit_cox_on_yf.R via source("code/update_L.R"). If update_L.R is not
#       sourced, compute_R_k() must be available from another source.
#
# DEPENDENCIES:
#   ebnm  (CRAN)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(l_{ik}) for a single factor k  [Cluster B: pure-genomics EBNM]
#'
#' Under the Cox-on-YF reformulation (eta = Z_F * beta_tilde), L appears only
#' in the genomics likelihood. The update is a VECTOR EBNM with n observations
#' driven entirely by the genomics residual R_k. No survival arguments.
#'
#' @param Tau         p-vector: feature-specific noise precision tau_j (> 0)
#' @param EF_k        p-vector: posterior mean of factor column k
#' @param EF2_k       p-vector: posterior SECOND MOMENT of factor column k
#'                    (E_q[f_{jk}^2] = Var + mean^2). Must satisfy EF2_k >= EF_k^2.
#' @param R_k         n x p matrix: partial residual R^{-k} (from compute_R_k)
#' @param prior_family character: EBNM prior family (default "point_exponential")
#' @param A_floor     numeric: minimum value for A_L (default 1e-10)
#'
#' @return Named list:
#'   $mean        -- n-vector: posterior mean E_q[l_{ik}]
#'   $second      -- n-vector: posterior 2nd moment E_q[l_{ik}^2]
#'   $sd          -- n-vector: posterior SD sqrt(Var_q(l_{ik}))
#'   $A           -- scalar: precision A_L (floored)
#'   $B           -- n-vector: signal B_L[i]
#'   $x           -- n-vector: EBNM pseudo-obs x_i = B_i/A
#'   $s           -- scalar: EBNM pseudo-noise s = 1/sqrt(A)
#'   $ebnm_result -- raw ebnm() return object
#'
#' @export
#' @family L_surv_YFB_update
#' @seealso \code{\link{update_F_surv_YFB_k}} for the dual-source F update
update_L_surv_YFB_k <- function(Tau, EF_k, EF2_k,
                                 R_k,
                                 prior_family = "point_exponential",
                                 A_floor      = 1e-10) {

  # ------------------------------------------------------------------
  # Precision A_L  (scalar — same for all subjects i)
  #
  #   A_L = sum_j(tau_j * E[f^2_{jk}])
  #
  # Uses E[f^2] (second moment) for error-in-variables correction.
  # ------------------------------------------------------------------
  A_gen <- sum(Tau * EF2_k)                                     # scalar
  A_L   <- max(A_gen, A_floor)                                  # scalar [floored]

  # ------------------------------------------------------------------
  # Signal B_L[i]  (n-vector)
  #
  #   B_L[i] = sum_j(tau_j * R^{-k}_{ij} * f_bar_{jk})
  #           = (R_k %*% (Tau * EF_k))[i]
  # ------------------------------------------------------------------
  B_L <- as.vector(R_k %*% (Tau * EF_k))                       # n-vector

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise
  # ------------------------------------------------------------------
  x_L <- B_L / A_L
  s_L <- 1.0 / sqrt(A_L)
  x_L[!is.finite(x_L)] <- 0

  res <- ebnm(x = x_L, s = s_L, prior_family = prior_family)

  l_mean   <- res$posterior$mean
  l_sd     <- res$posterior$sd
  l_second <- l_sd^2 + l_mean^2

  list(
    mean        = l_mean,
    second      = l_second,
    sd          = l_sd,
    A           = A_L,
    B           = B_L,
    x           = x_L,
    s           = s_L,
    ebnm_result = res
  )
}

# =============================================================================
# Convenience wrapper: full Gauss-Seidel loop over all K factors
# =============================================================================

#' Update q(L) for all K factors  [Cluster B: pure-genomics]
#'
#' Iterates over k = 1..K, computing the partial residual R_k and calling
#' update_L_surv_YFB_k() for each factor. Uses Gauss-Seidel ordering.
#'
#' Under Cluster B, no survival arguments are needed (w, z, EBeta, EBeta2
#' do not enter q(L); they flow through q(F) instead).
#'
#' @param Y           n x p matrix: observed genomics data
#' @param EL          n x K matrix: posterior means of loadings (will be updated)
#' @param EL2         n x K matrix: posterior second moments of loadings
#' @param EF          p x K matrix: posterior means of factors
#' @param EF2         p x K matrix: posterior second moments of factors
#' @param Tau         p-vector: feature-specific noise precision
#' @param prior_family character: EBNM prior family (default "point_exponential")
#' @param A_floor     numeric: precision floor (default 1e-10)
#'
#' @return Named list:
#'   $EL      -- n x K matrix of updated posterior means
#'   $EL2     -- n x K matrix of updated posterior second moments
#'   $details -- length-K list, each element is the full update_L_surv_YFB_k result
#' @export
#' @family L_surv_YFB_update
update_L_surv_YFB_all <- function(Y, EL, EL2, EF, EF2, Tau,
                                   prior_family = "point_exponential",
                                   A_floor      = 1e-10) {

  K        <- ncol(EL)
  EL_curr  <- EL
  EL2_curr <- EL2
  details  <- vector("list", K)

  for (k in seq_len(K)) {
    R_k <- compute_R_k(Y, EL_curr, EF, k)

    res_k <- update_L_surv_YFB_k(
      Tau          = Tau,
      EF_k         = EF[, k],
      EF2_k        = EF2[, k],
      R_k          = R_k,
      prior_family = prior_family,
      A_floor      = A_floor
    )

    EL_curr[, k]  <- res_k$mean
    EL2_curr[, k] <- res_k$second
    details[[k]]  <- res_k
  }

  list(
    EL      = EL_curr,
    EL2     = EL2_curr,
    details = details
  )
}
