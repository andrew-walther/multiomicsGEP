# =============================================================================
# code/update_F.R
#
# Modular q(f_{jk}) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE:
#   derivations/qF/qF_update_derivation.tex  (self-contained)
#   derivations/MF_UpdateDerivations/MF_V2_Companion.tex  Section 5
#   derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.tex
#
# MATHEMATICAL SUMMARY:
#   f_{jk} appears ONLY in the genomics likelihood (not survival).  The
#   coordinate-ascent update reduces to a p-dimensional EBNM problem:
#
#     A_F[j] = tau_j * sum_i E[l^2_{ik}]        [p-vector]
#     B_F[j] = tau_j * sum_i R^{-k}_{ij} * l_bar_{ik}   [p-vector]
#
#     EBNM(x = B_F/A_F, s = 1/sqrt(A_F))  <- VECTOR call (p observations)
#
#   IMPORTANT PROPERTY: tau_j CANCELS in x_j = B_F[j]/A_F[j] because both
#   A_F and B_F are proportional to tau_j.  However, tau_j does NOT cancel
#   in s_j = 1/sqrt(A_F[j]): features with higher precision get tighter
#   pseudo-noise scales and hence less shrinkage.
#
# DESIGN PATTERN:
#   Part of the modular update series:
#     update_beta.R  ->  update_L.R  ->  update_F.R  ->  update_tau.R
#
# DEPENDENCIES:
#   ebnm  (CRAN)
#   compute_R_k() from code/update_L.R (used by update_F_all)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(f_{jk}) for a single factor k
#'
#' Performs the EBNM-based coordinate-ascent update for the biological
#' factor column f_{.k} using ONLY the genomics likelihood (F does not
#' appear in the Cox survival term).
#'
#' The function takes pre-computed inputs and returns a self-contained
#' result list, making it independently testable.
#'
#' KEY PROPERTY: The EBNM pseudo-observation x_j = B_F[j]/A_F[j] is
#' INDEPENDENT of tau_j (the feature-specific precision cancels).  Only
#' the pseudo-noise s_j = 1/sqrt(A_F[j]) depends on tau_j: higher
#' precision features experience less shrinkage.
#'
#' @param Tau         p-vector: feature-specific noise precision tau_j (> 0)
#' @param EL_k        n-vector: posterior mean of loading column k (l_bar_{ik}).
#'                    NOTE: This should be the UPDATED EL[,k] from the L step
#'                    (Gauss-Seidel ordering).
#' @param EL2_k       n-vector: posterior SECOND MOMENT of loading column k
#'                    (E_q[l_{ik}^2] = Var + mean^2).  Must satisfy EL2_k >= EL_k^2.
#' @param R_k         n x p matrix: partial residual R^{-k} (from compute_R_k).
#'                    NOTE: R_k is still valid after the L update because it
#'                    excludes factor k — only EL[,k] changed, not EL[,k'].
#' @param prior_family character: EBNM prior family (default "point_normal")
#' @param A_floor     numeric: minimum value for each A_F[j] (default 1e-10)
#'
#' @return Named list:
#'   $mean        -- p-vector: posterior mean E_q[f_{jk}]
#'   $second      -- p-vector: posterior 2nd moment E_q[f_{jk}^2]
#'   $sd          -- p-vector: posterior SD sqrt(Var_q(f_{jk}))
#'   $A           -- p-vector: precision A_{jk} [floored]
#'   $B           -- p-vector: signal B_{jk}
#'   $x           -- p-vector: EBNM pseudo-obs x_j = B_j/A_j (tau-free)
#'   $s           -- p-vector: EBNM pseudo-noise s_j = 1/sqrt(A_j) (tau-dependent)
#'   $sum_EL2_k   -- scalar: sum_i E[l_{ik}^2] (diagnostic)
#'   $ebnm_result -- raw ebnm() return object
#'
#' @examples
#' library(ebnm)
#' set.seed(1); n <- 100; p <- 200
#' Tau    <- rep(2.0, p)
#' EL_k   <- rnorm(n);  EL2_k <- EL_k^2 + 0.1
#' R_k    <- matrix(rnorm(n * p), n, p)
#' res <- update_F_k(Tau, EL_k, EL2_k, R_k)
#' cat("Factor estimates (first 5):", round(res$mean[1:5], 3), "\n")
update_F_k <- function(Tau, EL_k, EL2_k, R_k,
                        prior_family = "point_normal",
                        A_floor      = 1e-10) {

  # ------------------------------------------------------------------
  # Precision A_{jk}  (p-vector)
  #
  #   A_F[j] = tau_j * sum_i E[l_{ik}^2]
  #
  # sum_i(EL2_k) is a SCALAR (total loading second moment across samples).
  # The error-in-variables correction enters through EL2_k (uses full
  # second moment, not squared mean): posterior uncertainty in L
  # appropriately inflates the effective noise for F_k.
  # ------------------------------------------------------------------
  sum_EL2_k <- sum(EL2_k)                             # scalar
  A_F       <- pmax(Tau * sum_EL2_k, A_floor)         # p-vector [A3]

  # ------------------------------------------------------------------
  # Signal B_{jk}  (p-vector)
  #
  #   B_F[j] = tau_j * sum_i R^{-k}_{ij} * l_bar_{ik}
  #          = tau_j * (t(R_k) %*% EL_k)[j]
  #
  # The efficient form t(R_k) %*% EL_k gives a p-vector directly.
  # ------------------------------------------------------------------
  B_F <- Tau * as.vector(t(R_k) %*% EL_k)             # p-vector

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise (p-vectors)
  #
  #   x_j = B_F[j] / A_F[j]
  #       = (tau_j * (t(R_k) %*% EL_k)[j]) / (tau_j * sum_EL2_k)
  #       = (t(R_k) %*% EL_k)[j] / sum_EL2_k
  #   -> tau_j CANCELS in x_j!
  #
  #   s_j = 1 / sqrt(A_F[j]) = 1 / sqrt(tau_j * sum_EL2_k)
  #   -> tau_j does NOT cancel in s_j
  # ------------------------------------------------------------------
  x_F <- B_F / A_F
  s_F <- 1.0 / sqrt(A_F)

  # ------------------------------------------------------------------
  # Solve the p-dimensional EBNM problem
  # Each feature j provides one observation x_j with noise s_j.
  # ------------------------------------------------------------------
  res <- ebnm(x = x_F, s = s_F, prior_family = prior_family)

  f_mean   <- res$posterior$mean                       # p-vector
  f_sd     <- res$posterior$sd                         # p-vector
  f_second <- f_sd^2 + f_mean^2                       # p-vector

  list(
    mean        = f_mean,
    second      = f_second,
    sd          = f_sd,
    A           = A_F,
    B           = B_F,
    x           = x_F,
    s           = s_F,
    sum_EL2_k   = sum_EL2_k,
    ebnm_result = res
  )
}

# =============================================================================
# Convenience wrapper: full Gauss-Seidel loop over all K factors
# =============================================================================

#' Update q(F) for all K factors (Gauss-Seidel CAVI loop)
#'
#' Iterates over k = 1..K, computing R_k and calling update_F_k() for each
#' factor.  Uses Gauss-Seidel ordering: once EF[,k] is updated, the new
#' values are used for subsequent factors k' > k.
#'
#' IMPORTANT: This wrapper assumes that EL and EL2 have already been updated
#' for the current iteration (e.g., by update_L_all).  In the full CAVI loop,
#' the F update for factor k uses the UPDATED EL[,k] from the L step.
#'
#' REQUIRES: compute_R_k() from code/update_L.R
#'
#' @param Y           n x p matrix: observed genomics data
#' @param EL          n x K matrix: posterior means of loadings (ALREADY UPDATED)
#' @param EL2         n x K matrix: posterior second moments of loadings
#' @param EF          p x K matrix: posterior means of factors (will be updated)
#' @param EF2         p x K matrix: posterior second moments of factors
#' @param Tau         p-vector: feature-specific noise precision
#' @param prior_family character: EBNM prior family (default "point_normal")
#' @param A_floor     numeric: precision floor (default 1e-10)
#'
#' @return Named list:
#'   $EF      -- p x K matrix of updated posterior means
#'   $EF2     -- p x K matrix of updated posterior second moments
#'   $details -- length-K list, each element is the full update_F_k result
update_F_all <- function(Y, EL, EL2, EF, EF2, Tau,
                          prior_family = "point_normal",
                          A_floor      = 1e-10) {

  p <- nrow(EF)
  K <- ncol(EF)
  EF_curr  <- EF                     # mutable copy for Gauss-Seidel
  EF2_curr <- EF2
  details  <- vector("list", K)

  for (k in seq_len(K)) {
    # Partial residual R_k (uses current EL and Gauss-Seidel EF)
    R_k <- compute_R_k(Y, EL, EF_curr, k)

    res_k <- update_F_k(
      Tau        = Tau,
      EL_k       = EL[, k],
      EL2_k      = EL2[, k],
      R_k        = R_k,
      prior_family = prior_family,
      A_floor    = A_floor
    )

    # Gauss-Seidel: update EF columns so next k uses fresh values
    EF_curr[, k]  <- res_k$mean
    EF2_curr[, k] <- res_k$second
    details[[k]]  <- res_k
  }

  list(
    EF      = EF_curr,
    EF2     = EF2_curr,
    details = details
  )
}
