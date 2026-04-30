# =============================================================================
# code/update_L.R
#
# Modular q(l_{ik}) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE (Cluster B):
#   derivations/cox_on_YF/qL_unsupervised_derivation.tex  (D1)
#   derivations/Modular_UpdateDerivations/  q(L) section (original, for comparison)
#
# MATHEMATICAL SUMMARY (Cluster B — Cox-on-YF reformulation):
#   Under eta_i = (y_i F) beta_tilde, the loading l_{ik} appears ONLY in the
#   genomics likelihood (Y ≈ LF').  The survival linear predictor uses the
#   observed projection Z_F = Y E[F], so L drops from the Cox term entirely.
#
#   The update reduces to a PURE-GENOMICS vector EBNM:
#
#     A_L = sum_j(tau_j * E[f^2_{jk}])   [scalar, same for all i]
#     B_L[i] = (R_k %*% (Tau * EF_k))[i] [n-vector]
#
#     EBNM(x = B_L/A_L, s = 1/sqrt(A_L))   <- VECTOR call (n observations)
#
#   Unlike the original dual-source update, there are no Cox weights w_i,
#   no EBeta_k, no EBeta2_k, and no alpha mixing — L is entirely driven by
#   genomics.  Survival signal enters the CAVI through q(F) instead (see
#   update_F.R and derivations/cox_on_YF/qF_supervised_derivation.tex).
#
# DESIGN PATTERN:
#   Part of the modular update series:
#     update_beta.R  ->  update_L.R  ->  update_F.R  ->  update_tau.R
#   Each exposes a _k() (single factor) and _all() (full loop) function.
#
# DEPENDENCIES:
#   ebnm  (CRAN)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Helper: partial residual matrix R^{-k}
# =============================================================================

#' Compute the partial residual matrix R^{-k}
#'
#' Removes all factors EXCEPT k from the data matrix Y, yielding the
#' residual that factor k should explain:
#'
#'   R_k = Y - EL %*% t(EF) + outer(EL[,k], EF[,k])
#'       = Y - sum_{k' != k} EL[,k'] %*% t(EF[,k'])
#'
#' This is an n x p matrix.  Used by both the L and F updates.
#'
#' @param Y   n x p matrix: observed genomics data
#' @param EL  n x K matrix: posterior means of all loadings
#' @param EF  p x K matrix: posterior means of all factors
#' @param k   integer: factor index to isolate (1-based)
#' @return n x p matrix: partial residual for factor k
#' @export
#' @family L_update
#' @seealso \code{\link{update_F_k}} which also uses R_k from this function
compute_R_k <- function(Y, EL, EF, k) {
  Y_hat <- EL %*% t(EF)                             # n x p: full reconstruction
  Y - Y_hat + outer(EL[, k], EF[, k])               # add back factor k
}

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(l_{ik}) for a single factor k  [Cluster B: pure-genomics EBNM]
#'
#' Under the Cox-on-YF reformulation (eta = Z_F * beta_tilde), L appears only
#' in the genomics likelihood.  The update is a VECTOR EBNM with n observations
#' driven entirely by the genomics residual R_k.
#'
#' @param Tau         p-vector: feature-specific noise precision tau_j (> 0)
#' @param EF_k        p-vector: posterior mean of factor column k (f_bar_{jk})
#' @param EF2_k       p-vector: posterior SECOND MOMENT of factor column k
#'                    (E_q[f_{jk}^2] = Var + mean^2).  Must satisfy EF2_k >= EF_k^2.
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
#' @family L_update
#' @examples
#' library(ebnm)
#' set.seed(1); n <- 100; p <- 50
#' Tau  <- rep(2.0, p)
#' EF_k <- rnorm(p);  EF2_k <- EF_k^2 + 0.1
#' R_k  <- matrix(rnorm(n * p), n, p)
#' res  <- update_L_k(Tau, EF_k, EF2_k, R_k)
#' cat("Loading column estimates (first 5):", round(res$mean[1:5], 3), "\n")
update_L_k <- function(Tau, EF_k, EF2_k,
                        R_k,
                        prior_family = "point_exponential",
                        A_floor      = 1e-10) {

  # ------------------------------------------------------------------
  # Precision A_L  (scalar — same for all subjects i)
  #
  #   A_L = sum_j(tau_j * E[f^2_{jk}])
  #
  # Uses E[f^2] (second moment) for error-in-variables correction:
  # posterior uncertainty in F inflates the effective noise for L,
  # preventing overfitting to uncertain factor estimates.
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
  #
  # A_L is a scalar, so x and s have the same shape as B_L:
  #   x[i] = B_L[i] / A_L   (n-vector)
  #   s    = 1/sqrt(A_L)     (scalar — same noise for all i)
  # ------------------------------------------------------------------
  x_L <- B_L / A_L
  s_L <- 1.0 / sqrt(A_L)
  x_L[!is.finite(x_L)] <- 0

  # ------------------------------------------------------------------
  # Solve the n-dimensional EBNM problem
  # ------------------------------------------------------------------
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
#' update_L_k() for each factor.  Uses Gauss-Seidel ordering: once EL[,k]
#' is updated, the new values are used for subsequent factors k' > k.
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
#'   $details -- length-K list, each element is the full update_L_k result
#' @export
#' @family L_update
update_L_all <- function(Y, EL, EL2, EF, EF2, Tau,
                          prior_family = "point_exponential",
                          A_floor      = 1e-10) {

  K        <- ncol(EL)
  EL_curr  <- EL
  EL2_curr <- EL2
  details  <- vector("list", K)

  for (k in seq_len(K)) {
    R_k <- compute_R_k(Y, EL_curr, EF, k)

    res_k <- update_L_k(
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
