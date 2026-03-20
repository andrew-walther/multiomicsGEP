# =============================================================================
# code/update_beta.R
#
# Modular q(beta_k) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE:
#   derivations/qB/qBeta_update_derivation.tex  (self-contained)
#   derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.tex
#   Section 6 (EBNM formulation for survival coefficients)
#
# MATHEMATICAL SUMMARY:
#   beta_k appears ONLY in the Cox survival likelihood.  The coordinate-ascent
#   update reduces to a 1D EBNM problem with:
#
#     A_k = sum_i W_{ii} * E_q[l_{ik}^2]     (precision; uses full 2nd moment)
#     B_k = sum_i W_{ii} * z_i^{-k} * l_bar_{ik}  (signal)
#     x_k = B_k / A_k,   s_k = 1 / sqrt(A_k)
#     (g_hat_{beta_k}, q_hat_{beta_k}) = EBNM(x_k, s_k)
#
#   The use of E_q[l^2] (not l_bar^2) is an error-in-variables correction:
#   posterior uncertainty in L appropriately inflates the effective noise for
#   beta_k, preventing overfitting to uncertain loadings.
#
# DESIGN PATTERN:
#   This module is the first of a series:
#     update_beta.R  ->  update_L.R  ->  update_F.R  ->  update_tau.R
#   Each exposes a _k() (single factor) and _all() (full loop) function.
#
# DEPENDENCIES:
#   ebnm  (CRAN)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Helper: partial working response z^{-k}_i
# =============================================================================

#' Compute the partial working response z^{-k}_i
#'
#' Removes the expected contribution of factor k from the working response z,
#' giving the residual signal that factor k's beta should explain.
#'
#' z^{-k}_i = z_i - sum_{k' != k} l_bar_{ik'} * beta_bar_{k'}
#'
#' This is equivalent to: z - (EL %*% EBeta) + EL[, k] * EBeta[k]
#' because sum_{k'} l_bar_{ik'} beta_k' - l_bar_{ik} beta_k
#'        = sum_{k' != k} l_bar_{ik'} beta_k'
#'
#' NOTE: z_no_k does NOT depend on l_{ik} or beta_k, only on k' != k.
#' This means the SAME z_no_k can be reused for both the L update and
#' the beta update of factor k (the L update changes EL[,k] but not EL[,k']
#' for k' != k).  See REVISED.tex Sec. 6 "z_no_k reuse rationale".
#'
#' @param z      n-vector: full working response z_i = eta_hat_i + u_i/W_{ii}
#' @param EL     n x K matrix: posterior means of loadings
#' @param EBeta  K-vector: current posterior means of survival coefficients
#' @param k      integer: factor index to exclude (1-based)
#' @return n-vector: partial working response for factor k
#' @export
#' @family beta_update
#' @seealso \code{\link{update_L_k}} which reuses z_no_k for the loading update
compute_z_no_k <- function(z, EL, EBeta, k) {
  eta_no_k <- as.vector(EL %*% EBeta) - EL[, k] * EBeta[k]
  z - eta_no_k
}

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(beta_k) for a single factor k
#'
#' Performs the EBNM-based coordinate-ascent update for the survival
#' coefficient beta_k given the Cox Taylor working quantities (z_no_k, w),
#' and the current posterior moments of the k-th loading column.
#'
#' The function is intentionally decoupled from the CAVI loop: it takes
#' pre-computed vectors and returns a self-contained result list.  This
#' makes it independently testable and reusable outside the full algorithm.
#'
#' @param w            n-vector: Cox neg-diagonal Hessian weights W_{ii} (> 0)
#' @param z_no_k       n-vector: partial working response z^{-k}_i.
#'                     Compute via compute_z_no_k() or inline as
#'                     z - (EL %*% EBeta) + EL[,k]*EBeta[k].
#' @param EL_k         n-vector: posterior mean of loadings for factor k
#'                     (l_bar_{ik} = E_q[l_{ik}])
#' @param EL2_k        n-vector: posterior SECOND MOMENT of loadings for k
#'                     (E_q[l_{ik}^2] = Var_q(l_{ik}) + l_bar_{ik}^2).
#'                     Must satisfy EL2_k >= EL_k^2.
#'                     [CRITICAL: use second moment, NOT squared mean]
#' @param prior_family character: EBNM prior family (default "point_normal").
#'                     "point_normal" promotes sparsity in beta.
#'                     Use "normal" for a purely Gaussian prior.
#' @param A_floor      numeric: minimum value for A_k to prevent 0-division.
#'                     Default 1e-10.  [A3 in V2.R]
#'
#' @return Named list:
#'   $mean        -- posterior mean E_q[beta_k]
#'   $second      -- posterior 2nd moment E_q[beta_k^2] = sd^2 + mean^2
#'   $sd          -- posterior standard deviation sqrt(Var_q(beta_k))
#'   $A           -- precision A_k = sum(w * EL2_k) [floored]
#'   $B           -- signal B_k = sum(w * z_no_k * EL_k)
#'   $x           -- EBNM pseudo-obs x_k = B_k/A_k
#'   $s           -- EBNM pseudo-noise s_k = 1/sqrt(A_k)
#'   $ebnm_result -- raw ebnm() return object (for diagnostics)
#'
#' @export
#' @family beta_update
#' @examples
#' library(ebnm)
#' set.seed(1); n <- 100
#' w      <- rep(2.0, n)
#' EL_k   <- rnorm(n)
#' EL2_k  <- EL_k^2 + 0.1       # second moment > squared mean
#' z_no_k <- EL_k * 1.5 + rnorm(n, sd = 0.5)
#' res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
#' cat("beta estimate:", round(res$mean, 3), "\n")
update_beta_k <- function(w, z_no_k, EL_k, EL2_k,
                          prior_family = "point_normal",
                          A_floor      = 1e-10) {

  # ------------------------------------------------------------------
  # Precision (A_k): error-in-variables correction
  #   sum_i W_{ii} * E_q[l_{ik}^2]
  # Using the full second moment (Var + mean^2) rather than the squared
  # mean inflates the effective noise for beta_k, preventing overfitting
  # to uncertain loadings (Companion.tex Sec. 6.2, Eq. A_k).
  # ------------------------------------------------------------------
  # Floor triggers when all weights are zero or all loadings are zero —
  # degenerate inputs that would otherwise cause division by zero in x_k and s_k.
  A_k <- max(sum(w * EL2_k), A_floor)

  # ------------------------------------------------------------------
  # Signal (B_k): weighted inner product of partial response and loading
  #   sum_i W_{ii} * z_i^{-k} * l_bar_{ik}
  # ------------------------------------------------------------------
  B_k <- sum(w * z_no_k * EL_k)

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise
  #   x_k = B_k / A_k  (normal likelihood mean)
  #   s_k = 1/sqrt(A_k) (normal likelihood sd)
  # ------------------------------------------------------------------
  x_k <- B_k / A_k
  s_k <- 1.0 / sqrt(A_k)

  # ------------------------------------------------------------------
  # Solve the 1D EBNM problem: EBNM(x_k, s_k)
  # Returns posterior (mean, sd) and fitted prior g_hat.
  # Note: the 1/sigma^2 prior precision term from the Gaussian prior
  # g(beta_k) ~ N(0, sigma^2) is subsumed into the EBNM prior
  # estimation — ebnm() learns the prior variance from data.
  # (REVISED.tex lines 911-914)
  # ------------------------------------------------------------------
  res <- ebnm(x = x_k, s = s_k, prior_family = prior_family)

  beta_mean   <- res$posterior$mean
  beta_sd     <- res$posterior$sd
  beta_second <- beta_sd^2 + beta_mean^2

  # Pack all diagnostics for caller / test inspection.
  # Returning A, B, x, s enables tests to verify math identities directly.
  list(
    mean        = beta_mean,
    second      = beta_second,
    sd          = beta_sd,
    A           = A_k,
    B           = B_k,
    x           = x_k,
    s           = s_k,
    ebnm_result = res
  )
}

# =============================================================================
# Convenience wrapper: full Gauss-Seidel loop over all K factors
# =============================================================================

#' Update q(beta) for all K factors (Gauss-Seidel CAVI loop)
#'
#' Iterates over k = 1..K, computing the partial working response z^{-k}_i
#' and calling update_beta_k() for each factor.  Uses Gauss-Seidel ordering:
#' once beta_k is updated, the new value is used when computing z^{-k'}_i
#' for subsequent factors k' > k.
#'
#' This is a "pure" function with respect to L and F: it does NOT modify
#' the loading or factor matrices.  It only updates beta.
#'
#' INTEGRATION NOTE:
#'   In the full CAVI loop (V2.R), z and w are recomputed at the start of
#'   each outer iteration (calc_cox_taylor).  The beta update uses the SAME
#'   z_no_k that was computed for the L update of that factor.  For
#'   standalone use of update_beta_all(), pass the current z and w directly.
#'
#' @param w            n-vector: Cox neg-diagonal Hessian weights W_{ii}
#' @param z            n-vector: full working response z_i = eta_hat + u/w
#' @param EL           n x K matrix: posterior means of all loadings
#' @param EL2          n x K matrix: posterior second moments of all loadings
#' @param EBeta        K-vector: current posterior means (warm start for
#'                     computing z_no_k; will be updated in-place within loop)
#' @param prior_family character: EBNM prior family (default "point_normal")
#' @param A_floor      numeric: precision floor (default 1e-10)
#'
#' @return Named list:
#'   $EBeta   -- K-vector of updated posterior means
#'   $EBeta2  -- K-vector of updated posterior second moments
#'   $details -- length-K list, each element is the full update_beta_k result
#' @export
#' @family beta_update
update_beta_all <- function(w, z, EL, EL2, EBeta,
                            prior_family = "point_normal",
                            A_floor      = 1e-10) {

  K          <- ncol(EL)
  # Mutable copy: Gauss-Seidel requires updating beta_k in-place so that
  # z_no_k for factor k' > k incorporates the freshly updated beta_k.
  EBeta_curr <- EBeta
  EBeta2_new <- numeric(K)
  details    <- vector("list", K)

  for (k in seq_len(K)) {
    # Partial working response (uses Gauss-Seidel EBeta_curr)
    z_no_k <- compute_z_no_k(z, EL, EBeta_curr, k)

    res_k <- update_beta_k(
      w          = w,
      z_no_k     = z_no_k,
      EL_k       = EL[, k],
      EL2_k      = EL2[, k],
      prior_family = prior_family,
      A_floor    = A_floor
    )

    # Gauss-Seidel: update EBeta_curr so next k uses the new value
    EBeta_curr[k]  <- res_k$mean
    EBeta2_new[k]  <- res_k$second
    details[[k]]   <- res_k
  }

  list(
    EBeta   = EBeta_curr,
    EBeta2  = EBeta2_new,
    details = details
  )
}
