# =============================================================================
# code/update_beta.R
#
# Modular q(beta_tilde_k) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE (Cluster B):
#   derivations/cox_on_YF/qBeta_YF_derivation.tex  (D3)
#   derivations/MF_UpdateDerivations/  q(beta) section (original, for comparison)
#
# MATHEMATICAL SUMMARY (Cluster B — Cox-on-YF reformulation):
#   Under eta_i = Z_F * beta_tilde where Z_F = Y * E[F] (observed projection
#   scores), beta_tilde_k appears ONLY in the Cox survival likelihood via Z_F.
#   The update reduces to a 1D EBNM problem:
#
#     A_k = alpha * sum_i w_i * ZF_{ik}^2      (ZF_ik^2 exact — no 2nd moment)
#     B_k = alpha * sum_i w_i * z_i^{-k} * ZF_{ik}
#     x_k = B_k / A_k   (alpha cancels in ratio)
#     s_k = 1 / sqrt(A_k)
#
#   Key difference from original q(beta): ZF_k = Y * E[f_k] is observed
#   (not latent), so no second-moment correction (EL2_k) is needed.
#   alpha controls shrinkage strength; alpha=1 gives the fully-scaled update.
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

#' Compute the partial working response z^{-k}_i  [Cluster B: ZF predictor]
#'
#' Removes the expected contribution of factor k from the working response z,
#' giving the residual signal that factor k's beta_tilde should explain.
#'
#' z^{-k}_i = z_i - sum_{k' != k} ZF_{ik'} * beta_tilde_{k'}
#'
#' This is equivalent to: z - (ZF %*% EBeta) + ZF[, k] * EBeta[k]
#' because the sum over k' != k is the full eta minus factor k.
#'
#' NOTE: z_no_k does NOT depend on ZF[,k] or EBeta[k], only on k' != k.
#' The SAME z_no_k can be reused for the F update of factor k (since
#' the F update changes EF[,k] which changes ZF[,k], but ZF is held fixed
#' within each outer CAVI iteration).
#'
#' @param z      n-vector: full working response z_i = eta_hat_i + u_i/w_i
#' @param ZF     n x K matrix: observed projection scores Y * E[F]
#' @param EBeta  K-vector: current posterior means of beta_tilde
#' @param k      integer: factor index to exclude (1-based)
#' @return n-vector: partial working response for factor k
#' @export
#' @family beta_update
#' @seealso \code{\link{update_F_k}} which uses z_no_k for the F update
compute_z_no_k <- function(z, ZF, EBeta, k) {
  eta_no_k <- as.vector(ZF %*% EBeta) - ZF[, k] * EBeta[k]
  z - eta_no_k
}

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(beta_tilde_k) for a single factor k  [Cluster B: ZF predictor]
#'
#' Performs the EBNM-based coordinate-ascent update for the reparameterised
#' survival coefficient beta_tilde_k given Cox Taylor working quantities
#' (z_no_k, w) and the k-th column of the observed projection scores ZF.
#'
#' ZF_k = Y * E[f_k] is treated as a fixed observed quantity per CAVI outer
#' iteration, so there is no second-moment correction — ZF_k^2 is exact.
#'
#' @param w            n-vector: Cox neg-diagonal Hessian weights w_i (> 0)
#' @param z_no_k       n-vector: partial working response z^{-k}_i.
#'                     Compute via compute_z_no_k() or inline as
#'                     z - (ZF %*% EBeta) + ZF[,k]*EBeta[k].
#' @param ZF_k         n-vector: k-th column of Y * E[F] (projection scores).
#'                     These are observed quantities — no posterior variance
#'                     correction is needed (unlike the original EL_k case).
#' @param prior_family character: EBNM prior family (default "point_normal").
#' @param alpha        numeric in [0, 1]: scaling weight on survival term.
#'                     A_k = alpha * sum(w * ZF_k^2).
#'                     alpha cancels in x_k = B_k/A_k; only s_k is affected.
#'                     Default 0.5.
#' @param A_floor      numeric: minimum value for A_k (default 1e-10).
#'
#' @return Named list:
#'   $mean        -- posterior mean E_q[beta_tilde_k]
#'   $second      -- posterior 2nd moment E_q[beta_tilde_k^2] = sd^2 + mean^2
#'   $sd          -- posterior SD sqrt(Var_q(beta_tilde_k))
#'   $A           -- precision A_k = alpha * sum(w * ZF_k^2) [floored]
#'   $B           -- signal B_k = alpha * sum(w * z_no_k * ZF_k)
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
#' ZF_k   <- rnorm(n)           # observed projection scores Y*E[f_k]
#' z_no_k <- ZF_k * 1.5 + rnorm(n, sd = 0.5)
#' res <- update_beta_k(w, z_no_k, ZF_k)
#' cat("beta_tilde estimate:", round(res$mean, 3), "\n")
update_beta_k <- function(w, z_no_k, ZF_k,
                          prior_family = "point_normal",
                          alpha        = 0.5,
                          A_floor      = 1e-10) {

  # ------------------------------------------------------------------
  # Precision (A_k): alpha * sum_i w_i * ZF_{ik}^2
  #
  # ZF_k is observed (not latent), so ZF_k^2 is exact — no second-moment
  # correction is needed.  This differs from the original beta update which
  # used E_q[l_{ik}^2] = Var_q + mean^2 to account for shrinkage uncertainty.
  # alpha controls the scale of shrinkage; alpha cancels in x_k = B_k/A_k.
  # ------------------------------------------------------------------
  A_k <- max(alpha * sum(w * ZF_k^2), A_floor)

  # ------------------------------------------------------------------
  # Signal (B_k): alpha * sum_i w_i * z_i^{-k} * ZF_{ik}
  #
  # Note: alpha cancels in x_k = B_k/A_k, so the EBNM pseudo-observation
  # is alpha-independent.  Only s_k = 1/sqrt(A_k) depends on alpha.
  # ------------------------------------------------------------------
  B_k <- alpha * sum(w * z_no_k * ZF_k)

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise
  #   x_k = B_k / A_k  (normal likelihood mean)
  #   s_k = 1/sqrt(A_k) (normal likelihood sd)
  # ------------------------------------------------------------------
  x_k <- B_k / A_k
  s_k <- 1.0 / sqrt(A_k)
  x_k[!is.finite(x_k)] <- 0
  s_k[!is.finite(s_k) | s_k <= 0] <- 1e5
  s_k <- pmax(s_k, 1e-8)

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

#' Update q(beta_tilde) for all K factors  [Cluster B: ZF predictor]
#'
#' Iterates over k = 1..K, computing the partial working response z^{-k}_i
#' and calling update_beta_k() for each factor.  Uses Gauss-Seidel ordering:
#' once beta_tilde_k is updated, the new value is used for z^{-k'}_i for
#' subsequent factors k' > k.
#'
#' ZF = Y * E[F] is pre-computed by the caller once per outer CAVI iteration
#' and held fixed across the k-loop.
#'
#' @param w            n-vector: Cox neg-diagonal Hessian weights w_i
#' @param z            n-vector: full working response z_i = eta_hat + u/w
#' @param ZF           n x K matrix: observed projection scores Y * E[F]
#' @param EBeta        K-vector: current posterior means (Gauss-Seidel start)
#' @param prior_family character: EBNM prior family (default "point_normal")
#' @param alpha        numeric in [0, 1]: survival mixing weight (default 0.5)
#' @param A_floor      numeric: precision floor (default 1e-10)
#'
#' @return Named list:
#'   $EBeta   -- K-vector of updated posterior means
#'   $EBeta2  -- K-vector of updated posterior second moments
#'   $details -- length-K list, each element is the full update_beta_k result
#' @export
#' @family beta_update
update_beta_all <- function(w, z, ZF, EBeta,
                            prior_family = "point_normal",
                            alpha        = 0.5,
                            A_floor      = 1e-10) {

  K          <- ncol(ZF)
  EBeta_curr <- EBeta
  EBeta2_new <- numeric(K)
  details    <- vector("list", K)

  for (k in seq_len(K)) {
    z_no_k <- compute_z_no_k(z, ZF, EBeta_curr, k)

    res_k <- update_beta_k(
      w          = w,
      z_no_k     = z_no_k,
      ZF_k       = ZF[, k],
      prior_family = prior_family,
      alpha      = alpha,
      A_floor    = A_floor
    )

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
