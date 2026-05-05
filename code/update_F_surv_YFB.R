# =============================================================================
# code/update_F_surv_YFB.R
#
# q(F) update for the Cox-on-YF model (Cluster B: eta = YFB)
#
# MODEL DIFFERENCE: Under eta = YFB, F appears in BOTH the genomics likelihood
# (Y ≈ LF') AND the Cox survival likelihood via ZF = YF. This is a DUAL-SOURCE
# update. Compare to update_F.R (Cluster A), which is PURE-GENOMICS only.
#
# ALPHA_F=0 DEFAULT: The dual-source F update with alpha>0 creates a positive-
# feedback instability (see DECISIONS.md 2026-04-30). With alpha_F=0, F is
# learned from pure genomics; EBeta becomes non-zero because A_beta = sum(w*ZF^2)
# is non-zero from SVD initialization regardless of EBeta. This is empirically
# equivalent to the EBMF warm-start result (unsupervised F + Cox on ZF scores).
#
# WHY ALPHA_F=0 RESOLVES EBeta=0:
#   1. EF initialized from SVD → non-zero from iter 1
#   2. ZF = Y·EF is non-zero from iter 1 (observed, not latent)
#   3. β update: A_β = Σᵢ w·ZF²_k ≠ 0 even before any CAVI iteration
#   4. No feedback loop: A_gen = Tau·sum(EL²_k) does not depend on EF directly
#   5. Train/test consistency preserved: prediction still uses Y_test·EF·β̃
#
# DERIVATION REFERENCE:
#   derivations/cox_on_YF/qF_supervised_derivation.tex  (dual-source derivation)
#   DECISIONS.md 2026-04-30  (alpha_F=0 rationale)
#
# MATHEMATICAL SUMMARY (dual-source, alpha_F > 0 case):
#   A_F[j] = (1-alpha) * tau_j * sum_i E[l^2_{ik}]           [genomics]
#           + alpha * E[beta_tilde_k^2] * sum_i w_i * y_ij^2  [survival]
#   B_F[j] = (1-alpha) * tau_j * (R_k' EL_k)[j]              [genomics]
#           + alpha * E[beta_tilde_k] * (Y'(w*z_no_k))[j]     [survival]
#   EBNM(x = B_F/A_F, s = 1/sqrt(A_F))   <- VECTOR (p observations)
#
# With alpha_F=0 (default): reduces to pure-genomics EBNM where tau cancels
# in x_j = B_gen/A_gen (same as Cluster A update_F.R but with explicit alpha).
#
# DESIGN PATTERN:
#   Part of the Cluster B modular update series:
#     update_beta.R  ->  update_L_surv_YFB.R  ->  update_F_surv_YFB.R  ->  update_tau.R
#
# DEPENDENCIES:
#   ebnm  (CRAN)
#   compute_R_k() from update_L.R (sourced by fit_cox_on_yf.R)
#   compute_z_no_k() from update_beta.R (sourced by fit_cox_on_yf.R)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(f_{jk}) for a single factor k  [Cluster B: dual-source, alpha_F=0 default]
#'
#' Under the Cox-on-YF reformulation, F appears in both the genomics likelihood
#' and the Cox survival likelihood via ZF = YF. The update is a dual-source
#' vector EBNM with p observations. With alpha=0 (default), reduces to
#' pure-genomics, preventing the positive-feedback instability.
#'
#' @param Tau         p-vector: feature-specific noise precision tau_j (> 0)
#' @param EL_k        n-vector: posterior mean of loading column k
#' @param EL2_k       n-vector: posterior SECOND MOMENT of loading column k.
#'                    Must satisfy EL2_k >= EL_k^2.
#' @param R_k         n x p matrix: partial residual R^{-k} (from compute_R_k)
#' @param EBeta_k     scalar: posterior mean of beta_tilde_k (default 0)
#' @param EBeta2_k    scalar: posterior 2nd moment of beta_tilde_k (default 0)
#' @param YtWY_diag   p-vector: diagonal of Y' diag(w) Y (default NULL → zeros)
#' @param YtWz_no_k   p-vector: Y' (w * z_no_k) (default NULL → zeros)
#' @param prior_family character: EBNM prior family (default "point_exponential")
#' @param alpha       numeric in [0, 1]: survival mixing weight.
#'                    Default 0 (pure genomics). Set alpha > 0 only with caution:
#'                    the dual-source F update with alpha > 0 has a positive-
#'                    feedback instability documented in DECISIONS.md 2026-04-30.
#' @param A_floor     numeric: minimum value for each A_F[j] (default 1e-10)
#'
#' @return Named list:
#'   $mean        -- p-vector: posterior mean E_q[f_{jk}]
#'   $second      -- p-vector: posterior 2nd moment E_q[f_{jk}^2]
#'   $sd          -- p-vector: posterior SD sqrt(Var_q(f_{jk}))
#'   $A           -- p-vector: combined precision (floored)
#'   $B           -- p-vector: combined signal
#'   $A_gen       -- p-vector: raw genomics precision
#'   $A_surv      -- p-vector: raw survival precision
#'   $B_gen       -- p-vector: raw genomics signal
#'   $B_surv      -- p-vector: raw survival signal
#'   $x           -- p-vector: EBNM pseudo-obs x_j = B_j/A_j
#'   $s           -- p-vector: EBNM pseudo-noise s_j = 1/sqrt(A_j)
#'   $sum_EL2_k   -- scalar: sum_i E[l_{ik}^2] (diagnostic)
#'   $ebnm_result -- raw ebnm() return object
#'
#' @export
#' @family F_surv_YFB_update
#' @seealso \code{\link{update_L_surv_YFB_k}} for the pure-genomics L update
update_F_surv_YFB_k <- function(Tau, EL_k, EL2_k, R_k,
                                 EBeta_k      = 0,
                                 EBeta2_k     = 0,
                                 YtWY_diag    = NULL,
                                 YtWz_no_k    = NULL,
                                 prior_family = "point_exponential",
                                 alpha        = 0,        # alpha_F=0: pure genomics
                                 A_floor      = 1e-10) {

  p <- length(Tau)

  if (is.null(YtWY_diag)) YtWY_diag <- rep(0, p)
  if (is.null(YtWz_no_k)) YtWz_no_k <- rep(0, p)

  # ------------------------------------------------------------------
  # Genomics precision and signal  (p-vectors)
  #
  #   A_gen[j] = tau_j * sum_i E[l_{ik}^2]
  #   B_gen[j] = tau_j * (R_k' EL_k)[j]
  #
  # With alpha=0, A_F = A_gen and B_F = B_gen. In the ratio x_j = B_gen/A_gen,
  # tau_j cancels, making the pseudo-observation tau-free (key Cluster A
  # property preserved here at alpha=0).
  # ------------------------------------------------------------------
  sum_EL2_k <- sum(EL2_k)
  A_gen     <- Tau * sum_EL2_k                             # p-vector
  B_gen     <- Tau * as.vector(t(R_k) %*% EL_k)           # p-vector

  # ------------------------------------------------------------------
  # Survival precision and signal  (p-vectors)
  #
  #   A_surv[j] = E[beta_tilde_k^2] * sum_i w_i * y_ij^2
  #   B_surv[j] = E[beta_tilde_k] * sum_i w_i * z_no_k_i * y_ij
  #
  # With alpha=0 (default), these are zero-weighted in the combination.
  # ------------------------------------------------------------------
  A_surv <- EBeta2_k * YtWY_diag                          # p-vector
  B_surv <- EBeta_k  * YtWz_no_k                          # p-vector

  # ------------------------------------------------------------------
  # Alpha-weighted combination and precision floor
  # ------------------------------------------------------------------
  A_F <- pmax((1 - alpha) * A_gen + alpha * A_surv, A_floor)
  B_F <- (1 - alpha) * B_gen + alpha * B_surv

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise (p-vectors)
  # ------------------------------------------------------------------
  x_F <- B_F / A_F
  s_F <- 1.0 / sqrt(A_F)
  x_F[!is.finite(x_F)] <- 0
  s_F[!is.finite(s_F) | s_F <= 0] <- 1e5
  s_F <- pmax(s_F, 1e-8)

  res <- ebnm(x = x_F, s = s_F, prior_family = prior_family)

  f_mean   <- res$posterior$mean
  f_sd     <- res$posterior$sd
  f_second <- f_sd^2 + f_mean^2

  list(
    mean        = f_mean,
    second      = f_second,
    sd          = f_sd,
    A           = A_F,
    B           = B_F,
    A_gen       = A_gen,
    A_surv      = A_surv,
    B_gen       = B_gen,
    B_surv      = B_surv,
    x           = x_F,
    s           = s_F,
    sum_EL2_k   = sum_EL2_k,
    ebnm_result = res
  )
}

# =============================================================================
# Convenience wrapper: full Gauss-Seidel loop over all K factors
# =============================================================================

#' Update q(F) for all K factors  [Cluster B: dual-source, alpha_F=0 default]
#'
#' Iterates over k = 1..K, computing R_k and calling update_F_surv_YFB_k() for
#' each factor. Uses Gauss-Seidel ordering. Survival arguments are optional.
#'
#' @param Y           n x p matrix: observed genomics data
#' @param EL          n x K matrix: posterior means of loadings (already updated)
#' @param EL2         n x K matrix: posterior second moments of loadings
#' @param EF          p x K matrix: posterior means of factors (will be updated)
#' @param EF2         p x K matrix: posterior second moments of factors
#' @param Tau         p-vector: feature-specific noise precision
#' @param w           n-vector: Cox weights (optional; NULL = no survival terms)
#' @param z           n-vector: full working response (optional)
#' @param EBeta       K-vector: posterior means of beta_tilde (optional)
#' @param EBeta2      K-vector: posterior 2nd moments of beta_tilde (optional)
#' @param prior_family character: EBNM prior family (default "point_exponential")
#' @param alpha       numeric in [0, 1]: survival mixing weight (default 0)
#' @param A_floor     numeric: precision floor (default 1e-10)
#'
#' @return Named list:
#'   $EF      -- p x K matrix of updated posterior means
#'   $EF2     -- p x K matrix of updated posterior second moments
#'   $details -- length-K list, each element is the full update_F_surv_YFB_k result
#' @export
#' @family F_surv_YFB_update
update_F_surv_YFB_all <- function(Y, EL, EL2, EF, EF2, Tau,
                                   w            = NULL,
                                   z            = NULL,
                                   EBeta        = NULL,
                                   EBeta2       = NULL,
                                   prior_family = "point_exponential",
                                   alpha        = 0,        # alpha_F=0 default
                                   A_floor      = 1e-10) {

  K        <- ncol(EF)
  EF_curr  <- EF
  EF2_curr <- EF2
  details  <- vector("list", K)

  use_survival <- !is.null(w) && !is.null(z) && !is.null(EBeta) && !is.null(EBeta2)

  if (use_survival && alpha > 0) {
    ZF_init   <- Y %*% EF
    YtWY_diag <- as.vector(t(Y^2) %*% w)
  }

  for (k in seq_len(K)) {
    R_k <- compute_R_k(Y, EL, EF_curr, k)

    if (use_survival && alpha > 0) {
      z_no_k    <- compute_z_no_k(z, ZF_init, EBeta, k)
      YtWz_no_k <- as.vector(t(Y) %*% (w * z_no_k))
    } else {
      YtWz_no_k <- NULL
    }

    res_k <- update_F_surv_YFB_k(
      Tau          = Tau,
      EL_k         = EL[, k],
      EL2_k        = EL2[, k],
      R_k          = R_k,
      EBeta_k      = if (use_survival) EBeta[k]  else 0,
      EBeta2_k     = if (use_survival) EBeta2[k] else 0,
      YtWY_diag    = if (use_survival && alpha > 0) YtWY_diag else NULL,
      YtWz_no_k    = YtWz_no_k,
      prior_family = prior_family,
      alpha        = alpha,
      A_floor      = A_floor
    )

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
