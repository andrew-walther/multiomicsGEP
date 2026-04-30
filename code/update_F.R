# =============================================================================
# code/update_F.R
#
# Modular q(f_{jk}) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE (Cluster B):
#   derivations/cox_on_YF/qF_supervised_derivation.tex  (D2)
#   derivations/MF_UpdateDerivations/  q(F) section (original, for comparison)
#
# MATHEMATICAL SUMMARY (Cluster B — Cox-on-YF reformulation):
#   Under eta_i = (Y_i F) beta_tilde, f_{jk} appears in BOTH the genomics
#   likelihood (Y ≈ LF') AND the Cox survival likelihood via ZF = YF.  The
#   update reduces to a dual-source p-dimensional EBNM:
#
#     A_F[j] = (1-alpha) * tau_j * sum_i E[l^2_{ik}]           [genomics]
#            + alpha * E[beta_tilde_k^2] * sum_i w_i * y_ij^2  [survival]
#
#     B_F[j] = (1-alpha) * tau_j * (R_k' EL_k)[j]              [genomics]
#            + alpha * E[beta_tilde_k] * (Y'(w*z_no_k))[j]     [survival]
#
#     EBNM(x = B_F/A_F, s = 1/sqrt(A_F))   <- VECTOR (p observations)
#
#   IMPORTANT: The genomics pseudo-observation x_j = B_gen/A_gen is still
#   tau-free (tau cancels in the ratio when alpha=0).  With alpha>0, the
#   full x_j = B_F/A_F is no longer tau-free.
#
#   Pre-computed inputs from the caller (fit_modular.R):
#     YtWY_diag[j] = sum_i w_i * y_ij^2 = diag(Y' diag(w) Y)[j]  [p-vec]
#     YtWz_no_k[j] = sum_i w_i * z_no_k_i * y_ij = (Y'(w*z_no_k))[j] [p-vec]
#   These are optional: when omitted (NULL), survival terms are zero.
#
# DESIGN PATTERN:
#   Part of the modular update series:
#     update_beta.R  ->  update_L.R  ->  update_F.R  ->  update_tau.R
#
# DEPENDENCIES:
#   ebnm  (CRAN)
#   compute_R_k() from code/update_L.R
#   compute_z_no_k() from code/update_beta.R (used by update_F_all when survival args provided)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Core: single-factor update
# =============================================================================

#' Update q(f_{jk}) for a single factor k  [Cluster B: dual-source EBNM]
#'
#' Under the Cox-on-YF reformulation, F appears in both the genomics likelihood
#' and the Cox survival likelihood.  The update is a dual-source vector EBNM
#' with p observations.  Survival terms are optional (default to zero).
#'
#' KEY PROPERTY (genomics-only limit): When survival args are omitted or zero,
#' the EBNM pseudo-observation x_j = B_gen/A_gen is tau-free (tau cancels in
#' the ratio).  With non-zero survival terms, x_j is no longer tau-free.
#'
#' @param Tau         p-vector: feature-specific noise precision tau_j (> 0)
#' @param EL_k        n-vector: posterior mean of loading column k
#' @param EL2_k       n-vector: posterior SECOND MOMENT of loading column k.
#'                    Must satisfy EL2_k >= EL_k^2.
#' @param R_k         n x p matrix: partial residual R^{-k} (from compute_R_k)
#' @param EBeta_k     scalar: posterior mean of beta_tilde_k (default 0)
#' @param EBeta2_k    scalar: posterior 2nd moment of beta_tilde_k (default 0)
#' @param YtWY_diag   p-vector: diagonal of Y' diag(w) Y, i.e. colSums(Y^2 * w)
#'                    (default NULL → treated as zeros → no survival contribution)
#' @param YtWz_no_k   p-vector: Y' (w * z_no_k), i.e. t(Y) %*% (w * z_no_k)
#'                    (default NULL → treated as zeros → no survival contribution)
#' @param prior_family character: EBNM prior family (default "point_exponential")
#' @param alpha       numeric in [0, 1]: survival mixing weight.
#'                    A_F = pmax((1-alpha)*A_gen + alpha*A_surv, floor)
#'                    B_F = (1-alpha)*B_gen + alpha*B_surv
#'                    alpha=0: pure genomics (original formula); alpha=1: pure survival.
#'                    Default 0.5.
#' @param A_floor     numeric: minimum value for each A_F[j] (default 1e-10)
#' @param normalize_AB logical: if TRUE, rescale A_surv and B_surv so that
#'                    mean(A_surv) matches mean(A_gen) before alpha-mixing.
#'                    Addresses the structural ~p/n scale imbalance where A_gen
#'                    (genomics precision, sums over n subjects) typically exceeds
#'                    A_surv (survival precision, = EBeta2 * YtWY_diag) by
#'                    orders of magnitude when EBeta is small.  Same mechanism as
#'                    the Cluster A normalize_AB fix applied to update_L_k().
#'                    Skipped automatically when mean(A_surv) <= 1e-12 (pure-zero
#'                    EBeta case) to avoid amplifying numerical noise.
#'                    Default FALSE for backward compatibility.
#'
#' @return Named list:
#'   $mean        -- p-vector: posterior mean E_q[f_{jk}]
#'   $second      -- p-vector: posterior 2nd moment E_q[f_{jk}^2]
#'   $sd          -- p-vector: posterior SD sqrt(Var_q(f_{jk}))
#'   $A           -- p-vector: combined precision (floored)
#'   $B           -- p-vector: combined signal
#'   $A_gen       -- p-vector: raw genomics precision (unfloored)
#'   $A_surv      -- p-vector: raw survival precision (unfloored)
#'   $B_gen       -- p-vector: raw genomics signal
#'   $B_surv      -- p-vector: raw survival signal
#'   $x           -- p-vector: EBNM pseudo-obs x_j = B_j/A_j
#'   $s           -- p-vector: EBNM pseudo-noise s_j = 1/sqrt(A_j)
#'   $sum_EL2_k   -- scalar: sum_i E[l_{ik}^2] (diagnostic)
#'   $ebnm_result -- raw ebnm() return object
#'
#' @export
#' @family F_update
#' @examples
#' library(ebnm)
#' set.seed(1); n <- 100; p <- 200
#' Tau  <- rep(2.0, p)
#' EL_k <- rnorm(n);  EL2_k <- EL_k^2 + 0.1
#' R_k  <- matrix(rnorm(n * p), n, p)
#' res  <- update_F_k(Tau, EL_k, EL2_k, R_k)
#' cat("Factor estimates (first 5):", round(res$mean[1:5], 3), "\n")
update_F_k <- function(Tau, EL_k, EL2_k, R_k,
                        EBeta_k      = 0,
                        EBeta2_k     = 0,
                        YtWY_diag    = NULL,
                        YtWz_no_k    = NULL,
                        prior_family = "point_exponential",
                        alpha        = 0.5,
                        A_floor      = 1e-10,
                        normalize_AB = FALSE) {

  p <- length(Tau)

  # Default: zero out survival contribution when pre-computed matrices not provided
  if (is.null(YtWY_diag)) YtWY_diag <- rep(0, p)
  if (is.null(YtWz_no_k)) YtWz_no_k <- rep(0, p)

  # ------------------------------------------------------------------
  # Genomics precision and signal  (p-vectors)
  #
  #   A_gen[j] = tau_j * sum_i E[l_{ik}^2]
  #   B_gen[j] = tau_j * (R_k' EL_k)[j]
  #
  # sum_EL2_k is a scalar that broadcasts across all p features.
  # EL2_k (second moment) corrects for posterior uncertainty in L
  # (error-in-variables), preventing overfitting.
  # ------------------------------------------------------------------
  sum_EL2_k <- sum(EL2_k)
  A_gen     <- Tau * sum_EL2_k                             # p-vector (raw)
  B_gen     <- Tau * as.vector(t(R_k) %*% EL_k)           # p-vector (raw)

  # ------------------------------------------------------------------
  # Survival precision and signal  (p-vectors)
  #
  #   A_surv[j] = E[beta_tilde_k^2] * sum_i w_i * y_ij^2
  #             = EBeta2_k * YtWY_diag[j]
  #   B_surv[j] = E[beta_tilde_k] * sum_i w_i * z_no_k_i * y_ij
  #             = EBeta_k * YtWz_no_k[j]
  #
  # YtWY_diag and YtWz_no_k are pre-computed by the caller to avoid
  # passing Y and w repeatedly.  When not provided, default to zero.
  # ------------------------------------------------------------------
  A_surv <- EBeta2_k * YtWY_diag                          # p-vector (raw)
  B_surv <- EBeta_k  * YtWz_no_k                          # p-vector (raw)

  # ------------------------------------------------------------------
  # Optional A/B rescaling  [Cluster B analogue of Cluster A Fix 4]
  #
  # A_gen = Tau * sum_EL2_k sums over n subjects; A_surv = EBeta2 * YtWY_diag
  # is bounded by EBeta2 (tiny at initialisation).  The structural ratio
  # A_gen / A_surv can exceed 10^3, making alpha effectively zero even when
  # set to 0.5.  Rescaling A_surv (and B_surv consistently) so that
  # mean(A_surv) approaches mean(A_gen) restores alpha as a true mixing knob.
  # The pseudo-observation direction B_surv/A_surv is unchanged by the
  # rescaling; only its weight relative to the genomics term changes.
  #
  # CAP at 100: Full normalisation (scale = m_gen/m_surv) can reach 10^3+
  # when EBeta is small but non-zero (warm-start ~ 0.01).  The corresponding
  # x_surv = EBeta/EBeta2 * (YtWz/YtWY) can then be O(100), swamping x_gen
  # and causing a catastrophic EF update that collapses the L decomposition
  # for the other factors.  Capping at 100 limits the survival-to-genomics
  # signal ratio to at most ~100×, keeping the update stable while still
  # providing meaningful survival guidance.
  #
  # Guard: skip when mean(A_surv) <= 1e-12 (EBeta = 0 exactly) or means are
  # non-finite (can arise from EBeta2 = 0 exactly and Inf in YtWY_diag).
  # ------------------------------------------------------------------
  A_surv_eff <- A_surv
  B_surv_eff <- B_surv
  if (normalize_AB) {
    m_surv <- mean(A_surv)
    m_gen  <- mean(A_gen)
    if (is.finite(m_surv) && is.finite(m_gen) && m_surv > 1e-12 && m_gen > 1e-12) {
      scale_surv <- min(m_gen / m_surv, 100)
      A_surv_eff <- A_surv * scale_surv
      B_surv_eff <- B_surv * scale_surv
    }
  }

  # ------------------------------------------------------------------
  # Alpha-weighted combination and precision floor
  # ------------------------------------------------------------------
  A_F <- pmax((1 - alpha) * A_gen + alpha * A_surv_eff, A_floor)
  B_F <- (1 - alpha) * B_gen + alpha * B_surv_eff

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise (p-vectors)
  # ------------------------------------------------------------------
  x_F <- B_F / A_F
  s_F <- 1.0 / sqrt(A_F)
  x_F[!is.finite(x_F)] <- 0
  s_F[!is.finite(s_F) | s_F <= 0] <- 1e5
  s_F <- pmax(s_F, 1e-8)

  # ------------------------------------------------------------------
  # Solve the p-dimensional EBNM problem
  # ------------------------------------------------------------------
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

#' Update q(F) for all K factors  [Cluster B: dual-source, optional survival]
#'
#' Iterates over k = 1..K, computing R_k and calling update_F_k() for each
#' factor.  Uses Gauss-Seidel ordering.  Survival terms are optional: if w,
#' z, EBeta, EBeta2 are provided, the dual-source update is used; otherwise
#' the pure-genomics fallback is used.
#'
#' When survival args are provided, ZF = Y * EF is computed once from the
#' initial EF (before any Gauss-Seidel updates) and held fixed as the
#' predictor for z_no_k computation throughout the k-loop.  YtWY_diag is
#' also pre-computed once per call.
#'
#' @param Y           n x p matrix: observed genomics data
#' @param EL          n x K matrix: posterior means of loadings (already updated)
#' @param EL2         n x K matrix: posterior second moments of loadings
#' @param EF          p x K matrix: posterior means of factors (will be updated)
#' @param EF2         p x K matrix: posterior second moments of factors
#' @param Tau         p-vector: feature-specific noise precision
#' @param w           n-vector: Cox weights (optional; NULL = no survival)
#' @param z           n-vector: full working response (optional)
#' @param EBeta       K-vector: posterior means of beta_tilde (optional)
#' @param EBeta2      K-vector: posterior 2nd moments of beta_tilde (optional)
#' @param prior_family character: EBNM prior family (default "point_exponential")
#' @param alpha       numeric in [0, 1]: survival mixing weight (default 0.5)
#' @param A_floor     numeric: precision floor (default 1e-10)
#' @param normalize_AB logical: passed through to update_F_k (default FALSE)
#'
#' @return Named list:
#'   $EF      -- p x K matrix of updated posterior means
#'   $EF2     -- p x K matrix of updated posterior second moments
#'   $details -- length-K list, each element is the full update_F_k result
#' @export
#' @family F_update
update_F_all <- function(Y, EL, EL2, EF, EF2, Tau,
                          w          = NULL,
                          z          = NULL,
                          EBeta      = NULL,
                          EBeta2     = NULL,
                          prior_family = "point_exponential",
                          alpha        = 0.5,
                          A_floor      = 1e-10,
                          normalize_AB = FALSE) {

  K        <- ncol(EF)
  EF_curr  <- EF
  EF2_curr <- EF2
  details  <- vector("list", K)

  use_survival <- !is.null(w) && !is.null(z) && !is.null(EBeta) && !is.null(EBeta2)

  if (use_survival) {
    # ZF from the initial EF (held fixed across the k-loop per CAVI contract)
    ZF_init   <- Y %*% EF
    # diag(Y' diag(w) Y): pre-compute once, reuse for all k
    YtWY_diag <- as.vector(t(Y^2) %*% w)
  }

  for (k in seq_len(K)) {
    R_k <- compute_R_k(Y, EL, EF_curr, k)

    if (use_survival) {
      z_no_k    <- compute_z_no_k(z, ZF_init, EBeta, k)
      YtWz_no_k <- as.vector(t(Y) %*% (w * z_no_k))
    } else {
      YtWz_no_k <- NULL
    }

    res_k <- update_F_k(
      Tau          = Tau,
      EL_k         = EL[, k],
      EL2_k        = EL2[, k],
      R_k          = R_k,
      EBeta_k      = if (use_survival) EBeta[k]  else 0,
      EBeta2_k     = if (use_survival) EBeta2[k] else 0,
      YtWY_diag    = if (use_survival) YtWY_diag else NULL,
      YtWz_no_k    = YtWz_no_k,
      prior_family = prior_family,
      alpha        = alpha,
      A_floor      = A_floor,
      normalize_AB = normalize_AB
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
