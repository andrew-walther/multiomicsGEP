# =============================================================================
# code/update_L.R
#
# Modular q(l_{ik}) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE:
#   derivations/qL/qL_update_derivation.tex  (self-contained)
#   derivations/MF_UpdateDerivations/MF_V2_Companion.tex  Section 4
#   derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.tex
#
# MATHEMATICAL SUMMARY:
#   l_{ik} appears in BOTH the genomics and survival likelihoods.  The
#   coordinate-ascent update reduces to an n-dimensional EBNM problem:
#
#     A_L[i] = sum_j(tau_j * E[f^2_{jk}])  +  W_ii * E[beta_k^2]
#              [scalar genomics term]           [n-vector survival term]
#
#     B_L[i] = sum_j(tau_j * R^{-k}_{ij} * f_bar_{jk})   (genomics)
#            + W_ii * z^{-k}_i * beta_bar_k                (survival)
#
#     EBNM(x = B_L/A_L, s = 1/sqrt(A_L))   <- VECTOR call (n observations)
#
#   Unlike the beta update (scalar EBNM), the L update is a VECTOR EBNM
#   because A_L and B_L vary across samples i (due to Cox weights W_ii).
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

#' Update q(l_{ik}) for a single factor k
#'
#' Performs the EBNM-based coordinate-ascent update for the patient loading
#' column l_{.k} given BOTH genomics and survival data.  Unlike the beta
#' update (scalar EBNM), this is a VECTOR EBNM with n observations.
#'
#' The function is intentionally decoupled from the CAVI loop: it takes
#' pre-computed matrices/vectors and returns a self-contained result list.
#'
#' @param Tau         p-vector: feature-specific noise precision tau_j (> 0)
#' @param EF_k        p-vector: posterior mean of factor column k (f_bar_{jk})
#' @param EF2_k       p-vector: posterior SECOND MOMENT of factor column k
#'                    (E_q[f_{jk}^2] = Var + mean^2).  Must satisfy EF2_k >= EF_k^2.
#' @param w           n-vector: Cox neg-diagonal Hessian weights W_{ii} (>= 0)
#' @param EBeta_k     scalar: posterior mean of survival coefficient beta_k
#' @param EBeta2_k    scalar: posterior second moment E_q[beta_k^2]
#' @param R_k         n x p matrix: partial residual R^{-k} (from compute_R_k)
#' @param z_no_k      n-vector: partial working response z^{-k}_i
#' @param prior_family character: EBNM prior family (default "point_normal")
#' @param alpha       numeric in [0, 1]: mixing weight on the survival term.
#'                    A_L = pmax((1-alpha)*A_gen + alpha*A_surv, A_floor)
#'                    B_L = (1-alpha)*B_gen + alpha*B_surv
#'                    alpha=0 gives pure-genomics EBNM (survival contribution zeroed);
#'                    alpha=1 gives pure-survival EBNM (genomics contribution zeroed).
#'                    Default 0.5. Note: $B_gen and $B_surv in the return list are
#'                    always the unweighted raw components (before alpha scaling),
#'                    so $B = (1-alpha)*$B_gen + alpha*$B_surv.
#' @param A_floor     numeric: minimum value for each A_L[i] (default 1e-10)
#'
#' @return Named list:
#'   $mean        -- n-vector: posterior mean E_q[l_{ik}]
#'   $second      -- n-vector: posterior 2nd moment E_q[l_{ik}^2]
#'   $sd          -- n-vector: posterior SD sqrt(Var_q(l_{ik}))
#'   $A           -- n-vector: precision A_{ik} [floored]
#'   $B           -- n-vector: total signal B_{ik} = (1-alpha)*B_gen + alpha*B_surv
#'   $B_gen       -- n-vector: raw (unweighted) genomics component of B
#'   $B_surv      -- n-vector: raw (unweighted) survival component of B
#'   $x           -- n-vector: EBNM pseudo-obs x_i = B_i/A_i
#'   $s           -- n-vector: EBNM pseudo-noise s_i = 1/sqrt(A_i)
#'   $ebnm_result -- raw ebnm() return object
#'
#' @export
#' @family L_update
#' @examples
#' library(ebnm)
#' set.seed(1); n <- 100; p <- 50
#' Tau     <- rep(2.0, p)
#' EF_k   <- rnorm(p);  EF2_k <- EF_k^2 + 0.1
#' w      <- rep(1.0, n)
#' EBeta_k <- 0.5;  EBeta2_k <- 0.5^2 + 0.05
#' R_k    <- matrix(rnorm(n * p), n, p)
#' z_no_k <- rnorm(n)
#' res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
#' cat("Loading column estimates (first 5):", round(res$mean[1:5], 3), "\n")
update_L_k <- function(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k,
                        R_k, z_no_k,
                        prior_family = "point_exponential",
                        alpha        = 0.5,
                        lambda       = 1.0,
                        A_floor      = 1e-10,
                        normalize_AB = FALSE) {

  # ------------------------------------------------------------------
  # Precision A_{ik}  (n-vector)
  #
  #   Genomics:  (1-alpha) * sum_j(tau_j * E[f^2_{jk}])  [scalar, same for all i]
  #   Survival:  alpha * W_{ii} * E[beta_k^2]             [n-vector, sample-specific]
  #
  # The genomics term uses E[f^2] (not f_bar^2) for error-in-variables
  # correction: posterior uncertainty in F inflates the effective noise
  # for L, preventing overfitting to uncertain factor estimates.
  # alpha controls the trade-off: alpha=0 is pure genomics, alpha=1 is pure survival.
  # ------------------------------------------------------------------
  # A_gen is a SCALAR (constant across all patients) because the genomics
  # precision sum_j(tau_j * EF2_jk) does not depend on sample i.
  A_gen  <- sum(Tau * EF2_k)                                     # scalar
  # A_surv is an n-VECTOR because Cox weights W_ii differ per patient.
  # This is why L requires a vector EBNM (unlike beta's scalar EBNM).
  A_surv <- lambda * w * EBeta2_k                                # λ-scaled n-vector

  # ------------------------------------------------------------------
  # Signal B_{ik}  (n-vector)
  #
  #   Genomics:  (1-alpha) * (R_k %*% (Tau * EF_k))[i]   [matrix-vector product]
  #   Survival:  alpha * W_{ii} * z^{-k}_i * beta_bar_k  [element-wise]
  #
  # The efficient form R_k %*% (Tau * EF_k) avoids an n x p sweep.
  # NOTE: $B_gen and $B_surv in the return list are the RAW (unweighted)
  # components — before alpha scaling — to preserve diagnostic utility.
  # The combined $B reflects the weighting: (1-alpha)*B_gen + alpha*B_surv.
  # ------------------------------------------------------------------
  B_gen  <- as.vector(R_k %*% (Tau * EF_k))                     # n-vector (raw)
  B_surv <- lambda * w * z_no_k * EBeta_k                        # λ-scaled n-vector (raw)

  # ------------------------------------------------------------------
  # Optional rescaling [Fix 4 of docs/beta_zero_fix_design.md §4.8]
  #
  # Without rescaling, A_gen ≈ sum_j(τ_j · E[f_jk²]) ~ p · scale, while
  # A_surv ~ O(1) per subject — a structural ~p× imbalance that makes the
  # survival term negligible regardless of EBeta. With normalize_AB = TRUE,
  # rescale A_surv (and B_surv consistently) so its typical magnitude
  # matches A_gen. This keeps the EBNM noise interpretation intact
  # (s_L = 1/sqrt(A_L) stays in original units, avoiding over-shrinkage)
  # while making α a meaningful fraction-of-influence knob between the
  # two sources.
  #
  # Departs from strict ELBO maximization — verify ELBO is monotone
  # empirically. When mean(A_surv) is at or below the floor (e.g., when
  # EBeta_k ≈ 0 at init), the rescale is skipped to avoid amplifying
  # numerical noise.
  # ------------------------------------------------------------------
  A_surv_eff <- A_surv
  B_surv_eff <- B_surv
  if (normalize_AB) {
    m_surv <- mean(A_surv)
    if (m_surv > 1e-12 && A_gen > 1e-12) {
      scale_surv <- A_gen / m_surv
      A_surv_eff <- A_surv * scale_surv
      B_surv_eff <- B_surv * scale_surv
    }
    # else: skip — survival term is effectively zero this iteration; let
    # the genomics term drive the update unmodified.
  }
  A_L <- pmax((1 - alpha) * A_gen + alpha * A_surv_eff, A_floor)    # n-vector [A3]
  B_L <- (1 - alpha) * B_gen + alpha * B_surv_eff                   # weighted combination

  # ------------------------------------------------------------------
  # EBNM pseudo-observation and noise (n-vectors)
  # ------------------------------------------------------------------
  x_L <- B_L / A_L
  s_L <- 1.0 / sqrt(A_L)
  x_L[!is.finite(x_L)] <- 0
  s_L[!is.finite(s_L) | s_L <= 0] <- 1e5
  s_L <- pmax(s_L, 1e-8)

  # ------------------------------------------------------------------
  # Solve the n-dimensional EBNM problem
  # Each sample i provides one observation x_i with noise s_i.
  # ------------------------------------------------------------------
  res <- ebnm(x = x_L, s = s_L, prior_family = prior_family)

  l_mean   <- res$posterior$mean                     # n-vector
  l_sd     <- res$posterior$sd                       # n-vector
  l_second <- l_sd^2 + l_mean^2                     # n-vector

  # Pack all diagnostics including both B_gen and B_surv components,
  # enabling tests and demos to inspect the dual-source contributions.
  list(
    mean        = l_mean,
    second      = l_second,
    sd          = l_sd,
    A           = A_L,
    B           = B_L,
    B_gen       = B_gen,
    B_surv      = B_surv,
    x           = x_L,
    s           = s_L,
    ebnm_result = res
  )
}

# =============================================================================
# Convenience wrapper: full Gauss-Seidel loop over all K factors
# =============================================================================

#' Update q(L) for all K factors (Gauss-Seidel CAVI loop)
#'
#' Iterates over k = 1..K, computing the partial residual R_k and partial
#' working response z_no_k, then calling update_L_k() for each factor.
#' Uses Gauss-Seidel ordering: once EL[,k] is updated, the new values
#' are used for subsequent factors k' > k.
#'
#' INTEGRATION NOTE:
#'   In the full CAVI loop (V2.R), the L update for factor k is immediately
#'   followed by the F update and beta update for the same k.  This _all
#'   wrapper updates ALL K loading columns before returning.  For standalone
#'   testing and demos; the CAVI loop uses update_L_k() directly.
#'
#' @param Y           n x p matrix: observed genomics data
#' @param EL          n x K matrix: posterior means of loadings (will be updated)
#' @param EL2         n x K matrix: posterior second moments of loadings
#' @param EF          p x K matrix: posterior means of factors
#' @param EF2         p x K matrix: posterior second moments of factors
#' @param Tau         p-vector: feature-specific noise precision
#' @param w           n-vector: Cox neg-diagonal Hessian weights
#' @param z           n-vector: full working response
#' @param EBeta       K-vector: posterior means of survival coefficients
#' @param EBeta2      K-vector: posterior second moments of survival coefficients
#' @param prior_family character: EBNM prior family (default "point_normal")
#' @param alpha       numeric in [0, 1]: survival mixing weight, passed to
#'                    update_L_k() for each factor (default 0.5)
#' @param A_floor     numeric: precision floor (default 1e-10)
#'
#' @return Named list:
#'   $EL      -- n x K matrix of updated posterior means
#'   $EL2     -- n x K matrix of updated posterior second moments
#'   $details -- length-K list, each element is the full update_L_k result
#' @export
#' @family L_update
update_L_all <- function(Y, EL, EL2, EF, EF2, Tau, w, z,
                          EBeta, EBeta2,
                          prior_family = "point_exponential",
                          alpha        = 0.5,
                          lambda       = 1.0,
                          A_floor      = 1e-10) {

  n <- nrow(EL)
  K <- ncol(EL)
  EL_curr  <- EL                   # mutable copy for Gauss-Seidel
  EL2_curr <- EL2
  details  <- vector("list", K)

  for (k in seq_len(K)) {
    # Partial residual R_k (uses current Gauss-Seidel EL)
    R_k <- compute_R_k(Y, EL_curr, EF, k)

    # Partial working response z^{-k} (uses current EBeta — not updated here).
    # Inlined rather than calling compute_z_no_k() from update_beta.R to avoid
    # forcing a source-order dependency between modules. This keeps update_L.R
    # self-contained for standalone use and testing.
    eta_no_k <- as.vector(EL_curr %*% EBeta) - EL_curr[, k] * EBeta[k]
    z_no_k   <- z - eta_no_k

    res_k <- update_L_k(
      Tau        = Tau,
      EF_k       = EF[, k],
      EF2_k      = EF2[, k],
      w          = w,
      EBeta_k    = EBeta[k],
      EBeta2_k   = EBeta2[k],
      R_k        = R_k,
      z_no_k     = z_no_k,
      prior_family = prior_family,
      alpha      = alpha,
      lambda     = lambda,
      A_floor    = A_floor
    )

    # Gauss-Seidel: update EL columns so next k uses fresh values
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
