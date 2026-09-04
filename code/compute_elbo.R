# ============================================================
# Script: compute_elbo.R
# Purpose: Full ELBO computation for the Supervised Bayesian MF model.
#          Provides the two missing terms beyond the genomics proxy:
#          (1) survival likelihood contribution and (2) KL divergences
#          for q(L), q(F), q(beta) extracted from EBNM outputs.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-01
# Dependencies: none (pure R arithmetic; sourced by fit_modular.R)
# ============================================================

# ==============================================================================
# Full ELBO Structure
# ==============================================================================
#
# The complete variational lower bound is:
#
#   ELBO = E_q[log p(Y | L, F, tau)]               -- genomics likelihood
#        + E_q[log PL(t, delta | L, beta)]          -- survival likelihood
#        + sum_k KL_L(k) + sum_k KL_F(k) + sum_k KL_beta(k)
#                                                    -- prior-posterior divergences
#
# where KL_X(k) = E_q[log g_X(theta_k)] - E_q[log q_X(theta_k)]  (NOTE: <= 0)
#
# Term 1 is already computed in update_tau.R as `elbo_proxy`.
# Terms 2-5 are computed here.
#
# Mathematical derivation:
#   See derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf,
#   Section "ELBO" (Equations 189-204).
#
# ==============================================================================

# ------------------------------------------------------------------------------
#' Compute KL Divergence Contribution from One EBNM Update
#'
#' In the CAVI framework each update q(theta_k) is solved as an EBNM problem
#' with pseudo-observations x_i = B_i / A_i and pseudo-noise s_i^2 = 1/A_i.
#' The ebnm package returns the marginal log-likelihood:
#'
#'   L_EBNM = sum_i log integral N(x_i; theta_i, 1/A_i) g(theta_i) dtheta_i
#'
#' This decomposes as:
#'
#'   L_EBNM = E_q[sum_i log N(x_i; theta_i, 1/A_i)]
#'          + E_q[log g(theta)] - E_q[log q(theta)]
#'
#' where the first term (expected pseudo log-likelihood) is computable from
#' returned posterior moments, so:
#'
#'   KL_k = E_q[log g(theta)] - E_q[log q(theta)]
#'        = L_EBNM - E_q[sum_i log N(x_i; theta_i, 1/A_i)]
#'
#' and:
#'
#'   E_q[log N(x_i; theta_i, 1/A_i)]
#'     = 0.5 * log(A_i / (2*pi)) - A_i/2 * (x_i^2 - 2*x_i*mu_i + E[theta_i^2])
#'
#' NOTE: KL_k <= 0 always.  When KL_k = 0 the posterior equals the prior
#' exactly (complete certainty or no data).
#'
#' Works for vectors (L_k: length n; F_k: length p) and scalars (beta_k: length 1).
#'
#' @param ebnm_log_lik  scalar: ebnm_result$log_likelihood from the update call
#' @param A             numeric vector (length n/p/1): EBNM precisions (= 1/s^2)
#' @param x             numeric vector: EBNM pseudo-observations (= B/A)
#' @param mean_q        numeric vector: posterior means E_q[theta_i]
#' @param second_q      numeric vector: posterior 2nd moments E_q[theta_i^2]
#'
#' @return scalar: KL contribution E_q[log g(theta)] - E_q[log q(theta)]  (<= 0)
#'
#' @examples
#' # Scalar beta_k case
#' compute_ebnm_kl(ebnm_log_lik = -0.5, A = 2, x = 1.0, mean_q = 0.8, second_q = 0.64)
#'
#' @seealso fit_modular.R (caller), derivations/MF_UpdateDerivations/
# ------------------------------------------------------------------------------
compute_ebnm_kl <- function(ebnm_log_lik, A, x, mean_q, second_q) {

  # E_q[sum_i log N(x_i; theta_i, 1/A_i)]
  # = sum_i [ 0.5*log(A_i/(2*pi)) - A_i/2*(x_i^2 - 2*x_i*mu_i + E[theta_i^2]) ]
  expected_pseudo_lik <- 0.5 * sum(log(A / (2 * pi))) -
                         0.5 * sum(A * (x^2 - 2 * x * mean_q + second_q))

  # KL = L_EBNM - E_q[pseudo lik]
  kl <- ebnm_log_lik - expected_pseudo_lik

  kl
}

# ------------------------------------------------------------------------------
#' Compute Survival Likelihood ELBO Contribution
#'
#' The Cox partial log-likelihood is non-conjugate, so the CAVI loop uses a
#' second-order Taylor expansion around the current posterior means eta_0:
#'
#'   log PL(eta) ~ log PL(eta_0)
#'                 + u' (eta - eta_0)          [linear: zero in expectation]
#'                 - (1/2) (eta - eta_0)' W (eta - eta_0)
#'
#' where u = score and W = diag(w) = negative diagonal Hessian (positive).
#'
#' Taking expectation under q(L, beta) (with eta_0 = E_q[eta]):
#'
#'   E_q[log PL] ~= logPL(eta_0) - (1/2) sum_i w_i * Var_q(eta_i)
#'
#' The linear term vanishes because u is evaluated at eta_0 = E_q[eta].
#'
#' Under the mean-field factorization (l_ik, beta_k independent across k and i):
#'
#'   Var_q(eta_i) = E_q[eta_i^2] - (E_q[eta_i])^2
#'                = sum_k E_q[l_ik^2] * E_q[beta_k^2]
#'                  + (E_q[eta_i])^2
#'                  - sum_k E_q[l_ik]^2 * E_q[beta_k]^2
#'                  - (E_q[eta_i])^2
#'                = sum_k EL2[i,k] * EBeta2[k] - sum_k EL[i,k]^2 * EBeta[k]^2
#'
#' (The (E_q[eta_i])^2 terms cancel; the residual is the sum of per-factor
#' variance contributions, i.e. Var_q(l_ik * beta_k) = EL2[i,k]*EBeta2[k]
#' - EL[i,k]^2 * EBeta[k]^2 summed over k.)
#'
#' When all posterior variances are zero (EL2 = EL^2, EBeta2 = EBeta^2) this
#' returns logPL exactly — the uncertainty correction vanishes.
#'
#' @param logPL   scalar: Cox partial log-likelihood at current posterior means,
#'                computed by calc_cox_taylor()
#' @param w       n-vector: Cox diagonal Hessian (negative, positive values)
#' @param EL      n x K matrix: posterior means E_q[l_{ik}]
#' @param EL2     n x K matrix: posterior 2nd moments E_q[l_{ik}^2]
#' @param EBeta   K-vector: posterior means E_q[beta_k]. When `cohort_idx` is
#'                supplied, a K x C matrix instead (one column per cohort;
#'                see `code/update_beta_cohort.R`).
#' @param EBeta2  K-vector (or K x C matrix, matching `EBeta`): posterior 2nd
#'                moments E_q[beta_k^2]
#' @param cohort_idx NULL (default, preserves the shared-beta path exactly),
#'                or an n-vector of integers in 1..C giving each patient's
#'                cohort index into `EBeta`/`EBeta2`'s columns. Each
#'                patient's own cohort column is used when forming
#'                `Var_q(eta_i)`, since `beta_k^{c(i)}` varies by patient.
#'
#' @return scalar: E_q[log PL(t, delta | L, beta)] (approximate)
#'
#' @examples
#' # Zero-variance case: should return logPL exactly
#' EL <- matrix(c(1, 2), 2, 1); EBeta <- 0.5
#' compute_survival_elbo(-3.0, w = c(0.1, 0.2), EL, EL^2, EBeta, EBeta^2)
#'
#' @seealso fit_modular.R (caller), derivations/MF_UpdateDerivations/
# ------------------------------------------------------------------------------
compute_survival_elbo <- function(logPL, w, EL, EL2, EBeta, EBeta2, cohort_idx = NULL) {

  if (is.null(cohort_idx)) {
    # Var_q(eta_i) = sum_k [EL2[i,k]*EBeta2[k] - EL[i,k]^2 * EBeta[k]^2]
    # sweep(EL2, 2, EBeta2, "*"):  each col k of EL2 multiplied by EBeta2[k]
    var_eta <- rowSums(sweep(EL2, 2, EBeta2, "*")) -
               rowSums(sweep(EL^2, 2, EBeta^2, "*"))
  } else {
    # Cohort-specific beta: each patient i uses column cohort_idx[i] of the
    # K x C matrices EBeta/EBeta2. t(EBeta2)[cohort_idx, ] broadcasts each
    # patient's own cohort's EBeta2 row (length K) into an n x K matrix in
    # one indexing operation (same pattern as compute_z_no_k_cohort()).
    var_eta <- rowSums(EL2 * t(EBeta2)[cohort_idx, , drop = FALSE]) -
               rowSums(EL^2 * t(EBeta^2)[cohort_idx, , drop = FALSE])
  }

  # Uncertainty correction: subtract (1/2) * sum_i w_i * Var_q(eta_i)
  logPL - 0.5 * sum(w * var_eta)
}

# ==============================================================================
# KL divergence for cohort f_c — Normal posterior vs Normal prior
# ==============================================================================

#' KL divergence for q(f_c) vs N(0, sigma^2) prior
#'
#' Computes the negative KL divergence (i.e. the ELBO contribution from the
#' cohort F rows) in closed form.
#'
#' KL[N(mu, 1/A) || N(0, sigma^2)] = (1/2) * (E[f^2]/sigma^2 + log(A*sigma^2) - 1)
#'
#' SIGN CONVENTION: Returns a NEGATIVE scalar, matching the compute_ebnm_kl()
#' convention (KL terms are subtracted from the ELBO by adding a negative value).
#' The caller uses:
#'   elbo_full <- elbo_full + compute_normal_kl(...)
#' which correctly penalises the ELBO.  NEVER write:
#'   elbo_full <- elbo_full - compute_normal_kl(...)
#' that double-negates and would increase the ELBO with stronger divergence.
#'
#' @param EF_cohort      p x (C-1) matrix: posterior means E[f_cj]
#' @param EF2_cohort     p x (C-1) matrix: posterior second moments E[f_cj^2]
#'                       Must satisfy EF2_cohort > EF_cohort^2 element-wise.
#' @param sigma_F_cohort scalar: prior SD (must match the value used in update_F_cohort_col)
#'
#' @return scalar <= 0: negative total KL summed over all genes and cohort columns
#'
#' @seealso update_F_cohort_all, fit_modular.R
compute_normal_kl <- function(EF_cohort, EF2_cohort, sigma_F_cohort) {
  sigma2  <- sigma_F_cohort^2
  # Posterior precision recovered from second moment and mean:  A = 1/Var = 1/(EF2 - EF^2)
  A       <- 1.0 / (EF2_cohort - EF_cohort^2)          # p x (C-1), all > 0
  kl_mat  <- 0.5 * (EF2_cohort / sigma2 + log(A * sigma2) - 1.0)  # p x (C-1), all >= 0
  -sum(kl_mat)    # negative: KL >= 0, so -KL <= 0
}
