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
#' @param EBeta   K-vector: posterior means E_q[beta_k]
#' @param EBeta2  K-vector: posterior 2nd moments E_q[beta_k^2]
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
compute_survival_elbo <- function(logPL, w, EL, EL2, EBeta, EBeta2) {

  # Var_q(eta_i) = sum_k [EL2[i,k]*EBeta2[k] - EL[i,k]^2 * EBeta[k]^2]
  # sweep(EL2, 2, EBeta2, "*"):  each col k of EL2 multiplied by EBeta2[k]
  var_eta <- rowSums(sweep(EL2, 2, EBeta2, "*")) -
             rowSums(sweep(EL^2, 2, EBeta^2, "*"))

  # Uncertainty correction: subtract (1/2) * sum_i w_i * Var_q(eta_i)
  logPL - 0.5 * sum(w * var_eta)
}
