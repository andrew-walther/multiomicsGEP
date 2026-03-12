# =============================================================================
# code/update_tau.R
#
# Modular q(tau_j) variational update for Supervised Bayesian MF.
#
# DERIVATION REFERENCE:
#   derivations/qTau/qTau_update_derivation.tex  (self-contained)
#   derivations/MF_UpdateDerivations/MF_V2_Companion.tex  Section 7
#   derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.tex
#
# MATHEMATICAL SUMMARY:
#   tau_j is the feature-specific noise precision.  Unlike L, F, and beta,
#   the tau update is a CLOSED-FORM MLE (not an EBNM problem):
#
#     R2_bar_{ij} = (Y_{ij} - sum_k l_bar_{ik} f_bar_{jk})^2
#                 + sum_k [ E[l^2_{ik}]*E[f^2_{jk}] - l_bar^2_{ik}*f_bar^2_{jk} ]
#                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#                   Variance correction term (prevents underestimating noise)
#
#     tau_hat_j = n / sum_i R2_bar_{ij}      (column-specific MLE)
#
#   The variance correction ensures that posterior uncertainty in L and F
#   is properly accounted for when estimating noise — without it, tau
#   would be systematically overestimated (noise underestimated).
#
# DESIGN PATTERN:
#   Part of the modular update series:
#     update_beta.R  ->  update_L.R  ->  update_F.R  ->  update_tau.R
#   Unlike the others, tau has NO per-k loop (all K factors at once)
#   and NO ebnm dependency.
#
# DEPENDENCIES:
#   None (base R only)
# =============================================================================

# =============================================================================
# Helper: variance correction term
# =============================================================================

#' Compute the posterior variance correction matrix
#'
#' The variance correction accounts for posterior uncertainty in L and F
#' when computing expected squared residuals:
#'
#'   Var_Term_{ij} = sum_k [ E[l^2_{ik}] * E[f^2_{jk}] - l_bar_{ik}^2 * f_bar_{jk}^2 ]
#'
#' Equivalently:
#'   Var_Term = (EL2 %*% t(EF2)) - ((EL)^2 %*% t((EF)^2))
#'
#' Each entry is >= 0 because E[X^2] >= (E[X])^2.  The Var_Term is zero
#' only when EL2 = EL^2 AND EF2 = EF^2 (zero posterior variance).
#'
#' @param EL   n x K matrix: posterior means of loadings
#' @param EL2  n x K matrix: posterior second moments of loadings
#' @param EF   p x K matrix: posterior means of factors
#' @param EF2  p x K matrix: posterior second moments of factors
#' @return n x p matrix: variance correction (non-negative)
compute_var_term <- function(EL, EL2, EF, EF2) {
  (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
}

# =============================================================================
# Helper: expected squared residual
# =============================================================================

#' Compute the expected squared residual matrix R2_bar
#'
#' The expected squared residual under the variational posterior:
#'
#'   R2_bar_{ij} = (Y_{ij} - sum_k l_bar_{ik} * f_bar_{jk})^2 + Var_Term_{ij}
#'
#' The first term is the squared residual at the posterior means.
#' The second term (variance correction) inflates this to account for
#' posterior uncertainty in L and F.
#'
#' @param Y    n x p matrix: observed genomics data
#' @param EL   n x K matrix: posterior means of loadings
#' @param EL2  n x K matrix: posterior second moments of loadings
#' @param EF   p x K matrix: posterior means of factors
#' @param EF2  p x K matrix: posterior second moments of factors
#' @return n x p matrix: expected squared residuals (non-negative)
compute_expected_residual_sq <- function(Y, EL, EL2, EF, EF2) {
  mean_resid_sq <- (Y - EL %*% t(EF))^2             # n x p: (Y - prediction)^2
  Var_Term      <- compute_var_term(EL, EL2, EF, EF2)  # n x p: variance correction
  mean_resid_sq + Var_Term
}

# =============================================================================
# Core: precision update
# =============================================================================

#' Update feature-specific noise precision tau_j (column-wise MLE)
#'
#' Computes the maximum-likelihood estimate of each tau_j from the
#' expected squared residuals under the variational posterior.  Unlike
#' the L, F, and beta updates, this is NOT an EBNM problem — it is a
#' closed-form MLE for each feature independently.
#'
#' Also computes the genomics ELBO proxy for convergence monitoring:
#'   ELBO_proxy = sum_j [ n/2 * log(tau_j) - tau_j/2 * sum_i R2_bar_{ij} ]
#'
#' @param Y          n x p matrix: observed genomics data
#' @param EL         n x K matrix: posterior means of loadings
#' @param EL2        n x K matrix: posterior second moments of loadings
#' @param EF         p x K matrix: posterior means of factors
#' @param EF2        p x K matrix: posterior second moments of factors
#' @param tau_floor  numeric: minimum denominator = n * tau_floor.
#'                   Prevents tau from becoming infinite when residuals are
#'                   near zero.  Default 1e-8 matches V2.R [A3].
#'
#' @return Named list:
#'   $Tau         -- p-vector: updated noise precision per feature
#'   $R2_bar      -- n x p matrix: expected squared residuals
#'   $Var_Term    -- n x p matrix: variance correction term
#'   $elbo_proxy  -- scalar: genomics ELBO proxy = E_q[log P(Y|L,F,tau)]
#'
#' @examples
#' set.seed(1); n <- 50; p <- 100; K <- 3
#' EL  <- matrix(rnorm(n * K), n, K)
#' EL2 <- EL^2 + 0.1
#' EF  <- matrix(rnorm(p * K), p, K)
#' EF2 <- EF^2 + 0.1
#' Y   <- EL %*% t(EF) + matrix(rnorm(n * p, sd = 0.5), n, p)
#' res <- update_tau(Y, EL, EL2, EF, EF2)
#' cat("Tau (first 5):", round(res$Tau[1:5], 2), "\n")
#' cat("ELBO proxy:", round(res$elbo_proxy, 1), "\n")
update_tau <- function(Y, EL, EL2, EF, EF2,
                        tau_floor = 1e-8) {

  n <- nrow(Y)

  # ------------------------------------------------------------------
  # Variance correction (n x p)
  # Accounts for posterior uncertainty: EL2*EF2 - EL^2*EF^2 >= 0
  # Without this, tau is systematically overestimated (noise underestimated).
  # See REVISED.tex correction R7 (sign error) and R8 (subscript fix).
  # ------------------------------------------------------------------
  Var_Term <- compute_var_term(EL, EL2, EF, EF2)

  # ------------------------------------------------------------------
  # Expected squared residual (n x p)
  # R2_bar = (Y - mean prediction)^2 + variance inflation
  # ------------------------------------------------------------------
  R2_bar <- (Y - EL %*% t(EF))^2 + Var_Term

  # ------------------------------------------------------------------
  # Column-specific MLE  (p-vector)
  #   tau_j = n / sum_i R2_bar_{ij}
  # Floor prevents division by zero when residuals are very small.
  # ------------------------------------------------------------------
  col_sums <- colSums(R2_bar)                         # p-vector
  Tau <- n / pmax(col_sums, n * tau_floor)            # p-vector [A3]

  # ------------------------------------------------------------------
  # ELBO proxy: genomics log-likelihood term
  #   E_q[log P(Y | L, F, tau)] = sum_j [n/2 log(tau_j) - tau_j/2 sum_i R2_bar_{ij}]
  # ------------------------------------------------------------------
  elbo_proxy <- sum(n / 2 * log(Tau) - Tau / 2 * col_sums)

  list(
    Tau        = Tau,
    R2_bar     = R2_bar,
    Var_Term   = Var_Term,
    elbo_proxy = elbo_proxy
  )
}
