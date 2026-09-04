# =============================================================================
# code/update_beta_cohort.R
#
# Cohort-specific survival coefficients beta_k^(c) for the Cox-on-YF (YFB)
# model: eta_i = sum_k (Y.F)_{ik} * beta_k^{c(i)}, one coefficient per
# (factor, cohort) pair instead of one shared beta_k per factor.
#
# RELATIONSHIP TO THE TWO EXISTING COHORT-AWARE EXTENSIONS (DECISIONS.md
# 2026-05-22, 2026-07-15) -- distinct, composable, NOT a replacement for
# either:
#   - `cohort_id` (code/update_F_cohort.R): fixed indicator columns
#     appended to L/F, absorbing per-gene PLATFORM offsets in the
#     GENOMICS reconstruction. By construction beta_cohort = 0 for those
#     columns -- cohort membership never reaches the survival term this
#     way. Unchanged by this file.
#   - `strata_id` (code/fit_cox_on_yf.R's calc_cox_taylor_yf): stratified
#     Cox partial likelihood, i.e. a cohort-specific BASELINE HAZARD, with
#     beta still shared across strata. Unchanged by this file.
#   - `beta_cohort_id` (this file): cohort-specific SURVIVAL COEFFICIENTS.
#     The only one of the three that lets a factor be prognostic in one
#     cohort and not another. In practice all three often share the same
#     underlying grouping vector (e.g. dataset_labels), but they are
#     different model components and conflating their arguments in
#     fit_cox_on_yf() would be a mistake -- hence a third, independent
#     argument rather than overloading cohort_id or strata_id.
#
# MATHEMATICAL SUMMARY:
#   Per factor k, per cohort c: A_k^(c) = alpha * sum_{i in c} w_i * ZF_{ik}^2
#                                B_k^(c) = alpha * sum_{i in c} w_i * z^{-k}_i * ZF_{ik}
#   These are EXACTLY update_beta_k()'s A_k/B_k formulas (code/update_beta.R),
#   restricted to the patient index subset for cohort c -- update_beta_k()
#   itself is not called and not edited; this file only reproduces its
#   A/B arithmetic per cohort.
#
#   PARTIAL POOLING (the point of this extension): rather than C
#   independent scalar ebnm() calls per factor (no sharing of information
#   across cohorts), the C pairs (x_k^{(1..C)}, s_k^{(1..C)}) are stacked
#   into one VECTOR ebnm() call per factor. The C cohort betas then share
#   one learned prior g_hat_k and shrink toward a common value when the
#   data do not strongly contradict that -- the "consensus factors are
#   more reproducible" idea, made operational. compute_ebnm_kl()
#   (code/compute_elbo.R) is length-generic and needs no change: kl_beta[k]
#   is still a single scalar per factor, summed over all C cohort
#   contributions inside compute_ebnm_kl()'s own sum().
#
# DESIGN PATTERN:
#   Sibling to update_beta_all() (code/update_beta.R), which
#   fit_cox_on_yf() does not call (it inlines the loop) -- so
#   update_beta_all() is untouched and this file adds a parallel path
#   used only when beta_cohort_id is supplied.
#
# DEPENDENCIES:
#   ebnm  (CRAN)
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# =============================================================================
# Helper: per-patient partial working response under cohort-specific beta
# =============================================================================

#' Compute the partial working response z^{-k}_i when beta is cohort-specific
#'
#' Cohort-specific analogue of compute_z_no_k() (code/update_beta.R), which
#' assumes one shared EBeta (a length-K vector) and cannot express a
#' per-patient lookup into a different cohort's coefficient column.
#'
#'   eta_i = sum_{k'} ZF_{i,k'} * EBeta_mat[k', cohort_idx_i]
#'   z^{-k}_i = z_i - eta_i + ZF_{i,k} * EBeta_mat[k, cohort_idx_i]
#'
#' `t(EBeta_mat)[cohort_idx, ]` broadcasts each patient's own cohort's
#' coefficient row across all K factors in one indexing operation (an
#' n x K matrix), so `rowSums(ZF * that)` gives eta_i for all patients at
#' once without an explicit per-patient loop.
#'
#' @param z           n-vector: full working response
#' @param ZF          n x K matrix: observed projection scores (Y . EF_norm)
#' @param EBeta_mat   K x C matrix: current posterior means, one column per cohort
#' @param cohort_idx  n-vector of integers in 1..C: each patient's cohort index
#' @param k           integer: factor index to exclude (1-based)
#' @return n-vector: partial working response for factor k
#' @export
#' @family beta_cohort_update
compute_z_no_k_cohort <- function(z, ZF, EBeta_mat, cohort_idx, k) {
  eta_full <- rowSums(ZF * t(EBeta_mat)[cohort_idx, , drop = FALSE])
  eta_no_k <- eta_full - ZF[, k] * EBeta_mat[k, cohort_idx]
  z - eta_no_k
}

# =============================================================================
# Core: single-factor update across all C cohorts, one vectorized ebnm() call
# =============================================================================

#' Update q(beta_k^{(1..C)}) for a single factor k across all cohorts
#'
#' Computes A_k^(c), B_k^(c) for each cohort c (update_beta_k()'s formulas,
#' restricted to cohort c's patient indices), then makes ONE vectorized
#' ebnm() call across the C (x, s) pairs so the C cohort coefficients
#' partially pool through a shared prior. See the file header for why this
#' is not C independent update_beta_k() calls.
#'
#' @param w           n-vector: Cox neg-diagonal Hessian weights W_{ii} (> 0)
#' @param z_no_k      n-vector: partial working response, from
#'                    compute_z_no_k_cohort()
#' @param ZF_k        n-vector: observed projection scores for factor k
#'                    (Y . EF_norm)[, k]
#' @param ZF2_k       n-vector: ZF_k^2 (ZF is observed, so its "second
#'                    moment" is its squared value -- same convention as
#'                    the non-cohort YFB beta update)
#' @param cohort_idx  n-vector of integers in 1..C: each patient's cohort index
#' @param C           integer: number of cohorts
#' @param prior_family character: EBNM prior family (default "normal",
#'                    matching fit_cox_on_yf()'s default prior_beta)
#' @param alpha       numeric in [0, 1]: survival mixing weight (default 0.5)
#' @param A_floor     numeric: minimum value for each A_k^(c) (default 1e-10)
#' @param survival_divisor numeric: divides each A_k^(c) and B_k^(c) (default 1)
#'
#' @return Named list:
#'   $mean    -- C-vector: posterior mean E_q[beta_k^(c)] per cohort
#'   $second  -- C-vector: posterior 2nd moment E_q[beta_k^(c)^2] per cohort
#'   $sd      -- C-vector: posterior SD per cohort
#'   $A       -- C-vector: precision A_k^(c) [floored]
#'   $B       -- C-vector: signal B_k^(c)
#'   $x       -- C-vector: EBNM pseudo-obs x_k^(c) = B_k^(c)/A_k^(c)
#'   $s       -- C-vector: EBNM pseudo-noise s_k^(c) = 1/sqrt(A_k^(c))
#'   $ebnm_result -- raw ebnm() return object (ONE call, length C)
#'
#' @export
#' @family beta_cohort_update
#' @seealso \code{\link{update_beta_k}} for the non-cohort, shared-beta update
update_beta_cohort_k <- function(w, z_no_k, ZF_k, ZF2_k, cohort_idx, C,
                                  prior_family = "normal",
                                  alpha        = 0.5,
                                  A_floor      = 1e-10,
                                  survival_divisor = 1) {

  A_vec <- numeric(C)
  B_vec <- numeric(C)
  for (c in seq_len(C)) {
    idx_c <- which(cohort_idx == c)
    # Identical formulas to update_beta_k()'s A_k/B_k, restricted to cohort c.
    A_vec[c] <- max(alpha * sum(w[idx_c] * ZF2_k[idx_c]) / survival_divisor, A_floor)
    B_vec[c] <- alpha * sum(w[idx_c] * z_no_k[idx_c] * ZF_k[idx_c]) / survival_divisor
  }

  x_vec <- B_vec / A_vec
  s_vec <- 1.0 / sqrt(A_vec)
  x_vec[!is.finite(x_vec)] <- 0
  s_vec[!is.finite(s_vec) | s_vec <= 0] <- 1e5
  s_vec <- pmax(s_vec, 1e-8)

  # ONE vectorized ebnm() call across all C cohorts -- partial pooling.
  res <- ebnm(x = x_vec, s = s_vec, prior_family = prior_family)

  beta_mean   <- res$posterior$mean
  beta_sd     <- res$posterior$sd
  beta_second <- beta_sd^2 + beta_mean^2

  list(
    mean        = beta_mean,
    second      = beta_second,
    sd          = beta_sd,
    A           = A_vec,
    B           = B_vec,
    x           = x_vec,
    s           = s_vec,
    ebnm_result = res
  )
}

# =============================================================================
# Convenience wrapper: full Gauss-Seidel loop over all K factors
# =============================================================================

#' Update q(beta) for all K factors, cohort-specific (Gauss-Seidel CAVI loop)
#'
#' Sibling to update_beta_all() (code/update_beta.R), which fit_cox_on_yf()
#' does not call (it inlines the loop) -- update_beta_all() is therefore
#' untouched by this addition.
#'
#' @param w           n-vector: Cox neg-diagonal Hessian weights
#' @param z           n-vector: full working response
#' @param ZF          n x K matrix: observed projection scores (Y . EF_norm)
#' @param EBeta_mat   K x C matrix: current posterior means (warm start;
#'                    updated in-place within the loop, Gauss-Seidel)
#' @param cohort_idx  n-vector of integers in 1..C
#' @param C           integer: number of cohorts
#' @param prior_family,alpha,A_floor,survival_divisor: forwarded to
#'   update_beta_cohort_k() for each factor
#'
#' @return Named list:
#'   $EBeta   -- K x C matrix of updated posterior means
#'   $EBeta2  -- K x C matrix of updated posterior second moments
#'   $details -- length-K list, each element the full update_beta_cohort_k() result
#' @export
#' @family beta_cohort_update
update_beta_cohort_all <- function(w, z, ZF, EBeta_mat, cohort_idx, C,
                                    prior_family = "normal",
                                    alpha        = 0.5,
                                    A_floor      = 1e-10,
                                    survival_divisor = 1) {

  K               <- ncol(ZF)
  EBeta_mat_curr  <- EBeta_mat
  EBeta2_mat_new  <- matrix(0, K, C)
  details         <- vector("list", K)

  for (k in seq_len(K)) {
    z_no_k <- compute_z_no_k_cohort(z, ZF, EBeta_mat_curr, cohort_idx, k)

    res_k <- update_beta_cohort_k(
      w = w, z_no_k = z_no_k, ZF_k = ZF[, k], ZF2_k = ZF[, k]^2,
      cohort_idx = cohort_idx, C = C,
      prior_family = prior_family, alpha = alpha,
      A_floor = A_floor, survival_divisor = survival_divisor
    )

    # Gauss-Seidel: update EBeta_mat_curr[k, ] so subsequent k' uses the new value.
    EBeta_mat_curr[k, ] <- res_k$mean
    EBeta2_mat_new[k, ] <- res_k$second
    details[[k]]        <- res_k
  }

  list(
    EBeta   = EBeta_mat_curr,
    EBeta2  = EBeta2_mat_new,
    details = details
  )
}

# =============================================================================
# External prediction fallback: pooled beta from pooled patient-sums
# =============================================================================

#' Pooled beta_k for external prediction (no cohort membership available)
#'
#' A new external cohort has no beta^(c) of its own. Rather than average
#' the C fitted cohort coefficients (rowMeans(EBeta_mat), WRONG with
#' unequal cohort sizes -- 144 vs 129 in this project's own training
#' data), this recomputes A_k/B_k from POOLED PATIENT-SUMS over all
#' training patients (as if beta were shared, i.e. exactly
#' update_beta_k()'s formula applied to the full training set at the
#' final converged ZF/w/z), then re-solves the scalar EBNM problem.
#' This is what "the shared part of the signal" means operationally --
#' see the roxygen on predict_cox_on_yf()'s cohort_id_test argument.
#'
#' @param w      n-vector: Cox neg-diagonal Hessian weights at convergence
#' @param z      n-vector: full working response at convergence
#' @param ZF     n x K matrix: observed projection scores at convergence
#' @param EBeta_pooled_init K-vector: Gauss-Seidel starting point (e.g. the
#'                rowMeans of the fitted cohort betas, or zeros)
#' @param prior_family,alpha,A_floor,survival_divisor: forwarded to
#'   update_beta_all() (code/update_beta.R) -- this IS exactly
#'   update_beta_all()'s computation (Gauss-Seidel over k, shared beta
#'   across all patients), applied at the converged w/z/ZF; no new
#'   arithmetic is introduced here, only the "pooled-sums, not
#'   rowMeans(EBeta)" framing from the roxygen above.
#'
#' @return K-vector: pooled posterior mean beta_k, one scalar per factor
#' @export
#' @family beta_cohort_update
#' @seealso \code{\link{update_beta_all}}
compute_pooled_beta <- function(w, z, ZF, EBeta_pooled_init,
                                 prior_family = "normal",
                                 alpha        = 0.5,
                                 A_floor      = 1e-10,
                                 survival_divisor = 1) {
  if (!exists("update_beta_all", mode = "function"))
    stop("update_beta_all() must be sourced (code/update_beta.R) before compute_pooled_beta().")

  res <- update_beta_all(w, z, ZF, ZF^2, EBeta_pooled_init,
                          prior_family = prior_family, alpha = alpha,
                          A_floor = A_floor, survival_divisor = survival_divisor)
  res$EBeta
}
