# ============================================================
# Script: code/compute_bic.R
# Purpose: Joint (genomics + survival) log-likelihood and BIC for a
#          fit_cox_on_yf() (YFB, eta = (YF)beta) model object.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-27
# Dependencies: code/fit_cox_on_yf.R (must be sourced first -- provides
#               calc_cox_taylor_yf()); no other dependency on compute_elbo.R
#               to avoid a source cycle (compute_elbo.R is sourced BY
#               fit_cox_on_yf.R).
# ============================================================

# ------------------------------------------------------------------------------
#' Joint log-likelihood and BIC for a fit_cox_on_yf() model
#'
#' Computes an ELBO-style joint log-likelihood, combining the genomics
#' reconstruction term already tracked during CAVI with a Cox partial
#' log-likelihood recomputed post-hoc from the fit's returned parameters, and
#' derives a BIC from it for comparing K_init across a sweep.
#'
#' \strong{This is an ELBO-style bound, not an exact marginal likelihood.}
#' The genomics term is the variational lower bound on
#' E_q\[log p(Y | L, F, tau)\] tracked at convergence
#' (`fit$history$elbo_proxy[fit$history$n_iter]`), not the true marginal
#' log p(Y). Treat cross-K comparisons as bound comparisons, not exact
#' likelihood-ratio comparisons.
#'
#' \strong{Genomics term}: read directly from
#' `fit$history$elbo_proxy[fit$history$n_iter]` (the value written by
#' `update_tau()` after the final factor sweep, from the same EL/EF/Tau that
#' are returned) rather than `tail(fit$history$elbo_proxy, 1)`, since a fit
#' that converges before `max_iter` leaves trailing NA entries in the history
#' vector. `update_tau()` omits the Gaussian normalizing constant
#' `-(n*p/2)*log(2*pi)`; `include_gaussian_const = TRUE` adds it back. This
#' constant does not depend on K, so it never changes the K ranking -- it
#' only makes the printed log-likelihood defensible as an actual Gaussian
#' log-density term.
#'
#' \strong{Survival term}: recomputed post-hoc, NOT read from the in-loop
#' `history$elbo_full` survival component. The in-loop `logPL` is evaluated
#' at the previous iterate's EF/EBeta and predates the post-loop
#' `sign_correction` beta flip, so it does not correspond to the returned
#' fit. Recomputing `ZF <- Y %*% sweep(fit$EF, 2, fit$EF_norms, "/")`,
#' `eta <- ZF %*% fit$EBeta`, then `calc_cox_taylor_yf(eta, time, status)`
#' matches exactly what `predict_cox_on_yf()` does at test time, so this
#' scores the model as actually deployed. `EF_norms` is floored at `1e-10`
#' inside `fit_cox_on_yf()` (not exactly zero), so an ARD-pruned all-zero EF
#' column divides cleanly to 0 rather than producing NaN.
#'
#' \strong{df}: `K_init * (n + p + 1)` -- i.e. proportional to the number of
#' factors the CAVI loop was *initialized* with (`length(fit$EBeta)`), not
#' the number of ARD-surviving active factors (K_eff). This is a deliberate
#' choice: K_eff is a post-hoc classification of the same fit, not a
#' separate model with fewer free parameters, so charging df by K_eff would
#' understate model complexity. `n + p + 1` per factor is (n loadings in
#' EL[,k]) + (p loadings in EF[,k]) + (1 survival coefficient EBeta\[k\]).
#'
#' @param fit A `fit_cox_on_yf()` return object. Must contain `EF`,
#'   `EF_norms`, `EBeta`, and `history$elbo_proxy`/`history$n_iter`.
#' @param Y n x p genomics matrix used to fit the model (same orientation as
#'   passed to `fit_cox_on_yf()`).
#' @param time n-vector of observed survival/censoring times.
#' @param status n-vector of event indicators (1 = event, 0 = censored).
#' @param strata NULL (default), or an n-vector of stratum labels, forwarded
#'   to `calc_cox_taylor_yf()` for a stratified partial likelihood.
#' @param include_gaussian_const Logical; add the Gaussian normalizing
#'   constant `-(n*p/2)*log(2*pi)` to the genomics term. Default TRUE.
#' @return A list with `loglik_genomics`, `loglik_survival`, `loglik_joint`,
#'   `df`, `bic`, `n`, `p`, `K_init`.
# ------------------------------------------------------------------------------
compute_joint_ll_bic <- function(fit, Y, time, status, strata = NULL,
                                  include_gaussian_const = TRUE) {

  if (!exists("calc_cox_taylor_yf")) {
    stop("calc_cox_taylor_yf() not found -- source code/fit_cox_on_yf.R before ",
         "code/compute_bic.R.")
  }

  required_fields <- c("EF", "EF_norms", "EBeta", "history")
  missing_fields <- setdiff(required_fields, names(fit))
  if (length(missing_fields) > 0) {
    stop("fit is missing required field(s): ", paste(missing_fields, collapse = ", "))
  }
  if (is.null(fit$history$elbo_proxy) || is.null(fit$history$n_iter)) {
    stop("fit$history must contain elbo_proxy and n_iter.")
  }

  n <- nrow(Y)
  p <- ncol(Y)

  if (length(time) != n || length(status) != n) {
    stop(sprintf("time/status length (%d/%d) must match nrow(Y) (%d).",
                 length(time), length(status), n))
  }
  if (nrow(fit$EF) != p) {
    stop(sprintf("fit$EF has %d rows, expected p = ncol(Y) = %d.", nrow(fit$EF), p))
  }
  K_init <- ncol(fit$EF)
  if (length(fit$EF_norms) != K_init) {
    stop(sprintf("fit$EF_norms has length %d, expected ncol(fit$EF) = %d.",
                 length(fit$EF_norms), K_init))
  }
  if (length(fit$EBeta) != K_init) {
    stop(sprintf("fit$EBeta has length %d, expected ncol(fit$EF) = %d.",
                 length(fit$EBeta), K_init))
  }

  # fit$EF_cohort columns contribute to elbo_proxy (via update_tau on the
  # cohort-augmented reconstruction) but are not charged in df above, since
  # df only counts the K_init genomics/survival factors. LL/BIC from a
  # cohort_id fit are therefore not comparable to a cohort_id = NULL fit at
  # the same K_init.
  if (!is.null(fit$EF_cohort)) {
    warning("fit has a non-NULL EF_cohort: elbo_proxy includes cohort columns ",
            "that df does not charge, so LL/BIC are not comparable to a ",
            "cohort_id = NULL fit at the same K_init.")
  }

  n_iter <- fit$history$n_iter
  if (length(fit$history$elbo_proxy) < n_iter) {
    stop(sprintf("fit$history$elbo_proxy has length %d, shorter than n_iter = %d.",
                 length(fit$history$elbo_proxy), n_iter))
  }

  # ---- Genomics term: E_q[log p(Y | L, F, tau)] --------------------------
  loglik_genomics <- fit$history$elbo_proxy[n_iter]
  if (include_gaussian_const) {
    loglik_genomics <- loglik_genomics - (n * p / 2) * log(2 * pi)
  }
  if (!is.finite(loglik_genomics)) {
    stop(sprintf("loglik_genomics is not finite (%.6g) at n_iter = %d.",
                 loglik_genomics, n_iter))
  }

  # ---- Survival term: recomputed post-hoc at the returned EF/EBeta -------
  ZF  <- Y %*% sweep(fit$EF, 2, fit$EF_norms, "/")
  eta <- as.vector(ZF %*% fit$EBeta)
  surv <- calc_cox_taylor_yf(eta, time, status, strata = strata)
  loglik_survival <- surv$logPL
  if (!is.finite(loglik_survival)) {
    stop(sprintf("loglik_survival is not finite (%.6g).", loglik_survival))
  }

  loglik_joint <- loglik_genomics + loglik_survival

  # ---- df and BIC: K_init-based, not K_eff-based --------------------------
  df  <- K_init * (n + p + 1)
  bic <- -2 * loglik_joint + log(n) * df

  list(
    loglik_genomics = loglik_genomics,
    loglik_survival = loglik_survival,
    loglik_joint    = loglik_joint,
    df              = df,
    bic             = bic,
    n               = n,
    p               = p,
    K_init          = K_init
  )
}
