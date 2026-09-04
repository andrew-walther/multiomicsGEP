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
#' column divides cleanly to 0 rather than producing NaN. Before scoring,
#' `eta` is independently re-oriented by a correct-direction
#' (`reverse = TRUE`) concordance check -- `fit$EBeta`'s sign, as returned
#' by `fit_cox_on_yf()`, is not a reliable hazard-direction indicator (see
#' the orientation-correction comment at the "Survival term" code block
#' below, and DECISIONS.md 2026-09-04).
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
#'
#' \strong{Three limitations of this BIC, all real (see DECISIONS.md
#' 2026-09-04 and code/compute_cv_loglik.R for a genuine held-out
#' alternative -- do not "cross-validate" this BIC itself; see below).}
#' \enumerate{
#'   \item The Cox term is a \emph{partial} likelihood with the baseline
#'     hazard profiled out. Its effective sample size is the event count
#'     (414 in the D4 training data), not `n` = 273 patients -- `log(n)`
#'     charges the wrong sample size to the survival half of the model.
#'   \item `df = K_init * (n + p + 1)` is charged by the starting K
#'     (`K_init`), not the ARD-surviving `K_eff` -- deliberate (see above),
#'     but it means BIC penalizes model complexity that most of the fit's
#'     factors have already been shrunk away from. Under a cohort-specific
#'     beta (`beta_cohort_id`, see code/update_beta_cohort.R), the `+1`
#'     term becomes `+C` (C cohorts), one survival coefficient per cohort
#'     per factor instead of one shared coefficient per factor.
#'   \item `log(n) * df` applies a single, `n`-scaled penalty against a
#'     genomics reconstruction term whose effective sample size is
#'     `n * p` (approx. 563,000 entries for D4), not `n`. The penalty is
#'     calibrated for the survival term's sample size, not the genomics
#'     term's.
#' }
#' \strong{Why there is no "CV-BIC."} BIC's log(n) penalty exists to
#' correct for in-sample optimism; a genuinely held-out likelihood
#' (code/compute_cv_loglik.R) already corrects for that same optimism
#' directly, and cross-validation is asymptotically equivalent to AIC, not
#' to BIC -- combining the two would penalize the same optimism twice
#' under two different justifications. The right comparison is: keep this
#' in-sample BIC with the three caveats above, and report
#' code/compute_cv_loglik.R's held-out log-likelihoods alongside it as a
#' separate, principled criterion -- not as a "cross-validated BIC."
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

  # Orientation correction (DECISIONS.md 2026-09-04): fit_cox_on_yf()'s own
  # post-loop sign-correction (fit_cox_on_yf.R, "Phase C") checks
  # concordance(Surv(time,status) ~ eta_final) WITHOUT reverse = TRUE.
  # survival::concordance()'s formula method defaults to assuming a larger
  # predictor means a LARGER (longer) response -- the opposite of a Cox risk
  # score, where a larger eta means a SHORTER survival time; the package's
  # own documentation calls this out explicitly as the one case requiring
  # reverse = TRUE. Without it, fit_cox_on_yf()'s check is inverted, so
  # fit$EBeta's returned sign is not a reliable hazard-direction indicator
  # (confirmed empirically: a fit with true training concordance ~0.90+ in
  # the correct direction reads as ~0.90+ under the UNCORRECTED check too,
  # so the "already fine, don't flip" branch fires when a flip was in fact
  # needed). This is harmless for every sign-invariant metric already
  # reported elsewhere in this project (external C-index via
  # oriented_cindex()'s max(c, 1-c); K_eff via classify_factors()'s
  # abs(EBeta) threshold; pathway-enrichment adverse/protective direction,
  # which by design uses each program's marginal association, not joint
  # beta sign -- DECISIONS.md 2026-06-16). It is NOT harmless here: Cox
  # partial log-likelihood is not sign-invariant, so scoring an
  # anti-oriented eta as-is understates the model's genuine fit. This
  # function orients eta itself (correct-direction concordance check, using
  # reverse = TRUE), independent of whatever fit_cox_on_yf() did, rather
  # than trusting fit$EBeta's sign.
  c_check <- tryCatch(
    as.numeric(survival::concordance(survival::Surv(time, status) ~ eta,
                                      reverse = TRUE)$concordance),
    error = function(e) NA_real_
  )
  if (is.finite(c_check) && c_check < 0.5) eta <- -eta

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
