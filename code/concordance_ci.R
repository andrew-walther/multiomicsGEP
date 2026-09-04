# ============================================================
# Script: concordance_ci.R
# Purpose: Bootstrap confidence intervals for Harrell's C-index (concordance),
#          and a paired-bootstrap CI for the difference in C-index between two
#          risk scores evaluated on the same patients.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-16
# Dependencies: survival
# ============================================================

suppressPackageStartupMessages(library(survival))

#' @title Concordance under a frozen, training-data-derived orientation
#'
#' @description Scores a risk vector as-is, under the convention "larger risk
#'   = shorter survival" (a Cox risk score), using
#'   \code{survival::concordance(..., reverse = TRUE)}. Unlike this project's
#'   older \code{oriented_cindex()} pattern (\code{max(c_raw, 1 - c_raw)}),
#'   this does NOT look at `time`/`status` to decide whether to flip the risk
#'   score -- the orientation must already be fixed, at fit time, from
#'   training data only (e.g. `fit_cox_on_yf()`'s Phase C sign correction,
#'   fixed 2026-09-04 -- see DECISIONS.md). A risk score that is genuinely
#'   anti-concordant on the data it is scored against will correctly report a
#'   value below 0.5 here, rather than having that below-chance result masked
#'   by taking the max with its complement.
#'
#' @param risk numeric vector of per-patient linear-predictor risk scores,
#'   already oriented ("higher = worse prognosis") by the caller.
#' @param time,status as in \code{bootstrap_concordance_ci()}.
#'
#' @return Numeric scalar: the concordance, NOT sign-adjusted.
#'
#' @family concordance-ci
frozen_reverse_cindex <- function(risk, time, status) {
  if (length(risk) != length(time) || length(time) != length(status))
    stop("risk, time, and status must all have the same length.")
  if (sd(risk) == 0) return(NA_real_)  # constant risk score: concordance undefined
  as.numeric(concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance)
}

#' @title Percentile-bootstrap confidence interval for a C-index
#'
#' @description Resamples patients with replacement and recomputes Harrell's
#'   concordance each time, returning a percentile CI. The risk score's
#'   orientation ("higher = worse prognosis") is fixed **once**, from the
#'   full sample, using this project's existing oriented-concordance
#'   convention (\code{max(c, 1-c)}; see \code{oriented_cindex()} in
#'   \code{results/benchmark_sim/run_desurv_comparison.R}); individual
#'   bootstrap replicates are NOT re-oriented. Re-orienting every replicate
#'   would force every resample's concordance to be $\ge 0.5$ by
#'   construction, which upward-biases the CI for a weak or null signal.
#'   Fixing the orientation once and letting replicates vary freely below
#'   0.5 is what gives an honest interval.
#'
#' @param risk numeric vector of per-patient linear-predictor risk scores.
#' @param time numeric vector of follow-up times (same length as `risk`).
#' @param status integer/logical event indicator, 1/TRUE = event (same length).
#' @param B integer number of bootstrap replicates (default 2000).
#' @param seed integer RNG seed, for reproducibility (default 1).
#' @param conf_level confidence level (default 0.95).
#' @param flip NULL (default) or a logical scalar. NULL preserves the
#'   original behavior: the flip is decided from the SAME `risk`/`time`/
#'   `status` being scored (this project's original `oriented_cindex()`-style
#'   convention; see the Description above for why this is fine for CI width
#'   but is a same-sample orientation decision). Pass a logical to instead use
#'   a FROZEN orientation established elsewhere (e.g. from training data
#'   only, per DECISIONS.md 2026-09-04) -- `estimate` is then the concordance
#'   of `risk` (or `-risk` if `flip = TRUE`) as-is, NOT adjusted toward 0.5,
#'   so a genuinely poor score can come back below 0.5.
#'
#' @return A list with `estimate` (oriented point-estimate C-index), `lower`,
#'   `upper` (percentile CI bounds), `se` (bootstrap standard deviation),
#'   `B`, and `flipped` (logical; whether the risk score was sign-flipped to
#'   establish the "higher = worse" orientation).
#'
#' @family concordance-ci
#' @seealso bootstrap_concordance_diff_ci, frozen_reverse_cindex
bootstrap_concordance_ci <- function(risk, time, status, B = 2000, seed = 1,
                                      conf_level = 0.95, flip = NULL) {
  n <- length(risk)
  if (length(time) != n || length(status) != n)
    stop("risk, time, and status must all have the same length.")
  if (anyNA(risk) || anyNA(time) || anyNA(status))
    stop("risk, time, and status must not contain NA values.")
  if (n < 10)
    stop("Too few observations (n < 10) for a bootstrap concordance CI.")
  if (sum(status) < 2)
    stop("Too few events (< 2) to compute a meaningful concordance CI.")

  # frozen_orientation: whether flip was supplied by the caller (TRUE/FALSE),
  # as opposed to being decided here from the same sample being scored (the
  # original, flip=NULL behavior). This also controls whether concordance()
  # is called with reverse=TRUE below: a caller-supplied flip is meant to
  # match frozen_reverse_cindex()'s "larger risk = shorter survival" (Cox)
  # convention, which requires reverse=TRUE; the legacy flip=NULL path uses
  # plain (non-reverse) concordance throughout, unchanged, for exact
  # backward compatibility.
  frozen_orientation <- !is.null(flip)
  if (!frozen_orientation) {
    c_full <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
    flip <- c_full < 0.5
    risk_use <- if (flip) -risk else risk
    estimate <- if (flip) 1 - c_full else c_full
  } else {
    if (!is.logical(flip) || length(flip) != 1 || is.na(flip))
      stop("flip must be NULL or a single non-NA logical.")
    risk_use <- if (flip) -risk else risk
    estimate <- as.numeric(concordance(Surv(time, status) ~ risk_use, reverse = TRUE)$concordance)
  }

  set.seed(seed)
  boot_c <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    # A resample can, by chance, contain fewer than 2 events; concordance() then
    # returns NaN silently rather than erroring. Fail loud here, at the point of
    # failure, with an actionable message, rather than let a downstream
    # quantile() error obscure the actual cause.
    if (sum(status[idx]) < 2)
      stop(sprintf(
        "Bootstrap replicate %d/%d has fewer than 2 events after resampling ",
        b, B), sprintf("(n=%d, seed=%d, overall events=%d). ", n, seed, sum(status)),
        "Concordance is undefined for this resample -- this can happen with a ",
        "low event-rate or small cohort; try a different seed or a larger B.")
    boot_c[b] <- as.numeric(
      concordance(Surv(time[idx], status[idx]) ~ risk_use[idx], reverse = frozen_orientation)$concordance
    )
  }

  alpha <- 1 - conf_level
  qs <- stats::quantile(boot_c, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

  list(estimate = estimate, lower = qs[1], upper = qs[2],
       se = stats::sd(boot_c), B = B, flipped = flip)
}

#' @title Paired-bootstrap confidence interval for a difference in C-index
#'
#' @description Compares two risk scores evaluated on the SAME patients (same
#'   `time`/`status`) -- e.g. the recommended model vs. a two-step baseline
#'   scored on the same held-out cohort. Each bootstrap replicate resamples
#'   patients ONCE and scores both risk vectors on that identical resample,
#'   preserving the pairing/correlation between the two models' errors on the
#'   same patients. This is the statistically correct way to test whether one
#'   model is significantly more concordant than another on shared data --
#'   comparing two independently-constructed CIs is not equivalent and can be
#'   overly conservative.
#'
#'   As in \code{bootstrap_concordance_ci()}, each risk score's orientation is
#'   fixed once from the full sample; a difference is only meaningful once
#'   both scores are on a "higher = worse" scale.
#'
#' @param risk_a,risk_b numeric vectors of risk scores from the two models
#'   being compared (same patients, same length).
#' @param time,status as in `bootstrap_concordance_ci()`.
#' @param B,seed,conf_level as in `bootstrap_concordance_ci()`.
#' @param flip_a,flip_b NULL (default) or logical scalars, analogous to
#'   `bootstrap_concordance_ci()`'s `flip`: NULL decides each score's flip
#'   from the same `time`/`status` being scored here (original behavior);
#'   a logical uses a pre-frozen orientation instead (e.g. each risk score is
#'   already correctly signed from its own model's training fit, so both are
#'   typically passed as `FALSE`).
#'
#' @return A list with `estimate` (point-estimate difference, C(risk_a) -
#'   C(risk_b)), `lower`, `upper` (percentile CI on the difference), `se`,
#'   `B`, and `significant` (logical; TRUE if the CI excludes 0).
#'
#' @family concordance-ci
#' @seealso bootstrap_concordance_ci, frozen_reverse_cindex
bootstrap_concordance_diff_ci <- function(risk_a, risk_b, time, status,
                                           B = 2000, seed = 1,
                                           conf_level = 0.95,
                                           flip_a = NULL, flip_b = NULL) {
  n <- length(risk_a)
  if (length(risk_b) != n || length(time) != n || length(status) != n)
    stop("risk_a, risk_b, time, and status must all have the same length.")
  if (anyNA(risk_a) || anyNA(risk_b) || anyNA(time) || anyNA(status))
    stop("risk_a, risk_b, time, and status must not contain NA values.")
  if (n < 10)
    stop("Too few observations (n < 10) for a paired bootstrap concordance CI.")
  if (sum(status) < 2)
    stop("Too few events (< 2) to compute a meaningful concordance CI.")

  # .score() preserves the ORIGINAL arithmetic (1 - c_full via the complement,
  # not a fresh concordance() call on the negated vector) when flip is NULL,
  # so existing flip_a/flip_b=NULL callers get bit-for-bit identical results
  # to before this parameter was added. A caller-supplied flip is scored with
  # reverse=TRUE, matching frozen_reverse_cindex()'s "larger risk = shorter
  # survival" (Cox) convention -- see bootstrap_concordance_ci() for why.
  .score <- function(risk, flip) {
    if (is.null(flip)) {
      c_full <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
      flip <- c_full < 0.5
      list(flip = flip, frozen = FALSE, c = if (flip) 1 - c_full else c_full)
    } else {
      if (!is.logical(flip) || length(flip) != 1 || is.na(flip))
        stop("flip_a/flip_b must be NULL or a single non-NA logical.")
      risk_use <- if (flip) -risk else risk
      list(flip = flip, frozen = TRUE,
           c = as.numeric(concordance(Surv(time, status) ~ risk_use, reverse = TRUE)$concordance))
    }
  }
  res_a <- .score(risk_a, flip_a); flip_a <- res_a$flip; frozen_a <- res_a$frozen
  res_b <- .score(risk_b, flip_b); flip_b <- res_b$flip; frozen_b <- res_b$frozen
  risk_a_use <- if (flip_a) -risk_a else risk_a
  risk_b_use <- if (flip_b) -risk_b else risk_b
  estimate <- res_a$c - res_b$c

  set.seed(seed)
  boot_diff <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    # See bootstrap_concordance_ci() for why this is checked explicitly: a
    # degenerate resample (< 2 events) would otherwise silently produce NaN.
    if (sum(status[idx]) < 2)
      stop(sprintf(
        "Bootstrap replicate %d/%d has fewer than 2 events after resampling ",
        b, B), sprintf("(n=%d, seed=%d, overall events=%d). ", n, seed, sum(status)),
        "Concordance is undefined for this resample -- this can happen with a ",
        "low event-rate or small cohort; try a different seed or a larger B.")
    ca <- as.numeric(concordance(Surv(time[idx], status[idx]) ~ risk_a_use[idx], reverse = frozen_a)$concordance)
    cb <- as.numeric(concordance(Surv(time[idx], status[idx]) ~ risk_b_use[idx], reverse = frozen_b)$concordance)
    boot_diff[b] <- ca - cb
  }

  alpha <- 1 - conf_level
  qs <- stats::quantile(boot_diff, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  significant <- (qs[1] > 0) || (qs[2] < 0)

  list(estimate = estimate, lower = qs[1], upper = qs[2],
       se = stats::sd(boot_diff), B = B, significant = significant)
}
