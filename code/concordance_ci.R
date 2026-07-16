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
#'
#' @return A list with `estimate` (oriented point-estimate C-index), `lower`,
#'   `upper` (percentile CI bounds), `se` (bootstrap standard deviation),
#'   `B`, and `flipped` (logical; whether the risk score was sign-flipped to
#'   establish the "higher = worse" orientation).
#'
#' @family concordance-ci
#' @seealso bootstrap_concordance_diff_ci
bootstrap_concordance_ci <- function(risk, time, status, B = 2000, seed = 1,
                                      conf_level = 0.95) {
  n <- length(risk)
  if (length(time) != n || length(status) != n)
    stop("risk, time, and status must all have the same length.")
  if (anyNA(risk) || anyNA(time) || anyNA(status))
    stop("risk, time, and status must not contain NA values.")
  if (n < 10)
    stop("Too few observations (n < 10) for a bootstrap concordance CI.")
  if (sum(status) < 2)
    stop("Too few events (< 2) to compute a meaningful concordance CI.")

  c_full <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  flip <- c_full < 0.5
  risk_use <- if (flip) -risk else risk
  estimate <- if (flip) 1 - c_full else c_full

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
      concordance(Surv(time[idx], status[idx]) ~ risk_use[idx])$concordance
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
#'
#' @return A list with `estimate` (point-estimate difference, C(risk_a) -
#'   C(risk_b)), `lower`, `upper` (percentile CI on the difference), `se`,
#'   `B`, and `significant` (logical; TRUE if the CI excludes 0).
#'
#' @family concordance-ci
#' @seealso bootstrap_concordance_ci
bootstrap_concordance_diff_ci <- function(risk_a, risk_b, time, status,
                                           B = 2000, seed = 1,
                                           conf_level = 0.95) {
  n <- length(risk_a)
  if (length(risk_b) != n || length(time) != n || length(status) != n)
    stop("risk_a, risk_b, time, and status must all have the same length.")
  if (anyNA(risk_a) || anyNA(risk_b) || anyNA(time) || anyNA(status))
    stop("risk_a, risk_b, time, and status must not contain NA values.")
  if (n < 10)
    stop("Too few observations (n < 10) for a paired bootstrap concordance CI.")
  if (sum(status) < 2)
    stop("Too few events (< 2) to compute a meaningful concordance CI.")

  c_full_a <- as.numeric(concordance(Surv(time, status) ~ risk_a)$concordance)
  c_full_b <- as.numeric(concordance(Surv(time, status) ~ risk_b)$concordance)
  flip_a <- c_full_a < 0.5
  flip_b <- c_full_b < 0.5
  risk_a_use <- if (flip_a) -risk_a else risk_a
  risk_b_use <- if (flip_b) -risk_b else risk_b
  estimate <- (if (flip_a) 1 - c_full_a else c_full_a) -
              (if (flip_b) 1 - c_full_b else c_full_b)

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
    ca <- as.numeric(concordance(Surv(time[idx], status[idx]) ~ risk_a_use[idx])$concordance)
    cb <- as.numeric(concordance(Surv(time[idx], status[idx]) ~ risk_b_use[idx])$concordance)
    boot_diff[b] <- ca - cb
  }

  alpha <- 1 - conf_level
  qs <- stats::quantile(boot_diff, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  significant <- (qs[1] > 0) || (qs[2] < 0)

  list(estimate = estimate, lower = qs[1], upper = qs[2],
       se = stats::sd(boot_diff), B = B, significant = significant)
}
