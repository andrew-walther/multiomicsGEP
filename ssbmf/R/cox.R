# Cox primitives in log space ----
#
# Breslow partial likelihood with strata:
#   l_PL(eta) = sum_{i: delta_i = 1} [ eta_i - log sum_{j in R_i} exp(eta_j) ],
#   R_i = { j : t_j >= t_i, s_j = s_i }.
# All risk-set sums are formed as cumulative log-sum-exp, so |eta| of order
# 1e4 stays finite.

#' Cumulative log-sum-exp, numerically stable
#' @param x numeric vector
#' @return numeric vector, element q = log sum_{r <= q} exp(x_r)
#' @keywords internal
#' @noRd
.logcumsumexp <- function(x) {
  if (!length(x)) return(x)
  M <- max(x)
  cs <- cumsum(exp(x - M))
  if (all(cs > 0) && is.finite(M)) return(log(cs) + M)
  # slow path: running maximum
  out <- numeric(length(x)); run <- -Inf
  for (q in seq_along(x)) {
    a <- x[q]
    run <- if (run == -Inf) a else max(run, a) + log1p(exp(-abs(run - a)))
    out[q] <- run
  }
  out
}

#' Pre-compute risk-set orderings for each stratum
#' @return list of strata, each with `idx` (sample indices sorted by
#'   decreasing time), `last_desc` (for each sorted position, the last position
#'   of its tie group in decreasing order) and `first_desc` (first position of
#'   its tie group).
#' @keywords internal
#' @noRd
.risk_index <- function(time, strata = NULL) {
  s <- if (is.null(strata)) rep(1L, length(time)) else as.integer(factor(strata))
  lapply(split(seq_along(time), s), function(idx) {
    o <- idx[order(-time[idx])]
    rl <- rle(time[o])
    ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1L
    list(idx = o, last_desc = rep(ends, rl$lengths), first_desc = rep(starts, rl$lengths))
  })
}

#' Log risk-set sums log sum_{j in R_i} exp(eta_j) (Breslow)
#'
#' @param eta numeric n-vector
#' @param time survival times
#' @param strata optional stratum labels
#' @param ri optional result of `.risk_index()`
#' @return numeric n-vector
#' @keywords internal
.logsumexp_risk <- function(eta, time, strata = NULL, ri = NULL) {
  if (is.null(ri)) ri <- .risk_index(time, strata)
  out <- numeric(length(eta))
  for (r in ri) {
    cs <- .logcumsumexp(eta[r$idx])
    out[r$idx] <- cs[r$last_desc]
  }
  out
}

#' Cox log partial likelihood (Breslow ties, optional strata)
#' @inheritParams .logsumexp_risk
#' @param status event indicator
#' @return scalar log partial likelihood
#' @keywords internal
.cox_loglik <- function(eta, time, status, strata = NULL, ri = NULL) {
  lse <- .logsumexp_risk(eta, time, strata, ri)
  sum(status * (eta - lse))
}

#' Moment-bound Cox term and its derivatives
#'
#' Evaluates the Gaussian-moment (Jensen) bound of the derivation, Section 3,
#'   S = sum_{i: delta_i = 1} [ eta_bar_i - log sum_{j in R_i} exp(a_j) ],
#'   a_j = eta_bar_j + v_j / 2,
#' and the quantities used by the updates: P_j = sum_{events i, j in R_i}
#' pi_ij with pi_ij = exp(a_j) / sum_{l in R_i} exp(a_l), and the residual
#' rho_j = delta_j - P_j. Then dS/d eta_bar = rho and dS/dv = -P/2.
#'
#' @param eta_bar,v predictor means and variances under q
#' @param status event indicator
#' @param ri result of `.risk_index()`
#' @return list(`S`, `P`, `rho`)
#' @keywords internal
.cox_moment <- function(eta_bar, v, status, ri) {
  a <- eta_bar + 0.5 * v
  n <- length(a)
  lse <- numeric(n); logP <- numeric(n)
  for (r in ri) {
    o <- r$idx
    cs <- .logcumsumexp(a[o])
    lse_o <- cs[r$last_desc]
    lse[o] <- lse_o
    # sum over events i with t_i <= t_j of exp(-lse_i): cumulative over
    # increasing time, i.e. reverse of the decreasing order, tie-inclusive.
    ev <- status[o] == 1
    term <- ifelse(ev, -lse_o, -Inf)
    rev_cs <- rev(.logcumsumexp(rev(term)))        # sum over positions q' >= q (times <= t_q)
    logP[o] <- a[o] + rev_cs[r$first_desc]
  }
  P <- exp(logP)
  list(S = sum(status * (eta_bar - lse)), P = P, rho = status - P)
}
