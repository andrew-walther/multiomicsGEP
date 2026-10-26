# ============================================================
# Script: verify_derivation.R
# Purpose: Numerical checks of ssbmf_derivation.qmd: (1) the Cox-term bounds
#          of Section 3 against Monte Carlo values of E_q log PL, and (2) the
#          stationarity of every coordinate update of Sections 4-5 by central
#          finite differences of the objective.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: survival; part 2 loads the ssbmf package with pkgload
# Usage: Rscript derivations/unified_ssbmf/verify_derivation.R [bounds|stationarity|all]
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args[1] else "all"

# Part 1: Cox-term bounds ----

#' Log of risk-set sums, log sum_{j in R_i} exp(a_j), Breslow ties, by stratum
#'
#' @param a numeric n-vector
#' @param time,strata survival times and stratum labels
#' @return numeric n-vector; element i is the log risk-set sum for R_i
lse_risk <- function(a, time, strata) {
  out <- numeric(length(a))
  for (s in unique(strata)) {
    idx <- which(strata == s)
    o <- idx[order(-time[idx])]                 # decreasing time
    cs <- numeric(length(o)); run <- -Inf
    for (q in seq_along(o)) { run <- max(run, a[o[q]]) + log1p(exp(-abs(run - a[o[q]]))); cs[q] <- run }
    # ties: every member of a tie group gets the sum through the last tied member
    tt <- time[o]
    last <- vapply(seq_along(o), function(q) max(which(tt == tt[q])), 1L)
    out[o] <- cs[last]
  }
  out
}

cox_ll <- function(eta, time, status, strata) sum(status * (eta - lse_risk(eta, time, strata)))

#' Draw one state from a mean-field q and return Monte Carlo and bound values
#'
#' F_jk ~ (1 - w) delta_0 + w N(m_jk, sd_jk^2) (spike-and-slab, signed like a
#' point-Laplace posterior) or its absolute value (nonnegative); beta_k ~ N(b_k, s_k^2).
check_bounds <- function(seed, n = 50, p = 200, K = 3, nonneg = FALSE, B = 20000, beta_scale = 0.4) {
  set.seed(seed)
  time <- round(rexp(n), 1) + 0.1                 # rounding creates ties
  status <- rbinom(n, 1, 0.7); strata <- rep(1:2, length.out = n)
  Yt <- scale(matrix(rnorm(n * p), n, p) %*% diag(runif(p, 0.5, 1.5)), scale = FALSE)
  w <- matrix(runif(p * K, 0.2, 0.9), p, K)
  m <- matrix(rnorm(p * K, 0, 0.05), p, K); sdF <- matrix(runif(p * K, 0.01, 0.05), p, K)
  b <- rnorm(K, 0, beta_scale); sb <- runif(K, 0.05, 0.3)
  # exact moments of the spike-and-slab (or folded) components
  if (!nonneg) {
    f1 <- w * m; f2 <- w * (m^2 + sdF^2)
  } else {
    # |N(m, sd^2)| moments by simulation (large) for the folded case
    draws <- array(abs(rnorm(p * K * 4000, m, sdF)), c(p, K, 4000))
    f1 <- w * apply(draws, 1:2, mean); f2 <- w * apply(draws^2, 1:2, mean)
  }
  Zbar <- Yt %*% f1; V <- (Yt^2) %*% (f2 - f1^2)
  eta_bar <- drop(Zbar %*% b)
  v <- drop((Zbar^2) %*% sb^2 + V %*% (b^2 + sb^2))
  # Monte Carlo E_q log PL
  mc <- numeric(B)
  for (r in seq_len(B)) {
    Fd <- matrix(rnorm(p * K, m, sdF), p, K); if (nonneg) Fd <- abs(Fd)
    Fd <- Fd * (matrix(runif(p * K), p, K) < w)
    bd <- rnorm(K, b, sb)
    mc[r] <- cox_ll(drop(Yt %*% Fd %*% bd), time, status, strata)
  }
  ev <- status == 1
  jensen <- sum(status * eta_bar) - sum(lse_risk(eta_bar + v / 2, time, strata)[ev])
  d_j <- vapply(seq_len(n), function(j) sum(ev & strata == strata[j] & time <= time[j]), 1)
  bohning <- cox_ll(eta_bar, time, status, strata) - 0.25 * sum(d_j * v)
  c(mc = mean(mc), mc_se = sd(mc) / sqrt(B), jensen = jensen, bohning = bohning,
    plugin = cox_ll(eta_bar, time, status, strata), sd_eta = sd(eta_bar), mean_sd_q_eta = mean(sqrt(v)))
}

run_bounds <- function() {
  cat("Part 1: Cox-term bounds versus Monte Carlo E_q log PL\n")
  grid <- expand.grid(seed = 1:6, nonneg = c(FALSE, TRUE), beta_scale = c(0.4, 1.5))
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(r) {
    g <- grid[r, ]
    cbind(g, t(check_bounds(g$seed, nonneg = g$nonneg, beta_scale = g$beta_scale, B = 4000)))
  }))
  res$gap_jensen <- res$mc - res$jensen          # >= 0 if (c) is a lower bound
  res$gap_bohning <- res$mc - res$bohning
  res$jensen_ok <- res$gap_jensen > -3 * res$mc_se
  res$bohning_ok <- res$gap_bohning > -3 * res$mc_se
  print(format(res[, c("seed", "nonneg", "beta_scale", "sd_eta", "mean_sd_q_eta", "mc", "mc_se",
                       "jensen", "gap_jensen", "bohning", "gap_bohning", "jensen_ok", "bohning_ok")],
               digits = 4), row.names = FALSE)
  cat(sprintf("\nJensen (c) below Monte Carlo (within 3 MC SE): %d of %d states; median gap %.3f nats\n",
              sum(res$jensen_ok), nrow(res), median(res$gap_jensen)))
  cat(sprintf("Bohning (b) below Monte Carlo: %d of %d states; median gap %.3f nats\n",
              sum(res$bohning_ok), nrow(res), median(res$gap_bohning)))
  invisible(res)
}

# Part 2: stationarity of the updates (requires the package) ----

run_stationarity <- function() {
  cat("\nPart 2: stationarity of each coordinate update by central finite differences\n")
  pkgload::load_all("ssbmf", quiet = TRUE, export_all = TRUE)   # run from the repository root
  ok <- ssbmf:::.verify_stationarity(seed = 1, verbose = TRUE)
  cat(sprintf("Maximum absolute finite-difference gradient after updates: %.2e per posterior SD of each entry (threshold 1e-6)\n",
              max(ok$max_grad)))
  invisible(ok)
}

if (what %in% c("bounds", "all")) bounds <- run_bounds()
if (what %in% c("stationarity", "all")) stat <- run_stationarity()
