# =============================================================================
# tests/test_stratified_cox.R
# Tests for the stratified Cox partial likelihood (Item 3: study-specific
# baseline hazard, the 6/18 lab-meeting "+ strata(study)" feedback).
#
# Covers:
#   Section 1-3  calc_cox_taylor_yf(..., strata=)  (Cluster B / YFB helper)
#   Section 4    error handling (strata length mismatch)
#   Section 5    calc_cox_taylor(..., strata=)     (Cluster A / LB helper)
#   Section 6    fit_cox_on_yf(..., strata_id=)    (YFB integration)
#   Section 7    fit_supervised_mf_modular(..., strata_id=)  (LB integration)
#
# Stratified partial likelihood forms Breslow risk sets *within* each stratum;
# the baseline hazard still cancels per-stratum, so no parametric h0 is
# introduced (DECISIONS.md, Item 3).
#
# calc_cox_taylor, calc_cox_taylor_yf, fit_supervised_mf_modular, and
# fit_cox_on_yf are all sourced by run_tests.R before this file loads.
# =============================================================================

suppressPackageStartupMessages(library(survival))

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------
local({
  set.seed(303)
  n <- 60L
  .sc <<- list(
    n      = n,
    x      = rnorm(n),
    time   = rexp(n, 0.1),          # continuous => ties have measure zero
    status = rbinom(n, 1L, 0.7),
    g3     = sample(c("A", "B", "C"), n, replace = TRUE),
    beta   = 0.7
  )
})

# =============================================================================
# Section 1: single-stratum reduction (CORRECTNESS ANCHOR)
#   strata=rep(1,n) must reproduce the unstratified result bit-for-bit.
# =============================================================================

run_test("StratCox-T1: calc_cox_taylor_yf single-stratum == unstratified (u,w,logPL)", {
  eta    <- .sc$x * .sc$beta
  base   <- calc_cox_taylor_yf(eta, .sc$time, .sc$status)
  strat1 <- calc_cox_taylor_yf(eta, .sc$time, .sc$status, strata = rep(1L, .sc$n))
  assert_near(strat1$logPL, base$logPL, tol = 1e-12, msg = "logPL must be identical")
  assert_near(max(abs(strat1$u - base$u)), 0, tol = 1e-12, msg = "u must be identical")
  assert_near(max(abs(strat1$w - base$w)), 0, tol = 1e-12, msg = "w must be identical")
})

# =============================================================================
# Section 2: additive decomposition over strata
#   Stratified logPL = sum of per-stratum logPL; u,w scatter back by index.
# =============================================================================

run_test("StratCox-T2: calc_cox_taylor_yf stratified == manual per-stratum decomposition", {
  eta   <- .sc$x * .sc$beta
  g     <- .sc$g3
  strat <- calc_cox_taylor_yf(eta, .sc$time, .sc$status, strata = g)

  u_ref <- numeric(.sc$n); w_ref <- numeric(.sc$n); logpl_ref <- 0
  for (lev in unique(g)) {
    idx <- which(g == lev)
    r <- calc_cox_taylor_yf(eta[idx], .sc$time[idx], .sc$status[idx])
    u_ref[idx] <- r$u
    w_ref[idx] <- r$w
    logpl_ref  <- logpl_ref + r$logPL
  }
  assert_near(strat$logPL, logpl_ref, tol = 1e-12, msg = "logPL = sum over strata")
  assert_near(max(abs(strat$u - u_ref)), 0, tol = 1e-12, msg = "u matches decomposition")
  assert_near(max(abs(strat$w - w_ref)), 0, tol = 1e-12, msg = "w matches decomposition")
})

# =============================================================================
# Section 3: independent cross-check vs survival::coxph stratified Breslow PL
#   Evaluate the partial log-likelihood at a FIXED coefficient (iter.max=0).
# =============================================================================

run_test("StratCox-T3: calc_cox_taylor_yf logPL matches coxph + strata() (Breslow)", {
  eta <- .sc$x * .sc$beta
  res <- calc_cox_taylor_yf(eta, .sc$time, .sc$status, strata = .sc$g3)
  df  <- data.frame(time = .sc$time, status = .sc$status, x = .sc$x, g = .sc$g3)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x + strata(g), data = df,
          init = .sc$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(res$logPL, tail(fit$loglik, 1), tol = 1e-6,
              msg = "stratified logPL must match coxph reference")
})

run_test("StratCox-T3b: unstratified calc_cox_taylor_yf matches coxph (Breslow) — Breslow sanity", {
  eta <- .sc$x * .sc$beta
  res <- calc_cox_taylor_yf(eta, .sc$time, .sc$status)   # strata NULL
  df  <- data.frame(time = .sc$time, status = .sc$status, x = .sc$x)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x, data = df,
          init = .sc$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(res$logPL, tail(fit$loglik, 1), tol = 1e-6,
              msg = "helper must be Breslow, matching coxph")
})

# =============================================================================
# Section 3b: tied-event-time Breslow correctness (fixed 2026-09-04 — the
# unstratified fixture above uses rexp() continuous times, so ties have
# measure zero and never exercised this code path; real training data has
# tied event times, so a dedicated tied-time fixture is required)
# =============================================================================

local({
  set.seed(404)
  n <- 12L
  .tie <<- list(
    n      = n,
    x      = rnorm(n),
    # Deliberately tied times: three distinct values, several rows each.
    time   = rep(c(1, 2, 3, 4), each = 3L),
    status = c(1,0,1, 1,1,0, 0,1,1, 1,0,1),
    beta   = 0.6
  )
})

run_test("StratCox-T3c: calc_cox_taylor_yf logPL matches coxph (Breslow) with tied times", {
  eta <- .tie$x * .tie$beta
  res <- calc_cox_taylor_yf(eta, .tie$time, .tie$status)
  df  <- data.frame(time = .tie$time, status = .tie$status, x = .tie$x)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x, data = df,
          init = .tie$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(res$logPL, tail(fit$loglik, 1), tol = 1e-6,
              msg = "tied-time logPL must match coxph Breslow reference")
})

run_test("StratCox-T3d: calc_cox_taylor_yf u == coxph martingale residuals with tied times", {
  eta <- .tie$x * .tie$beta
  res <- calc_cox_taylor_yf(eta, .tie$time, .tie$status)
  df  <- data.frame(time = .tie$time, status = .tie$status, x = .tie$x)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x, data = df,
          init = .tie$beta, iter.max = 0, ties = "breslow")
  )
  resid_ref <- residuals(fit, type = "martingale")
  assert_near(max(abs(res$u - resid_ref)), 0, tol = 1e-6,
              msg = "tied-time u must match coxph martingale residuals")
})

run_test("StratCox-T3e: permuting tied rows leaves logPL, u, w unchanged (order invariance)", {
  eta <- .tie$x * .tie$beta
  res_base <- calc_cox_taylor_yf(eta, .tie$time, .tie$status)

  # Permute only within each tied-time block, then apply the same permutation
  # to eta so (eta, time, status) triples stay attached to the same "subject".
  # perm[j] = which original subject now sits at row j.
  set.seed(1)
  blocks <- split(seq_len(.tie$n), .tie$time)
  shuffled <- lapply(blocks, function(idx) if (length(idx) > 1) sample(idx) else idx)
  perm <- integer(.tie$n)
  perm[unlist(blocks)] <- unlist(shuffled)

  res_perm <- calc_cox_taylor_yf(eta[perm], .tie$time[perm], .tie$status[perm])
  assert_near(res_perm$logPL, res_base$logPL, tol = 1e-10,
              msg = "logPL must be invariant to permuting tied rows")
  # res_perm's row j corresponds to original subject perm[j]; scatter back to
  # original subject identity (row i) before comparing to res_base.
  inv_perm <- order(perm)
  assert_near(max(abs(res_perm$u[inv_perm] - res_base$u)), 0, tol = 1e-10,
              msg = "u must be invariant (per-subject) to permuting tied rows")
  assert_near(max(abs(res_perm$w[inv_perm] - res_base$w)), 0, tol = 1e-10,
              msg = "w must be invariant (per-subject) to permuting tied rows")
})

run_test("StratCox-T3f: calc_cox_taylor_yf stratified with ties within a stratum matches coxph", {
  eta <- .tie$x * .tie$beta
  g   <- rep(c("A", "B"), length.out = .tie$n)
  res <- calc_cox_taylor_yf(eta, .tie$time, .tie$status, strata = g)
  df  <- data.frame(time = .tie$time, status = .tie$status, x = .tie$x, g = g)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x + strata(g), data = df,
          init = .tie$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(res$logPL, tail(fit$loglik, 1), tol = 1e-6,
              msg = "stratified tied-time logPL must match coxph reference")
})

# =============================================================================
# Section 4: error handling
# =============================================================================

run_test("StratCox-T4: strata length mismatch raises an error", {
  eta <- .sc$x * .sc$beta
  got_error <- tryCatch({
    calc_cox_taylor_yf(eta, .sc$time, .sc$status, strata = c(1L, 2L, 3L))
    FALSE
  }, error = function(e) TRUE)
  assert_true(got_error, msg = "mismatched strata length must error")
})

# =============================================================================
# Section 5: LB helper mirror (calc_cox_taylor)
# =============================================================================

run_test("StratCox-T5a: calc_cox_taylor single-stratum == unstratified", {
  eta    <- .sc$x * .sc$beta
  base   <- calc_cox_taylor(eta, .sc$time, .sc$status)
  strat1 <- calc_cox_taylor(eta, .sc$time, .sc$status, strata = rep(1L, .sc$n))
  assert_near(strat1$logPL, base$logPL, tol = 1e-12, msg = "logPL identical")
  assert_near(max(abs(strat1$u - base$u)), 0, tol = 1e-12, msg = "u identical")
  assert_near(max(abs(strat1$w - base$w)), 0, tol = 1e-12, msg = "w identical")
})

run_test("StratCox-T5b: calc_cox_taylor stratified logPL matches coxph + strata()", {
  eta <- .sc$x * .sc$beta
  res <- calc_cox_taylor(eta, .sc$time, .sc$status, strata = .sc$g3)
  df  <- data.frame(time = .sc$time, status = .sc$status, x = .sc$x, g = .sc$g3)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x + strata(g), data = df,
          init = .sc$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(res$logPL, tail(fit$loglik, 1), tol = 1e-6,
              msg = "LB stratified logPL must match coxph")
})

# =============================================================================
# Section 6: fit_cox_on_yf integration (YFB)
# =============================================================================

run_test("StratCox-T6: fit_cox_on_yf strata_id=rep(1,n) == strata_id=NULL (EBeta anchor)", {
  n <- 50L; p <- 20L; K <- 2L
  set.seed(11); Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1L, 0.7)

  set.seed(11)
  f0 <- suppressMessages(fit_cox_on_yf(Y, time, status, K = K, max_iter = 5L,
                                       verbose = FALSE))
  set.seed(11)
  f1 <- suppressMessages(fit_cox_on_yf(Y, time, status, K = K, max_iter = 5L,
                                       verbose = FALSE, strata_id = rep(1L, n)))
  assert_near(max(abs(f1$EBeta - f0$EBeta)), 0, tol = 1e-8,
              msg = "single-stratum fit must reduce to the current model")
})

run_test("StratCox-T7: fit_cox_on_yf with 2 strata runs, finite EBeta of length K", {
  n <- 50L; p <- 20L; K <- 2L
  set.seed(12); Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1L, 0.7)
  g <- rep(c("A", "B"), each = n / 2)
  fit <- suppressMessages(fit_cox_on_yf(Y, time, status, K = K, max_iter = 5L,
                                        verbose = FALSE, strata_id = g))
  assert_length(fit$EBeta, K)
  assert_finite(fit$EBeta)
})

# =============================================================================
# Section 7: fit_supervised_mf_modular integration (LB)
# =============================================================================

run_test("StratCox-T8: fit_supervised_mf_modular strata_id=rep(1,n) == NULL (EBeta anchor)", {
  n <- 50L; p <- 20L; K <- 2L
  set.seed(21); Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1L, 0.7)

  set.seed(21)
  f0 <- suppressMessages(fit_supervised_mf_modular(Y, time, status, K = K,
                                                   max_iter = 5L, verbose = FALSE))
  set.seed(21)
  f1 <- suppressMessages(fit_supervised_mf_modular(Y, time, status, K = K,
                                                   max_iter = 5L, verbose = FALSE,
                                                   strata_id = rep(1L, n)))
  assert_near(max(abs(f1$EBeta - f0$EBeta)), 0, tol = 1e-8,
              msg = "single-stratum LB fit must reduce to the current model")
})

run_test("StratCox-T9: fit_supervised_mf_modular with 2 strata runs, finite EBeta length K", {
  n <- 50L; p <- 20L; K <- 2L
  set.seed(22); Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1L, 0.7)
  g <- rep(c("A", "B"), each = n / 2)
  fit <- suppressMessages(fit_supervised_mf_modular(Y, time, status, K = K,
                                                    max_iter = 5L, verbose = FALSE,
                                                    strata_id = g))
  assert_length(fit$EBeta, K)
  assert_finite(fit$EBeta)
})

# =============================================================================
# Section 8: NA handling (fail-loud)
#   as.factor() silently drops NA from levels, so an NA-labelled sample would be
#   excluded from every risk set and keep u=0,w=0 -> 0/0=NaN downstream. Must
#   error instead of silently corrupting the fit.
# =============================================================================

run_test("StratCox-T10a: calc_cox_taylor_yf errors on NA in strata", {
  eta  <- .sc$x * .sc$beta
  g_na <- .sc$g3; g_na[1] <- NA
  got  <- tryCatch({ calc_cox_taylor_yf(eta, .sc$time, .sc$status, strata = g_na); FALSE },
                   error = function(e) TRUE)
  assert_true(got, msg = "NA strata must error, not silently drop samples")
})

run_test("StratCox-T10b: calc_cox_taylor errors on NA in strata", {
  eta  <- .sc$x * .sc$beta
  g_na <- .sc$g3; g_na[1] <- NA
  got  <- tryCatch({ calc_cox_taylor(eta, .sc$time, .sc$status, strata = g_na); FALSE },
                   error = function(e) TRUE)
  assert_true(got, msg = "NA strata must error, not silently drop samples")
})

run_test("StratCox-T11a: fit_cox_on_yf errors on NA in strata_id", {
  n <- 40L; p <- 15L; K <- 2L
  set.seed(5); Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1L, 0.7)
  g <- rep(c("A", "B"), each = n / 2); g[3] <- NA
  got <- tryCatch({
    suppressMessages(fit_cox_on_yf(Y, time, status, K = K, max_iter = 2L,
                                   verbose = FALSE, strata_id = g)); FALSE
  }, error = function(e) TRUE)
  assert_true(got, msg = "NA strata_id must error before fitting")
})

run_test("StratCox-T11b: fit_supervised_mf_modular errors on NA in strata_id", {
  n <- 40L; p <- 15L; K <- 2L
  set.seed(6); Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1L, 0.7)
  g <- rep(c("A", "B"), each = n / 2); g[3] <- NA
  got <- tryCatch({
    suppressMessages(fit_supervised_mf_modular(Y, time, status, K = K, max_iter = 2L,
                                               verbose = FALSE, strata_id = g)); FALSE
  }, error = function(e) TRUE)
  assert_true(got, msg = "NA strata_id must error before fitting")
})

# =============================================================================
# Section 9: martingale-residual oracle for the stratified score u
#   The Cox score u_i = delta_i - theta_i*H_i is exactly the martingale residual,
#   so survival::residuals(type="martingale") is an independent oracle for u
#   (the coxph cross-checks above validate logPL only).
# =============================================================================

run_test("StratCox-T12a: calc_cox_taylor_yf u == coxph martingale residuals (stratified)", {
  eta <- .sc$x * .sc$beta
  res <- calc_cox_taylor_yf(eta, .sc$time, .sc$status, strata = .sc$g3)
  df  <- data.frame(time = .sc$time, status = .sc$status, x = .sc$x, g = .sc$g3)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x + strata(g), data = df,
          init = .sc$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(max(abs(res$u - residuals(fit, type = "martingale"))), 0, tol = 1e-8,
              msg = "stratified score u must equal coxph martingale residuals")
})

run_test("StratCox-T12b: calc_cox_taylor u == coxph martingale residuals (stratified)", {
  eta <- .sc$x * .sc$beta
  res <- calc_cox_taylor(eta, .sc$time, .sc$status, strata = .sc$g3)
  df  <- data.frame(time = .sc$time, status = .sc$status, x = .sc$x, g = .sc$g3)
  fit <- suppressWarnings(
    coxph(Surv(time, status) ~ x + strata(g), data = df,
          init = .sc$beta, iter.max = 0, ties = "breslow")
  )
  assert_near(max(abs(res$u - residuals(fit, type = "martingale"))), 0, tol = 1e-8,
              msg = "stratified score u must equal coxph martingale residuals")
})

# =============================================================================
# Section 13: fit_cox_on_yf() Phase C orientation fix (review finding, Step 2,
# fixed 2026-09-04 -- see DECISIONS.md). Phase C's concordance() call now
# passes reverse=TRUE, the correct semantics for a Cox risk score. This is
# the single, frozen, training-data-only orientation decision for the fit --
# no downstream evaluator may re-derive it from the data it later scores.
# =============================================================================

run_test("StratCox-T13: fit_cox_on_yf sign_correction=TRUE guarantees correct-direction training concordance >= 0.5 across seeds", {
  # This is exactly the property Phase C exists to guarantee. Before the
  # 2026-09-04 fix, the concordance() check omitted reverse=TRUE and was
  # therefore inverted -- a fit could converge to a genuinely anti-hazard
  # orientation (true reverse-direction concordance well below 0.5) while
  # Phase C's own (uncorrected) check read "already fine, don't flip." This
  # test would fail under the pre-fix code for at least one of these seeds.
  n <- 60L; p <- 15L; K <- 3L
  set.seed(202)
  L_true <- matrix(rexp(n * K, 1), n, K)
  F_true <- matrix(rexp(p * K, 1), p, K)
  beta_true <- c(1.5, -1.2, 0.8)
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.3), n, p)
  eta_true <- L_true %*% beta_true
  time   <- rexp(n, rate = exp(0.5 * scale(eta_true)[, 1]))
  status <- rbinom(n, 1L, 0.8)

  for (s in c(1, 2, 3, 4, 5)) {
    set.seed(s)
    fit <- fit_cox_on_yf(Y, time, status, K = K, max_iter = 40, verbose = FALSE,
                          sign_correction = TRUE)
    ZF  <- Y %*% sweep(fit$EF, 2, fit$EF_norms, "/")
    eta_final <- as.vector(ZF %*% fit$EBeta)
    c_true <- as.numeric(concordance(Surv(time, status) ~ eta_final, reverse = TRUE)$concordance)
    assert_true(c_true >= 0.5 - 1e-8,
                msg = sprintf("seed %d: correct-direction training concordance below 0.5 (%.4f) -- Phase C failed to flip", s, c_true))
  }
})
