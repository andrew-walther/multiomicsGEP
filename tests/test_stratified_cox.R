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
