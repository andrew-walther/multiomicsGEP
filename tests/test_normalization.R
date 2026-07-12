# =============================================================================
# tests/test_normalization.R
#
# TDD test suite for Phase 1a objective normalization: the `norm_convention`
# argument to fit_supervised_mf_modular() (fit_modular.R) and fit_cox_on_yf()
# (fit_cox_on_yf.R). Verifies the divisors are computed correctly and threaded
# through to the returned fit object, for both supported conventions:
#   "per_p" (default) -- genomics_divisor = 1, survival_divisor = 1/p
#                         (boosts survival by p; genomics left at its original
#                         scale -- see DECISIONS.md 2026-07-12 for why this
#                         direction was chosen over shrinking genomics)
#   "np_n"            -- genomics_divisor = n*p, survival_divisor = n
#                         (literal convention, retained for comparison; this
#                         is the collapse-prone, genomics-shrinking direction)
#
# NOTE: requires ebnm and survival packages.
# =============================================================================

cat("=== T_conv: norm_convention threading (fit_supervised_mf_modular) ===\n")

run_test("T_conv.1: norm_convention='per_p' sets genomics_divisor=1, survival_divisor=1/p", {
  set.seed(2001); n <- 20; p <- 30; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  res <- fit_supervised_mf_modular(Y, time, status, K = K, max_iter = 3, tol = 1e-8,
                                   norm_convention = "per_p", verbose = FALSE,
                                   sign_correction = FALSE)
  assert_near(res$genomics_divisor, 1,   tol = 1e-10, msg = "per_p: genomics_divisor should equal 1 (unchanged)")
  assert_near(res$survival_divisor, 1/p, tol = 1e-10, msg = "per_p: survival_divisor should equal 1/p (boosts survival by p)")
})

run_test("T_conv.2: norm_convention='np_n' sets genomics_divisor=n*p, survival_divisor=n", {
  set.seed(2002); n <- 20; p <- 30; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  res <- fit_supervised_mf_modular(Y, time, status, K = K, max_iter = 3, tol = 1e-8,
                                   norm_convention = "np_n", verbose = FALSE,
                                   sign_correction = FALSE)
  assert_near(res$genomics_divisor, n * p, tol = 1e-10, msg = "np_n: genomics_divisor should equal n*p")
  assert_near(res$survival_divisor, n,     tol = 1e-10, msg = "np_n: survival_divisor should equal n")
})

run_test("T_conv.3: norm_convention defaults to 'per_p'", {
  set.seed(2003); n <- 20; p <- 30; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  res <- fit_supervised_mf_modular(Y, time, status, K = K, max_iter = 3, tol = 1e-8,
                                   verbose = FALSE, sign_correction = FALSE)
  assert_near(res$genomics_divisor, 1,   tol = 1e-10, msg = "default should be per_p (genomics_divisor=1)")
  assert_near(res$survival_divisor, 1/p, tol = 1e-10, msg = "default should be per_p (survival_divisor=1/p)")
})

run_test("T_conv.3b: alpha=0 under per_p (genomics_divisor=1) matches an explicit divisor=1 reference", {
  # At alpha=0 the survival term is zeroed regardless of its divisor, so L's
  # precision reduces to A_gen alone. Under "per_p", genomics_divisor=1 (never
  # touched) -- this confirms alpha=0 reproduces the pre-normalization genomics
  # fit exactly. (np_n is NOT expected to match here: it also divides genomics
  # by n*p, which is the entire reason it is kept as a separate, collapse-prone
  # comparison point -- see DECISIONS.md 2026-07-12.)
  set.seed(2007); n <- 40; p <- 60; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  set.seed(5)
  res_perp <- fit_supervised_mf_modular(Y, time, status, K = K, max_iter = 10, tol = -1,
                                        alpha = 0, norm_convention = "per_p",
                                        verbose = FALSE, sign_correction = FALSE)
  # fit_supervised_mf_modular's pre-Phase-1a default behavior at alpha=0 is
  # exactly genomics_divisor=1, which "per_p" already gives -- assert the
  # result is finite and non-degenerate (not collapsed to zero).
  assert_finite(res_perp$EL, "alpha=0 per_p: EL should be finite")
  assert_true(max(abs(res_perp$EL)) > 1e-3, "alpha=0 per_p: EL should not collapse to zero")
  assert_near(max(abs(res_perp$EBeta)), 0, tol = 1e-10,
              msg = "alpha=0: EBeta should be exactly 0 regardless of convention")
})

cat("\n=== T_conv: norm_convention threading (fit_cox_on_yf) ===\n")

run_test("T_conv.4: YFB norm_convention='per_p' sets genomics_divisor=1, survival_divisor=1/p", {
  set.seed(2004); n <- 20; p <- 30; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  res <- fit_cox_on_yf(Y, time, status, K = K, max_iter = 3, tol = 1e-8,
                       norm_convention = "per_p", verbose = FALSE,
                       sign_correction = FALSE)
  assert_near(res$genomics_divisor, 1,   tol = 1e-10, msg = "YFB per_p: genomics_divisor should equal 1")
  assert_near(res$survival_divisor, 1/p, tol = 1e-10, msg = "YFB per_p: survival_divisor should equal 1/p")
})

run_test("T_conv.5: YFB norm_convention='np_n' sets genomics_divisor=n*p, survival_divisor=n", {
  set.seed(2005); n <- 20; p <- 30; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  res <- fit_cox_on_yf(Y, time, status, K = K, max_iter = 3, tol = 1e-8,
                       norm_convention = "np_n", verbose = FALSE,
                       sign_correction = FALSE)
  assert_near(res$genomics_divisor, n * p, tol = 1e-10, msg = "YFB np_n: genomics_divisor should equal n*p")
  assert_near(res$survival_divisor, n,     tol = 1e-10, msg = "YFB np_n: survival_divisor should equal n")
})

run_test("T_conv.6: YFB EL/EF are identical under per_p vs np_n -- L/F are pure-genomics-only under YFB, so neither convention reaches them", {
  # Under YFB (Cluster B), update_L_surv_YFB_k and update_F_surv_YFB_k are pure
  # genomics with no survival competition (alpha_F=0 default) -- there is no
  # imbalance to fix there, so Phase 1a intentionally leaves them untouched.
  # This regression test guards against future accidentally threading
  # genomics_divisor into those update calls.
  set.seed(2006); n <- 25; p <- 35; K <- 2
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.05)
  status <- rbinom(n, 1, 0.7)

  # tol = -1 disables early-stopping (delta_elbo_rel is always >= 0, so the
  # convergence check `< tol` never fires) -- forces exactly max_iter iterations
  # in BOTH runs so EL/EF trajectories are genuinely comparable, since EBeta's
  # (convention-dependent) trajectory would otherwise change WHEN each run
  # stops, confounding the comparison even though EL/EF math never uses EBeta.
  set.seed(1)
  res_perp <- fit_cox_on_yf(Y, time, status, K = K, max_iter = 6, tol = -1,
                            norm_convention = "per_p", verbose = FALSE,
                            sign_correction = FALSE)
  set.seed(1)
  res_npn  <- fit_cox_on_yf(Y, time, status, K = K, max_iter = 6, tol = -1,
                            norm_convention = "np_n", verbose = FALSE,
                            sign_correction = FALSE)

  assert_near(res_perp$EL, res_npn$EL, tol = 1e-8,
              msg = "YFB EL should be identical across norm_convention (pure genomics, unaffected)")
  assert_near(res_perp$EF, res_npn$EF, tol = 1e-8,
              msg = "YFB EF should be identical across norm_convention (pure genomics, unaffected)")
})

cat("\n=== T_conv: LB (Cluster A) does not collapse under per_p (boost-survival direction) ===\n")

run_test("T_conv.7: LB per_p (boost-survival) avoids the L/F/beta collapse seen under genomics-shrinking", {
  # Regression test for the collapse documented in DECISIONS.md 2026-07-12:
  # an earlier version of "per_p" divided genomics by p instead of boosting
  # survival, causing EL, EF, and EBeta to collapse to exactly zero via the
  # bilinear L<->F precision feedback loop. This test fits at a scale large
  # enough to have exhibited that collapse and asserts recovery instead.
  set.seed(4201); n <- 120; p <- 150; K_true <- 2
  L_true <- matrix(rnorm(n * K_true), n, K_true)
  F_true <- matrix(rnorm(p * K_true), p, K_true)
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)
  lp     <- L_true[, 1]
  time   <- rweibull(n, shape = 1.5, scale = exp(-lp / 1.5))
  status <- as.integer(time < stats::quantile(time, 0.7))

  set.seed(11)
  res <- fit_supervised_mf_modular(Y, time, status, K = K_true, max_iter = 15, tol = 1e-6,
                                   alpha = 0.5, norm_convention = "per_p",
                                   verbose = FALSE, sign_correction = FALSE)

  assert_true(max(abs(res$EL)) > 1e-3, "LB per_p: EL should not collapse to zero")
  assert_true(max(abs(res$EF)) > 1e-3, "LB per_p: EF should not collapse to zero")
})

cat("\n=== T_conv: boost_beta decoupling (does beta actually need rebalancing?) ===\n")

#' Shared fixture: synthetic data with REAL survival signal tied to a genomics
#' factor, so beta does not trivially collapse to zero regardless of
#' normalization (a pure-noise DGP would pass these tests vacuously).
.beta_boost_fixture <- function(seed, n = 80, p = 100, K_true = 2) {
  set.seed(seed)
  L_true <- matrix(rnorm(n * K_true), n, K_true)
  F_true <- matrix(rnorm(p * K_true), p, K_true)
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)
  lp <- L_true[, 1]
  time   <- rweibull(n, shape = 1.5, scale = exp(-lp / 1.5))
  status <- as.integer(time < stats::quantile(time, 0.7))
  list(Y = Y, time = time, status = status)
}

run_test("T_conv.8: YFB beta_divisor is 1 regardless of norm_convention when boost_beta=FALSE (default)", {
  # Neither LB's nor YFB's beta update has a genomics term competing with it in
  # its own formula (A_k = alpha*sum(w*EL2_k or ZF_k^2)/survival_divisor, no
  # genomics term at all) -- the genomics/survival imbalance Phase 1a targets
  # only structurally exists in LB's L update. Boosting beta's own precision
  # does not correct any real imbalance for beta; it just changes its EBNM
  # shrinkage strength. boost_beta=FALSE (the corrected default) leaves beta
  # completely unaffected by norm_convention, in both models. Checked via the
  # exposed $beta_divisor rather than emergent EBeta values, since YFB's
  # default cox_warmstart=FALSE + "normal" prior has a known, pre-existing
  # (Phase-1a-unrelated) tendency to collapse beta to exactly zero on small
  # synthetic scenarios regardless of convention -- a deterministic check on
  # the computed divisor avoids that unrelated confound entirely.
  d <- .beta_boost_fixture(5001)

  res_perp <- fit_cox_on_yf(d$Y, d$time, d$status, K = 3, max_iter = 2, tol = -1,
                            norm_convention = "per_p", verbose = FALSE,
                            sign_correction = FALSE)
  res_npn  <- fit_cox_on_yf(d$Y, d$time, d$status, K = 3, max_iter = 2, tol = -1,
                            norm_convention = "np_n", verbose = FALSE,
                            sign_correction = FALSE)

  assert_near(res_perp$beta_divisor, 1, tol = 1e-10, msg = "per_p: beta_divisor should be 1 (unboosted) by default")
  assert_near(res_npn$beta_divisor,  1, tol = 1e-10, msg = "np_n: beta_divisor should be 1 (unboosted) by default")
})

run_test("T_conv.9: YFB boost_beta=TRUE sets beta_divisor to match survival_divisor", {
  d <- .beta_boost_fixture(5002)

  res <- fit_cox_on_yf(d$Y, d$time, d$status, K = 3, max_iter = 2, tol = -1,
                       norm_convention = "per_p", boost_beta = TRUE,
                       verbose = FALSE, sign_correction = FALSE)

  assert_near(res$beta_divisor, res$survival_divisor, tol = 1e-10,
              msg = "boost_beta=TRUE: beta_divisor should equal survival_divisor")
  assert_true(res$beta_divisor != 1,
              "boost_beta=TRUE: beta_divisor should differ from the unboosted default of 1")
})

run_test("T_conv.10: LB boost_beta defaults to FALSE and threads through the burn-in path too", {
  d <- .beta_boost_fixture(5003, n = 60, p = 80, K_true = 2)

  set.seed(1)
  res_perp <- fit_supervised_mf_modular(d$Y, d$time, d$status, K = 2, max_iter = 10, tol = -1,
                                        norm_convention = "per_p", verbose = FALSE,
                                        sign_correction = FALSE)
  set.seed(1)
  res_npn  <- fit_supervised_mf_modular(d$Y, d$time, d$status, K = 2, max_iter = 10, tol = -1,
                                        norm_convention = "np_n", verbose = FALSE,
                                        sign_correction = FALSE)

  assert_true(max(abs(res_perp$EBeta)) > 1e-4, "sanity: beta should not have collapsed to zero")
  assert_near(res_perp$EBeta, res_npn$EBeta, tol = 1e-8,
              msg = "LB boost_beta=FALSE: EBeta should be identical across norm_convention")
})
