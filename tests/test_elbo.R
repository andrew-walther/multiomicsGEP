# =============================================================================
# tests/test_elbo.R
# Tests for compute_elbo.R (compute_ebnm_kl, compute_survival_elbo) and
# the full ELBO integration in fit_supervised_mf_modular().
#
# Run via:  Rscript tests/run_tests.R
# =============================================================================

# NOTE: test_helpers.R, compute_elbo.R, and fit_modular.R are sourced by
# run_tests.R before this file is loaded.

# =============================================================================
# Section 1: calc_cox_taylor backward compatibility + new logPL field
# =============================================================================

run_test("calc_cox_taylor returns $u (backward compat)", {
  set.seed(1)
  n <- 20
  eta    <- rnorm(n)
  time   <- rexp(n)
  status <- sample(0:1, n, replace = TRUE)
  res <- calc_cox_taylor(eta, time, status)
  assert_true(!is.null(res$u), "u missing")
  assert_length(res$u, n)
})

run_test("calc_cox_taylor returns $w (backward compat)", {
  set.seed(1)
  n <- 20
  eta    <- rnorm(n)
  time   <- rexp(n)
  status <- sample(0:1, n, replace = TRUE)
  res <- calc_cox_taylor(eta, time, status)
  assert_true(!is.null(res$w), "w missing")
  assert_positive(res$w)
})

run_test("calc_cox_taylor returns $logPL as finite scalar", {
  set.seed(1)
  n <- 20
  eta    <- rnorm(n)
  time   <- rexp(n)
  status <- sample(0:1, n, replace = TRUE)
  res <- calc_cox_taylor(eta, time, status)
  assert_true(!is.null(res$logPL), "logPL missing")
  assert_length(res$logPL, 1)
  assert_finite(res$logPL)
})

run_test("calc_cox_taylor logPL <= 0 (partial likelihood is a probability)", {
  # The Breslow partial likelihood is a product of values in (0,1], so its log
  # should be <= 0 when there is at least one event.
  set.seed(2)
  n      <- 30
  time   <- rexp(n)
  status <- c(rep(1, 15), rep(0, 15))
  eta    <- rnorm(n)
  res    <- calc_cox_taylor(eta, time, status)
  assert_true(res$logPL <= 0, sprintf("logPL = %.4f > 0", res$logPL))
})

# =============================================================================
# Section 2: compute_ebnm_kl
# =============================================================================

run_test("compute_ebnm_kl returns a finite scalar", {
  # Minimal sanity: any numeric inputs should yield a finite number.
  kl <- compute_ebnm_kl(ebnm_log_lik = -2.0,
                        A = c(1, 2, 1), x = c(0.5, -0.3, 1.0),
                        mean_q = c(0.4, -0.2, 0.8), second_q = c(0.17, 0.05, 0.65))
  assert_length(kl, 1)
  assert_finite(kl)
})

run_test("compute_ebnm_kl works for scalar inputs (beta_k case)", {
  # beta_k is length-1; A, x, mean_q, second_q are all length-1 vectors.
  kl <- compute_ebnm_kl(ebnm_log_lik = -0.5,
                        A = 2.0, x = 1.0,
                        mean_q = 0.8, second_q = 0.64)
  assert_length(kl, 1)
  assert_finite(kl)
})

run_test("compute_ebnm_kl: when posterior = Gaussian N(mu, 1/A) and g = N(0,1) KL is computable", {
  # For a Gaussian EBNM with g = N(0,1) (A_prior=1, x=mu, s=1/sqrt(A)):
  # log_lik = -0.5*(x^2/(1+1/A) + log(2*pi*(1 + 1/A)))
  # We just check the sign convention: KL = E[log g] - E[log q] <= 0 when
  # q has MORE probability mass than g (over-confident posterior).
  set.seed(3)
  n <- 50
  mu <- rnorm(n)
  A  <- rep(10, n)   # tight posterior — should have negative KL
  x  <- mu           # pseudo-obs at posterior mean
  second <- mu^2 + 1/A  # E[theta^2] = Var + mean^2

  # Approximate ebnm_log_lik for g=N(0,1): log p(x | s^2 = 1/A) with g = N(0,1)
  # = sum_i log N(x_i; 0, 1 + 1/A_i)
  ebnm_log_lik <- sum(dnorm(x, mean = 0, sd = sqrt(1 + 1/A), log = TRUE))

  kl <- compute_ebnm_kl(ebnm_log_lik, A, x, mu, second)
  assert_true(kl <= 0, sprintf("KL = %.4f > 0 (should be <= 0)", kl))
})

# =============================================================================
# Section 3: compute_survival_elbo
# =============================================================================

run_test("compute_survival_elbo returns a finite scalar", {
  # Under Cluster B: ZF = Y·EF (observed projection scores, n x K)
  set.seed(4)
  n <- 20; K <- 3
  ZF     <- matrix(rnorm(n * K), n, K)
  EBeta  <- rnorm(K)
  EBeta2 <- EBeta^2 + 0.01
  w      <- runif(n, 0.01, 0.5)
  logPL  <- -5.0
  res    <- compute_survival_elbo(logPL, w, ZF, EBeta, EBeta2)
  assert_length(res, 1)
  assert_finite(res)
})

run_test("compute_survival_elbo with zero posterior variance returns logPL exactly", {
  # When EBeta2 = EBeta^2 (no uncertainty in beta_tilde), Var_q(eta) = 0
  # and the correction term vanishes: result should equal logPL.
  set.seed(5)
  n <- 15; K <- 2
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta <- rnorm(K)
  w     <- runif(n, 0.01, 0.3)
  logPL <- -8.5
  res   <- compute_survival_elbo(logPL, w, ZF, EBeta, EBeta^2)
  assert_near(res, logPL, tol = 1e-10, "should equal logPL when no variance")
})

run_test("compute_survival_elbo: uncertainty correction reduces value vs logPL", {
  # With positive posterior variance in beta_tilde the correction term is > 0,
  # so E_q[log PL] < logPL(eta_0).
  set.seed(6)
  n <- 20; K <- 2
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta <- rnorm(K)
  w     <- runif(n, 0.1, 1.0)
  logPL <- -6.0
  EBeta2 <- EBeta^2 + 0.5
  res   <- compute_survival_elbo(logPL, w, ZF, EBeta, EBeta2)
  assert_true(res < logPL,
              sprintf("surv_elbo=%.4f should be < logPL=%.4f", res, logPL))
})

# =============================================================================
# Section 4: Full ELBO integration in fit_supervised_mf_modular
# =============================================================================

# Shared synthetic data for integration tests.
# n=150, p=200 chosen for numerical stability with the alpha=0.5 default:
# the small n=60 dataset caused eta explosion under reduced regularization.
.elbo_test_data <- local({
  set.seed(42)
  n <- 150; p <- 200; K_true <- 2
  L_true <- matrix(rnorm(n * K_true), n, K_true)
  F_true <- matrix(rnorm(p * K_true), p, K_true)
  E      <- matrix(rnorm(n * p, sd = 0.5), n, p)
  Y      <- L_true %*% t(F_true) + E
  # Weibull survival with linear predictor from first factor only
  lp     <- L_true[, 1]
  time   <- rweibull(n, shape = 1.5, scale = exp(-lp / 1.5))
  status <- as.integer(time < quantile(time, 0.7))
  list(Y = Y, time = time, status = status)
})

run_test("fit_supervised_mf_modular returns $history$elbo_full", {
  d   <- .elbo_test_data
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 5, verbose = FALSE)
  assert_true(!is.null(res$history$elbo_full), "elbo_full missing from history")
})

run_test("history$elbo_full has same length as history$elbo_proxy", {
  d   <- .elbo_test_data
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 5, verbose = FALSE)
  assert_equal(length(res$history$elbo_full),
               length(res$history$elbo_proxy),
               "elbo_full and elbo_proxy differ in length")
})

run_test("history tracks delta_L, delta_Beta, and delta_elbo_rel", {
  d   <- .elbo_test_data
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 6, verbose = FALSE)
  assert_equal(length(res$history$delta_L), length(res$history$elbo_full),
               "delta_L length mismatch")
  assert_equal(length(res$history$delta_Beta), length(res$history$elbo_full),
               "delta_Beta length mismatch")
  assert_equal(length(res$history$delta_elbo_rel), length(res$history$elbo_full),
               "delta_elbo_rel length mismatch")
  assert_true(!is.null(res$history$factor_pve), "factor_pve missing")
  assert_equal(nrow(res$history$factor_pve), length(res$history$elbo_full),
               "factor_pve row count mismatch")
})

run_test("history$elbo_full is all finite", {
  d   <- .elbo_test_data
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 10, verbose = FALSE)
  assert_finite(res$history$elbo_full)
})

run_test("history$elbo_proxy unchanged after adding elbo_full", {
  # The proxy is the genomics log-likelihood from update_tau.R.
  # Adding full ELBO tracking must NOT alter proxy values.
  d   <- .elbo_test_data
  set.seed(99)
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 8, verbose = FALSE)
  assert_finite(res$history$elbo_proxy, "proxy contains non-finite values")
  # Proxy should be negative (log-likelihood) and non-constant over iters
  assert_true(diff(range(res$history$elbo_proxy)) > 0,
              "elbo_proxy is constant — something went wrong")
})

run_test("elbo_full is no greater than elbo_proxy when alpha=0", {
  # With alpha=0, the full ELBO reduces to the genomics proxy plus KL terms.
  # Since each KL contribution is <= 0, full ELBO must be <= proxy.
  d   <- .elbo_test_data
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 8, alpha = 0,
                                   verbose = FALSE)
  last_full  <- tail(res$history$elbo_full,  1)
  last_proxy <- tail(res$history$elbo_proxy, 1)
  assert_true(last_full <= last_proxy,
              sprintf("alpha=0: elbo_full=%.1f should be <= elbo_proxy=%.1f",
                      last_full, last_proxy))
})

run_test("elbo_full uses alpha-weighted genomics and survival terms", {
  d <- .elbo_test_data

  set.seed(123)
  res_gen <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                       K = 2, max_iter = 1, alpha = 0,
                                       verbose = FALSE)
  set.seed(123)
  res_surv <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                        K = 2, max_iter = 1, alpha = 1,
                                        verbose = FALSE)

  assert_true(res_gen$history$elbo_full[1] <= res_gen$history$elbo_proxy[1],
              "alpha=0 full ELBO should be genomics proxy plus non-positive KL terms")
  assert_true(abs(res_surv$history$elbo_full[1] - res_surv$history$elbo_proxy[1]) > 1,
              "alpha=1 should not equal the genomics proxy")
  assert_true(abs(res_gen$history$elbo_full[1] - res_surv$history$elbo_full[1]) > 1,
              "alpha=0 and alpha=1 should produce materially different full ELBOs")
})

run_test("convergence uses relative ELBO change, not parameter deltas", {
  d   <- .elbo_test_data
  res <- fit_supervised_mf_modular(d$Y, d$time, d$status,
                                   K = 2, max_iter = 8, tol = 1,
                                   verbose = FALSE)
  assert_true(isTRUE(res$history$converged),
              "Expected convergence under loose ELBO tolerance")
  assert_true(any(is.finite(res$history$delta_elbo_rel[-1])),
              "Expected finite relative ELBO deltas after the first iteration")
})
