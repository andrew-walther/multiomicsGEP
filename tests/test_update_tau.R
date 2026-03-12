# =============================================================================
# tests/test_update_tau.R
#
# TDD test suite for the tau (precision) CAVI update.
# Written BEFORE implementation (red phase).
#
# Tests cover:
#   T1: Mathematical identity checks (Var_Term, R2_bar, Tau)
#   T2: Perfect reconstruction (zero residual)
#   T3: Known noise recovery (heteroscedastic, sample size effect)
#   T4: Variance correction (Var_Term inflates R2_bar)
#   T5: ELBO proxy (formula, finiteness, optimality)
#   T6: Dimension and shape
#   T7: Numerical stability (zeros, extremes, floor)
#   T8: Pipeline interaction (valid p-vector, downstream use)
#   T9: Consistency with V2.R inline code
#
# NOTE: update_tau.R has no ebnm dependency
# DO NOT source test_helpers.R or update_tau.R here - done by run_tests.R
# =============================================================================

cat("=== T1: Mathematical Identity Checks ===\n")

run_test("T1.1: Var_Term = EL2 %*% t(EF2) - EL^2 %*% t(EF^2)", {
  set.seed(101); n <- 4; p <- 6; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0, 0.5), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0, 0.5), p, K)

  result <- compute_var_term(EL, EL2, EF, EF2)
  expected <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
  assert_near(result, expected, tol = 1e-12, msg = "Var_Term formula mismatch")
})

run_test("T1.2: Var_Term >= 0 element-wise (within floating-point tolerance)", {
  set.seed(102); n <- 4; p <- 6; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.5), n, K)  # EL2 > EL^2
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.5), p, K)  # EF2 > EF^2

  result <- compute_var_term(EL, EL2, EF, EF2)
  assert_true(all(result >= -1e-12),
              sprintf("Var_Term has negative entries: min = %.3e", min(result)))
})

run_test("T1.3: R2_bar = (Y - EL*EF')^2 + Var_Term", {
  set.seed(103); n <- 4; p <- 6; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0, 0.3), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0, 0.3), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  result <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
  var_term <- compute_var_term(EL, EL2, EF, EF2)
  expected <- (Y - EL %*% t(EF))^2 + var_term
  assert_near(result, expected, tol = 1e-12, msg = "R2_bar formula mismatch")
})

run_test("T1.4: R2_bar >= 0 element-wise", {
  set.seed(104); n <- 4; p <- 6; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.5), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.5), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  result <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
  assert_true(all(result >= -1e-12),
              sprintf("R2_bar has negative entries: min = %.3e", min(result)))
})

run_test("T1.5: Tau[j] = n / colSums(R2_bar)[j] per-feature", {
  set.seed(105); n <- 10; p <- 8; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.3), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.3), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  R2_bar <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
  expected_Tau <- n / pmax(colSums(R2_bar), n * 1e-8)
  assert_near(res$Tau, expected_Tau, tol = 1e-12, msg = "Tau per-feature mismatch")
})

cat("\n=== T2: Perfect Reconstruction ===\n")

run_test("T2.1: Var_Term = 0 when EL2 = EL^2 AND EF2 = EF^2", {
  set.seed(201); n <- 5; p <- 4; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2                                # zero posterior variance
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2                                # zero posterior variance

  result <- compute_var_term(EL, EL2, EF, EF2)
  assert_near(result, matrix(0, n, p), tol = 1e-12,
              msg = "Var_Term should be zero when no posterior variance")
})

run_test("T2.2: R2_bar = 0 when Y = EL*EF' exactly and EL2=EL^2, EF2=EF^2", {
  set.seed(202); n <- 5; p <- 4; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2
  Y    <- EL %*% t(EF)                        # perfect reconstruction

  result <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
  assert_near(result, matrix(0, n, p), tol = 1e-12,
              msg = "R2_bar should be zero for perfect reconstruction")
})

run_test("T2.3: Tau hits ceiling (1/tau_floor) when R2_bar = 0", {
  set.seed(203); n <- 5; p <- 4; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2
  Y    <- EL %*% t(EF)

  tau_floor <- 1e-8
  res <- update_tau(Y, EL, EL2, EF, EF2, tau_floor = tau_floor)
  # When R2_bar = 0, colSums(R2_bar) = 0, pmax(0, n*tau_floor) = n*tau_floor
  # Tau = n / (n * tau_floor) = 1/tau_floor
  expected_ceiling <- 1 / tau_floor
  assert_near(res$Tau, rep(expected_ceiling, p), tol = 1e-4,
              msg = "Tau should hit ceiling 1/tau_floor for zero residual")
})

cat("\n=== T3: Known Noise Recovery ===\n")

run_test("T3.1: Tau close to true 1/sigma_j^2 for known noise (n=200)", {
  set.seed(301); n <- 200; p <- 50; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2                                # no posterior variance for clean test
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2
  sigma_true <- rep(0.5, p)                   # homoscedastic noise
  E    <- matrix(rnorm(n * p, sd = rep(sigma_true, each = n)), n, p)
  Y    <- EL %*% t(EF) + E

  res <- update_tau(Y, EL, EL2, EF, EF2)
  true_tau <- 1 / sigma_true^2
  # Allow 50% relative error for n=200
  rel_error <- abs(res$Tau - true_tau) / true_tau
  mean_rel_error <- mean(rel_error)
  assert_true(mean_rel_error < 0.50,
              sprintf("Mean relative Tau error = %.3f (expected < 0.50)", mean_rel_error))
})

run_test("T3.2: Heteroscedastic: high-precision group has higher Tau", {
  set.seed(302); n <- 200; p <- 100; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2
  sigma <- c(rep(1.0, 50), rep(0.3, 50))     # group 1: sigma=1, group 2: sigma=0.3
  E    <- matrix(rnorm(n * p, sd = rep(sigma, each = n)), n, p)
  Y    <- EL %*% t(EF) + E

  res <- update_tau(Y, EL, EL2, EF, EF2)
  mean_tau_group1 <- mean(res$Tau[1:50])
  mean_tau_group2 <- mean(res$Tau[51:100])
  assert_true(mean_tau_group2 > mean_tau_group1,
              sprintf("Group 2 (sigma=0.3) should have higher Tau: group1=%.2f, group2=%.2f",
                      mean_tau_group1, mean_tau_group2))
})

run_test("T3.3: More samples -> better Tau estimation (lower relative error)", {
  set.seed(303); p <- 30; K <- 2
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2
  sigma_true <- rep(0.8, p)
  true_tau <- 1 / sigma_true^2

  # Small sample: n=50
  n_small <- 50
  EL_small  <- matrix(rnorm(n_small * K), n_small, K)
  EL2_small <- EL_small^2
  E_small <- matrix(rnorm(n_small * p, sd = rep(sigma_true, each = n_small)), n_small, p)
  Y_small <- EL_small %*% t(EF) + E_small
  res_small <- update_tau(Y_small, EL_small, EL2_small, EF, EF2)
  err_small <- mean(abs(res_small$Tau - true_tau) / true_tau)

  # Large sample: n=500
  n_large <- 500
  EL_large  <- matrix(rnorm(n_large * K), n_large, K)
  EL2_large <- EL_large^2
  E_large <- matrix(rnorm(n_large * p, sd = rep(sigma_true, each = n_large)), n_large, p)
  Y_large <- EL_large %*% t(EF) + E_large
  res_large <- update_tau(Y_large, EL_large, EL2_large, EF, EF2)
  err_large <- mean(abs(res_large$Tau - true_tau) / true_tau)

  assert_true(err_large < err_small,
              sprintf("Larger n should give lower error: n=500 err=%.3f, n=50 err=%.3f",
                      err_large, err_small))
})

cat("\n=== T4: Variance Correction ===\n")

run_test("T4.1: Posterior variance inflates R2_bar above naive (Y-EL*EF')^2", {
  set.seed(401); n <- 10; p <- 8; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.1, 0.5), n, K)  # posterior variance > 0
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.1, 0.5), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  R2_bar <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
  naive  <- (Y - EL %*% t(EF))^2

  # R2_bar = naive + Var_Term, and Var_Term > 0 when posterior variance > 0
  diff <- R2_bar - naive
  assert_true(all(diff >= -1e-12),
              sprintf("R2_bar should be >= naive; min diff = %.3e", min(diff)))
  assert_true(sum(diff) > 0,
              "Var_Term should make R2_bar strictly larger than naive in total")
})

run_test("T4.2: Ignoring Var_Term (EL2=EL^2, EF2=EF^2) overestimates Tau", {
  set.seed(402); n <- 50; p <- 20; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.1, 0.5), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.1, 0.5), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  # Correct Tau (with posterior variance)
  res_correct <- update_tau(Y, EL, EL2, EF, EF2)

  # Incorrect Tau (ignoring posterior variance)
  res_naive <- update_tau(Y, EL, EL^2, EF, EF^2)

  # Ignoring Var_Term means smaller R2_bar → larger Tau (overestimate precision)
  assert_true(mean(res_naive$Tau) > mean(res_correct$Tau),
              sprintf("Naive Tau (%.3f) should be larger than correct (%.3f)",
                      mean(res_naive$Tau), mean(res_correct$Tau)))
})

run_test("T4.3: More posterior uncertainty -> larger Var_Term", {
  set.seed(403); n <- 10; p <- 8; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EF   <- matrix(rnorm(p * K), p, K)

  # Low posterior variance
  EL2_low <- EL^2 + 0.01
  EF2_low <- EF^2 + 0.01
  vt_low  <- compute_var_term(EL, EL2_low, EF, EF2_low)

  # High posterior variance
  EL2_high <- EL^2 + 1.0
  EF2_high <- EF^2 + 1.0
  vt_high  <- compute_var_term(EL, EL2_high, EF, EF2_high)

  assert_true(sum(vt_high) > sum(vt_low),
              sprintf("Higher uncertainty should give larger Var_Term: low=%.3f, high=%.3f",
                      sum(vt_low), sum(vt_high)))
})

cat("\n=== T5: ELBO Proxy ===\n")

run_test("T5.1: ELBO formula: sum(n/2 * log(Tau) - Tau/2 * colSums(R2_bar))", {
  set.seed(501); n <- 20; p <- 10; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.3), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.3), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  R2_bar <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
  expected_elbo <- sum(n / 2 * log(res$Tau) - res$Tau / 2 * colSums(R2_bar))
  assert_near(res$elbo_proxy, expected_elbo, tol = 1e-10,
              msg = "ELBO proxy formula mismatch")
})

run_test("T5.2: ELBO proxy is finite for valid inputs", {
  set.seed(502); n <- 20; p <- 10; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.3), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.3), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  assert_finite(res$elbo_proxy, msg = "ELBO proxy should be finite")
})

run_test("T5.3: ELBO at MLE Tau > ELBO at perturbed Tau", {
  set.seed(503); n <- 50; p <- 15; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.2), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.2), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  R2_bar <- res$R2_bar
  Tau_mle <- res$Tau
  elbo_mle <- res$elbo_proxy

  # Perturb Tau by multiplying by a random factor
  set.seed(504)
  Tau_perturbed <- Tau_mle * exp(rnorm(p, sd = 0.5))  # random perturbation
  elbo_perturbed <- sum(n / 2 * log(Tau_perturbed) - Tau_perturbed / 2 * colSums(R2_bar))

  assert_true(elbo_mle >= elbo_perturbed - 1e-6,
              sprintf("MLE ELBO (%.2f) should be >= perturbed ELBO (%.2f)",
                      elbo_mle, elbo_perturbed))
})

cat("\n=== T6: Dimension and Shape ===\n")

run_test("T6.1: Var_Term dimensions n x p", {
  set.seed(601); n <- 7; p <- 11; K <- 3
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.1
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.1

  result <- compute_var_term(EL, EL2, EF, EF2)
  assert_true(nrow(result) == n && ncol(result) == p,
              sprintf("Expected %dx%d, got %dx%d", n, p, nrow(result), ncol(result)))
})

run_test("T6.2: Tau is p-vector", {
  set.seed(602); n <- 7; p <- 11; K <- 3
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.1
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.1
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  assert_length(res$Tau, p, "Tau should be a p-vector")
})

cat("\n=== T7: Numerical Stability ===\n")

run_test("T7.1: Y = 0 matrix -> Tau is valid (not NaN/Inf)", {
  n <- 5; p <- 4; K <- 2
  EL   <- matrix(1, n, K)
  EL2  <- EL^2 + 0.1
  EF   <- matrix(1, p, K)
  EF2  <- EF^2 + 0.1
  Y    <- matrix(0, n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  assert_finite(res$Tau, msg = "Tau should be finite for Y=0")
  assert_positive(res$Tau, msg = "Tau should be positive for Y=0")
})

run_test("T7.2: n=1, p=1, K=1 -> minimal dimensions work", {
  n <- 1; p <- 1; K <- 1
  EL   <- matrix(0.5, n, K)
  EL2  <- matrix(0.5^2 + 0.01, n, K)
  EF   <- matrix(1.0, p, K)
  EF2  <- matrix(1.0^2 + 0.01, p, K)
  Y    <- matrix(0.8, n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  assert_length(res$Tau, 1, "Tau should have length 1")
  assert_finite(res$Tau, msg = "Tau should be finite for minimal input")
  assert_positive(res$Tau, msg = "Tau should be positive for minimal input")
  assert_finite(res$elbo_proxy, msg = "ELBO proxy should be finite for minimal input")
})

run_test("T7.3: EL with entries ~1e3 -> no overflow in matrix products", {
  set.seed(703); n <- 10; p <- 5; K <- 2
  EL   <- matrix(rnorm(n * K, mean = 0, sd = 1e3), n, K)
  EL2  <- EL^2 + 1e4      # large second moments
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.1
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  assert_finite(res$Tau, msg = "Tau should be finite with large EL entries")
  assert_positive(res$Tau, msg = "Tau should be positive with large EL entries")
  assert_finite(res$elbo_proxy, msg = "ELBO proxy should be finite with large EL entries")
})

run_test("T7.4: tau_floor prevents zero/infinite Tau", {
  n <- 5; p <- 4; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2
  Y    <- EL %*% t(EF)                        # perfect fit → colSums(R2_bar) ~ 0

  tau_floor <- 1e-6
  res <- update_tau(Y, EL, EL2, EF, EF2, tau_floor = tau_floor)
  # Tau should be capped at 1/tau_floor = 1e6 (not Inf)
  assert_true(all(is.finite(res$Tau)),
              "Tau should be finite even when R2_bar is zero")
  assert_true(all(res$Tau <= 1 / tau_floor + 1),
              sprintf("Tau should not exceed 1/tau_floor; max Tau = %.2e", max(res$Tau)))
})

cat("\n=== T8: Pipeline Interaction ===\n")

run_test("T8.1: Tau from update_tau is a valid p-vector of positive values", {
  set.seed(801); n <- 30; p <- 20; K <- 3
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.3), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.3), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  assert_length(res$Tau, p, "Tau must have length p")
  assert_positive(res$Tau, "Tau must be strictly positive")
  assert_finite(res$Tau, "Tau must be finite")
})

run_test("T8.2: Tau can be used in update_L_k's A_L computation without error", {
  # In V2.R, A_L[k] = sum(Tau * EF2[,k])  (element-wise product then sum)
  # Tau (p-vec) * EF2[,k] (p-vec) → p-vec → sum → scalar
  set.seed(802); n <- 20; p <- 15; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.1
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.1
  Y    <- matrix(rnorm(n * p), n, p)

  res <- update_tau(Y, EL, EL2, EF, EF2)
  # Simulate downstream A_L computation for each factor
  for (k in 1:K) {
    A_L_k <- sum(res$Tau * EF2[, k])
    assert_finite(A_L_k, sprintf("A_L[%d] should be finite", k))
    assert_true(A_L_k > 0, sprintf("A_L[%d] should be positive (got %.3e)", k, A_L_k))
  }
})

cat("\n=== T9: V2.R Consistency ===\n")

run_test("T9.1: compute_var_term matches V2.R line 374 exactly", {
  set.seed(777); n <- 50; p <- 30; K <- 3
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.5), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.5), p, K)

  # V2.R line 374 inline
  Var_Term_v2 <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))

  # Modular function
  Var_Term_mod <- compute_var_term(EL, EL2, EF, EF2)

  assert_near(Var_Term_mod, Var_Term_v2, tol = 1e-12,
              msg = "compute_var_term does not match V2.R line 374")
})

run_test("T9.2: update_tau matches V2.R lines 374-376 and 379 exactly", {
  set.seed(777); n <- 50; p <- 30; K <- 3
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + matrix(runif(n * K, 0.01, 0.5), n, K)
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + matrix(runif(p * K, 0.01, 0.5), p, K)
  Y    <- matrix(rnorm(n * p), n, p)

  # V2.R inline code (lines 374-376, 379):
  Var_Term_v2 <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
  R2_bar_v2   <- (Y - EL %*% t(EF))^2 + Var_Term_v2
  Tau_v2      <- n / pmax(colSums(R2_bar_v2), n * 1e-8)
  elbo_v2     <- sum(n / 2 * log(Tau_v2) - Tau_v2 / 2 * colSums(R2_bar_v2))

  # Modular function
  res_mod <- update_tau(Y, EL, EL2, EF, EF2)

  assert_near(res_mod$Var_Term, Var_Term_v2, tol = 1e-12,
              msg = "Var_Term mismatch vs V2.R")
  assert_near(res_mod$R2_bar, R2_bar_v2, tol = 1e-12,
              msg = "R2_bar mismatch vs V2.R")
  assert_near(res_mod$Tau, Tau_v2, tol = 1e-12,
              msg = "Tau mismatch vs V2.R")
  assert_near(res_mod$elbo_proxy, elbo_v2, tol = 1e-10,
              msg = "elbo_proxy mismatch vs V2.R")
})
