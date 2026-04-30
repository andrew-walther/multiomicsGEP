# =============================================================================
# tests/test_update_F.R
#
# TDD test suite for the q(f_j) CAVI update (factor matrix F).
# Written BEFORE implementation (red phase).
#
# Tests cover:
#   T1: Mathematical identity checks (A_F, B_F, x, s, tau cancellation)
#   T2: WLS / non-informative limit
#   T3: Known signal recovery (K=1)
#   T4: Multi-factor recovery (K=5) via update_F_all
#   T5: Null factor shrinkage
#   T6: Tau differential shrinkage
#   T7: Numerical stability (degenerate inputs, extremes)
#   T8: R_k and Gauss-Seidel ordering
#   T9: Consistency with V2.R inline code (lines 334-340)
#
# NOTE: requires ebnm package.  Install with:
#   install.packages("ebnm")
#
# FUNCTION SIGNATURES (from code/update_F.R):
#   update_F_k(Tau, EL_k, EL2_k, R_k, ...)
#     -> list with mean(p-vec), second(p-vec), sd(p-vec),
#        A(p-vec), B(p-vec), x(p-vec), s(p-vec),
#        sum_EL2_k(scalar), ebnm_result
#
#   update_F_all(Y, EL, EL2, EF, EF2, Tau, ...)
#     -> list with EF(p x K), EF2(p x K), details
#
#   compute_R_k(Y, EL, EF, k) -- from update_L.R
#     -> n x p residual matrix with factor k removed
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# DO NOT source test_helpers.R or update_F.R here - done by run_tests.R

# ---- helper: compute R_k locally (no dependency on update_L.R) ----
.compute_R_k_local <- function(Y, EL, EF, k) {
  # R_k = Y - sum_{k' != k} EL[,k'] %*% t(EF[,k'])
  # Equivalent to: Y - EL %*% t(EF) + EL[,k] %*% t(EF[,k])
  Y - EL %*% t(EF) + EL[, k] %o% EF[, k]
}

cat("=== T1: Mathematical Identity Checks ===\n")

run_test("T1.1: A_F[j] = (1-alpha)*Tau[j]*sum(EL2_k) with default alpha=0.5", {
  set.seed(1); n <- 5; p <- 8
  Tau   <- abs(rnorm(p)) + 0.1                   # p-vector
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + runif(n, 0, 0.2)            # second moments
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  sum_EL2_k <- sum(EL2_k)
  expected_A <- pmax(0.5 * Tau * sum_EL2_k, 1e-10)     # p-vector with alpha=0.5
  assert_near(res$A, expected_A, msg = "A_F[j] mismatch")
})

run_test("T1.2: B_F[j] = (1-alpha)*Tau[j]*(t(R_k)%*%EL_k)[j] with default alpha=0.5", {
  set.seed(2); n <- 5; p <- 8
  Tau   <- abs(rnorm(p)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  expected_B <- 0.5 * Tau * as.vector(t(R_k) %*% EL_k)
  assert_near(res$B, expected_B, msg = "B_F[j] mismatch")
})

run_test("T1.3: x is a p-vector (dimension check)", {
  set.seed(3); n <- 5; p <- 8
  Tau   <- abs(rnorm(p)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  assert_length(res$x, p, "x should be p-vector")
})

run_test("T1.4: tau CANCELS in x_j = B_j/A_j (Tau=1 vs Tau=1000 identical)", {
  # x_j = Tau[j] * (t(R_k) %*% EL_k)[j] / (Tau[j] * sum(EL2_k))
  #      = (t(R_k) %*% EL_k)[j] / sum(EL2_k)
  # The Tau[j] cancels in the ratio, so x should be identical regardless of Tau
  set.seed(4); n <- 10; p <- 15
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  Tau_low  <- rep(1, p)
  Tau_high <- rep(1000, p)

  res_low  <- update_F_k(Tau_low,  EL_k, EL2_k, R_k)
  res_high <- update_F_k(Tau_high, EL_k, EL2_k, R_k)

  assert_near(res_low$x, res_high$x, tol = 1e-12,
              msg = "x_j should be identical for Tau=1 vs Tau=1000")
})

run_test("T1.5: tau does NOT cancel in s_j (Tau=1 vs Tau=1000 differ)", {
  # s_j = 1/sqrt(A_j) = 1/sqrt(Tau[j]*sum(EL2_k))
  # s depends on Tau, so s with Tau=1 should differ from Tau=1000
  set.seed(5); n <- 10; p <- 15
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  Tau_low  <- rep(1, p)
  Tau_high <- rep(1000, p)

  res_low  <- update_F_k(Tau_low,  EL_k, EL2_k, R_k)
  res_high <- update_F_k(Tau_high, EL_k, EL2_k, R_k)

  # s_j at Tau=1000 should be much smaller than at Tau=1
  assert_true(all(res_high$s < res_low$s),
              "s_j with Tau=1000 should be smaller than with Tau=1")
  # Verify they are not identical
  assert_true(max(abs(res_low$s - res_high$s)) > 0.01,
              "s_j should differ between Tau=1 and Tau=1000")
})

run_test("T1.6: second moment = sd^2 + mean^2 (element-wise)", {
  set.seed(6); n <- 8; p <- 12
  Tau   <- abs(rnorm(p)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  expected_second <- res$sd^2 + res$mean^2
  assert_near(res$second, expected_second, msg = "second moment identity failed")
})

cat("\n=== T2: WLS / Non-Informative Limit ===\n")

run_test("T2.1: EL2_k = EL_k^2 (zero variance) -> x_j = OLS form", {
  # When EL2_k = EL_k^2 (no posterior variance in L), the EBNM pseudo-obs
  # x_j = (t(R_k) %*% EL_k)[j] / sum(EL_k^2) which is the OLS form
  set.seed(10); n <- 50; p <- 20
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2                 # zero posterior variance
  Tau   <- abs(rnorm(p)) + 0.5
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  ols_x <- as.vector(t(R_k) %*% EL_k) / sum(EL_k^2)
  assert_near(res$x, ols_x, tol = 1e-10,
              msg = "x_j should equal OLS form when EL2=EL^2")
})

run_test("T2.2: large Tau -> posterior mean close to x_j (less shrinkage)", {
  # With high precision (large Tau), the EBNM observation has very low noise
  # s_j = 1/sqrt(Tau*sum_EL2), so the posterior mean should closely match x_j.
  # F_true is non-negative (abs(rnorm)) so x_j > 0; point_exponential prior
  # is permissive for positive pseudo-obs, so the high-SNR limit holds.
  set.seed(11); n <- 50; p <- 30
  EL_k  <- rnorm(n, sd = 1)
  EL2_k <- EL_k^2 + 1e-4         # near-zero variance
  F_true <- abs(rnorm(p, sd = 2))  # non-negative true factor (matches point_exp prior)
  R_k   <- EL_k %o% F_true + matrix(rnorm(n * p, sd = 0.01), n, p)
  Tau   <- rep(1e6, p)            # extremely high precision

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  # With very low noise and positive x_j, posterior mean should be close to x_j
  diffs <- abs(res$mean - res$x)
  assert_true(max(diffs) < 1.0,
              sprintf("max|mean-x| = %.3f, expected close for large Tau", max(diffs)))
})

cat("\n=== T3: Known Signal Recovery K=1 ===\n")

run_test("T3.1: dense factor column recovery (cor > 0.7)", {
  set.seed(20); n <- 100; p <- 50
  F_true <- rnorm(p, sd = 1.5)           # dense factor column
  EL_k   <- rnorm(n)
  EL2_k  <- EL_k^2 + 0.01
  # Y = EL_k %o% F_true + noise
  R_k    <- EL_k %o% F_true + matrix(rnorm(n * p, sd = 0.5), n, p)
  Tau    <- rep(4.0, p)                   # 1/0.5^2 = 4

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  correlation <- cor(res$mean, F_true)
  assert_true(correlation > 0.7,
              sprintf("Dense factor recovery too low: cor = %.3f", correlation))
})

run_test("T3.2: sparse factor (50/500 nonzero) -- recovery of active set", {
  set.seed(21); n <- 100; p <- 500
  F_true <- rep(0, p)
  active <- sample(p, 50)
  F_true[active] <- rnorm(50, sd = 3.0)
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.01
  R_k   <- EL_k %o% F_true + matrix(rnorm(n * p, sd = 0.3), n, p)
  Tau   <- rep(1 / 0.3^2, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  # The top-50 absolute values in estimated F should overlap well with active set
  est_top <- order(abs(res$mean), decreasing = TRUE)[1:50]
  overlap <- length(intersect(est_top, active))
  assert_true(overlap >= 25,
              sprintf("Sparse active set overlap = %d/50, expected >= 25", overlap))
})

run_test("T3.3: higher Tau -> better recovery (compare cor at Tau=1 vs Tau=10)", {
  set.seed(22); n <- 80; p <- 60
  F_true <- rnorm(p, sd = 2)
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.01
  noise <- matrix(rnorm(n * p, sd = 1.0), n, p)
  R_k   <- EL_k %o% F_true + noise

  Tau_low  <- rep(1, p)
  Tau_high <- rep(10, p)

  res_low  <- update_F_k(Tau_low,  EL_k, EL2_k, R_k)
  res_high <- update_F_k(Tau_high, EL_k, EL2_k, R_k)

  cor_low  <- cor(res_low$mean,  F_true)
  cor_high <- cor(res_high$mean, F_true)
  assert_true(cor_high >= cor_low - 0.05,
              sprintf("Higher Tau should give better recovery: cor_low=%.3f, cor_high=%.3f",
                      cor_low, cor_high))
})

cat("\n=== T4: Multi-Factor K=5 via update_F_all ===\n")

run_test("T4.1: output dimensions p x K", {
  set.seed(30); n <- 40; p <- 60; K <- 5
  Y    <- matrix(rnorm(n * p), n, p)
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.05
  Tau  <- abs(rnorm(p)) + 0.1

  res <- update_F_all(Y, EL, EL2, EF, EF2, Tau)
  assert_true(nrow(res$EF) == p && ncol(res$EF) == K,
              sprintf("EF dimension wrong: got %d x %d, expected %d x %d",
                      nrow(res$EF), ncol(res$EF), p, K))
  assert_true(nrow(res$EF2) == p && ncol(res$EF2) == K,
              sprintf("EF2 dimension wrong: got %d x %d, expected %d x %d",
                      nrow(res$EF2), ncol(res$EF2), p, K))
  assert_length(res$details, K, "details length")
})

run_test("T4.2: all second moments >= squared means", {
  set.seed(31); n <- 40; p <- 60; K <- 5
  Y    <- matrix(rnorm(n * p), n, p)
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.05
  Tau  <- abs(rnorm(p)) + 0.1

  res <- update_F_all(Y, EL, EL2, EF, EF2, Tau)
  diff <- res$EF2 - res$EF^2
  assert_true(all(diff >= -1e-10),
              sprintf("EF2 < EF^2 somewhere: min diff = %.4e", min(diff)))
})

run_test("T4.3: all outputs finite", {
  set.seed(32); n <- 40; p <- 60; K <- 5
  Y    <- matrix(rnorm(n * p), n, p)
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.05
  Tau  <- abs(rnorm(p)) + 0.1

  res <- update_F_all(Y, EL, EL2, EF, EF2, Tau)
  assert_finite(as.vector(res$EF),  "EF contains non-finite")
  assert_finite(as.vector(res$EF2), "EF2 contains non-finite")
})

cat("\n=== T5: Null Factor Shrinkage ===\n")

run_test("T5.1: noise R_k -> factor column near zero", {
  set.seed(40); n <- 100; p <- 80
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.05
  R_k   <- matrix(rnorm(n * p, sd = 1.0), n, p)  # pure noise, no signal in F
  Tau   <- rep(1.0, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  # Point-normal prior should shrink toward zero when no signal
  mean_abs <- mean(abs(res$mean))
  assert_true(mean_abs < 1.0,
              sprintf("Null factor not shrunk: mean|EF| = %.3f", mean_abs))
})

run_test("T5.2: zero EL_k -> sum_EL2_k=0 -> A_F hits floor -> B_F=0 -> EF=0", {
  n <- 30; p <- 20
  EL_k  <- rep(0, n)
  EL2_k <- rep(0, n)               # zero loadings -> no information
  R_k   <- matrix(rnorm(n * p), n, p)
  Tau   <- rep(1.0, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  assert_near(res$sum_EL2_k, 0, tol = 1e-15, msg = "sum_EL2_k should be 0")
  assert_near(res$B, rep(0, p), tol = 1e-12, msg = "B_F should be 0 when EL_k=0")
  assert_near(res$mean, rep(0, p), tol = 1e-6, msg = "EF should be 0 when EL_k=0")
})

cat("\n=== T6: Tau Differential Shrinkage ===\n")

run_test("T6.1: high-tau features get less shrinkage than low-tau features", {
  # For the same true signal strength, features with higher Tau (lower noise)
  # should have posterior means closer to the truth (less shrinkage)
  set.seed(50); n <- 100; p <- 40
  F_true <- rep(2.0, p)           # uniform true factor
  EL_k  <- rnorm(n, sd = 1)
  EL2_k <- EL_k^2 + 0.01
  R_k   <- EL_k %o% F_true + matrix(rnorm(n * p, sd = 1.0), n, p)

  # Half features have high Tau, half have low Tau
  Tau <- rep(0.5, p)
  Tau[1:(p/2)] <- 50.0            # high precision for first half

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  # High-Tau features should be closer to x_j (less shrinkage)
  shrinkage_high <- mean(abs(res$mean[1:(p/2)] - res$x[1:(p/2)]))
  shrinkage_low  <- mean(abs(res$mean[(p/2+1):p] - res$x[(p/2+1):p]))
  # Note: x_j is the same regardless of Tau (T1.4 shows Tau cancels in x),
  # but s_j is smaller for high Tau, meaning less EBNM shrinkage
  # So high-Tau features have |mean - x| <= low-Tau features (on average)
  assert_true(shrinkage_high <= shrinkage_low + 0.2,
              sprintf("High-tau shrinkage=%.3f should be <= low-tau=%.3f",
                      shrinkage_high, shrinkage_low))
})

run_test("T6.2: larger EL2 (more loading uncertainty) -> more shrinkage", {
  # More uncertainty in loadings -> larger sum_EL2_k -> larger A_F -> smaller s_j
  # But also B_F stays the same (uses EL_k, not EL2_k), so x_j gets smaller
  # Net effect: more shrinkage because |x_j| is smaller relative to s_j
  set.seed(51); n <- 100; p <- 30
  EL_k   <- rnorm(n)
  R_k    <- matrix(rnorm(n * p, sd = 0.5), n, p)
  # Add signal: R_k ~ EL_k %o% F_true + noise
  F_true <- rnorm(p, sd = 2)
  R_k    <- EL_k %o% F_true + R_k
  Tau    <- rep(4.0, p)

  # Low loading variance
  EL2_low  <- EL_k^2 + 0.01
  res_low  <- update_F_k(Tau, EL_k, EL2_low, R_k)

  # High loading variance (same mean, much larger second moment)
  EL2_high <- EL_k^2 + 5.0
  res_high <- update_F_k(Tau, EL_k, EL2_high, R_k)

  # x_j should be smaller in absolute value for high EL2 (denominator is larger)
  assert_true(mean(abs(res_high$x)) < mean(abs(res_low$x)),
              "Higher EL2 should yield smaller |x_j|")
})

run_test("T6.3: uniform Tau -> all s_j equal (homoscedastic limit)", {
  set.seed(52); n <- 20; p <- 25
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)
  Tau   <- rep(3.0, p)            # uniform precision

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  # All s_j should be identical since Tau is uniform and sum_EL2_k is a scalar
  assert_near(res$s, rep(res$s[1], p), tol = 1e-15,
              msg = "s_j should all be equal when Tau is uniform")
})

cat("\n=== T7: Numerical Stability ===\n")

run_test("T7.1: all-zero Tau -> A_F hits floor, no NaN/Inf", {
  n <- 15; p <- 10
  Tau   <- rep(0, p)
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  assert_near(res$A, rep(1e-10, p), msg = "A should equal floor when Tau=0")
  assert_finite(c(res$mean, res$second, res$sd, res$x, res$s))
})

run_test("T7.2: extreme Tau=1e10 -> no overflow", {
  set.seed(60); n <- 15; p <- 10
  Tau   <- rep(1e10, p)
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.01
  R_k   <- matrix(rnorm(n * p, sd = 0.001), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B))
})

run_test("T7.3: p=1, n=1 -> works correctly with default alpha=0.5", {
  Tau   <- 2.0
  EL_k  <- 1.5
  EL2_k <- 1.5^2 + 0.1
  R_k   <- matrix(3.0, nrow = 1, ncol = 1)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k)
  expected_sum_EL2 <- sum(EL2_k)
  expected_A <- pmax(0.5 * Tau * expected_sum_EL2, 1e-10)
  expected_B <- 0.5 * Tau * as.vector(t(R_k) %*% EL_k)
  assert_near(res$A, expected_A, msg = "A mismatch for p=1,n=1")
  assert_near(res$B, expected_B, msg = "B mismatch for p=1,n=1")
  assert_length(res$mean, 1, "mean should have length 1")
  assert_finite(c(res$mean, res$second))
})

run_test("T7.4: custom A_floor is respected", {
  n <- 10; p <- 5
  Tau   <- rep(0, p)
  EL_k  <- rep(0, n)
  EL2_k <- rep(0, n)
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k, A_floor = 1.0)
  assert_near(res$A, rep(1.0, p), msg = "A should equal custom floor")
  assert_near(res$s, rep(1 / sqrt(1.0), p), msg = "s_j should use custom floor")
})

cat("\n=== T8: R_k and Gauss-Seidel ===\n")

run_test("T8.1: R_k uses correct EL_k (verify B_F computation) with default alpha=0.5", {
  # Verify that update_F_k computes B_F using the correct residual
  set.seed(70); n <- 20; p <- 15; K <- 3
  Y   <- matrix(rnorm(n * p), n, p)
  EL  <- matrix(rnorm(n * K), n, K)
  EF  <- matrix(rnorm(p * K), p, K)
  Tau <- abs(rnorm(p)) + 0.1
  k <- 2

  # Compute R_k manually
  R_k <- .compute_R_k_local(Y, EL, EF, k)
  EL2_k <- EL[, k]^2 + 0.05

  res <- update_F_k(Tau, EL[, k], EL2_k, R_k)

  # Verify B_F matches manual computation with this R_k (default alpha=0.5)
  expected_B <- 0.5 * Tau * as.vector(t(R_k) %*% EL[, k])
  assert_near(res$B, expected_B, tol = 1e-12,
              msg = "B_F should match manual R_k computation with alpha=0.5")
})

run_test("T8.2: update_F_all K columns all finite with valid second moments", {
  set.seed(71); n <- 30; p <- 40; K <- 4
  Y    <- matrix(rnorm(n * p), n, p)
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EF   <- matrix(rnorm(p * K), p, K)
  EF2  <- EF^2 + 0.05
  Tau  <- abs(rnorm(p)) + 0.1

  res <- update_F_all(Y, EL, EL2, EF, EF2, Tau)
  assert_finite(as.vector(res$EF),  "EF not finite")
  assert_finite(as.vector(res$EF2), "EF2 not finite")
  # All second moments >= squared means
  diff <- res$EF2 - res$EF^2
  assert_true(all(diff >= -1e-10),
              sprintf("Variance negative somewhere: min diff = %.4e", min(diff)))
})

cat("\n=== T9: V2.R Consistency ===\n")

run_test("T9.1: update_F_k with explicit alpha reproduces alpha-weighted reference exactly", {
  # V2.R used the unscaled formula (alpha=0 in the new parameterisation).
  # We verify with alpha=0.3: compute reference from scratch and confirm the
  # modular function matches it. Note: x_j = B/A is alpha-independent (cancels).
  set.seed(888); n <- 50; p <- 100
  Tau   <- abs(rnorm(p)) + 0.5
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + runif(n, 0.01, 0.1)  # second moments
  R_k   <- matrix(rnorm(n * p), n, p)
  test_alpha <- 0.3

  # --- alpha-weighted reference ---
  sum_EL2_k_ref <- sum(EL2_k)                                          # scalar
  A_F_ref <- pmax((1 - test_alpha) * Tau * sum_EL2_k_ref, 1e-10)       # p-vector
  B_F_ref <- (1 - test_alpha) * Tau * as.vector(t(R_k) %*% EL_k)       # p-vector
  res_ref <- ebnm(x = B_F_ref / A_F_ref, s = 1 / sqrt(A_F_ref),
                  prior_family = "point_normal")
  EF_ref  <- res_ref$posterior$mean
  EF2_ref <- res_ref$posterior$sd^2 + res_ref$posterior$mean^2

  # --- Modular function with same alpha ---
  res_mod <- update_F_k(Tau, EL_k, EL2_k, R_k, prior_family = "point_normal",
                        alpha = test_alpha)

  assert_near(res_mod$A,         A_F_ref,         tol = 1e-12, msg = "A_F mismatch vs reference")
  assert_near(res_mod$B,         B_F_ref,         tol = 1e-12, msg = "B_F mismatch vs reference")
  assert_near(res_mod$mean,      EF_ref,          tol = 1e-12, msg = "posterior mean mismatch vs reference")
  assert_near(res_mod$second,    EF2_ref,         tol = 1e-12, msg = "second moment mismatch vs reference")
  assert_near(res_mod$sum_EL2_k, sum_EL2_k_ref,   tol = 1e-15, msg = "sum_EL2_k mismatch vs reference")
})

cat("\n=== T_alpha: Alpha Mixing Parameter Edge Cases ===\n")

run_test("T_alpha.1: alpha=1 -> A_F and B_F both zero (floor activates, F zeroed)", {
  # alpha=1 means (1-alpha)=0, so A_F = pmax(0, floor) = floor and B_F = 0.
  # The EBNM gets x=B/A=0 and the posterior mean should be ~0 (F is killed).
  set.seed(801); n <- 20; p <- 15
  Tau   <- abs(rnorm(p)) + 0.5
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k, alpha = 1)
  assert_near(res$B, rep(0, p), tol = 1e-12, msg = "B_F should be 0 when alpha=1")
  assert_near(res$A, rep(1e-10, p), tol = 1e-12, msg = "A_F should equal floor when alpha=1")
  assert_near(res$mean, rep(0, p), tol = 1e-6, msg = "EF should be ~0 when alpha=1")
})

run_test("T_alpha.2: alpha=0 -> gives original unscaled formula (A = Tau*sum_EL2)", {
  # alpha=0 means (1-alpha)=1, so A_F = pmax(Tau*sum_EL2_k, floor) and
  # B_F = Tau*(t(R_k)%*%EL_k) — identical to V2.R.
  set.seed(802); n <- 20; p <- 15
  Tau   <- abs(rnorm(p)) + 0.5
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k, alpha = 0)
  expected_A <- pmax(Tau * sum(EL2_k), 1e-10)
  expected_B <- Tau * as.vector(t(R_k) %*% EL_k)
  assert_near(res$A, expected_A, tol = 1e-12, msg = "alpha=0: A_F should equal unscaled formula")
  assert_near(res$B, expected_B, tol = 1e-12, msg = "alpha=0: B_F should equal unscaled formula")
})

cat("\n=== T_surv: Dual-Source Survival Terms (Cluster B) ===\n")

run_test("T_surv.1: A_surv[j] = EBeta2_k * YtWY_diag[j]", {
  # Under Cluster B, the survival precision for feature j is E[beta_tilde_k^2]
  # times the pre-computed diagonal of Y'diag(w)Y.
  set.seed(901); n <- 30; p <- 20
  Tau       <- rep(1.0, p)
  EL_k      <- rnorm(n)
  EL2_k     <- EL_k^2 + 0.1
  R_k       <- matrix(rnorm(n * p), n, p)
  EBeta2_k  <- 0.8
  YtWY_diag <- abs(rnorm(p)) + 0.5    # simulated diag(Y'diag(w)Y)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k,
                    EBeta2_k  = EBeta2_k,
                    YtWY_diag = YtWY_diag)
  expected_A_surv <- EBeta2_k * YtWY_diag
  assert_near(res$A_surv, expected_A_surv, tol = 1e-12,
              msg = "A_surv[j] != EBeta2_k * YtWY_diag[j]")
})

run_test("T_surv.2: B_surv[j] = EBeta_k * YtWz_no_k[j]", {
  # The survival signal for feature j is E[beta_tilde_k] times the
  # pre-computed Y'(w * z_no_k).
  set.seed(902); n <- 30; p <- 20
  Tau       <- rep(1.0, p)
  EL_k      <- rnorm(n)
  EL2_k     <- EL_k^2 + 0.1
  R_k       <- matrix(rnorm(n * p), n, p)
  EBeta_k   <- 1.3
  YtWz_no_k <- rnorm(p)

  res <- update_F_k(Tau, EL_k, EL2_k, R_k,
                    EBeta_k   = EBeta_k,
                    YtWz_no_k = YtWz_no_k)
  expected_B_surv <- EBeta_k * YtWz_no_k
  assert_near(res$B_surv, expected_B_surv, tol = 1e-12,
              msg = "B_surv[j] != EBeta_k * YtWz_no_k[j]")
})

run_test("T_surv.3: alpha=0 with non-zero survival args -> A_F = A_gen (survival suppressed)", {
  # When alpha=0, the (1-alpha) weight on genomics is 1 and alpha weight on
  # survival is 0, so survival terms drop out regardless of EBeta2_k / YtWY_diag.
  set.seed(903); n <- 25; p <- 15
  Tau       <- abs(rnorm(p)) + 0.5
  EL_k      <- rnorm(n)
  EL2_k     <- EL_k^2 + 0.1
  R_k       <- matrix(rnorm(n * p), n, p)
  EBeta2_k  <- 5.0
  YtWY_diag <- abs(rnorm(p)) + 1.0   # non-zero, would inflate A if alpha > 0

  res_surv <- update_F_k(Tau, EL_k, EL2_k, R_k,
                          EBeta2_k  = EBeta2_k,
                          YtWY_diag = YtWY_diag,
                          alpha = 0)
  res_base <- update_F_k(Tau, EL_k, EL2_k, R_k, alpha = 0)

  assert_near(res_surv$A, res_base$A, tol = 1e-12,
              msg = "alpha=0: A_F should equal A_gen regardless of survival args")
  assert_near(res_surv$B, res_base$B, tol = 1e-12,
              msg = "alpha=0: B_F should equal B_gen regardless of survival args")
})

run_test("T_surv.4: alpha>0 with EBeta2_k>0 -> A_F > (1-alpha)*A_gen (survival inflates A)", {
  # With alpha > 0 and a non-zero EBeta2_k, the survival term adds positive
  # mass to A_F, so A_F should exceed the genomics-only contribution.
  set.seed(904); n <- 30; p <- 20
  Tau       <- rep(2.0, p)
  EL_k      <- rnorm(n)
  EL2_k     <- EL_k^2 + 0.1
  R_k       <- matrix(rnorm(n * p), n, p)
  EBeta2_k  <- 2.0
  YtWY_diag <- rep(1.0, p)   # each feature contributes equally

  res <- update_F_k(Tau, EL_k, EL2_k, R_k,
                    EBeta2_k  = EBeta2_k,
                    YtWY_diag = YtWY_diag,
                    alpha = 0.5)

  # Genomics-only A (alpha = 0 baseline)
  res_gen <- update_F_k(Tau, EL_k, EL2_k, R_k, alpha = 0)
  A_gen_contrib <- 0.5 * res_gen$A   # (1-alpha) * A_gen

  # A_F should be strictly greater because A_surv > 0
  assert_true(all(res$A > A_gen_contrib - 1e-10),
              "alpha=0.5, EBeta2>0: A_F should exceed (1-alpha)*A_gen")
  assert_true(mean(res$A) > mean(A_gen_contrib),
              "mean A_F should exceed mean (1-alpha)*A_gen with survival")
})

run_test("T_surv.5: YtWY_diag=0 -> A_surv=0 -> A_F = (1-alpha)*A_gen", {
  # When there is no survival weight on any feature (YtWY_diag=0), the
  # survival precision term is zero and A_F reduces to (1-alpha)*A_gen.
  set.seed(905); n <- 20; p <- 12
  Tau       <- rep(1.5, p)
  EL_k      <- rnorm(n)
  EL2_k     <- EL_k^2 + 0.1
  R_k       <- matrix(rnorm(n * p), n, p)
  EBeta2_k  <- 3.0
  alpha     <- 0.5

  res <- update_F_k(Tau, EL_k, EL2_k, R_k,
                    EBeta2_k  = EBeta2_k,
                    YtWY_diag = rep(0, p),   # zero survival weights
                    alpha = alpha)

  expected_A <- pmax((1 - alpha) * Tau * sum(EL2_k), 1e-10)
  assert_near(res$A, expected_A, tol = 1e-12,
              msg = "YtWY_diag=0: A_F should equal (1-alpha)*A_gen")
  assert_near(res$A_surv, rep(0, p), tol = 1e-15,
              msg = "A_surv should be zero when YtWY_diag=0")
})

run_test("T_surv.6: EBeta_k=0 -> B_surv=0 -> B_F = (1-alpha)*B_gen", {
  # When the posterior mean of beta_tilde_k is zero, the survival signal
  # B_surv vanishes and B_F is driven entirely by genomics.
  set.seed(906); n <- 20; p <- 12
  Tau       <- abs(rnorm(p)) + 0.5
  EL_k      <- rnorm(n)
  EL2_k     <- EL_k^2 + 0.1
  R_k       <- matrix(rnorm(n * p), n, p)
  alpha     <- 0.5

  res <- update_F_k(Tau, EL_k, EL2_k, R_k,
                    EBeta_k   = 0,
                    YtWz_no_k = rnorm(p),   # non-zero, but zeroed by EBeta_k=0
                    alpha = alpha)

  expected_B <- (1 - alpha) * Tau * as.vector(t(R_k) %*% EL_k)
  assert_near(res$B, expected_B, tol = 1e-12,
              msg = "EBeta_k=0: B_F should equal (1-alpha)*B_gen")
  assert_near(res$B_surv, rep(0, p), tol = 1e-15,
              msg = "B_surv should be zero when EBeta_k=0")
})
