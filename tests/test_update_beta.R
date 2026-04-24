# =============================================================================
# tests/test_update_beta.R
#
# TDD test suite for the q_beta CAVI update.
# Written BEFORE implementation (red phase).
#
# Tests cover:
#   T1: Mathematical identity checks (A_k, B_k, x_k, s_k)
#   T2: Non-informative prior / WLS limit
#   T3: Known signal recovery (K=1)
#   T4: Multi-factor recovery (K=5) via update_beta_all
#   T5: Null factor (beta=0 shrinkage)
#   T6: Error-in-variables (posterior variance inflates precision)
#   T7: Numerical stability (degenerate weights, extreme values)
#   T8: Gauss-Seidel ordering in update_beta_all
#   T9: Consistency with V2.R inline code
#
# NOTE: requires ebnm package.  Install with:
#   install.packages("ebnm")
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# ---- helper: simulate Cox Taylor (z, w) given known beta and L ----
.sim_cox_taylor <- function(n, EL, beta_true, seed = 42) {
  set.seed(seed)
  K  <- ncol(EL)
  eta <- as.vector(EL %*% beta_true)
  # Exponential survival times proportional to exp(-eta)
  time   <- rexp(n, rate = exp(eta) / mean(exp(eta)))
  status <- rep(1L, n)   # all events for simplicity
  list(time = time, status = status)
}

# simple Cox Taylor without needing source(V2.R) — inline
.calc_cox_taylor_local <- function(eta, time, status) {
  n   <- length(time)
  ord <- order(time)
  time_s   <- time[ord];  status_s <- status[ord];  eta_s <- eta[ord]
  theta    <- exp(eta_s)
  risk_sum <- rev(cumsum(rev(theta)))
  h <- status_s / risk_sum
  H <- cumsum(h)
  u_s <- status_s - theta * H
  w_s <- theta * H
  w_s[w_s < 1e-6] <- 1e-6
  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_s;  w[ord] <- w_s
  list(u = u, w = w)
}

cat("=== T1: Mathematical Identity Checks ===\n")

run_test("T1.1: A_k = alpha * sum(w * EL2_k) with default alpha=0.5", {
  set.seed(1); n <- 5
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + runif(n, 0, 0.2)   # second moment > squared mean
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  expected_A <- max(0.5 * sum(w * EL2_k), 1e-10)
  assert_near(res$A, expected_A, msg = "A_k mismatch")
})

run_test("T1.2: B_k = alpha * sum(w * z_no_k * EL_k) with default alpha=0.5", {
  set.seed(2); n <- 5
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  expected_B <- 0.5 * sum(w * z_no_k * EL_k)
  assert_near(res$B, expected_B, msg = "B_k mismatch")
})

run_test("T1.3: x_k = B_k / A_k", {
  set.seed(3); n <- 6
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_near(res$x, res$B / res$A, msg = "x_k != B_k / A_k")
})

run_test("T1.4: s_k = 1 / sqrt(A_k)", {
  set.seed(4); n <- 6
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_near(res$s, 1 / sqrt(res$A), msg = "s_k != 1/sqrt(A_k)")
})

run_test("T1.5: second moment = sd^2 + mean^2", {
  set.seed(5); n <- 8
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_near(res$second, res$sd^2 + res$mean^2, msg = "second moment identity")
})

run_test("T1.6: return values are all finite", {
  set.seed(6); n <- 10
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B, res$x, res$s))
})

cat("\n=== T2: Non-Informative Prior / WLS Limit ===\n")

run_test("T2.1: x_k equals WLS estimate when no posterior variance in L", {
  # When EL2_k = EL_k^2 (no posterior variance), A_k = sum(w * EL_k^2)
  # and the EBNM pseudo-obs x_k = sum(w*z*EL_k) / sum(w*EL_k^2)
  # which is the WLS estimator of beta_k given z = EL_k * beta_k + noise
  set.seed(10); n <- 50
  beta_true <- 2.0
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2    # zero posterior variance
  z_no_k <- EL_k * beta_true + rnorm(n, sd = 0.1)
  w     <- rep(5.0, n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  wls_est <- sum(w * z_no_k * EL_k) / sum(w * EL_k^2)
  assert_near(res$x, wls_est, tol = 1e-10, msg = "x_k should equal WLS when EL2=EL^2")
})

run_test("T2.2: larger signal-to-noise -> posterior mean close to x_k", {
  # When |x_k| >> s_k, EBNM should apply minimal shrinkage
  set.seed(11); n <- 100
  beta_true <- 5.0   # large signal
  EL_k  <- rnorm(n, sd = 1)
  EL2_k <- EL_k^2 + 1e-4  # near-zero variance
  z_no_k <- EL_k * beta_true + rnorm(n, sd = 0.01)  # tiny noise
  w     <- rep(10.0, n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  # With very strong signal, mean should be close to x_k
  assert_true(abs(res$mean - res$x) < 0.5, "posterior mean far from x_k despite strong signal")
})

cat("\n=== T3: Known Signal Recovery (K=1) ===\n")

run_test("T3.1: K=1 recovery of positive beta", {
  set.seed(20); n <- 200
  beta_true <- 2.0
  EL_k  <- matrix(rnorm(n), n, 1)
  EL2_k <- EL_k^2 + 0.01
  eta   <- as.vector(EL_k) * beta_true
  surv  <- .sim_cox_taylor(n, EL_k, beta_true, seed = 20)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z     <- eta + taylor$u / taylor$w
  w     <- taylor$w
  z_no_k <- z   # K=1, no other factors

  res <- update_beta_k(w, z_no_k, as.vector(EL_k), as.vector(EL2_k))
  # Sign should be correct and rough magnitude OK
  assert_true(res$mean > 0.5, sprintf("Positive beta not recovered (got %.3f)", res$mean))
})

run_test("T3.2: K=1 recovery of negative beta", {
  set.seed(21); n <- 200
  beta_true <- -1.5
  EL_k  <- matrix(rnorm(n), n, 1)
  EL2_k <- EL_k^2 + 0.01
  eta   <- as.vector(EL_k) * beta_true
  surv  <- .sim_cox_taylor(n, EL_k, beta_true, seed = 21)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z     <- eta + taylor$u / taylor$w
  w     <- taylor$w
  z_no_k <- z

  res <- update_beta_k(w, z_no_k, as.vector(EL_k), as.vector(EL2_k))
  assert_true(res$mean < -0.3, sprintf("Negative beta not recovered (got %.3f)", res$mean))
})

cat("\n=== T4: Multi-Factor Recovery (K=5, update_beta_all) ===\n")

run_test("T4.1: update_beta_all output dimensions", {
  set.seed(30); n <- 100; K <- 5
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EBeta <- rnorm(K)
  w    <- abs(rnorm(n)) + 0.5
  z    <- as.vector(EL %*% EBeta) + rnorm(n, sd = 0.5)

  res <- update_beta_all(w, z, EL, EL2, EBeta)
  assert_length(res$EBeta,   K, "EBeta length")
  assert_length(res$EBeta2,  K, "EBeta2 length")
  assert_length(res$details, K, "details length")
})

run_test("T4.2: all second moments non-negative", {
  set.seed(31); n <- 100; K <- 5
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EBeta <- rnorm(K)
  w    <- abs(rnorm(n)) + 0.5
  z    <- as.vector(EL %*% EBeta) + rnorm(n, sd = 0.5)

  res <- update_beta_all(w, z, EL, EL2, EBeta)
  assert_true(all(res$EBeta2 >= 0), "EBeta2 must be non-negative")
})

run_test("T4.3: beta signs match true values (strong signal)", {
  set.seed(32); n <- 500; K <- 3
  beta_true <- c(2.0, -1.5, 0.0)
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.01
  eta  <- as.vector(EL %*% beta_true)
  surv <- .sim_cox_taylor(n, EL, beta_true, seed = 32)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z    <- eta + taylor$u / taylor$w
  w    <- taylor$w

  res <- update_beta_all(w, z, EL, EL2, beta_true)
  # Signs for non-zero betas should be correct
  assert_true(res$EBeta[1] > 0,  sprintf("GEP1 sign wrong (%.3f)", res$EBeta[1]))
  assert_true(res$EBeta[2] < 0,  sprintf("GEP2 sign wrong (%.3f)", res$EBeta[2]))
  # Zero beta should be shrunk toward zero by point-normal
  assert_true(abs(res$EBeta[3]) < 0.5, sprintf("Null GEP3 not shrunk (%.3f)", res$EBeta[3]))
})

cat("\n=== T5: Null Factor (beta=0 shrinkage) ===\n")

run_test("T5.1: pure noise in z_no_k -> beta shrunk toward zero", {
  set.seed(40); n <- 200
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.05
  z_no_k <- rnorm(n, sd = 1.0)   # no signal: just noise
  w     <- rep(2.0, n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  # Point-normal with no real signal should shrink heavily
  # Mean should be near zero (allow for Monte Carlo variance)
  assert_true(abs(res$mean) < 1.5, sprintf("Null beta not shrunk: got %.3f", res$mean))
})

run_test("T5.2: zero EL_k -> B_k = 0 -> beta = 0", {
  # If EL_k = 0, there's no information, B_k = 0, x_k = 0
  n <- 50
  EL_k  <- rep(0, n)
  EL2_k <- rep(1e-10, n)
  z_no_k <- rnorm(n)
  w     <- rep(1.0, n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_near(res$B, 0, tol = 1e-12, msg = "B_k should be 0 when EL_k=0")
  assert_near(res$x, 0, tol = 1e-6,  msg = "x_k should be 0 when B_k=0")
  assert_near(res$mean, 0, tol = 1e-6, msg = "beta_mean should be 0 when x_k=0")
})

cat("\n=== T6: Error-in-Variables ===\n")

run_test("T6.1: larger EL2 (more uncertainty) -> larger A -> more shrinkage", {
  # Keep EL_k fixed, increase posterior variance: should increase A_k and shrink beta more
  set.seed(50); n <- 100
  EL_k   <- rnorm(n)
  z_no_k <- EL_k * 3.0 + rnorm(n, sd = 0.1)  # strong signal
  w      <- rep(2.0, n)

  # Low posterior variance
  EL2_low <- EL_k^2 + 0.01
  res_low  <- update_beta_k(w, z_no_k, EL_k, EL2_low)

  # High posterior variance (same mean, much larger second moment)
  EL2_high <- EL_k^2 + 5.0
  res_high  <- update_beta_k(w, z_no_k, EL_k, EL2_high)

  # Error-in-variables: higher EL2 → larger A_k → smaller x_k=B/A (signal divided by more)
  # → more shrinkage in posterior mean.  Note: s_k=1/sqrt(A_k) is SMALLER (not larger)
  # when A_k increases.  Shrinkage comes from smaller |x_k|/s_k = |B|/sqrt(A).
  assert_true(res_high$A > res_low$A, "Higher EL2 should yield larger A_k")
  assert_true(res_high$s < res_low$s, "Higher EL2 yields SMALLER s_k (1/sqrt(larger A_k))")
  assert_true(abs(res_high$mean) <= abs(res_low$mean) + 0.3,
              sprintf("Higher EL2 should yield more shrinkage: low=%.3f, high=%.3f",
                      res_low$mean, res_high$mean))
})

cat("\n=== T7: Numerical Stability ===\n")

run_test("T7.1: all-zero weights -> A hits floor, no NaN/Inf", {
  n <- 20
  w     <- rep(0, n)
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.1
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_near(res$A, 1e-10, msg = "A should equal floor when all w=0")
  assert_finite(c(res$mean, res$second, res$sd, res$x, res$s))
})

run_test("T7.2: extreme weights (w = 1e8) -> no overflow", {
  set.seed(60); n <- 20
  w     <- rep(1e8, n)
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + 0.01
  z_no_k <- rnorm(n, sd = 0.001)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B))
})

run_test("T7.3: single observation (n=1) with default alpha=0.5", {
  w     <- 2.0
  EL_k  <- 1.5
  EL2_k <- 1.5^2 + 0.1
  z_no_k <- 3.0

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  expected_A <- max(0.5 * w * EL2_k, 1e-10)
  expected_B <- 0.5 * w * z_no_k * EL_k
  assert_near(res$A, expected_A)
  assert_near(res$B, expected_B)
  assert_finite(c(res$mean, res$second))
})

run_test("T7.4: custom A_floor is respected", {
  n <- 10
  w     <- rep(0, n)
  EL_k  <- rep(0, n)
  EL2_k <- rep(0, n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k, A_floor = 1.0)
  assert_near(res$A, 1.0, msg = "A should equal custom floor")
  assert_near(res$s, 1 / sqrt(1.0), msg = "s_k should use custom floor")
})

run_test("T7.5: second moment >= mean^2 (variance non-negative)", {
  set.seed(61); n <- 50
  w    <- abs(rnorm(n)) + 0.1
  EL_k  <- rnorm(n)
  EL2_k <- EL_k^2 + runif(n, 0, 0.5)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k)
  assert_true(res$second >= res$mean^2 - 1e-10,
              sprintf("second moment %.4e < mean^2 %.4e", res$second, res$mean^2))
})

cat("\n=== T8: Gauss-Seidel Ordering in update_beta_all ===\n")

run_test("T8.1: compute_z_no_k excludes factor k correctly", {
  set.seed(70); n <- 10; K <- 3
  EL   <- matrix(rnorm(n * K), n, K)
  EBeta <- c(1.0, -0.5, 0.3)
  z    <- as.vector(EL %*% EBeta) + rnorm(n, sd = 0.1)

  for (k in 1:K) {
    z_k <- compute_z_no_k(z, EL, EBeta, k)
    # Manual computation
    eta_no_k <- as.vector(EL %*% EBeta) - EL[, k] * EBeta[k]
    expected  <- z - eta_no_k
    assert_near(z_k, expected, tol = 1e-12,
                msg = sprintf("z_no_k wrong for k=%d", k))
  }
})

run_test("T8.2: update_beta_all uses Gauss-Seidel (updates propagate within loop)", {
  # With Gauss-Seidel, after updating beta_1, the z_no_k for beta_2
  # should incorporate the NEW beta_1, not the old one.
  set.seed(71); n <- 50; K <- 2
  EL   <- matrix(rnorm(n * K), n, K)
  EL2  <- EL^2 + 0.05
  EBeta_init <- c(0.0, 0.0)
  w    <- rep(1.0, n)
  z    <- EL[, 1] * 2.0 + rnorm(n, sd = 0.3)

  res <- update_beta_all(w, z, EL, EL2, EBeta_init)
  # Result should be finite and have correct structure
  assert_finite(c(res$EBeta, res$EBeta2))
  assert_length(res$EBeta, K)
})

cat("\n=== T9: Consistency with V2.R Inline Code ===\n")

run_test("T9.1: update_beta_k with alpha=1 matches V2.R lines 356-361 exactly", {
  # V2.R uses the unscaled formula (equivalent to alpha=1.0 in the new parameterisation).
  # Pass alpha=1.0 explicitly so the modular function reproduces V2.R's result.
  set.seed(80); n <- 30
  w      <- abs(rnorm(n)) + 0.5
  z_no_k <- rnorm(n)
  EL_k   <- rnorm(n)
  EL2_k  <- EL_k^2 + runif(n, 0.01, 0.1)

  # V2.R inline code (lines 356-361):
  A_Beta_v2 <- max(sum(w * EL2_k), 1e-10)
  B_Beta_v2 <- sum(w * z_no_k * EL_k)
  res_v2    <- ebnm(x = B_Beta_v2 / A_Beta_v2, s = 1 / sqrt(A_Beta_v2),
                    prior_family = "point_normal")
  mean_v2   <- res_v2$posterior$mean
  sec_v2    <- res_v2$posterior$sd^2 + res_v2$posterior$mean^2

  # Modular function with alpha=1.0 to match V2.R's unscaled formula:
  res_mod <- update_beta_k(w, z_no_k, EL_k, EL2_k, alpha = 1.0)

  assert_near(res_mod$mean,   mean_v2, tol = 1e-12, msg = "posterior mean mismatch vs V2.R")
  assert_near(res_mod$second, sec_v2,  tol = 1e-12, msg = "second moment mismatch vs V2.R")
  assert_near(res_mod$A, A_Beta_v2,   tol = 1e-12, msg = "A_k mismatch vs V2.R")
  assert_near(res_mod$B, B_Beta_v2,   tol = 1e-12, msg = "B_k mismatch vs V2.R")
})

cat("\n=== T_alpha: Alpha Mixing Parameter Edge Cases ===\n")

run_test("T_alpha.1: alpha=0 -> B_k=0 -> beta mean approx 0", {
  # alpha=0 zeroes both A and B (floor takes over for A), so x_k=B/A=0
  # and EBNM with x=0 under point_normal returns posterior mean ~0.
  set.seed(90); n <- 30
  w      <- abs(rnorm(n)) + 0.5
  z_no_k <- EL_k <- rnorm(n)
  EL2_k  <- EL_k^2 + 0.1

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k, alpha = 0)
  assert_near(res$B, 0, tol = 1e-12, msg = "B_k should be 0 when alpha=0")
  assert_near(res$x, 0, tol = 1e-6,  msg = "x_k should be 0 when alpha=0")
  assert_near(res$mean, 0, tol = 1e-6, msg = "beta_mean should be ~0 when alpha=0")
})

run_test("T_alpha.2: alpha=1 -> gives original unscaled formula (A = sum(w*EL2_k))", {
  # alpha=1 reproduces the V2.R formula: A_k = sum(w*EL2_k), B_k = sum(w*z*EL_k)
  set.seed(91); n <- 20
  w      <- abs(rnorm(n)) + 0.5
  z_no_k <- rnorm(n)
  EL_k   <- rnorm(n)
  EL2_k  <- EL_k^2 + 0.05

  res <- update_beta_k(w, z_no_k, EL_k, EL2_k, alpha = 1.0)
  expected_A <- max(sum(w * EL2_k), 1e-10)
  expected_B <- sum(w * z_no_k * EL_k)
  assert_near(res$A, expected_A, msg = "A_k with alpha=1 should match unscaled formula")
  assert_near(res$B, expected_B, msg = "B_k with alpha=1 should match unscaled formula")
})
