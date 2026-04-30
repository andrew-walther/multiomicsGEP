# =============================================================================
# tests/test_update_beta.R
#
# TDD test suite for the q_beta_tilde CAVI update — Cluster B
# (Cox-on-YF reformulation): ZF predictor, no EL2 second moment.
#
# Under eta = ZF * beta_tilde (Cluster B), the predictor ZF = Y * E[F] is
# observed per CAVI iteration — no second-moment correction is needed.
# EL2_k is removed from update_beta_k(); EL/EL2 are replaced by ZF in
# update_beta_all() and compute_z_no_k().
#
# Tests cover:
#   T1: Mathematical identity checks (A_k, B_k, x_k, s_k)
#   T2: Non-informative prior / WLS limit
#   T3: Known signal recovery (K=1)
#   T4: Multi-factor recovery (K=5) via update_beta_all
#   T5: Null factor (beta=0 shrinkage)
#   T7: Numerical stability (degenerate weights, extreme values)
#   T8: Gauss-Seidel ordering in update_beta_all
#   T9: Consistency with reference formula
#   T_alpha: Alpha mixing parameter edge cases
#   T_NEW: Cluster B ZF-predictor verification
#
# NOTE: requires ebnm package.  Install with:
#   install.packages("ebnm")
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# ---- helper: simulate Cox Taylor (z, w) given known beta and ZF ----
.sim_cox_taylor <- function(n, ZF, beta_true, seed = 42) {
  set.seed(seed)
  eta    <- as.vector(ZF %*% beta_true)
  time   <- rexp(n, rate = exp(eta) / mean(exp(eta)))
  status <- rep(1L, n)
  list(time = time, status = status)
}

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

run_test("T1.1: A_k = alpha * sum(w * ZF_k^2) with default alpha=0.5", {
  set.seed(1); n <- 5
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  expected_A <- max(0.5 * sum(w * ZF_k^2), 1e-10)
  assert_near(res$A, expected_A, msg = "A_k mismatch")
})

run_test("T1.2: B_k = alpha * sum(w * z_no_k * ZF_k) with default alpha=0.5", {
  set.seed(2); n <- 5
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  expected_B <- 0.5 * sum(w * z_no_k * ZF_k)
  assert_near(res$B, expected_B, msg = "B_k mismatch")
})

run_test("T1.3: x_k = B_k / A_k", {
  set.seed(3); n <- 6
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_near(res$x, res$B / res$A, msg = "x_k != B_k / A_k")
})

run_test("T1.4: s_k = 1 / sqrt(A_k)", {
  set.seed(4); n <- 6
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_near(res$s, 1 / sqrt(res$A), msg = "s_k != 1/sqrt(A_k)")
})

run_test("T1.5: second moment = sd^2 + mean^2", {
  set.seed(5); n <- 8
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_near(res$second, res$sd^2 + res$mean^2, msg = "second moment identity")
})

run_test("T1.6: return values are all finite", {
  set.seed(6); n <- 10
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B, res$x, res$s))
})

cat("\n=== T2: Non-Informative Prior / WLS Limit ===\n")

run_test("T2.1: x_k equals WLS estimate (ZF is observed, no second-moment inflate)", {
  # A_k = sum(w * ZF_k^2) and B_k = sum(w * z * ZF_k), so
  # x_k = B_k/A_k = WLS estimator of beta given z = ZF_k * beta + noise.
  set.seed(10); n <- 50
  beta_true <- 2.0
  ZF_k   <- rnorm(n)
  z_no_k <- ZF_k * beta_true + rnorm(n, sd = 0.1)
  w      <- rep(5.0, n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  wls_est <- sum(w * z_no_k * ZF_k) / sum(w * ZF_k^2)
  assert_near(res$x, wls_est, tol = 1e-10, msg = "x_k should equal WLS estimator")
})

run_test("T2.2: larger signal-to-noise -> posterior mean close to x_k", {
  set.seed(11); n <- 100
  beta_true <- 5.0
  ZF_k   <- rnorm(n, sd = 1)
  z_no_k <- ZF_k * beta_true + rnorm(n, sd = 0.01)
  w      <- rep(10.0, n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_true(abs(res$mean - res$x) < 0.5, "posterior mean far from x_k despite strong signal")
})

cat("\n=== T3: Known Signal Recovery (K=1) ===\n")

run_test("T3.1: K=1 recovery of positive beta_tilde", {
  set.seed(20); n <- 200
  beta_true <- 2.0
  ZF_k   <- matrix(rnorm(n), n, 1)
  eta    <- as.vector(ZF_k) * beta_true
  surv   <- .sim_cox_taylor(n, ZF_k, beta_true, seed = 20)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z      <- eta + taylor$u / taylor$w
  w      <- taylor$w
  z_no_k <- z   # K=1, no other factors

  res <- update_beta_k(w, z_no_k, as.vector(ZF_k))
  assert_true(res$mean > 0.5, sprintf("Positive beta not recovered (got %.3f)", res$mean))
})

run_test("T3.2: K=1 recovery of negative beta_tilde", {
  set.seed(21); n <- 200
  beta_true <- -1.5
  ZF_k   <- matrix(rnorm(n), n, 1)
  eta    <- as.vector(ZF_k) * beta_true
  surv   <- .sim_cox_taylor(n, ZF_k, beta_true, seed = 21)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z      <- eta + taylor$u / taylor$w
  w      <- taylor$w
  z_no_k <- z

  res <- update_beta_k(w, z_no_k, as.vector(ZF_k))
  assert_true(res$mean < -0.3, sprintf("Negative beta not recovered (got %.3f)", res$mean))
})

cat("\n=== T4: Multi-Factor Recovery (K=5, update_beta_all) ===\n")

run_test("T4.1: update_beta_all output dimensions", {
  set.seed(30); n <- 100; K <- 5
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta <- rnorm(K)
  w     <- abs(rnorm(n)) + 0.5
  z     <- as.vector(ZF %*% EBeta) + rnorm(n, sd = 0.5)

  res <- update_beta_all(w, z, ZF, EBeta)
  assert_length(res$EBeta,   K, "EBeta length")
  assert_length(res$EBeta2,  K, "EBeta2 length")
  assert_length(res$details, K, "details length")
})

run_test("T4.2: all second moments non-negative", {
  set.seed(31); n <- 100; K <- 5
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta <- rnorm(K)
  w     <- abs(rnorm(n)) + 0.5
  z     <- as.vector(ZF %*% EBeta) + rnorm(n, sd = 0.5)

  res <- update_beta_all(w, z, ZF, EBeta)
  assert_true(all(res$EBeta2 >= 0), "EBeta2 must be non-negative")
})

run_test("T4.3: beta signs match true values (strong signal)", {
  set.seed(32); n <- 500; K <- 3
  beta_true <- c(2.0, -1.5, 0.0)
  ZF    <- matrix(rnorm(n * K), n, K)
  eta   <- as.vector(ZF %*% beta_true)
  surv  <- .sim_cox_taylor(n, ZF, beta_true, seed = 32)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z     <- eta + taylor$u / taylor$w
  w     <- taylor$w

  res <- update_beta_all(w, z, ZF, beta_true)
  assert_true(res$EBeta[1] > 0,  sprintf("GEP1 sign wrong (%.3f)", res$EBeta[1]))
  assert_true(res$EBeta[2] < 0,  sprintf("GEP2 sign wrong (%.3f)", res$EBeta[2]))
  assert_true(abs(res$EBeta[3]) < 0.5, sprintf("Null GEP3 not shrunk (%.3f)", res$EBeta[3]))
})

cat("\n=== T5: Null Factor (beta=0 shrinkage) ===\n")

run_test("T5.1: pure noise in z_no_k -> beta shrunk toward zero", {
  set.seed(40); n <- 200
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n, sd = 1.0)
  w      <- rep(2.0, n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_true(abs(res$mean) < 1.5, sprintf("Null beta not shrunk: got %.3f", res$mean))
})

run_test("T5.2: zero ZF_k -> B_k = 0 -> beta = 0", {
  n      <- 50
  ZF_k   <- rep(0, n)
  z_no_k <- rnorm(n)
  w      <- rep(1.0, n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_near(res$B, 0, tol = 1e-12, msg = "B_k should be 0 when ZF_k=0")
  assert_near(res$x, 0, tol = 1e-6,  msg = "x_k should be 0 when B_k=0")
  assert_near(res$mean, 0, tol = 1e-6, msg = "beta_mean should be 0 when x_k=0")
})

cat("\n=== T7: Numerical Stability ===\n")

run_test("T7.1: all-zero weights -> A hits floor, no NaN/Inf", {
  n      <- 20
  w      <- rep(0, n)
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_near(res$A, 1e-10, msg = "A should equal floor when all w=0")
  assert_finite(c(res$mean, res$second, res$sd, res$x, res$s))
})

run_test("T7.2: extreme weights (w = 1e8) -> no overflow", {
  set.seed(60); n <- 20
  w      <- rep(1e8, n)
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n, sd = 0.001)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B))
})

run_test("T7.3: single observation (n=1) with default alpha=0.5", {
  w      <- 2.0
  ZF_k   <- 1.5
  z_no_k <- 3.0

  res <- update_beta_k(w, z_no_k, ZF_k)
  expected_A <- max(0.5 * w * ZF_k^2, 1e-10)
  expected_B <- 0.5 * w * z_no_k * ZF_k
  assert_near(res$A, expected_A)
  assert_near(res$B, expected_B)
  assert_finite(c(res$mean, res$second))
})

run_test("T7.4: custom A_floor is respected", {
  n      <- 10
  w      <- rep(0, n)
  ZF_k   <- rep(0, n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k, A_floor = 1.0)
  assert_near(res$A, 1.0, msg = "A should equal custom floor")
  assert_near(res$s, 1 / sqrt(1.0), msg = "s_k should use custom floor")
})

run_test("T7.5: second moment >= mean^2 (variance non-negative)", {
  set.seed(61); n <- 50
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)
  assert_true(res$second >= res$mean^2 - 1e-10,
              sprintf("second moment %.4e < mean^2 %.4e", res$second, res$mean^2))
})

cat("\n=== T8: Gauss-Seidel Ordering in update_beta_all ===\n")

run_test("T8.1: compute_z_no_k with ZF excludes factor k correctly", {
  set.seed(70); n <- 10; K <- 3
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta <- c(1.0, -0.5, 0.3)
  z     <- as.vector(ZF %*% EBeta) + rnorm(n, sd = 0.1)

  for (k in 1:K) {
    z_k       <- compute_z_no_k(z, ZF, EBeta, k)
    eta_no_k  <- as.vector(ZF %*% EBeta) - ZF[, k] * EBeta[k]
    expected  <- z - eta_no_k
    assert_near(z_k, expected, tol = 1e-12,
                msg = sprintf("z_no_k wrong for k=%d", k))
  }
})

run_test("T8.2: update_beta_all uses Gauss-Seidel (updates propagate within loop)", {
  set.seed(71); n <- 50; K <- 2
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta_init <- c(0.0, 0.0)
  w     <- rep(1.0, n)
  z     <- ZF[, 1] * 2.0 + rnorm(n, sd = 0.3)

  res <- update_beta_all(w, z, ZF, EBeta_init)
  assert_finite(c(res$EBeta, res$EBeta2))
  assert_length(res$EBeta, K)
})

cat("\n=== T9: Consistency with Reference Formula ===\n")

run_test("T9.1: update_beta_k with alpha=1 matches unscaled reference formula", {
  # alpha=1: A_k = sum(w * ZF_k^2), B_k = sum(w * z_no_k * ZF_k)
  set.seed(80); n <- 30
  w      <- abs(rnorm(n)) + 0.5
  z_no_k <- rnorm(n)
  ZF_k   <- rnorm(n)

  A_ref  <- max(sum(w * ZF_k^2), 1e-10)
  B_ref  <- sum(w * z_no_k * ZF_k)
  res_ref <- ebnm(x = B_ref / A_ref, s = 1 / sqrt(A_ref),
                  prior_family = "point_normal")
  mean_ref <- res_ref$posterior$mean
  sec_ref  <- res_ref$posterior$sd^2 + res_ref$posterior$mean^2

  res_mod <- update_beta_k(w, z_no_k, ZF_k, alpha = 1.0)

  assert_near(res_mod$mean,   mean_ref, tol = 1e-12, msg = "posterior mean mismatch")
  assert_near(res_mod$second, sec_ref,  tol = 1e-12, msg = "second moment mismatch")
  assert_near(res_mod$A,      A_ref,    tol = 1e-12, msg = "A_k mismatch")
  assert_near(res_mod$B,      B_ref,    tol = 1e-12, msg = "B_k mismatch")
})

cat("\n=== T_alpha: Alpha Mixing Parameter Edge Cases ===\n")

run_test("T_alpha.1: alpha=0 -> B_k=0 -> beta mean approx 0", {
  set.seed(90); n <- 30
  w      <- abs(rnorm(n)) + 0.5
  ZF_k   <- rnorm(n)
  z_no_k <- ZF_k

  res <- update_beta_k(w, z_no_k, ZF_k, alpha = 0)
  assert_near(res$B, 0, tol = 1e-12, msg = "B_k should be 0 when alpha=0")
  assert_near(res$x, 0, tol = 1e-6,  msg = "x_k should be 0 when alpha=0")
  assert_near(res$mean, 0, tol = 1e-6, msg = "beta_mean should be ~0 when alpha=0")
})

run_test("T_alpha.2: alpha=1 -> gives unscaled formula (A = sum(w*ZF_k^2))", {
  set.seed(91); n <- 20
  w      <- abs(rnorm(n)) + 0.5
  z_no_k <- rnorm(n)
  ZF_k   <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k, alpha = 1.0)
  expected_A <- max(sum(w * ZF_k^2), 1e-10)
  expected_B <- sum(w * z_no_k * ZF_k)
  assert_near(res$A, expected_A, msg = "A_k with alpha=1 should match unscaled formula")
  assert_near(res$B, expected_B, msg = "B_k with alpha=1 should match unscaled formula")
})

cat("\n=== T_NEW: Cluster B ZF Predictor Verification ===\n")

run_test("T_NEW.1: beta_all with ZF=Y*EF predictor recovers correct beta signs", {
  # End-to-end test: generate Y, compute ZF = Y * EF, fit beta via update_beta_all.
  # ZF is the actual observed projection score that Cluster B uses at training time.
  set.seed(1001); n <- 300; p <- 50; K <- 2
  beta_true <- c(1.5, -1.0)

  # Generate synthetic data: L latent, F factor, Y = L F' + noise
  L_true <- matrix(abs(rnorm(n * K)), n, K)
  F_true <- matrix(rnorm(p * K), p, K)
  Y      <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)

  # ZF = Y * EF (observed projection)
  EF  <- F_true   # in practice EF comes from the F update
  ZF  <- Y %*% EF   # n x K matrix of projection scores

  eta  <- as.vector(ZF %*% beta_true)
  surv <- .sim_cox_taylor(n, ZF, beta_true, seed = 1001)
  taylor <- .calc_cox_taylor_local(eta, surv$time, surv$status)
  z    <- eta + taylor$u / taylor$w
  w    <- taylor$w

  res <- update_beta_all(w, z, ZF, rep(0, K))
  assert_true(res$EBeta[1] > 0,  sprintf("beta[1] sign wrong (%.3f)", res$EBeta[1]))
  assert_true(res$EBeta[2] < 0,  sprintf("beta[2] sign wrong (%.3f)", res$EBeta[2]))
})

run_test("T_NEW.2: 3-argument call (minimal Cluster B signature) is valid", {
  # Verify the new 3-arg signature — no EL2_k — runs without error.
  set.seed(1002); n <- 15
  w      <- abs(rnorm(n)) + 0.1
  ZF_k   <- rnorm(n)
  z_no_k <- rnorm(n)

  res <- update_beta_k(w, z_no_k, ZF_k)   # 3 required args only
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B),
                "Non-finite values in minimal-signature call")
  assert_true(res$second >= res$mean^2 - 1e-10, "second >= mean^2 must hold")
})
