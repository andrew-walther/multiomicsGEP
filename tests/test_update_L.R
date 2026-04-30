# =============================================================================
# tests/test_update_L.R
#
# TDD test suite for the q_L CAVI update (patient loadings) — Cluster B
# (Cox-on-YF reformulation): pure-genomics EBNM.
#
# Under eta = Z_F * beta_tilde (Cluster B), L appears only in the genomics
# likelihood.  All survival arguments (w, EBeta, EBeta2, z_no_k, alpha,
# lambda, normalize_AB) have been removed from update_L_k() and update_L_all().
#
# Tests cover:
#   T1: Mathematical identity checks (A, B, x, s, second)
#   T3: Known signal recovery (K=1)
#   T4: Multi-factor recovery (K=5) via update_L_all
#   T5: Null factor shrinkage
#   T6: Error-in-variables (EF2 inflates precision A)
#   T7: Numerical stability (degenerate inputs, extreme values)
#   T8: Gauss-Seidel ordering and compute_R_k
#   T9: Consistency with pure-genomics formula reference
#   T_NEW: Cluster B simplified-signature verification
#
# NOTE: requires ebnm package.  Install with:
#   install.packages("ebnm")
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# DO NOT source test_helpers.R or update_L.R here - done by run_tests.R

cat("=== T1: Mathematical Identity Checks ===\n")

run_test("T1.1: A_L = max(sum(Tau*EF2_k), A_floor) — scalar", {
  set.seed(101); n <- 5; p <- 3
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + runif(p, 0, 0.2)
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)

  expected_A <- max(sum(Tau * EF2_k), 1e-10)
  assert_near(res$A, expected_A, tol = 1e-10, msg = "A_L mismatch")
})

run_test("T1.2: B_L[i] = (R_k %*% (Tau * EF_k))[i]", {
  set.seed(102); n <- 5; p <- 3
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)

  expected_B <- as.vector(R_k %*% (Tau * EF_k))
  assert_near(res$B, expected_B, tol = 1e-10, msg = "B_L mismatch")
})

run_test("T1.3: x is n-vector", {
  set.seed(103); n <- 7; p <- 4
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_length(res$x, n, "x should be n-vector")
})

run_test("T1.4: s = 1/sqrt(A_L)", {
  set.seed(104); n <- 6; p <- 3
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_near(res$s, 1 / sqrt(res$A), tol = 1e-10, msg = "s != 1/sqrt(A)")
})

run_test("T1.5: second moment = sd^2 + mean^2 (element-wise for all n)", {
  set.seed(105); n <- 8; p <- 5
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_near(res$second, res$sd^2 + res$mean^2, tol = 1e-10,
              msg = "second moment identity violated")
})

cat("\n=== T3: Known Signal Recovery K=1 ===\n")

run_test("T3.1: recovery of positive loading column (cor > 0.7)", {
  set.seed(301); n <- 100; p <- 50
  L_true <- abs(rnorm(n)) + 0.5       # positive loadings
  F_true <- rnorm(p)
  Tau    <- rep(5.0, p)
  Y      <- outer(L_true, F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)

  EL_init <- matrix(rnorm(n), n, 1)
  EF      <- matrix(F_true + rnorm(p, sd = 0.1), p, 1)
  EF2     <- EF^2 + 0.01

  # K=1: adding back factor 1 equals Y itself
  R_k <- Y - outer(EL_init[,1], EF[,1]) + outer(EL_init[,1], EF[,1])

  res <- update_L_k(Tau, EF[,1], EF2[,1], R_k)
  cc <- cor(res$mean, L_true)
  assert_true(cc > 0.7, sprintf("Loading recovery too weak: cor=%.3f", cc))
})

cat("\n=== T4: Multi-Factor K=5 via update_L_all ===\n")

run_test("T4.1: output dimensions n x K", {
  set.seed(401); n <- 30; p <- 20; K <- 5
  Y   <- matrix(rnorm(n * p), n, p)
  EL  <- matrix(rnorm(n * K), n, K)
  EL2 <- EL^2 + 0.05
  EF  <- matrix(rnorm(p * K), p, K)
  EF2 <- EF^2 + 0.05
  Tau <- abs(rnorm(p)) + 0.1

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau)
  assert_true(nrow(res$EL) == n && ncol(res$EL) == K,
              sprintf("EL dims: expected %dx%d, got %dx%d", n, K, nrow(res$EL), ncol(res$EL)))
  assert_true(nrow(res$EL2) == n && ncol(res$EL2) == K,
              sprintf("EL2 dims: expected %dx%d, got %dx%d", n, K, nrow(res$EL2), ncol(res$EL2)))
  assert_length(res$details, K, "details length")
})

run_test("T4.2: all second moments >= squared means (element-wise)", {
  set.seed(402); n <- 30; p <- 20; K <- 5
  Y   <- matrix(rnorm(n * p), n, p)
  EL  <- matrix(rnorm(n * K), n, K)
  EL2 <- EL^2 + 0.05
  EF  <- matrix(rnorm(p * K), p, K)
  EF2 <- EF^2 + 0.05
  Tau <- abs(rnorm(p)) + 0.1

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau)
  diff <- res$EL2 - res$EL^2
  assert_true(all(diff >= -1e-10),
              sprintf("Second moment < squared mean; min diff = %.4e", min(diff)))
})

run_test("T4.3: all outputs finite", {
  set.seed(403); n <- 30; p <- 20; K <- 5
  Y   <- matrix(rnorm(n * p), n, p)
  EL  <- matrix(rnorm(n * K), n, K)
  EL2 <- EL^2 + 0.05
  EF  <- matrix(rnorm(p * K), p, K)
  EF2 <- EF^2 + 0.05
  Tau <- abs(rnorm(p)) + 0.1

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau)
  assert_finite(res$EL,  "EL has non-finite values")
  assert_finite(res$EL2, "EL2 has non-finite values")
})

cat("\n=== T5: Null Factor Shrinkage ===\n")

run_test("T5.1: pure noise R_k -> loading column near zero", {
  set.seed(501); n <- 100; p <- 30
  Tau   <- rep(1.0, p)
  EF_k  <- rnorm(p, sd = 0.01)       # near-zero factor
  EF2_k <- EF_k^2 + 0.001
  R_k   <- matrix(rnorm(n * p, sd = 1), n, p)  # pure noise

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_true(mean(abs(res$mean)) < 1.0,
              sprintf("Null factor not shrunk: mean(|L|) = %.3f", mean(abs(res$mean))))
})

run_test("T5.2: EF_k=0, EF2_k=0 -> A=floor, B=0, EL[,k]=~0", {
  n <- 20; p <- 10
  Tau   <- rep(1.0, p)
  EF_k  <- rep(0, p)
  EF2_k <- rep(0, p)
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_near(res$A, 1e-10, tol = 1e-12, msg = "A should equal floor")
  assert_near(res$B, rep(0, n), tol = 1e-10, msg = "B should be 0")
  assert_near(res$mean, rep(0, n), tol = 1e-6, msg = "EL[,k] should be 0 when B=0")
})

cat("\n=== T6: Error-in-Variables (EF2 inflates A) ===\n")

run_test("T6.1: larger EF2 -> larger A_gen -> smaller s (more effective shrinkage)", {
  set.seed(601); n <- 50; p <- 20
  Tau  <- rep(2.0, p)
  EF_k <- rnorm(p)
  R_k  <- matrix(rnorm(n * p), n, p)

  EF2_low  <- EF_k^2 + 0.01
  res_low  <- update_L_k(Tau, EF_k, EF2_low, R_k)

  EF2_high <- EF_k^2 + 5.0
  res_high <- update_L_k(Tau, EF_k, EF2_high, R_k)

  # Higher EF2 -> larger A_gen -> larger A
  assert_true(res_high$A >= res_low$A - 1e-10,
              "Higher EF2 should yield larger A")
  # Larger A -> smaller s
  assert_true(res_high$s <= res_low$s + 1e-10,
              "Higher EF2 yields SMALLER s (1/sqrt(larger A))")
})

cat("\n=== T7: Numerical Stability ===\n")

run_test("T7.1: all-zero Tau -> A hits floor, no NaN/Inf", {
  n <- 15; p <- 8
  Tau   <- rep(0, p)
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_near(res$A, 1e-10, tol = 1e-12, msg = "A should equal floor when Tau=0")
  assert_finite(c(res$mean, res$second, res$sd, res$x, res$s),
                "Non-finite values with zero Tau")
})

run_test("T7.2: extreme Tau=1e8 with p=10 -> no overflow", {
  set.seed(702); n <- 10; p <- 10
  Tau   <- rep(1e8, p)
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.01
  R_k   <- matrix(rnorm(n * p, sd = 0.001), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B),
                "Non-finite values with extreme Tau")
})

run_test("T7.3: n=1, p=1 -> works correctly", {
  Tau   <- 3.0
  EF_k  <- 1.5
  EF2_k <- 1.5^2 + 0.1
  R_k   <- matrix(0.7, 1, 1)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)

  expected_A <- max(sum(Tau * EF2_k), 1e-10)
  expected_B <- as.vector(R_k %*% (Tau * EF_k))

  assert_near(res$A, expected_A, tol = 1e-10, msg = "A wrong for n=1, p=1")
  assert_near(res$B, expected_B, tol = 1e-10, msg = "B wrong for n=1, p=1")
  assert_finite(c(res$mean, res$second), "Non-finite values for n=1, p=1")
})

run_test("T7.4: custom A_floor is respected", {
  n <- 10; p <- 5
  Tau   <- rep(0, p)
  EF_k  <- rep(0, p)
  EF2_k <- rep(0, p)
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k, A_floor = 1.0)
  assert_near(res$A, 1.0, tol = 1e-10, msg = "A should equal custom floor")
  assert_near(res$s, 1 / sqrt(1.0), tol = 1e-10, msg = "s should use custom floor")
})

cat("\n=== T8: Gauss-Seidel and R_k ===\n")

run_test("T8.1: compute_R_k excludes factor k correctly", {
  set.seed(801); n <- 8; p <- 6; K <- 3
  EL <- matrix(rnorm(n * K), n, K)
  EF <- matrix(rnorm(p * K), p, K)
  Y  <- EL %*% t(EF) + matrix(rnorm(n * p, sd = 0.1), n, p)

  for (k in 1:K) {
    R_k      <- compute_R_k(Y, EL, EF, k)
    Y_hat    <- EL %*% t(EF)
    expected <- Y - Y_hat + outer(EL[, k], EF[, k])
    assert_near(R_k, expected, tol = 1e-12,
                msg = sprintf("compute_R_k wrong for k=%d", k))
  }
})

run_test("T8.2: compute_R_k has correct dimensions n x p", {
  set.seed(802); n <- 12; p <- 7; K <- 4
  Y  <- matrix(rnorm(n * p), n, p)
  EL <- matrix(rnorm(n * K), n, K)
  EF <- matrix(rnorm(p * K), p, K)

  R_k <- compute_R_k(Y, EL, EF, k = 2)
  assert_true(nrow(R_k) == n && ncol(R_k) == p,
              sprintf("R_k dims: expected %dx%d, got %dx%d", n, p, nrow(R_k), ncol(R_k)))
})

run_test("T8.3: update_L_all uses Gauss-Seidel (EL[,1] changes before k=2's R_k)", {
  set.seed(803); n <- 30; p <- 15; K <- 2
  Y   <- matrix(rnorm(n * p), n, p)
  EL  <- matrix(rnorm(n * K), n, K)
  EL2 <- EL^2 + 0.05
  EF  <- matrix(rnorm(p * K), p, K)
  EF2 <- EF^2 + 0.05
  Tau <- abs(rnorm(p)) + 0.1

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau)

  # Gauss-Seidel: EL[,1] is updated before k=2's R_k is computed.
  # Recompute k=2 using the OLD EL[,1] (Jacobi-style) and verify B differs.
  R_k2_old <- compute_R_k(Y, EL, EF, k = 2)      # old EL[,1]
  B_old     <- as.vector(R_k2_old %*% (Tau * EF[, 2]))
  B_actual  <- res$details[[2]]$B

  EL1_changed <- !all(abs(res$EL[, 1] - EL[, 1]) < 1e-15)
  if (EL1_changed) {
    assert_true(max(abs(B_actual - B_old)) > 1e-10,
                "update_L_all should use Gauss-Seidel: k=2's B should differ from Jacobi")
  }
  assert_finite(c(res$EL, res$EL2), "Gauss-Seidel outputs non-finite")
})

cat("\n=== T9: Pure-Genomics Formula Reference Consistency ===\n")

run_test("T9.1: update_L_k matches manually computed pure-genomics EBNM exactly", {
  set.seed(999); n <- 50; p <- 30
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + runif(p, 0.01, 0.1)
  R_k   <- matrix(rnorm(n * p), n, p)

  # Manual reference: A = max(sum(Tau*EF2_k), floor), B = R_k %*% (Tau*EF_k)
  A_ref      <- max(sum(Tau * EF2_k), 1e-10)
  B_ref      <- as.vector(R_k %*% (Tau * EF_k))
  ebnm_ref   <- ebnm(x = B_ref / A_ref, s = 1 / sqrt(A_ref),
                      prior_family = "point_normal")
  mean_ref   <- ebnm_ref$posterior$mean
  second_ref <- ebnm_ref$posterior$sd^2 + ebnm_ref$posterior$mean^2

  res <- update_L_k(Tau, EF_k, EF2_k, R_k, prior_family = "point_normal")

  assert_near(res$A,      A_ref,      tol = 1e-12, msg = "A mismatch vs reference")
  assert_near(res$B,      B_ref,      tol = 1e-12, msg = "B mismatch vs reference")
  assert_near(res$mean,   mean_ref,   tol = 1e-12, msg = "posterior mean mismatch")
  assert_near(res$second, second_ref, tol = 1e-12, msg = "second moment mismatch")
})

run_test("T9.2: compute_R_k matches V2.R lines 290-291 exactly", {
  set.seed(998); n <- 40; p <- 25; K <- 4
  Y  <- matrix(rnorm(n * p), n, p)
  EL <- matrix(rnorm(n * K), n, K)
  EF <- matrix(rnorm(p * K), p, K)

  for (k in 1:K) {
    Y_hat_v2 <- EL %*% t(EF)
    R_k_v2   <- Y - Y_hat_v2 + outer(EL[, k], EF[, k])
    R_k_mod  <- compute_R_k(Y, EL, EF, k)
    assert_near(R_k_mod, R_k_v2, tol = 1e-12,
                msg = sprintf("compute_R_k mismatch vs V2.R for k=%d", k))
  }
})

cat("\n=== T_NEW: Cluster B Simplified Signature ===\n")

run_test("T_NEW.1: A is scalar (length 1), not n-vector", {
  # Under Cluster B, A_L = sum_j(tau_j * E[f^2_jk]) is the same for all
  # patients i (no Cox weights W_ii). The return value $A must be length 1.
  set.seed(1001); n <- 20; p <- 10
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)
  assert_length(res$A, 1L, "A should be scalar (length 1) under Cluster B")
  assert_length(res$s, 1L, "s should be scalar (length 1) under Cluster B")
})

run_test("T_NEW.2: 4-argument minimal call is valid and returns correct structure", {
  # Verify the new 4-arg signature — no w, EBeta, EBeta2, z_no_k —
  # runs without error and produces valid EBNM output.
  set.seed(1002); n <- 15; p <- 8
  Tau   <- abs(rnorm(p)) + 0.1
  EF_k  <- rnorm(p)
  EF2_k <- EF_k^2 + 0.1
  R_k   <- matrix(rnorm(n * p), n, p)

  res <- update_L_k(Tau, EF_k, EF2_k, R_k)   # 4 required args only
  assert_length(res$mean,   n, "mean should be n-vector")
  assert_length(res$second, n, "second should be n-vector")
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B),
                "Non-finite values in minimal-signature call")
  assert_true(all(res$second >= res$mean^2 - 1e-10),
              "second >= mean^2 must hold")
})
