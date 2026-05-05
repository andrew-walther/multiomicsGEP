# =============================================================================
# tests/test_update_L.R
#
# TDD test suite for the q_L CAVI update (patient loadings).
# Written BEFORE implementation (red phase).
#
# Tests cover:
#   T1: Mathematical identity checks (A, B, B_gen, B_surv, x, s, second)
#   T2: Genomics-only limit (w=0 or EBeta=0)
#   T3: Known signal recovery (K=1)
#   T4: Multi-factor recovery (K=5) via update_L_all
#   T5: Null factor shrinkage
#   T6: Error-in-variables (EF2/EBeta2 inflate precision)
#   T7: Numerical stability (degenerate inputs, extreme values)
#   T8: Gauss-Seidel ordering and compute_R_k
#   T9: Consistency with V2.R inline code
#
# NOTE: requires ebnm package.  Install with:
#   install.packages("ebnm")
# =============================================================================

suppressPackageStartupMessages(library(ebnm))

# DO NOT source test_helpers.R or update_L.R here - done by run_tests.R

cat("=== T1: Mathematical Identity Checks ===\n")

run_test("T1.1: A_L[i] = (1-alpha)*sum(Tau*EF2_k) + alpha*w[i]*EBeta2_k with default alpha=0.5", {
  set.seed(101); n <- 5; p <- 3
  Tau     <- abs(rnorm(p)) + 0.1       # p-vector
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + runif(p, 0, 0.2)
  w       <- abs(rnorm(n)) + 0.1       # n-vector
  EBeta_k  <- 1.3
  EBeta2_k <- EBeta_k^2 + 0.05
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)

  # Default alpha=0.5: A_L[i] = 0.5*sum(Tau*EF2_k) + 0.5*w[i]*EBeta2_k, pmax(..., 1e-10)
  gen_scalar <- sum(Tau * EF2_k)
  expected_A <- pmax(0.5 * gen_scalar + 0.5 * w * EBeta2_k, 1e-10)
  assert_near(res$A, expected_A, tol = 1e-10, msg = "A_L mismatch")
})

run_test("T1.2: B_L_gen[i] = (R_k %*% (Tau * EF_k))[i]", {
  set.seed(102); n <- 5; p <- 3
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.1
  EBeta_k  <- 0.5
  EBeta2_k <- EBeta_k^2 + 0.05
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)

  expected_B_gen <- as.vector(R_k %*% (Tau * EF_k))
  assert_near(res$B_gen, expected_B_gen, tol = 1e-10, msg = "B_gen mismatch")
})

run_test("T1.3: B_L_surv[i] = w[i] * z_no_k[i] * EBeta_k", {
  set.seed(103); n <- 5; p <- 3
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.1
  EBeta_k  <- -0.7
  EBeta2_k <- EBeta_k^2 + 0.05
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)

  expected_B_surv <- w * z_no_k * EBeta_k
  assert_near(res$B_surv, expected_B_surv, tol = 1e-10, msg = "B_surv mismatch")
})

run_test("T1.4: x is n-vector", {
  set.seed(104); n <- 7; p <- 4
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.1
  EBeta_k  <- 1.0
  EBeta2_k <- EBeta_k^2 + 0.1
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_length(res$x, n, "x should be n-vector")
})

run_test("T1.5: s[i] = 1/sqrt(A_L[i]) for all i", {
  set.seed(105); n <- 6; p <- 3
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.1
  EBeta_k  <- 0.8
  EBeta2_k <- EBeta_k^2 + 0.05
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_near(res$s, 1 / sqrt(res$A), tol = 1e-10, msg = "s != 1/sqrt(A)")
})

run_test("T1.6: second moment = sd^2 + mean^2 (element-wise for all n)", {
  set.seed(106); n <- 8; p <- 5
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.1
  EBeta_k  <- 1.2
  EBeta2_k <- EBeta_k^2 + 0.1
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_near(res$second, res$sd^2 + res$mean^2, tol = 1e-10,
              msg = "second moment identity violated")
})

cat("\n=== T2: Genomics-Only Limit ===\n")

run_test("T2.1: w=0 -> B_surv=0, combined B is (1-alpha)*B_gen", {
  # With w=0, the raw survival signal B_surv = w*z*EBeta = 0.
  # The combined B = (1-alpha)*B_gen + alpha*B_surv = (1-alpha)*B_gen.
  # With default alpha=0.5: B = 0.5 * B_gen.
  set.seed(201); n <- 10; p <- 5
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- rep(0, n)                # no survival contribution
  EBeta_k  <- 2.0
  EBeta2_k <- EBeta_k^2 + 0.1
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  # B_surv (raw, unweighted) must be 0 when w=0
  assert_near(res$B_surv, rep(0, n), tol = 1e-12, msg = "B_surv should be 0 when w=0")
  # B = (1-0.5)*B_gen + 0.5*0 = 0.5*B_gen with default alpha=0.5
  expected_B <- 0.5 * res$B_gen
  assert_near(res$B, expected_B, tol = 1e-12, msg = "B should equal (1-alpha)*B_gen when w=0")
})

run_test("T2.2: EBeta_k=0 -> B_surv=0", {
  set.seed(202); n <- 10; p <- 5
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.1
  EBeta_k  <- 0                       # zero survival coefficient
  EBeta2_k <- 0.05                    # could still have variance
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_near(res$B_surv, rep(0, n), tol = 1e-12, msg = "B_surv should be 0 when EBeta_k=0")
})

run_test("T2.3: w=0 -> A_L is constant across samples (scalar broadcast) with default alpha=0.5", {
  set.seed(203); n <- 8; p <- 4
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- rep(0, n)
  EBeta_k  <- 1.0
  EBeta2_k <- EBeta_k^2 + 0.1
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  # When w=0, A_L = (1-0.5)*sum(Tau*EF2_k) + 0 = 0.5*sum(Tau*EF2_k) for all i (constant)
  expected_A_scalar <- 0.5 * sum(Tau * EF2_k)
  assert_near(res$A, rep(pmax(expected_A_scalar, 1e-10), n), tol = 1e-10,
              msg = "A_L should be constant across samples when w=0")
})

cat("\n=== T3: Known Signal Recovery K=1 ===\n")

run_test("T3.1: recovery of positive loading column (cor > 0.7)", {
  set.seed(301); n <- 100; p <- 50
  L_true <- abs(rnorm(n)) + 0.5       # positive loadings
  F_true <- rnorm(p)
  Tau    <- rep(5.0, p)
  Y      <- outer(L_true, F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)

  # Initialise EL and EF near truth
  EL_init <- matrix(rnorm(n), n, 1)
  EF      <- matrix(F_true + rnorm(p, sd = 0.1), p, 1)
  EF2     <- EF^2 + 0.01

  # Compute R_k for K=1: R_k = Y (since we add back factor 1)
  R_k <- Y - outer(EL_init[,1], EF[,1]) + outer(EL_init[,1], EF[,1])

  # Survival terms (mild)
  w        <- rep(0.1, n)
  EBeta_k  <- 0.5
  EBeta2_k <- EBeta_k^2 + 0.01
  z_no_k   <- L_true * EBeta_k + rnorm(n, sd = 0.3)

  res <- update_L_k(Tau, EF[,1], EF2[,1], w, EBeta_k, EBeta2_k, R_k, z_no_k)
  cc <- cor(res$mean, L_true)
  assert_true(cc > 0.7, sprintf("Loading recovery too weak: cor=%.3f", cc))
})

run_test("T3.2: survival signal helps loading recovery", {
  set.seed(302); n <- 100; p <- 20
  L_true <- rnorm(n, mean = 2, sd = 1)
  F_true <- rnorm(p, sd = 0.3)         # weak genomics signal
  Tau    <- rep(1.0, p)
  Y      <- outer(L_true, F_true) + matrix(rnorm(n * p, sd = 1), n, p)

  EF      <- matrix(F_true, p, 1)
  EF2     <- EF^2 + 0.01
  R_k     <- Y  # K=1

  # Without survival (w=0)
  res_no_surv <- update_L_k(Tau, EF[,1], EF2[,1],
                             w = rep(0, n), EBeta_k = 0, EBeta2_k = 0,
                             R_k = R_k, z_no_k = rep(0, n))

  # With survival signal (w > 0, beta strong)
  w        <- rep(2.0, n)
  EBeta_k  <- 1.5
  EBeta2_k <- EBeta_k^2 + 0.01
  z_no_k   <- L_true * EBeta_k + rnorm(n, sd = 0.5)

  res_with_surv <- update_L_k(Tau, EF[,1], EF2[,1],
                               w, EBeta_k, EBeta2_k, R_k, z_no_k)

  cor_no   <- cor(res_no_surv$mean, L_true)
  cor_with <- cor(res_with_surv$mean, L_true)
  assert_true(cor_with > cor_no - 0.05,
              sprintf("Survival should help: cor_surv=%.3f, cor_nosurv=%.3f",
                      cor_with, cor_no))
})

cat("\n=== T4: Multi-Factor K=5 via update_L_all ===\n")

run_test("T4.1: output dimensions n x K", {
  set.seed(401); n <- 30; p <- 20; K <- 5
  Y      <- matrix(rnorm(n * p), n, p)
  EL     <- matrix(rnorm(n * K), n, K)
  EL2    <- EL^2 + 0.05
  EF     <- matrix(rnorm(p * K), p, K)
  EF2    <- EF^2 + 0.05
  Tau    <- abs(rnorm(p)) + 0.1
  w      <- abs(rnorm(n)) + 0.1
  z      <- rnorm(n)
  EBeta  <- rnorm(K)
  EBeta2 <- EBeta^2 + 0.05

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau, w, z, EBeta, EBeta2)
  assert_true(nrow(res$EL) == n && ncol(res$EL) == K,
              sprintf("EL dims: expected %dx%d, got %dx%d", n, K, nrow(res$EL), ncol(res$EL)))
  assert_true(nrow(res$EL2) == n && ncol(res$EL2) == K,
              sprintf("EL2 dims: expected %dx%d, got %dx%d", n, K, nrow(res$EL2), ncol(res$EL2)))
  assert_length(res$details, K, "details length")
})

run_test("T4.2: all second moments >= squared means (element-wise)", {
  set.seed(402); n <- 30; p <- 20; K <- 5
  Y      <- matrix(rnorm(n * p), n, p)
  EL     <- matrix(rnorm(n * K), n, K)
  EL2    <- EL^2 + 0.05
  EF     <- matrix(rnorm(p * K), p, K)
  EF2    <- EF^2 + 0.05
  Tau    <- abs(rnorm(p)) + 0.1
  w      <- abs(rnorm(n)) + 0.1
  z      <- rnorm(n)
  EBeta  <- rnorm(K)
  EBeta2 <- EBeta^2 + 0.05

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau, w, z, EBeta, EBeta2)
  diff <- res$EL2 - res$EL^2
  assert_true(all(diff >= -1e-10),
              sprintf("Second moment < squared mean; min diff = %.4e", min(diff)))
})

run_test("T4.3: all outputs finite", {
  set.seed(403); n <- 30; p <- 20; K <- 5
  Y      <- matrix(rnorm(n * p), n, p)
  EL     <- matrix(rnorm(n * K), n, K)
  EL2    <- EL^2 + 0.05
  EF     <- matrix(rnorm(p * K), p, K)
  EF2    <- EF^2 + 0.05
  Tau    <- abs(rnorm(p)) + 0.1
  w      <- abs(rnorm(n)) + 0.1
  z      <- rnorm(n)
  EBeta  <- rnorm(K)
  EBeta2 <- EBeta^2 + 0.05

  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau, w, z, EBeta, EBeta2)
  assert_finite(res$EL, "EL has non-finite values")
  assert_finite(res$EL2, "EL2 has non-finite values")
})

cat("\n=== T5: Null Factor Shrinkage ===\n")

run_test("T5.1: pure noise R_k and z_no_k -> loading column near zero", {
  set.seed(501); n <- 100; p <- 30
  Tau     <- rep(1.0, p)
  EF_k    <- rnorm(p, sd = 0.01)       # near-zero factor
  EF2_k   <- EF_k^2 + 0.001
  w       <- rep(0.5, n)
  EBeta_k  <- 0.01
  EBeta2_k <- EBeta_k^2 + 0.001
  R_k     <- matrix(rnorm(n * p, sd = 1), n, p)  # pure noise
  z_no_k  <- rnorm(n, sd = 1)                     # pure noise

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  # Point-normal should shrink most of these toward zero
  assert_true(mean(abs(res$mean)) < 1.0,
              sprintf("Null factor not shrunk: mean(|L|) = %.3f", mean(abs(res$mean))))
})

run_test("T5.2: zero EF_k and w=0 -> A=floor, B=0, EL[,k]=0", {
  n <- 20; p <- 10
  Tau     <- rep(1.0, p)
  EF_k    <- rep(0, p)                 # zero factor
  EF2_k   <- rep(0, p)
  w       <- rep(0, n)                 # no survival
  EBeta_k  <- 0
  EBeta2_k <- 0
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_near(res$A, rep(1e-10, n), tol = 1e-12, msg = "A should equal floor")
  assert_near(res$B, rep(0, n), tol = 1e-10, msg = "B should be 0")
  assert_near(res$mean, rep(0, n), tol = 1e-6, msg = "EL[,k] should be 0 when B=0")
})

cat("\n=== T6: Error-in-Variables ===\n")

run_test("T6.1: larger EF2 -> larger A_gen -> more shrinkage", {
  set.seed(601); n <- 50; p <- 20
  Tau     <- rep(2.0, p)
  EF_k    <- rnorm(p)
  w       <- rep(1.0, n)
  EBeta_k  <- 0.5
  EBeta2_k <- EBeta_k^2 + 0.01
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  # Low posterior variance in F
  EF2_low <- EF_k^2 + 0.01
  res_low <- update_L_k(Tau, EF_k, EF2_low, w, EBeta_k, EBeta2_k, R_k, z_no_k)

  # High posterior variance in F (same mean, much larger second moment)
  EF2_high <- EF_k^2 + 5.0
  res_high <- update_L_k(Tau, EF_k, EF2_high, w, EBeta_k, EBeta2_k, R_k, z_no_k)

  # Higher EF2 -> larger A_gen part -> larger A overall
  assert_true(all(res_high$A >= res_low$A - 1e-10),
              "Higher EF2 should yield larger A")
  # Larger A -> smaller s
  assert_true(all(res_high$s <= res_low$s + 1e-10),
              "Higher EF2 yields SMALLER s (1/sqrt(larger A))")
})

run_test("T6.2: larger EBeta2 -> larger A_surv -> more shrinkage", {
  set.seed(602); n <- 50; p <- 20
  Tau     <- rep(2.0, p)
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- rep(2.0, n)
  EBeta_k  <- 1.0
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  # Low EBeta2 (near EBeta_k^2)
  EBeta2_low <- EBeta_k^2 + 0.01
  res_low <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_low, R_k, z_no_k)

  # High EBeta2 (large uncertainty in beta)
  EBeta2_high <- EBeta_k^2 + 10.0
  res_high <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_high, R_k, z_no_k)

  # Higher EBeta2 -> larger survival contribution to A
  assert_true(all(res_high$A >= res_low$A - 1e-10),
              "Higher EBeta2 should yield larger A")
})

run_test("T6.3: B decomposition: B = (1-alpha)*B_gen + alpha*B_surv exactly", {
  # Verify the weighted decomposition identity holds for any alpha.
  # When genomics and survival signals reinforce, both B_gen and B_surv
  # should have the same sign on average; the weighted sum B reflects both.
  set.seed(603); n <- 50; p <- 20
  L_true  <- abs(rnorm(n)) + 1.0       # positive loadings
  F_true  <- rnorm(p)
  Tau     <- rep(2.0, p)
  EF_k    <- F_true
  EF2_k   <- EF_k^2 + 0.01
  R_k     <- outer(L_true, F_true) + matrix(rnorm(n * p, sd = 0.3), n, p)

  # Survival signal aligned with L_true
  EBeta_k  <- 1.0
  EBeta2_k <- EBeta_k^2 + 0.01
  z_no_k   <- L_true * EBeta_k + rnorm(n, sd = 0.2)
  w        <- rep(2.0, n)

  test_alpha <- 0.4
  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k, alpha = test_alpha)

  # The weighted decomposition identity must hold exactly:
  # B = (1 - alpha)*B_gen + alpha*B_surv (raw unweighted components)
  expected_B <- (1 - test_alpha) * res$B_gen + test_alpha * res$B_surv
  assert_near(res$B, expected_B, tol = 1e-12,
              msg = "B decomposition: B != (1-alpha)*B_gen + alpha*B_surv")

  # With signals reinforcing and alpha=0.4, the combined x should have
  # the same sign as both B_gen and B_surv on average (directional test).
  # Check that the signs of x and B_gen agree for most samples.
  sign_agree <- mean(sign(res$x) == sign(res$B_gen))
  assert_true(sign_agree >= 0.6,
              sprintf("Sign agreement of x and B_gen = %.2f, expected >= 0.6 when signals reinforce",
                      sign_agree))
})

cat("\n=== T7: Numerical Stability ===\n")

run_test("T7.1: all-zero Tau and w -> A hits floor, no NaN/Inf", {
  n <- 15; p <- 8
  Tau     <- rep(0, p)
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- rep(0, n)
  EBeta_k  <- 1.0
  EBeta2_k <- EBeta_k^2 + 0.1
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_near(res$A, rep(1e-10, n), tol = 1e-12, msg = "A should equal floor when Tau=0, w=0")
  assert_finite(c(res$mean, res$second, res$sd, res$x, res$s),
                "Non-finite values with zero Tau and w")
})

run_test("T7.2: extreme Tau=1e8 with p=10 -> no overflow", {
  set.seed(702); n <- 10; p <- 10
  Tau     <- rep(1e8, p)
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.01
  w       <- rep(1.0, n)
  EBeta_k  <- 0.5
  EBeta2_k <- EBeta_k^2 + 0.01
  R_k     <- matrix(rnorm(n * p, sd = 0.001), n, p)
  z_no_k  <- rnorm(n, sd = 0.001)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)
  assert_finite(c(res$mean, res$second, res$sd, res$A, res$B),
                "Non-finite values with extreme Tau")
})

run_test("T7.3: n=1, p=1 -> works correctly with default alpha=0.5", {
  Tau     <- 3.0
  EF_k    <- 1.5
  EF2_k   <- 1.5^2 + 0.1
  w       <- 2.0
  EBeta_k  <- 0.8
  EBeta2_k <- 0.8^2 + 0.05
  R_k     <- matrix(0.7, 1, 1)
  z_no_k  <- 1.2

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)

  expected_A <- pmax(0.5 * sum(Tau * EF2_k) + 0.5 * w * EBeta2_k, 1e-10)
  expected_B_gen  <- as.vector(R_k %*% (Tau * EF_k))
  expected_B_surv <- w * z_no_k * EBeta_k
  expected_B <- 0.5 * expected_B_gen + 0.5 * expected_B_surv

  assert_near(res$A, expected_A, tol = 1e-10, msg = "A wrong for n=1, p=1")
  assert_near(res$B, expected_B, tol = 1e-10, msg = "B wrong for n=1, p=1")
  assert_finite(c(res$mean, res$second), "Non-finite values for n=1, p=1")
})

run_test("T7.4: custom A_floor is respected", {
  n <- 10; p <- 5
  Tau     <- rep(0, p)
  EF_k    <- rep(0, p)
  EF2_k   <- rep(0, p)
  w       <- rep(0, n)
  EBeta_k  <- 0
  EBeta2_k <- 0
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k,
                    A_floor = 1.0)
  assert_near(res$A, rep(1.0, n), tol = 1e-10, msg = "A should equal custom floor")
  assert_near(res$s, rep(1 / sqrt(1.0), n), tol = 1e-10, msg = "s should use custom floor")
})

cat("\n=== T8: Gauss-Seidel and R_k ===\n")

run_test("T8.1: compute_R_k excludes factor k correctly", {
  set.seed(801); n <- 8; p <- 6; K <- 3
  EL <- matrix(rnorm(n * K), n, K)
  EF <- matrix(rnorm(p * K), p, K)

  # Fake Y = EL %*% t(EF) + noise
  Y <- EL %*% t(EF) + matrix(rnorm(n * p, sd = 0.1), n, p)

  for (k in 1:K) {
    R_k <- compute_R_k(Y, EL, EF, k)

    # Manual: R_k = Y - Y_hat + outer(EL[,k], EF[,k])
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
  Y      <- matrix(rnorm(n * p), n, p)
  EL     <- matrix(rnorm(n * K), n, K)
  EL2    <- EL^2 + 0.05
  EF     <- matrix(rnorm(p * K), p, K)
  EF2    <- EF^2 + 0.05
  Tau    <- abs(rnorm(p)) + 0.1
  w      <- abs(rnorm(n)) + 0.1
  z      <- rnorm(n)
  EBeta  <- c(1.0, -0.5)
  EBeta2 <- EBeta^2 + 0.05

  # Run update_L_all
  res <- update_L_all(Y, EL, EL2, EF, EF2, Tau, w, z, EBeta, EBeta2)

  # The Gauss-Seidel property means EL[,1] is updated before k=2's R_k is computed.
  # We can verify by checking that the k=2 detail's A and B are different from
  # what they would be if the OLD EL[,1] were used.

  # Recompute k=2 using the OLD EL[,1] (Jacobi-style)
  R_k2_old <- compute_R_k(Y, EL, EF, k = 2)
  B_gen_old <- as.vector(R_k2_old %*% (Tau * EF[, 2]))

  # The actual B_gen from k=2 should differ (Gauss-Seidel used updated EL[,1])
  B_gen_actual <- res$details[[2]]$B_gen

  # They should be different unless EL[,1] didn't change (astronomically unlikely)
  EL1_changed <- !all(abs(res$EL[, 1] - EL[, 1]) < 1e-15)
  if (EL1_changed) {
    assert_true(max(abs(B_gen_actual - B_gen_old)) > 1e-10,
                "update_L_all should use Gauss-Seidel: k=2's B_gen should differ from Jacobi")
  }
  # Structure check regardless
  assert_finite(c(res$EL, res$EL2), "Gauss-Seidel outputs non-finite")
})

cat("\n=== T9: V2.R Consistency ===\n")

run_test("T9.1: update_L_k with explicit alpha reproduces alpha-weighted reference exactly", {
  # V2.R used alpha=1 (unweighted sum). With the new parameterisation, V2.R is NOT
  # directly reproducible via alpha=1 because alpha=1 means pure survival (no genomics).
  # Instead, we verify with alpha=0.7: compute the reference from scratch and confirm
  # the modular function matches it exactly.
  set.seed(999); n <- 50; p <- 30
  Tau      <- abs(rnorm(p)) + 0.1
  EF_k     <- rnorm(p)
  EF2_k    <- EF_k^2 + runif(p, 0.01, 0.1)
  w        <- abs(rnorm(n)) + 0.5
  EBeta_k  <- 1.3
  EBeta2_k <- EBeta_k^2 + runif(1, 0.01, 0.1)
  R_k      <- matrix(rnorm(n * p), n, p)
  z_no_k   <- rnorm(n)
  test_alpha <- 0.7

  # ---------- alpha-weighted reference ----------
  # $B_gen and $B_surv are the RAW (unweighted) components.
  # $B is the weighted combination: (1-alpha)*B_gen + alpha*B_surv.
  A_L_ref <- pmax((1 - test_alpha) * sum(Tau * EF2_k) + test_alpha * w * EBeta2_k, 1e-10)
  B_L_gen_ref  <- as.vector(R_k %*% (Tau * EF_k))       # raw
  B_L_surv_ref <- w * z_no_k * EBeta_k                  # raw
  B_L_ref <- (1 - test_alpha) * B_L_gen_ref + test_alpha * B_L_surv_ref

  res_ref <- ebnm(x = B_L_ref / A_L_ref, s = 1 / sqrt(A_L_ref),
                  prior_family = "point_normal")
  mean_ref   <- res_ref$posterior$mean
  second_ref <- res_ref$posterior$sd^2 + res_ref$posterior$mean^2

  # ---------- Modular function with same alpha ----------
  res_mod <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k,
                        prior_family = "point_normal", alpha = test_alpha)

  assert_near(res_mod$A,      A_L_ref,       tol = 1e-12, msg = "A mismatch vs reference")
  assert_near(res_mod$B,      B_L_ref,       tol = 1e-12, msg = "B mismatch vs reference")
  # $B_gen and $B_surv are raw unweighted components
  assert_near(res_mod$B_gen,  B_L_gen_ref,   tol = 1e-12, msg = "B_gen mismatch (should be raw/unweighted)")
  assert_near(res_mod$B_surv, B_L_surv_ref,  tol = 1e-12, msg = "B_surv mismatch (should be raw/unweighted)")
  # Verify the weighted decomposition: B = (1-alpha)*B_gen + alpha*B_surv
  assert_near(res_mod$B,
              (1 - test_alpha) * res_mod$B_gen + test_alpha * res_mod$B_surv,
              tol = 1e-12, msg = "B should equal (1-alpha)*B_gen + alpha*B_surv")
  assert_near(res_mod$mean,   mean_ref,   tol = 1e-12, msg = "posterior mean mismatch vs reference")
  assert_near(res_mod$second, second_ref, tol = 1e-12, msg = "second moment mismatch vs reference")
})

run_test("T9.2: compute_R_k matches V2.R lines 290-291 exactly", {
  set.seed(998); n <- 40; p <- 25; K <- 4
  Y  <- matrix(rnorm(n * p), n, p)
  EL <- matrix(rnorm(n * K), n, K)
  EF <- matrix(rnorm(p * K), p, K)

  for (k in 1:K) {
    # V2.R inline (lines 290-291)
    Y_hat_v2 <- EL %*% t(EF)
    R_k_v2   <- Y - Y_hat_v2 + outer(EL[, k], EF[, k])

    # Modular function
    R_k_mod <- compute_R_k(Y, EL, EF, k)

    assert_near(R_k_mod, R_k_v2, tol = 1e-12,
                msg = sprintf("compute_R_k mismatch vs V2.R for k=%d", k))
  }
})

cat("\n=== T_alpha: Alpha Mixing Parameter Edge Cases ===\n")

run_test("T_alpha.1: alpha=0 -> pure genomics (A_L=A_gen, B_L=B_gen)", {
  # alpha=0 means no survival contribution: A = (1-0)*A_gen + 0*A_surv = A_gen
  # and B = (1-0)*B_gen + 0*B_surv = B_gen (raw, unweighted).
  set.seed(901); n <- 20; p <- 10
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.5
  EBeta_k  <- 1.5
  EBeta2_k <- EBeta_k^2 + 0.05
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k, alpha = 0)

  A_gen_expected <- pmax(rep(sum(Tau * EF2_k), n), 1e-10)
  B_gen_expected <- as.vector(R_k %*% (Tau * EF_k))

  assert_near(res$A, A_gen_expected, tol = 1e-10, msg = "alpha=0: A_L should equal A_gen")
  assert_near(res$B, B_gen_expected, tol = 1e-10, msg = "alpha=0: B_L should equal B_gen")
  # $B_gen is the raw component and should equal $B when alpha=0
  assert_near(res$B_gen, B_gen_expected, tol = 1e-10, msg = "B_gen should be unweighted")
})

run_test("T_alpha.2: alpha=1 -> pure survival (A_L=A_surv, B_L=B_surv)", {
  # alpha=1 means no genomics contribution: A = 0*A_gen + 1*A_surv = A_surv
  # and B = 0*B_gen + 1*B_surv = B_surv.
  set.seed(902); n <- 20; p <- 10
  Tau     <- abs(rnorm(p)) + 0.1
  EF_k    <- rnorm(p)
  EF2_k   <- EF_k^2 + 0.1
  w       <- abs(rnorm(n)) + 0.5
  EBeta_k  <- 1.5
  EBeta2_k <- EBeta_k^2 + 0.05
  R_k     <- matrix(rnorm(n * p), n, p)
  z_no_k  <- rnorm(n)

  res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k, alpha = 1)

  A_surv_expected <- pmax(w * EBeta2_k, 1e-10)
  B_surv_expected <- w * z_no_k * EBeta_k

  assert_near(res$A, A_surv_expected, tol = 1e-10, msg = "alpha=1: A_L should equal A_surv")
  assert_near(res$B, B_surv_expected, tol = 1e-10, msg = "alpha=1: B_L should equal B_surv")
  # $B_surv is the raw component and should equal $B when alpha=1
  assert_near(res$B_surv, B_surv_expected, tol = 1e-10, msg = "B_surv should be unweighted")
})
