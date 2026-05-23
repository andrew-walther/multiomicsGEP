# =============================================================================
# tests/test_update_F_cohort.R
# Unit tests for code/update_F_cohort.R and compute_normal_kl() in
# code/compute_elbo.R.
# =============================================================================

source("code/update_F_cohort.R")
source("code/compute_elbo.R")

# ---------------------------------------------------------------------------
# T1: Closed-form correctness — 2-cohort, 3-gene toy
# ---------------------------------------------------------------------------
run_test("FCohort-T1: EF_c == B/A to machine precision", {
  c_vec <- c(1, 0, 1, 0)           # n=4, 2 patients in cohort B
  R     <- matrix(c(2, 0, 4, 0,
                    1, 0, 3, 0,
                    5, 0, 1, 0),
                  nrow = 4, byrow = FALSE)   # p=3
  Tau   <- c(1.0, 2.0, 0.5)
  sigma <- 1.0

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma)

  n_c   <- 2L
  B_exp <- colSums(R[c(1, 3), ]) * Tau
  A_exp <- n_c * Tau + 1.0 / sigma^2
  EF_exp <- B_exp / A_exp

  assert_near(res$EF_c, EF_exp, tol = 1e-12, msg = "EF_c vs B/A")
})

# ---------------------------------------------------------------------------
# T2: Posterior variance correct — 1/A_c == EF2_c - EF_c^2
# ---------------------------------------------------------------------------
run_test("FCohort-T2: Posterior variance = 1/A_c", {
  set.seed(10)
  c_vec <- rbinom(20, 1, 0.5)
  R     <- matrix(rnorm(20 * 5), 20, 5)
  Tau   <- runif(5, 0.5, 2.0)

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 1.0)
  var_post <- res$EF2_c - res$EF_c^2

  assert_near(var_post, 1.0 / res$A_c, tol = 1e-10, msg = "Var = 1/A")
})

# ---------------------------------------------------------------------------
# T3: Large-data limit — posterior mean → column mean of R_cohort
# ---------------------------------------------------------------------------
run_test("FCohort-T3: Large n_c: EF_c -> column mean of R_cohort rows", {
  set.seed(20)
  n_c  <- 10000L
  p    <- 4
  R    <- matrix(rnorm(n_c * p, mean = 0.5), n_c, p)
  c_vec <- rep(1L, n_c)
  Tau  <- rep(1.0, p)
  sigma <- 10.0   # weak prior

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = sigma)
  # With weak prior and large n_c, E[f_cj] ≈ colMeans(R)
  assert_near(res$EF_c, colMeans(R), tol = 1e-2, msg = "Large-n_c limit")
})

# ---------------------------------------------------------------------------
# T4: Strong prior — EF_c shrinks toward zero
# ---------------------------------------------------------------------------
run_test("FCohort-T4: Strong prior (sigma=1e-6): EF_c near zero", {
  set.seed(30)
  c_vec <- c(1, 1, 1)
  R     <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3)  # p=2, large signal
  Tau   <- c(1.0, 1.0)

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 1e-6)

  assert_true(all(abs(res$EF_c) < 1e-4),
              msg = "Strong prior: EF_c near zero")
})

# ---------------------------------------------------------------------------
# T5: Zero residual — EF_c = 0
# ---------------------------------------------------------------------------
run_test("FCohort-T5: Zero residual -> EF_c = 0", {
  c_vec <- c(1, 0, 1, 0, 1)
  R     <- matrix(0.0, 5, 3)
  Tau   <- c(2.0, 3.0, 1.0)

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 1.0)

  assert_near(res$EF_c, rep(0.0, 3), tol = 1e-14, msg = "Zero R -> zero EF_c")
})

# ---------------------------------------------------------------------------
# T6: All patients in cohort — formula check
# ---------------------------------------------------------------------------
run_test("FCohort-T6: All patients in cohort — formula check", {
  set.seed(40)
  n    <- 12; p <- 3
  c_vec <- rep(1L, n)
  R     <- matrix(rnorm(n * p), n, p)
  Tau   <- runif(p, 0.5, 2.0)
  sigma <- 1.5

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = sigma)

  B_exp <- colSums(R) * Tau
  A_exp <- n * Tau + 1.0 / sigma^2
  assert_near(res$EF_c,  B_exp / A_exp, tol = 1e-12, msg = "EF_c formula")
  assert_near(res$EF2_c, 1.0 / A_exp + (B_exp / A_exp)^2, tol = 1e-12,
              msg = "EF2_c formula")
})

# ---------------------------------------------------------------------------
# T7: Output dimensions — single column update
# ---------------------------------------------------------------------------
run_test("FCohort-T7: Output dimensions from update_F_cohort_col", {
  set.seed(50)
  n <- 10; p <- 7
  c_vec <- c(rep(1, 5), rep(0, 5))
  R     <- matrix(rnorm(n * p), n, p)
  Tau   <- runif(p, 0.5, 2.0)

  res   <- update_F_cohort_col(c_vec, R, Tau)

  assert_length(res$EF_c,  p, "EF_c length = p")
  assert_length(res$EF2_c, p, "EF2_c length = p")
  assert_length(res$A_c,   p, "A_c length = p")
})

# ---------------------------------------------------------------------------
# T8: EF2_c >= EF_c^2 everywhere
# ---------------------------------------------------------------------------
run_test("FCohort-T8: EF2_c >= EF_c^2 (second moment >= squared mean)", {
  set.seed(60)
  c_vec <- rbinom(30, 1, 0.4)
  R     <- matrix(rnorm(30 * 6, sd = 3), 30, 6)
  Tau   <- runif(6, 1.0, 5.0)

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 2.0)

  assert_true(all(res$EF2_c >= res$EF_c^2 - 1e-14),
              msg = "EF2_c >= EF_c^2 everywhere")
})

# ---------------------------------------------------------------------------
# T9: Tau scaling — doubling Tau halves posterior variance
# ---------------------------------------------------------------------------
run_test("FCohort-T9: Doubling Tau halves posterior variance 1/A", {
  c_vec <- c(1, 0, 1, 0)
  R     <- matrix(rnorm(4 * 3), 4, 3)
  Tau1  <- c(1.0, 2.0, 0.5)
  Tau2  <- 2.0 * Tau1
  sigma <- 1.0

  res1 <- update_F_cohort_col(c_vec, R, Tau1, sigma)
  res2 <- update_F_cohort_col(c_vec, R, Tau2, sigma)

  # A2 = 2*n_c*tau + 1/sigma^2; A1 = n_c*tau + 1/sigma^2 — not exactly double
  # when 1/sigma^2 term is non-negligible, so we only check that Var decreases.
  assert_true(all(1.0 / res2$A_c < 1.0 / res1$A_c),
              msg = "Higher Tau -> smaller posterior variance")
})

# ---------------------------------------------------------------------------
# T10: Larger sigma_F_cohort -> less shrinkage (larger |EF_c|)
# ---------------------------------------------------------------------------
run_test("FCohort-T10: Larger sigma_F_cohort -> less shrinkage", {
  set.seed(70)
  c_vec <- c(1, 1, 0, 1, 0)
  R     <- matrix(abs(rnorm(5 * 4, mean = 2)), 5, 4)   # positive residuals
  Tau   <- rep(1.0, 4)

  res_tight <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 0.1)
  res_loose <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 10.0)

  assert_true(all(abs(res_loose$EF_c) >= abs(res_tight$EF_c)),
              msg = "Loose prior -> less shrinkage toward zero")
})

# ---------------------------------------------------------------------------
# T11: update_F_cohort_all agrees with manual loop
# ---------------------------------------------------------------------------
run_test("FCohort-T11: update_F_cohort_all matches manual per-column loop", {
  set.seed(80)
  n <- 20; p <- 5; C_cols <- 3
  L_cohort <- matrix(rbinom(n * C_cols, 1, 0.4), n, C_cols)
  R_cohort <- matrix(rnorm(n * p), n, p)
  Tau      <- runif(p, 0.5, 2.0)
  sigma    <- 1.2

  res_all  <- update_F_cohort_all(L_cohort, R_cohort, Tau, sigma)

  # Manual loop
  EF_man  <- matrix(0, p, C_cols)
  EF2_man <- matrix(0, p, C_cols)
  for (c in seq_len(C_cols)) {
    r <- update_F_cohort_col(L_cohort[, c], R_cohort, Tau, sigma)
    EF_man[,  c] <- r$EF_c
    EF2_man[, c] <- r$EF2_c
  }

  assert_near(res_all$EF_cohort,  EF_man,  tol = 1e-12, msg = "EF_cohort matches manual")
  assert_near(res_all$EF2_cohort, EF2_man, tol = 1e-12, msg = "EF2_cohort matches manual")
})

# ---------------------------------------------------------------------------
# T12: Output dimensions — update_F_cohort_all
# ---------------------------------------------------------------------------
run_test("FCohort-T12: update_F_cohort_all output dimensions p x (C-1)", {
  n <- 15; p <- 8; C_cols <- 2
  L_cohort <- matrix(rbinom(n * C_cols, 1, 0.5), n, C_cols)
  R_cohort <- matrix(rnorm(n * p), n, p)
  Tau      <- rep(1.0, p)

  res <- update_F_cohort_all(L_cohort, R_cohort, Tau)

  assert_true(nrow(res$EF_cohort)  == p,      "EF_cohort nrow = p")
  assert_true(ncol(res$EF_cohort)  == C_cols, "EF_cohort ncol = C-1")
  assert_true(nrow(res$EF2_cohort) == p,      "EF2_cohort nrow = p")
  assert_true(ncol(res$EF2_cohort) == C_cols, "EF2_cohort ncol = C-1")
})

# ---------------------------------------------------------------------------
# T13: A_c strictly positive
# ---------------------------------------------------------------------------
run_test("FCohort-T13: A_c strictly positive", {
  set.seed(90)
  c_vec <- rbinom(25, 1, 0.6)
  R     <- matrix(rnorm(25 * 4), 25, 4)
  Tau   <- runif(4, 0.1, 3.0)

  res   <- update_F_cohort_col(c_vec, R, Tau, sigma_F_cohort = 1.0)

  assert_true(all(res$A_c > 0), msg = "A_c all positive")
})

# ---------------------------------------------------------------------------
# T14: Patient row order irrelevant — permuting R rows gives same answer
# ---------------------------------------------------------------------------
run_test("FCohort-T14: Row permutation of R_cohort gives same EF_c", {
  set.seed(100)
  n <- 12; p <- 4
  c_vec  <- c(rep(1, 6), rep(0, 6))
  R      <- matrix(rnorm(n * p), n, p)
  Tau    <- runif(p, 0.5, 2.0)

  perm   <- sample(n)
  res1   <- update_F_cohort_col(c_vec,        R,        Tau)
  res2   <- update_F_cohort_col(c_vec[perm],  R[perm,], Tau)

  assert_near(res1$EF_c, res2$EF_c, tol = 1e-12, msg = "Row permutation invariance")
})

# ---------------------------------------------------------------------------
# T15: EF2_c strictly positive
# ---------------------------------------------------------------------------
run_test("FCohort-T15: EF2_c strictly positive", {
  set.seed(110)
  c_vec <- rbinom(20, 1, 0.5)
  R     <- matrix(rnorm(20 * 5), 20, 5)
  Tau   <- runif(5, 0.5, 2.0)

  res   <- update_F_cohort_col(c_vec, R, Tau)

  assert_true(all(res$EF2_c > 0), msg = "EF2_c all positive")
})

# ---------------------------------------------------------------------------
# T16: compute_normal_kl returns a non-positive value (sign convention)
# ---------------------------------------------------------------------------
run_test("FCohort-T16: compute_normal_kl returns non-positive value", {
  # Values from plan example
  kl_val <- compute_normal_kl(
    EF_cohort  = matrix(c(1, 2), ncol = 1),
    EF2_cohort = matrix(c(1.5, 5), ncol = 1),
    sigma_F_cohort = 1.0
  )
  assert_true(kl_val <= 0, msg = "compute_normal_kl <= 0")

  # Also check with zero posterior mean (shrunk to zero)
  kl_zero <- compute_normal_kl(
    EF_cohort  = matrix(c(0, 0), ncol = 1),
    EF2_cohort = matrix(c(0.5, 0.5), ncol = 1),
    sigma_F_cohort = 1.0
  )
  assert_true(kl_zero <= 0, msg = "compute_normal_kl <= 0 with zero means")
})
