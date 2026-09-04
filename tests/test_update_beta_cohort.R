# ============================================================
# Script:  test_update_beta_cohort.R
# Purpose: Tests for code/update_beta_cohort.R -- cohort-specific survival
#          coefficients beta_k^(c) (update_beta_cohort_k/_all,
#          compute_z_no_k_cohort, compute_pooled_beta).
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_update_beta_cohort.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages(library(ebnm))
if (!exists("update_beta_k")) source("code/update_beta.R")
if (!exists("update_beta_cohort_k")) source("code/update_beta_cohort.R")

cat("\n--- test_update_beta_cohort.R ---\n")

set.seed(4041)
n <- 120; K <- 4

# ============================================================
# compute_z_no_k_cohort() ----
# ============================================================

run_test("BC-T1: compute_z_no_k_cohort with C=1 matches compute_z_no_k", {
  ZF    <- matrix(rnorm(n * K), n, K)
  EBeta <- rnorm(K)
  z     <- rnorm(n)
  cohort_idx <- rep(1L, n)
  EBeta_mat  <- matrix(EBeta, K, 1)

  for (k in 1:K) {
    a <- compute_z_no_k_cohort(z, ZF, EBeta_mat, cohort_idx, k)
    b <- compute_z_no_k(z, ZF, EBeta, k)
    assert_near(a, b, tol = 1e-10, sprintf("mismatch at k=%d", k))
  }
})

run_test("BC-T2: compute_z_no_k_cohort correctly looks up each patient's own cohort column", {
  ZF <- matrix(1, 4, 2)  # 4 patients, 2 factors, all ZF=1 for simplicity
  cohort_idx <- c(1L, 2L, 1L, 2L)
  # rows = factor (1,2), cols = cohort (1,2): factor1 = c(10,100), factor2 = c(20,200)
  EBeta_mat  <- rbind(c(10, 100), c(20, 200))
  z <- rep(0, 4)
  z_no_1 <- compute_z_no_k_cohort(z, ZF, EBeta_mat, cohort_idx, k = 1)
  # eta_full[i] = EBeta_mat[1,c(i)] + EBeta_mat[2,c(i)]; eta_no_1 = eta_full - EBeta_mat[1,c(i)] = EBeta_mat[2,c(i)]
  expected_eta_no_1 <- c(20, 200, 20, 200)
  assert_near(z_no_1, -expected_eta_no_1, tol = 1e-10, "z_no_k_cohort did not use each patient's own cohort column")
})

# ============================================================
# update_beta_cohort_k() ----
# ============================================================

run_test("BC-T3: update_beta_cohort_k with C=1 matches update_beta_k exactly", {
  w      <- runif(n, 0.5, 2)
  ZF_k   <- rnorm(n)
  z_no_k <- ZF_k * 1.3 + rnorm(n, sd = 0.4)
  cohort_idx <- rep(1L, n)

  ref <- update_beta_k(w, z_no_k, ZF_k, ZF_k^2, prior_family = "normal", alpha = 0.5)
  res <- update_beta_cohort_k(w, z_no_k, ZF_k, ZF_k^2, cohort_idx, C = 1,
                               prior_family = "normal", alpha = 0.5)

  assert_near(res$mean, ref$mean, tol = 1e-8, "mean mismatch at C=1")
  assert_near(res$A, ref$A, tol = 1e-8, "A mismatch at C=1")
  assert_near(res$B, ref$B, tol = 1e-8, "B mismatch at C=1")
})

run_test("BC-T4: update_beta_cohort_k recovers genuinely different signs across cohorts", {
  set.seed(55)
  n_c <- 300
  cohort_idx <- c(rep(1L, n_c), rep(2L, n_c))
  w <- rep(1, 2 * n_c)
  ZF_k <- rnorm(2 * n_c)
  # cohort 1: positive signal; cohort 2: negative signal (strong, low noise)
  true_beta <- c(2.0, -2.0)
  z_no_k <- ZF_k * true_beta[cohort_idx] + rnorm(2 * n_c, sd = 0.1)

  res <- update_beta_cohort_k(w, z_no_k, ZF_k, ZF_k^2, cohort_idx, C = 2,
                               prior_family = "normal", alpha = 1.0)
  assert_true(res$mean[1] > 1.0, "cohort 1 should recover a strong positive coefficient")
  assert_true(res$mean[2] < -1.0, "cohort 2 should recover a strong negative coefficient")
})

run_test("BC-T5: partial pooling -- identical true signal across cohorts shrinks estimates toward each other", {
  set.seed(56)
  n_c <- 15  # small per-cohort n -- noisy independent OLS-style estimates would disagree more
  cohort_idx <- c(rep(1L, n_c), rep(2L, n_c), rep(3L, n_c))
  C <- 3
  w <- rep(1, n_c * C)
  ZF_k <- rnorm(n_c * C)
  true_beta_shared <- 1.5
  z_no_k <- ZF_k * true_beta_shared + rnorm(n_c * C, sd = 1.5)  # noisy, small n per cohort

  res_pooled <- update_beta_cohort_k(w, z_no_k, ZF_k, ZF_k^2, cohort_idx, C = C,
                                      prior_family = "normal", alpha = 1.0)

  # Independent per-cohort OLS-style estimate (equivalent to update_beta_k on each cohort alone)
  indep <- numeric(C)
  for (c in seq_len(C)) {
    idx <- which(cohort_idx == c)
    r <- update_beta_k(w[idx], z_no_k[idx], ZF_k[idx], ZF_k[idx]^2,
                        prior_family = "normal", alpha = 1.0)
    indep[c] <- r$mean
  }

  spread_pooled <- max(res_pooled$mean) - min(res_pooled$mean)
  spread_indep  <- max(indep) - min(indep)
  assert_true(spread_pooled <= spread_indep,
              "pooled (shared-prior) cohort estimates should not spread out MORE than independent per-cohort estimates")
})

# ============================================================
# update_beta_cohort_all() ----
# ============================================================

run_test("BC-T6: update_beta_cohort_all returns K x C matrices", {
  C <- 3
  cohort_idx <- sample(1:C, n, replace = TRUE)
  ZF <- matrix(rnorm(n * K), n, K)
  w  <- runif(n, 0.5, 2)
  z  <- rnorm(n)
  EBeta_init <- matrix(0, K, C)

  res <- update_beta_cohort_all(w, z, ZF, EBeta_init, cohort_idx, C,
                                 prior_family = "normal", alpha = 0.5)
  assert_equal(dim(res$EBeta), as.integer(c(K, C)), "EBeta shape mismatch")
  assert_equal(dim(res$EBeta2), as.integer(c(K, C)), "EBeta2 shape mismatch")
  assert_finite(res$EBeta, "EBeta must be finite")
  assert_true(all(res$EBeta2 >= res$EBeta^2 - 1e-8), "EBeta2 must be >= EBeta^2 (valid second moment)")
})

run_test("BC-T7: update_beta_cohort_all with C=1 matches update_beta_all exactly", {
  cohort_idx <- rep(1L, n)
  ZF <- matrix(rnorm(n * K), n, K)
  w  <- runif(n, 0.5, 2)
  z  <- rnorm(n)
  EBeta_init_vec <- rep(0, K)
  EBeta_init_mat <- matrix(EBeta_init_vec, K, 1)

  ref <- update_beta_all(w, z, ZF, ZF^2, EBeta_init_vec, prior_family = "normal", alpha = 0.5)
  res <- update_beta_cohort_all(w, z, ZF, EBeta_init_mat, cohort_idx, C = 1,
                                 prior_family = "normal", alpha = 0.5)

  assert_near(as.vector(res$EBeta), ref$EBeta, tol = 1e-8, "EBeta mismatch at C=1")
  assert_near(as.vector(res$EBeta2), ref$EBeta2, tol = 1e-8, "EBeta2 mismatch at C=1")
})

# ============================================================
# compute_pooled_beta() ----
# ============================================================

run_test("BC-T8: compute_pooled_beta equals update_beta_all applied to the same data", {
  ZF <- matrix(rnorm(n * K), n, K)
  w  <- runif(n, 0.5, 2)
  z  <- rnorm(n)
  init <- rep(0, K)

  pooled <- compute_pooled_beta(w, z, ZF, init, prior_family = "normal", alpha = 0.5)
  ref    <- update_beta_all(w, z, ZF, ZF^2, init, prior_family = "normal", alpha = 0.5)$EBeta

  assert_near(pooled, ref, tol = 1e-10, "compute_pooled_beta should exactly reproduce update_beta_all's EBeta")
})

run_test("BC-T9: compute_pooled_beta is NOT the same as rowMeans when cohort sizes are unequal", {
  set.seed(57)
  n1 <- 200; n2 <- 40  # deliberately unequal, like this project's own 144 vs 129 (exaggerated here)
  cohort_idx <- c(rep(1L, n1), rep(2L, n2))
  C <- 2
  ZF_k <- rnorm(n1 + n2)
  w <- rep(1, n1 + n2)
  true_beta <- c(1.0, 3.0)  # different magnitude per cohort
  z_no_k <- ZF_k * true_beta[cohort_idx] + rnorm(n1 + n2, sd = 0.3)
  ZF <- matrix(ZF_k, ncol = 1)
  z  <- z_no_k  # single factor, so z_no_k(k=1) == z itself in this construction

  fit <- update_beta_cohort_all(w, z, ZF, matrix(0, 1, C), cohort_idx, C,
                                 prior_family = "normal", alpha = 1.0)
  naive_rowmean <- mean(fit$EBeta[1, ])
  pooled        <- compute_pooled_beta(w, z, ZF, 0, prior_family = "normal", alpha = 1.0)

  # rowMeans would weight both cohorts equally (50/50); the pooled-patient-sum
  # computation weights by patient count (n1 >> n2), so it should sit closer
  # to the larger cohort's own coefficient than the naive rowMean does.
  dist_pooled_to_c1   <- abs(pooled[1] - fit$EBeta[1, 1])
  dist_rowmean_to_c1  <- abs(naive_rowmean - fit$EBeta[1, 1])
  assert_true(dist_pooled_to_c1 < dist_rowmean_to_c1,
              "pooled-patient-sum beta should lie closer to the larger cohort's estimate than a naive rowMean does")
})

report_results("test_update_beta_cohort.R")
