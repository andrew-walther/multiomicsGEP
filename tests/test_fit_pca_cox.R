# =============================================================================
# tests/test_fit_pca_cox.R
#
# TDD test suite for fit_pca_cox() / predict_pca_cox()
# (results/multi_cohort_sim/fit_pca_cox.R) -- the PCA+Cox two-step baseline
# used by Phase 2's joint-vs-2-step value-add simulation
# (docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md).
#
# NOTE: requires the survival package.
# =============================================================================

suppressPackageStartupMessages(library(survival))

cat("=== T1: fit_pca_cox() structure ===\n")

run_test("T1.1: EF has dimensions p x K", {
  set.seed(1); n <- 60; p <- 40; K <- 3
  Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1, 0.7)

  fit <- fit_pca_cox(Y, time, status, K = K)
  assert_true(nrow(fit$EF) == p && ncol(fit$EF) == K,
              sprintf("EF dims: expected %dx%d, got %dx%d", p, K, nrow(fit$EF), ncol(fit$EF)))
})

run_test("T1.2: cox_coef has length K", {
  set.seed(2); n <- 60; p <- 40; K <- 3
  Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1, 0.7)

  fit <- fit_pca_cox(Y, time, status, K = K)
  assert_length(fit$cox_coef, K, "cox_coef should be length K")
  assert_finite(fit$cox_coef, "cox_coef should be finite")
})

run_test("T1.3: center is a p-vector matching colMeans(Y)", {
  set.seed(3); n <- 50; p <- 20; K <- 2
  Y <- matrix(rnorm(n * p, mean = 5), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1, 0.7)

  fit <- fit_pca_cox(Y, time, status, K = K)
  assert_near(fit$center, colMeans(Y), tol = 1e-8, msg = "center should equal training colMeans")
})

cat("\n=== T2: predict_pca_cox() consistency ===\n")

run_test("T2.1: predicting on the training data reproduces prcomp scores exactly", {
  set.seed(4); n <- 50; p <- 30; K <- 3
  Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1, 0.7)

  fit  <- fit_pca_cox(Y, time, status, K = K)
  pred <- predict_pca_cox(fit, Y)

  pca_ref <- prcomp(Y, rank. = K, center = TRUE, scale. = FALSE)
  assert_near(pred$scores, pca_ref$x, tol = 1e-8,
              msg = "predict_pca_cox on training data should reproduce prcomp$x exactly")
})

run_test("T2.2: risk_scores = scores %*% cox_coef exactly", {
  set.seed(5); n <- 40; p <- 25; K <- 2
  Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1, 0.7)

  fit  <- fit_pca_cox(Y, time, status, K = K)
  pred <- predict_pca_cox(fit, Y)

  assert_near(pred$risk_scores, as.vector(pred$scores %*% fit$cox_coef), tol = 1e-10,
              msg = "risk_scores should equal scores %*% cox_coef")
})

run_test("T2.3: held-out projection uses the TRAINING center, not the test set's own mean", {
  set.seed(6); n_tr <- 50; n_te <- 20; p <- 20; K <- 2
  Y_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  time_tr <- rexp(n_tr, 0.1); status_tr <- rbinom(n_tr, 1, 0.7)
  Y_te <- matrix(rnorm(n_te * p, mean = 3), n_te, p)  # deliberately shifted mean

  fit  <- fit_pca_cox(Y_tr, time_tr, status_tr, K = K)
  pred <- predict_pca_cox(fit, Y_te)

  expected_scores <- sweep(Y_te, 2, fit$center, "-") %*% fit$EF
  assert_near(pred$scores, expected_scores, tol = 1e-8,
              msg = "held-out scores should use the training center, not recentre on test data")
})

cat("\n=== T3: Signal recovery ===\n")

run_test("T3.1: recovers a known prognostic direction with positive concordance", {
  set.seed(7); n <- 200; p <- 100; K_true <- 2
  L_true <- matrix(rnorm(n * K_true), n, K_true)
  F_true <- matrix(0, p, K_true)
  F_true[1:20, 1] <- rnorm(20, 0, 3)   # first 20 genes load on factor 1 (the prognostic one)
  F_true[21:40, 2] <- rnorm(20, 0, 3)  # next 20 genes load on factor 2 (irrelevant)
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)

  beta_true <- c(1.5, 0)  # only factor 1 is prognostic
  eta <- as.vector(L_true %*% beta_true)
  time   <- rweibull(n, shape = 1.5, scale = exp(-eta / 1.5))
  status <- as.integer(time < stats::quantile(time, 0.7))

  fit  <- fit_pca_cox(Y, time, status, K = 2)
  pred <- predict_pca_cox(fit, Y)

  c_idx <- as.numeric(survival::concordance(survival::Surv(time, status) ~ pred$risk_scores)$concordance)
  assert_true(max(c_idx, 1 - c_idx) > 0.65,
              sprintf("PCA+Cox should recover a real prognostic signal; oriented C=%.3f", max(c_idx, 1 - c_idx)))
})

cat("\n=== T4: Edge cases ===\n")

run_test("T4.1: K=1 works", {
  set.seed(8); n <- 30; p <- 15
  Y <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1); status <- rbinom(n, 1, 0.7)

  fit <- fit_pca_cox(Y, time, status, K = 1)
  assert_true(ncol(fit$EF) == 1, "K=1: EF should have 1 column")
  assert_length(fit$cox_coef, 1, "K=1: cox_coef should be length 1")
})

run_test("T4.2: rank-deficient Cox fit (perfectly separated / degenerate) returns finite zero coefficients, not an error", {
  set.seed(9); n <- 20; p <- 10; K <- 3
  Y <- matrix(rnorm(n * p), n, p)
  # Degenerate: all events at the same time, no variation for Cox to fit
  time   <- rep(1, n)
  status <- rep(1, n)

  fit <- tryCatch(fit_pca_cox(Y, time, status, K = K), error = function(e) NULL)
  assert_true(!is.null(fit), "fit_pca_cox should not error on a degenerate Cox fit")
  assert_finite(fit$cox_coef, "cox_coef should be finite (zeroed) rather than NA on a degenerate fit")
})
