test_that("weight_rule = 'none' reproduces the default fit exactly", {
  x <- sim_dd(n = 80, p = c(a = 40), K = 2, seed = 21)
  f0 <- ssbmf_fit(x$dd, 2, ssbmf_control(max_sweeps = 30))
  f1 <- ssbmf_fit(x$dd, 2, ssbmf_control(max_sweeps = 30, weight_rule = "none"))
  expect_identical(f0$F, f1$F); expect_identical(f0$beta, f1$beta)
})

test_that("effective-feature weight is near 1 for independent and below 1 for block-correlated features", {
  set.seed(1)
  R_ind <- matrix(rnorm(150 * 1000), 150, 1000)
  expect_gt(ssbmf:::.effective_features(scale(R_ind, scale = FALSE), 149) / 1000, 0.9)
  G <- matrix(rnorm(150 * 50), 150, 50)
  R_blk <- sqrt(0.5) * matrix(rnorm(150 * 1000), 150, 1000) + sqrt(0.5) * G[, rep(1:50, each = 20)]
  expect_lt(ssbmf:::.effective_features(scale(R_blk, scale = FALSE), 149) / 1000, 0.3)
})

test_that("the effective-feature rule keeps the fit monotone", {
  x <- sim_dd(n = 120, p = c(a = 300), K = 3, seed = 22, feature_correlation = "block")
  fit <- ssbmf_fit(x$dd, 3, ssbmf_control(weight_rule = "effective_features"))
  expect_equal(fit$step_decreases, 0L)
  expect_lt(fit$weights$block[["a"]], 1)
})

test_that("ssbmf_cv() can cross-validate the survival weight", {
  x <- sim_dd(n = 80, p = c(a = 30), K = 2, beta = c(1, 0), seed = 23)
  cv <- ssbmf_cv(x$dd, K_grid = 2, nfolds = 3, w_grid = c(1, 100), control = ssbmf_control(max_sweeps = 50))
  expect_identical(nrow(cv$cv), 2L)
  expect_true(cv$best_w_surv %in% c(1, 100))
  expect_identical(cv$fit$weights$surv, cv$best_w_surv)
})
