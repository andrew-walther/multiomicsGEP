x <- sim_dd(n = 150, p = c(expression = 60, methylation = 80), K = 3, beta = c(1.2, 0, 0), seed = 11)
fit <- ssbmf_fit(x$dd, K = 3)
val <- x$sim$validation

test_that("prediction aligns new data by feature name", {
  r0 <- predict(fit, newdata = val$X)
  perm <- lapply(val$X, function(A) A[sample(nrow(A)), ])
  expect_equal(predict(fit, newdata = perm), r0)
  extra <- lapply(val$X, function(A) rbind(A, extra_feature = rnorm(ncol(A))))
  expect_equal(predict(fit, newdata = extra), r0)
  miss <- val$X; miss$expression <- miss$expression[-5, ]
  expect_error(predict(fit, newdata = miss), rownames(val$X$expression)[5])
  tr <- lapply(val$X, t)
  expect_equal(predict(fit, newdata = tr, orientation = "samples_by_features"), r0)
})

test_that("predict() without newdata reproduces the fitted predictor", {
  expect_equal(predict(fit), fit$eta, tolerance = 1e-10)
  expect_equal(predict(fit, newdata = x$dd), fit$eta, tolerance = 1e-10)
})

test_that("risk has the frozen direction: higher risk, shorter survival", {
  r <- predict(fit, newdata = val$X, type = "risk")
  cc <- survival::concordance(survival::Surv(val$time, val$status) ~ r, reverse = TRUE)$concordance
  expect_gt(cc, 0.6)
})

test_that("coef() has one row per program with ordered intervals", {
  cf <- coef(fit)
  expect_identical(nrow(cf), 3L)
  expect_true(all(cf$lower <= cf$mean & cf$mean <= cf$upper))
  expect_true(all(cf$lfsr >= 0 & cf$lfsr <= 1))
})

test_that("classify labels exactly the planted prognostic program survival-active", {
  roles <- ssbmf_classify(fit)
  Fh <- rbind(fit$F$expression$mean, fit$F$methylation$mean)
  Ft <- rbind(x$sim$truth$F$expression, x$sim$truth$F$methylation)
  k1 <- which.max(abs(cor(Fh, Ft[, 1])))
  expect_identical(roles$role[k1], "survival_active")
  expect_identical(sum(roles$role == "survival_active"), 1L)
})

test_that("canonical scaling reproduces predict() and the factorization", {
  Fc <- ssbmf_loadings(fit); Lc <- ssbmf_activities(fit); cf <- coef(fit)
  eta <- Reduce(`+`, lapply(names(Fc), function(m) sweep(x$dd$Y[[m]], 2, fit$center[[m]]) %*% Fc[[m]])) %*% cf$mean
  expect_equal(unname(drop(eta)), unname(fit$eta), tolerance = 1e-10)
  expect_equal(Lc %*% t(Fc$expression), fit$L$mean %*% t(fit$F$expression$mean), tolerance = 1e-10)
  sc <- predict(fit, type = "scores")
  expect_equal(unname(apply(sc %*% diag(1 / ssbmf:::.canon_scale(fit)), 2, sd)[fit$score_sd > 0]),
               rep(1, sum(fit$score_sd > 0)), tolerance = 1e-10)
})

test_that("top features and summary print", {
  tf <- ssbmf_top_features(fit, n = 5)
  expect_identical(nrow(tf), 3L * 5L * 2L)
  expect_output(print(fit), "ssbmf_fit")
  expect_output(summary(fit), "Coefficients")
})
