test_that("ssbmf_cv() returns per-K summaries, applies the one-SE rule and refits", {
  x <- sim_dd(n = 100, p = c(a = 40), K = 2, beta = c(1, 0), seed = 12)
  cv <- ssbmf_cv(x$dd, K_grid = c(1, 2, 4), nfolds = 3, control = ssbmf_control(max_sweeps = 100))
  expect_identical(cv$cv$K, c(1, 2, 4))
  expect_true(all(cv$fold_C$C >= 0 & cv$fold_C$C <= 1))
  i <- which.max(cv$cv$mean_C)
  ok <- cv$cv$K[cv$cv$mean_C >= cv$cv$mean_C[i] - cv$cv$se_C[i]]
  expect_identical(cv$best, min(ok))
  expect_s3_class(cv$fit, "ssbmf_fit")
  expect_identical(length(cv$fit$beta$mean), as.integer(cv$best))
  # folds are balanced on events
  expect_lt(diff(range(tapply(x$dd$status, cv$folds, sum))), 2)
})

test_that("held-out C uses strata when given", {
  set.seed(1); r <- rnorm(40); t <- rexp(40); s <- rbinom(40, 1, 0.7); g <- rep(1:2, 20)
  ref <- survival::concordance(survival::Surv(t, s) ~ r + strata(g), reverse = TRUE)$concordance
  expect_equal(ssbmf:::.cindex(r, t, s, factor(g)), unname(ref))
})
