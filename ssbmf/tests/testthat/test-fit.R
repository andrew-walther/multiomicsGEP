test_that("two-block fit is monotone and converges", {
  x <- sim_dd(n = 120, p = c(expression = 50, methylation = 50), K = 3, beta = c(0.7, 0.35, 0), seed = 1)
  fit <- ssbmf_fit(x$dd, K = 3)
  expect_equal(fit$step_decreases, 0L)
  expect_true(fit$converged)
  expect_true(all(diff(fit$objective$trace$total) >= -1e-8 * abs(fit$objective$final[["total"]])))
})

test_that("with w_surv = 0 the fit agrees with a flashier backfit from the same start", {
  skip_if_not_installed("flashier")
  x <- sim_dd(n = 100, p = c(a = 80), K = 2, seed = 5)
  ctl <- ssbmf_control(w_surv = 0, tol = 1e-10, max_sweeps = 3000)
  fit <- ssbmf_fit(x$dd, K = 2, control = ctl)
  st0 <- ssbmf:::.init_state(ssbmf:::.prep_data(x$dd, ctl), 2)
  fl <- flashier::flash_init(x$dd$Y$a, var_type = 2) |> flashier::flash_set_verbose(0) |>
    flashier::flash_add_intercept(rowwise = TRUE) |>
    flashier::flash_factors_init(init = list(st0$EL, st0$EF[[1]]),
                                 ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_laplace)) |>
    flashier::flash_backfit(maxiter = 3000, tol = 1e-10)
  ours <- sweep(fit$L$mean %*% t(fit$F$a$mean), 2, fit$mu$a, "+")
  theirs <- fitted(fl)
  rel <- sqrt(sum((ours - theirs)^2)) / sqrt(sum(scale(theirs, scale = FALSE)^2))
  # flashier treats the intercept as a factor with a prior and posterior
  # variance, ssbmf as a point estimate; the measured difference is 0.0034
  expect_lt(rel, 0.01)
})

test_that("fit is scale-stable as p grows", {
  sim <- ssbmf_simulate(n = 100, p = c(a = 5000), K = 2, feature_correlation = "block", seed = 6)
  keep <- 1:500
  d_small <- ssbmf_data(list(a = sim$X$a[keep, ]), sim$time, sim$status)
  d_big <- ssbmf_data(sim$X, sim$time, sim$status)
  ctl <- ssbmf_control(max_sweeps = 200)
  e_small <- sd(ssbmf_fit(d_small, 2, ctl)$eta); e_big <- sd(ssbmf_fit(d_big, 2, ctl)$eta)
  expect_lt(max(e_small, e_big) / min(e_small, e_big), 2)
})

test_that("matrix input and one-element list input give identical fits", {
  sim <- ssbmf_simulate(n = 60, p = c(block1 = 40), K = 2, seed = 7)
  f1 <- ssbmf_fit(ssbmf_data(sim$X$block1, sim$time, sim$status), 2, ssbmf_control(max_sweeps = 20))
  f2 <- ssbmf_fit(ssbmf_data(list(block1 = sim$X$block1), sim$time, sim$status), 2, ssbmf_control(max_sweeps = 20))
  expect_identical(f1$F, f2$F); expect_identical(f1$beta, f2$beta)
})

test_that("non-ssbmf_data input is refused", {
  expect_error(ssbmf_fit(matrix(rnorm(20), 4), K = 2), "ssbmf_data")
})
