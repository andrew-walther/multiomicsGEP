test_that("dimensions and orientation", {
  sim <- ssbmf_simulate(n = 50, p = c(e = 40, m = 70), K = 3, n_validation = 20, seed = 1)
  expect_identical(dim(sim$X$e), c(40L, 50L)); expect_identical(dim(sim$X$m), c(70L, 50L))
  expect_identical(dim(sim$validation$X$m), c(70L, 20L))
  expect_identical(colnames(sim$X$e), colnames(sim$X$m))
  expect_length(sim$time, 50)
})

test_that("censoring is about 0.30", {
  for (seed in 1:3) {
    sim <- ssbmf_simulate(n = 200, p = c(e = 40), K = 2, seed = seed)
    expect_lt(abs(1 - mean(sim$status) - 0.30), 0.05)
  }
})

test_that("the low-variance prognostic program explains under 2% of variance", {
  sim <- ssbmf_simulate(n = 200, p = c(e = 300, m = 400), K = 3, low_variance_prognostic = TRUE, seed = 2)
  expect_true(all(sim$truth$pve[3, ] < 0.02))
  expect_identical(sim$truth$beta, c(0, 0, 1))
})

test_that("same seed gives identical output", {
  a <- ssbmf_simulate(n = 30, p = c(e = 20), K = 2, seed = 9)
  b <- ssbmf_simulate(n = 30, p = c(e = 20), K = 2, seed = 9)
  expect_identical(a, b)
})
