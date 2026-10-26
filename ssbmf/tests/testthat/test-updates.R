test_that("each update alone does not lower the objective (20 random states)", {
  worst <- 0
  for (seed in 1:20) {
    x <- sim_dd(n = 40, p = c(a = 30, b = 60), K = 3, seed = seed)
    s <- state_after(x$dd, K = 3, sweeps = seed %% 3)
    st <- s$st; d <- s$d
    ups <- list(function(st) ssbmf:::.update_mu(st, d), function(st) ssbmf:::.update_tau(st, d),
                function(st) ssbmf:::.update_L(st, d, 2), function(st) ssbmf:::.update_F(st, d, 1, 1),
                function(st) ssbmf:::.update_F(st, d, 2, 3), function(st) ssbmf:::.update_beta(st, d))
    for (u in ups) {
      new <- u(st)
      # recompute from scratch to make sure the cached objective is exact
      tot <- ssbmf:::.objective(new, d)$obj$total
      expect_equal(tot, new$obj$total, tolerance = 1e-10)
      worst <- min(worst, (tot - st$obj$total) / abs(st$obj$total))
    }
  }
  expect_gte(worst, -1e-8)
})

test_that("updates are stationary: finite-difference gradients vanish", {
  res <- ssbmf:::.verify_stationarity(seed = 2, verbose = FALSE)
  expect_lt(max(res$max_grad), 1e-5)
})

test_that("with w_surv = 0 the F update is the unsupervised EBMF update", {
  x <- sim_dd(n = 40, p = c(a = 30), K = 2, seed = 3)
  s <- state_after(x$dd, K = 2, sweeps = 2, control = ssbmf_control(w_surv = 0))
  st <- s$st; d <- s$d
  new <- ssbmf:::.update_F(st, d, 1, 1)
  tau <- st$tau[[1]]
  aR <- tau * sum(st$EL2[, 1])
  hR <- tau * (drop(crossprod(st$R[[1]], st$EL[, 1])) + st$EF[[1]][, 1] * sum(st$EL[, 1]^2))
  ref <- ebnm::ebnm(hR / aR, 1 / sqrt(aR), prior_family = "point_laplace", g_init = st$gF[[1]][[1]])
  expect_equal(unname(new$EF[[1]][, 1]), ref$posterior$mean, tolerance = 1e-6)
  expect_identical(new$b, numeric(2))
})

test_that("beta update leaves the zero state", {
  x <- sim_dd(n = 150, p = c(a = 60), K = 2, beta = c(1.5, 0), seed = 4)
  s <- state_after(x$dd, K = 2, sweeps = 0)
  st <- s$st; d <- s$d
  expect_identical(st$b, numeric(2)); expect_identical(st$s2, numeric(2))
  for (i in 1:5) {
    st <- ssbmf:::.update_mu(st, d); st <- ssbmf:::.update_tau(st, d)
    for (k in 1:2) st <- ssbmf:::.update_L(st, d, k)
    for (k in 1:2) st <- ssbmf:::.update_F(st, d, 1, k)
    st <- ssbmf:::.update_beta(st, d)
  }
  # canonical scale: per SD of the program score; the planted program is the one matching T1
  zsd <- apply(st$Zbar, 2, sd)
  k1 <- which.max(abs(cor(st$EF[[1]], x$sim$truth$F$a[, 1])))
  expect_gt(abs(st$b[k1] * zsd[k1]), 0.5)
})
