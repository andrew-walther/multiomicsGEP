test_that("KL equals the closed form for the normal family", {
  set.seed(3); x <- rnorm(50, 0, 2); s <- runif(50, 0.5, 1.5)
  f <- ssbmf:::.ebnm_fit(x, s, "normal")
  sig2 <- f$g$sd[1]^2
  v <- sig2 * s^2 / (sig2 + s^2); m <- sig2 * x / (sig2 + s^2)
  kl <- sum(0.5 * (v / sig2 + m^2 / sig2 - 1 - log(v / sig2)))
  expect_equal(f$kl, kl, tolerance = 1e-8)
  expect_equal(f$mean, m, tolerance = 1e-8)
})

test_that("KL equals the closed form for the point-normal family", {
  set.seed(4); x <- c(rnorm(30), rnorm(20, 3)); s <- rep(1, 50)
  f <- ssbmf:::.ebnm_fit(x, s, "point_normal")
  pi0 <- f$g$pi[1]; sig2 <- f$g$sd[2]^2
  # posterior: (1 - w) delta_0 + w N(m, v)
  lw1 <- log(1 - pi0) + dnorm(x, 0, sqrt(sig2 + s^2), log = TRUE)
  lw0 <- log(pi0) + dnorm(x, 0, s, log = TRUE)
  w <- 1 / (1 + exp(lw0 - lw1))
  v <- sig2 * s^2 / (sig2 + s^2); m <- sig2 * x / (sig2 + s^2)
  klN <- 0.5 * (v / sig2 + m^2 / sig2 - 1 - log(v / sig2))
  kl <- sum(ifelse(w < 1, (1 - w) * log((1 - w) / pi0), 0) + w * log(w / (1 - pi0)) + w * klN)
  expect_equal(f$kl, kl, tolerance = 1e-6)
})

test_that("KL is nonnegative and second moments dominate squared means for every family", {
  set.seed(5)
  for (fam in c("point_exponential", "point_laplace", "normal", "point_normal")) {
    x <- c(rexp(40), rnorm(20)); s <- runif(60, 0.3, 1)
    f <- ssbmf:::.ebnm_fit(x, s, fam)
    expect_gte(f$kl, -1e-8)
    expect_true(all(f$mean2 >= f$mean^2 - 1e-12))
    if (fam == "point_exponential") expect_true(all(f$mean >= 0))
  }
})

test_that("entries with no information get the prior as posterior", {
  f <- ssbmf:::.ebnm_fit(c(1, 2, NaN), c(1, 1, Inf), "normal")
  expect_true(all(is.finite(f$mean)))
  expect_equal(f$mean[3], 0, tolerance = 1e-6)
})
