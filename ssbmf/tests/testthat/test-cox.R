test_that("cox loglik matches coxph (Breslow, ties, strata)", {
  set.seed(1); n <- 60; x <- rnorm(n); eta <- 0.7 * x
  time <- round(rexp(n), 1) + 0.1; status <- rbinom(n, 1, 0.7); s <- rep(1:2, 30)
  ref <- survival::coxph(survival::Surv(time, status) ~ x + strata(s), ties = "breslow",
                         init = 0.7, control = survival::coxph.control(iter.max = 0))$loglik[1]
  expect_equal(ssbmf:::.cox_loglik(eta, time, status, s), ref, tolerance = 1e-8)
  ref0 <- survival::coxph(survival::Surv(time, status) ~ x, ties = "breslow",
                          init = 0.7, control = survival::coxph.control(iter.max = 0))$loglik[1]
  expect_equal(ssbmf:::.cox_loglik(eta, time, status), ref0, tolerance = 1e-8)
})

test_that("log-space risk sums are finite for |eta| up to 1e4", {
  expect_true(all(is.finite(ssbmf:::.logsumexp_risk(c(1e4, -1e4, 0, 5e3), 4:1, NULL))))
  expect_true(all(is.finite(ssbmf:::.logsumexp_risk(c(-1e4, 1e4, 0, 5e3), 4:1, NULL))))
})

test_that("moment Cox term reduces to the partial likelihood and has the stated derivatives", {
  set.seed(2); n <- 40
  time <- round(rexp(n), 1) + 0.1; status <- rbinom(n, 1, 0.6); s <- rep(1:2, 20)
  ri <- ssbmf:::.risk_index(time, s)
  eta <- rnorm(n); v <- runif(n, 0, 0.5)
  cm0 <- ssbmf:::.cox_moment(eta, numeric(n), status, ri)
  expect_equal(cm0$S, ssbmf:::.cox_loglik(eta, time, status, s), tolerance = 1e-10)
  cm <- ssbmf:::.cox_moment(eta, v, status, ri)
  h <- 1e-6
  g_eta <- vapply(1:n, function(j) { e1 <- e2 <- eta; e1[j] <- e1[j] + h; e2[j] <- e2[j] - h
    (ssbmf:::.cox_moment(e1, v, status, ri)$S - ssbmf:::.cox_moment(e2, v, status, ri)$S) / (2 * h) }, 0)
  g_v <- vapply(1:n, function(j) { v1 <- v2 <- v; v1[j] <- v1[j] + h; v2[j] <- v2[j] - h
    (ssbmf:::.cox_moment(eta, v1, status, ri)$S - ssbmf:::.cox_moment(eta, v2, status, ri)$S) / (2 * h) }, 0)
  expect_equal(cm$rho, g_eta, tolerance = 1e-6)
  expect_equal(-cm$P / 2, g_v, tolerance = 1e-6)
})
