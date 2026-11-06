# ============================================================
# Script: test_survival_block_weight.R
# Purpose: Tests for attempt S-13: the survival-block-weighted objective (untempered beta) and
#          the joint Gaussian q(beta) with per-program type-II ML prior variances.
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: ebnm, survival; code/update_beta.R, code/fit_ssbmf_yfb.R,
#               code/multimodal_yfb_updates.R
# ============================================================
# Derivation: derivations/SSBMF-YFB-model-derivations/components/survival_block_weight/2026-11-05-survival-block-weight-derivation.qmd

## Small canonical fixture ----

#' Small z-scored dataset with one prognostic program
#' @param seed RNG seed. @param n,p Dimensions. @param b True log-HR per SD of the score.
#' @return list(Y, time, status).
#' @examples
#' # sbw_fixture(1)
sbw_fixture <- function(seed, n = 80, p = 40, b = 1.2) {
  set.seed(seed)
  L <- matrix(rnorm(n * 3), n); F <- matrix(rnorm(p * 3) * rbinom(p * 3, 1, .4), p)
  Y <- scale(L %*% t(F) + matrix(rnorm(n * p), n))
  eta <- b * as.vector(scale(L[, 2]))
  t <- rexp(n, .1 * exp(eta)); cz <- rexp(n, .05)
  list(Y = Y, time = pmin(t, cz), status = as.numeric(t <= cz))
}

#' Canonical fit on the fixture
#' @param d Fixture. @param ... Passed to fit_ssbmf_yfb().
#' @return Fit.
#' @examples
#' # sbw_fit(sbw_fixture(1), w_S = 1)
sbw_fit <- function(d, ...) {
  set.seed(7)
  args <- modifyList(list(K = 3, max_iter = 25, tol = 1e-8, prior_LF = "point_laplace",
                          prior_beta = "normal", canonicalize = TRUE, verbose = FALSE), list(...))
  do.call(fit_ssbmf_yfb, c(list(d$Y, d$time, d$status), args))
}

## Joint beta update ----

run_test("S13: joint beta with K = 1 equals the mean-field Normal ebnm posterior", local({
  set.seed(3); n <- 50
  Z <- matrix(rnorm(n), n); VZ <- matrix(.05, n, 1); w <- runif(n, .2, 1)
  z <- 0.8 * Z[, 1] + rnorm(n, sd = 1.5)
  j <- update_beta_joint_ard(Z, VZ, w, z)
  mf <- update_beta_k(w, z, Z[, 1], Z[, 1]^2 + VZ[, 1], prior_family = "normal", alpha = 1)
  assert_near(j$mean, mf$mean, 1e-6)
  assert_near(diag(j$M), mf$second, 1e-6)
}))

run_test("S13: joint beta with orthogonal scores equals K separate mean-field updates", local({
  set.seed(4); n <- 60
  Q <- qr.Q(qr(matrix(rnorm(n * 3), n))) * sqrt(n)
  w <- rep(1, n); VZ <- matrix(0, n, 3)
  z <- Q %*% c(.6, 0, -.3) + rnorm(n)
  j <- update_beta_joint_ard(Q, VZ, w, as.vector(z))
  for (k in 1:3) {
    zk <- as.vector(z - Q[, -k] %*% j$mean[-k])
    mf <- update_beta_k(w, zk, Q[, k], Q[, k]^2, prior_family = "normal", alpha = 1)
    assert_near(j$mean[k], mf$mean, 1e-5, sprintf("program %d", k))
  }
}))

run_test("S13: collinear scores widen the joint marginal SD relative to mean-field", local({
  set.seed(5); n <- 100
  a <- rnorm(n); Z <- cbind(a, a + rnorm(n, sd = .3)); VZ <- matrix(0, n, 2); w <- rep(1, n)
  z <- 1 * a + rnorm(n)
  j <- update_beta_joint_ard(Z, VZ, w, z)
  inc <- j$sigma2 > 0
  assert_true(any(inc), "at least one program included")
  # mean-field conditional SD for an included program: 1 / sqrt(G_kk + 1/sigma2_k)
  G <- crossprod(Z * w, Z)
  k <- which(inc)[1]
  sd_mf <- 1 / sqrt(G[k, k] + 1 / j$sigma2[k])
  if (all(inc)) assert_true(sqrt(j$Sigma[k, k]) > sd_mf * 1.5, "joint SD inflated by collinearity")
  else assert_near(sqrt(j$Sigma[k, k]), sd_mf, 1e-8)
}))

run_test("S13: a pure-noise program gets prior variance 0 and beta exactly 0", local({
  set.seed(6); n <- 200
  Z <- cbind(rnorm(n), rnorm(n)); VZ <- matrix(0, n, 2); w <- rep(1, n)
  z <- 1.5 * Z[, 1] + rnorm(n)
  # remove any chance correlation of z with program 2 so its quality factor is 0
  Z[, 2] <- residuals(lm(Z[, 2] ~ Z[, 1] + z))
  j <- update_beta_joint_ard(Z, VZ, w, z)
  assert_equal(j$sigma2[2], 0); assert_equal(j$mean[2], 0); assert_equal(j$M[2, 2], 0)
  assert_true(j$sigma2[1] > 0 && abs(j$mean[1] - 1.5) < .3)
}))

run_test("S13: the joint beta step does not decrease its block objective", local({
  set.seed(8); n <- 80
  Z <- matrix(rnorm(n * 4), n) %*% matrix(c(1, .5, 0, 0, 0, 1, .4, 0, 0, 0, 1, .3, 0, 0, 0, 1), 4)
  VZ <- matrix(.02, n, 4); w <- runif(n, .3, 1); z <- Z %*% c(.5, 0, -.4, 0) + rnorm(n)
  j <- update_beta_joint_ard(Z, VZ, w, as.vector(z))
  # objective of the block at (m, Sigma, sigma2): E working loglik - KL
  obj <- function(m, S, s2) {
    G <- crossprod(Z * w, Z) + diag(colSums(w * VZ)); b <- as.vector(crossprod(Z, w * z))
    sum(b * m) - .5 * sum(G * (tcrossprod(m) + S)) - beta_joint_kl(m, S, s2)
  }
  o_star <- obj(j$mean, j$Sigma, j$sigma2)
  for (s2 in list(rep(.1, 4), rep(1, 4), j$sigma2 * 2 + .01)) {
    G <- crossprod(Z * w, Z) + diag(colSums(w * VZ)); b <- as.vector(crossprod(Z, w * z))
    S <- solve(G + diag(1 / s2)); m <- as.vector(S %*% b)
    assert_true(obj(m, S, s2) <= o_star + 1e-8)
  }
}))

## F update with joint beta moments ----

run_test("S13: F update with mean-field moments passed as beta-weighted signal is unchanged", local({
  set.seed(9); d <- sbw_fixture(9, n = 40, p = 20)
  Y <- d$Y; n <- nrow(Y); p <- ncol(Y)
  EL <- rnorm(n); EL2 <- EL^2 + .1; EF <- rnorm(p) * .2; EF2 <- EF^2 + .01
  R <- Y - outer(EL, EF) + outer(EL, EF); w <- runif(n, .2, 1); zno <- rnorm(n)
  Z <- as.vector(Y %*% EF); VZ <- as.vector(Y^2 %*% (EF2 - EF^2))
  prior <- multimodal_yfb_fit_signed_prior(rnorm(p), rep(1, p), "point_laplace")
  Tau <- rep(1.3, p); eb <- .7; eb2 <- .6
  a <- multimodal_yfb_update_F_mk(Y, Tau, EL, EL2, R, w, w * zno, Z, VZ, EF, EF2, eb, eb2, prior,
                                  prior_update = "ebnm")
  b <- multimodal_yfb_update_F_mk(Y, Tau, EL, EL2, R, w, eb * w * zno, Z, VZ, EF, EF2, eb, eb2, prior,
                                  prior_update = "ebnm", h_beta_weighted = TRUE)
  assert_near(a$EF, b$EF, 1e-12); assert_near(a$EF2, b$EF2, 1e-12)
}))

## Fitter: survival-block weight ----

run_test("S13: survival_block_weight at w_S = 1 reproduces the default fit", local({
  d <- sbw_fixture(11)
  a <- sbw_fit(d, w_S = 1); b <- sbw_fit(d, w_S = 1, survival_block_weight = TRUE)
  assert_near(a$EF, b$EF, 1e-10); assert_near(a$EBeta, b$EBeta, 1e-10)
  assert_near(a$history$elbo_full, b$history$elbo_full, 1e-8)
}))

run_test("S13: the survival-block weight leaves beta untempered (larger variance than A0)", local({
  d <- sbw_fixture(12)
  u <- sbw_fit(d, w_S = 20, survival_block_weight = TRUE)
  t0 <- sbw_fit(d, w_S = 20)
  vu <- u$EBeta2 - u$EBeta^2; vt <- t0$EBeta2 - t0$EBeta^2
  live <- abs(u$EBeta) > 1e-3 & abs(t0$EBeta) > 1e-3   # programs neither fit shrinks to 0
  assert_true(any(live), "some program live in both fits")
  # A0 tempers the beta likelihood by 20, so its posterior variance is about 1/20 as large
  assert_true(median(vu[live] / vt[live]) > 5,
              sprintf("variance ratio %.2f", median(vu[live] / vt[live])))
  # and the untempered fit calls no more programs survival-active than the tempered one
  act <- function(f) sum(abs(f$EBeta) >= 1.96 * sqrt(pmax(f$EBeta2 - f$EBeta^2, 0)) &
                         f$EBeta2 - f$EBeta^2 > 0)
  assert_true(act(u) <= act(t0), sprintf("active %d vs %d", act(u), act(t0)))
}))

run_test("S13: the monitored objective rises under the survival-block weight (w_S = 10)", local({
  d <- sbw_fixture(13)
  f <- sbw_fit(d, w_S = 10, survival_block_weight = TRUE)
  e <- f$history$elbo_full[seq_len(f$history$n_iter)]
  dec <- diff(e) < -1e-6 * abs(e[-1])
  assert_true(sum(dec) <= 1, sprintf("%d decreases", sum(dec)))
}))

run_test("S13: joint beta fit runs, rises, and stores a valid covariance", local({
  d <- sbw_fixture(14)
  f <- sbw_fit(d, w_S = 5, survival_block_weight = TRUE, beta_posterior = "joint")
  e <- f$history$elbo_full[seq_len(f$history$n_iter)]
  assert_true(sum(diff(e) < -1e-6 * abs(e[-1])) <= 1)
  S <- f$EBeta_cov
  assert_true(is.matrix(S) && all(dim(S) == 3) && isSymmetric(S, tol = 1e-8))
  assert_true(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
  assert_near(f$EBeta2, f$EBeta^2 + diag(S), 1e-10)
}))

run_test("S13: joint beta with diagonal-like scores is close to the mean-field fit", local({
  d <- sbw_fixture(15)
  a <- sbw_fit(d, w_S = 1, survival_block_weight = TRUE)
  b <- sbw_fit(d, w_S = 1, survival_block_weight = TRUE, beta_posterior = "joint")
  # risk scores agree closely (programs differ only through beta's posterior correlation)
  ra <- as.vector(d$Y %*% a$EF %*% a$EBeta); rb <- as.vector(d$Y %*% b$EF %*% b$EBeta)
  assert_true(cor(ra, rb) > .95, sprintf("risk correlation %.3f", cor(ra, rb)))
}))

run_test("S13: invalid option combinations fail loudly", local({
  d <- sbw_fixture(16)
  assert_true(inherits(try(fit_ssbmf_yfb(d$Y, d$time, d$status, K = 2, max_iter = 2,
    canonicalize = FALSE, survival_block_weight = TRUE, verbose = FALSE), silent = TRUE), "try-error"))
  assert_true(inherits(try(sbw_fit(d, w_S = 1, beta_posterior = "joint", prior_beta = "point_normal"),
    silent = TRUE), "try-error"))
}))
