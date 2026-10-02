# ============================================================
# Script: test_multimodal_yfb_pruning.R
# Purpose: Tests for the multimodal YFB parsimony framework: coordinate KL,
#          residual-precision models, ELBO nullcheck pruning, and intercept.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: ebnm; code/multimodal_yfb_updates.R, code/fit_multimodal_yfb.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("multimodal_yfb_coordinate_kl")) source("code/multimodal_yfb_updates.R")
if (!exists("fit_multimodal_yfb")) source("code/fit_multimodal_yfb.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")

run_test("MMYFB-Prune-T1: coordinate KL identity matches the closed-form Normal KL", {
  # Normal prior N(0, v0), pseudo-observation x with variance s2:
  # posterior N(m, v) with v = 1/(1/s2 + 1/v0), m = v x / s2
  x <- c(-2, 0.3, 1.7); s2 <- c(0.5, 1, 2); v0 <- 1.5
  v <- 1 / (1 / s2 + 1 / v0); m <- v * x / s2
  log_ml <- dnorm(x, 0, sqrt(s2 + v0), log = TRUE)
  closed <- sum(0.5 * ((v + m^2) / v0 - 1 - log(v / v0)))
  assert_near(multimodal_yfb_coordinate_kl(x, s2, m, v + m^2, log_ml), closed, tol = 1e-10)
  # the signed-posterior log_ml feeds the identity consistently
  post <- multimodal_yfb_signed_posterior(1 / s2[1], x[1] / s2[1], list(family = "normal", variance = v0))
  assert_near(post$log_ml, log_ml[1], tol = 1e-12)
  # no information: KL = 0
  assert_equal(multimodal_yfb_coordinate_kl(NA, Inf, 0, 1, NA), 0)
})

run_test("MMYFB-Prune-T2: point-Laplace KL is nonnegative and zero when the prior is a point mass", {
  prior <- list(family = "point_laplace", pi = 0.3, rate = 2)
  kl <- vapply(c(-3, -0.2, 0, 0.5, 4), function(x) {
    p <- multimodal_yfb_signed_posterior(4, 4 * x, prior)
    multimodal_yfb_coordinate_kl(p$x, p$s2, p$mean, p$second, p$log_ml)
  }, numeric(1))
  assert_true(all(kl > -1e-10))
  pm <- multimodal_yfb_signed_posterior(4, 2, list(family = "point_laplace", point_mass = TRUE))
  assert_equal(multimodal_yfb_coordinate_kl(pm$x, pm$s2, pm$mean, pm$second, pm$log_ml), 0)
})

run_test("MMYFB-Prune-T3: residual-precision models return the stated estimates", {
  r2 <- list(expression = c(10, 20, 40), methylation = c(5, 5))
  f <- multimodal_yfb_tau_from_residuals(r2, 10, "feature")
  assert_near(f$Tau$expression, 10 / c(10, 20, 40))
  m <- multimodal_yfb_tau_from_residuals(r2, 10, "modality")
  assert_near(m$Tau$expression, rep(30 / 70, 3))
  eb <- multimodal_yfb_tau_from_residuals(list(expression = c(r2$expression, 1e-12)), 10, "feature_eb")
  # the Gamma prior bounds a near-zero residual's precision: E[tau] <= (a + n/2) / b
  h <- eb$hyper$expression
  assert_true(max(eb$Tau$expression) <= (h[["a"]] + 5) / h[["b"]] + 1e-8)
})

run_test("MMYFB-Prune-T4: nullcheck gives zero ELBO change for an all-zero factor", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 60, n_validation = 20,
                                        p_expression = 12, p_methylation = 12, seed = 41)
  fit <- fit_multimodal_yfb(data$training$Y, data$training$time, data$training$event, K = 2,
                            control = list(max_outer = 15L))
  # append an exactly-zero third factor with zero KL
  fit$EL <- cbind(fit$EL, 0); fit$EL2 <- cbind(fit$EL2, 0)
  fit$EF <- lapply(fit$EF, function(x) cbind(x, 0)); fit$EF2 <- lapply(fit$EF2, function(x) cbind(x, 0))
  fit$EBeta <- c(fit$EBeta, 0); fit$EBeta2 <- c(fit$EBeta2, 0)
  fit$kl <- list(L = c(fit$kl$L, 0), F = lapply(fit$kl$F, function(x) c(x, 0)), beta = c(fit$kl$beta, 0))
  delta <- multimodal_yfb_nullcheck(fit, data$training$Y, data$training$time, data$training$event)
  assert_near(delta[3], 0, tol = 1e-8)
})

run_test("MMYFB-Prune-T5: pruning removes spare factors at moderate signal-to-noise", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 120, n_validation = 60,
                                        p_expression = 30, p_methylation = 30, seed = 42,
                                        truncate = FALSE, noise_scale = 0.25)
  fit <- fit_multimodal_yfb(data$training$Y, data$training$time, data$training$event, K = 6,
                            control = list(intercept = TRUE, prior_update = "ebnm", prune = TRUE,
                              prior_F = list(expression = "point_laplace", methylation = "point_laplace")))
  # true K = 3; pruning must remove at least two of the three spare factors
  assert_true(ncol(fit$EL) <= 4L)
  assert_true(NROW(fit$diagnostics$pruning) >= 2L)
})

run_test("MMYFB-Prune-T6: with an intercept, a constant feature shift does not change the fit", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 60, n_validation = 30,
                                        p_expression = 12, p_methylation = 12, seed = 43,
                                        truncate = FALSE, noise_scale = 0.5)
  ctl <- list(intercept = TRUE, max_outer = 20L,
              prior_F = list(expression = "point_laplace", methylation = "point_laplace"))
  shift <- function(Y) Map(function(x, s) x + s, Y, list(5, 2))
  a <- fit_multimodal_yfb(data$training$Y, data$training$time, data$training$event, 2, ctl)
  b <- fit_multimodal_yfb(shift(data$training$Y), data$training$time, data$training$event, 2, ctl)
  assert_near(a$EF$expression, b$EF$expression, tol = 1e-6)
  assert_near(b$mu$expression - a$mu$expression, rep(5, 12), tol = 1e-6)
})

run_test("MMYFB-Prune-T7: the compiled loading sweep reproduces the R sweep for every prior family", {
  if (!multimodal_yfb_use_cpp()) {
    cat("  (compiled sweep unavailable; skipping the comparison)\n")
  } else {
    set.seed(51)
    n <- 40; p <- 25
    Y <- matrix(rnorm(n * p), n, p); R <- matrix(rnorm(n * p), n, p)
    EL <- rexp(n)
    args <- list(Y_m = Y, Tau_m = rgamma(p, 2, 2), EL_k = EL, EL2_k = EL^2 + rexp(n),
                 R_mk = R, w = runif(n), h_minus_k = rnorm(n), EZ_k = rnorm(n),
                 VZ_k = runif(n), EF_k = EF0 <- rnorm(p, sd = 0.3), EF2_k = EF0^2 + 0.1,
                 EBeta_k = 0.4, EBeta2_k = 0.3)
    priors <- list(list(pi = 0.4, rate = 1.5, point_mass = FALSE),
                   list(family = "point_laplace", pi = 0.4, rate = 1.5, point_mass = FALSE),
                   list(family = "normal", variance = 0.8, point_mass = FALSE),
                   list(family = "point_laplace", point_mass = TRUE))
    for (pr in priors) {
      a <- do.call(multimodal_yfb_update_F_mk, c(args, list(prior = pr)))
      old <- options(multimodal_yfb.engine = "R")
      b <- do.call(multimodal_yfb_update_F_mk, c(args, list(prior = pr)))
      options(old)
      for (nm in c("EF", "EF2", "EZ", "VZ", "kl")) assert_near(a[[nm]], b[[nm]], tol = 1e-10)
      da <- as.matrix(a$details); db <- as.matrix(b$details)
      # log p(x) is NA by design for point-mass / zero-information coordinates
      assert_true(identical(is.na(da), is.na(db)))
      assert_near(da[!is.na(da)], db[!is.na(db)], tol = 1e-10)
    }
  }
})
