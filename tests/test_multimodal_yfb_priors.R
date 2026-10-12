# ============================================================
# Script: test_multimodal_yfb_priors.R
# Purpose: Tests for the signed (point-Laplace, Normal) loading priors in the
#          multimodal YFB model.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: ebnm; code/multimodal_yfb_updates.R, code/fit_multimodal_yfb.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("multimodal_yfb_signed_posterior")) source("code/multimodal_yfb_updates.R")
if (!exists("fit_multimodal_yfb")) source("code/fit_multimodal_yfb.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")

run_test("MMYFB-Prior-T1: point-Laplace posterior matches ebnm with the prior fixed", {
  x <- c(-6, -1.2, -0.1, 0, 0.4, 2.5, 9)
  s <- c(0.5, 1, 2, 1, 0.3, 1.5, 0.8)
  prior <- list(family = "point_laplace", pi = 0.3, rate = 1.7)
  g <- ebnm::laplacemix(pi = c(0.7, 0.3), mean = c(0, 0), scale = c(0, 1 / 1.7))
  ref <- ebnm::ebnm_point_laplace(x, s, g_init = g, fix_g = TRUE)$posterior
  ours <- lapply(seq_along(x), function(i) {
    A <- 1 / s[i]^2
    multimodal_yfb_signed_posterior(A, x[i] * A, prior)
  })
  mean_ours <- vapply(ours, `[[`, numeric(1), "mean")
  sd_ours <- sqrt(vapply(ours, `[[`, numeric(1), "second") - mean_ours^2)
  assert_near(mean_ours, ref$mean, tol = 1e-6)
  assert_near(sd_ours, ref$sd, tol = 1e-6)
})

run_test("MMYFB-Prior-T2: Normal posterior matches ebnm with the prior fixed", {
  x <- c(-3, 0.2, 4); s <- c(1, 0.5, 2)
  g <- ashr::normalmix(pi = 1, mean = 0, sd = 1.5)
  ref <- ebnm::ebnm_normal(x, s, g_init = g, fix_g = TRUE)$posterior
  ours <- lapply(seq_along(x), function(i) {
    A <- 1 / s[i]^2
    multimodal_yfb_signed_posterior(A, x[i] * A, list(family = "normal", variance = 1.5^2))
  })
  mean_ours <- vapply(ours, `[[`, numeric(1), "mean")
  assert_near(mean_ours, ref$mean, tol = 1e-10)
  assert_near(sqrt(vapply(ours, `[[`, numeric(1), "second") - mean_ours^2), ref$sd, tol = 1e-10)
})

run_test("MMYFB-Prior-T3: empirical-Bayes prior fit recovers the generating prior", {
  set.seed(11)
  n <- 20000
  theta <- ifelse(runif(n) < 0.3, rexp(n, 2) * sample(c(-1, 1), n, TRUE), 0)
  x <- theta + rnorm(n, sd = 0.2)
  laplace <- multimodal_yfb_fit_signed_prior(x, rep(0.04, n), "point_laplace")
  assert_near(laplace$pi, 0.3, tol = 0.03)
  assert_near(laplace$rate, 2, tol = 0.2)
  theta_n <- rnorm(n, sd = 1.3)
  normal <- multimodal_yfb_fit_signed_prior(theta_n + rnorm(n, sd = 0.5), rep(0.25, n), "normal")
  assert_near(normal$variance, 1.3^2, tol = 0.1)
  # no informative pseudo-observations: point mass
  empty <- multimodal_yfb_fit_signed_prior(c(NA, 1), c(Inf, Inf), "point_laplace")
  assert_true(empty$point_mass)
})

run_test("MMYFB-Prior-T4: a zero-information coordinate returns the prior moments", {
  lap <- multimodal_yfb_signed_posterior(0, 0, list(family = "point_laplace", pi = 0.4, rate = 2))
  assert_equal(lap$mean, 0)
  assert_near(lap$second, 2 * 0.4 / 4)
  nor <- multimodal_yfb_signed_posterior(0, 0, list(family = "normal", variance = 3))
  assert_near(nor$second, 3)
})

run_test("MMYFB-Prior-T5: canonicalization rescales Laplace rates and Normal variances", {
  Y <- list(expression = matrix(c(1, 3), 2, 1), methylation = matrix(c(2, 0), 2, 1))
  EF <- list(expression = matrix(2, 1, 1), methylation = matrix(1, 1, 1))
  EF2 <- lapply(EF, function(x) x^2 + 0.1)
  prior_F <- list(expression = list(list(family = "point_laplace", pi = 0.5, rate = 1, point_mass = FALSE)),
                  methylation = list(list(family = "normal", variance = 4, point_mass = FALSE)))
  out <- multimodal_yfb_canonicalize_factors(
    Y, matrix(1, 2, 1), matrix(1.1, 2, 1), EF, EF2, 1, 1.1,
    list(list(pi = 1, rate = 1, point_mass = FALSE)), prior_F)
  c_k <- out$scale[1]
  assert_near(out$prior_F$expression[[1]]$rate, c_k)
  assert_near(out$prior_F$methylation[[1]]$variance, 4 / c_k^2)
})

run_test("MMYFB-Prior-T6: fitter accepts signed loading priors and returns signed loadings", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 60,
                                        n_validation = 30, p_expression = 15,
                                        p_methylation = 15, seed = 31)
  fit <- fit_multimodal_yfb(data$training$Y, data$training$time, data$training$event,
                            K = 3, control = list(max_outer = 30L,
                              prior_F = list(expression = "point_laplace",
                                             methylation = "normal")))
  assert_equal(fit$prior_F$expression[[1]]$family, "point_laplace")
  assert_equal(fit$prior_F$methylation[[1]]$family, "normal")
  assert_true(all(is.finite(unlist(fit$EF))))
  assert_true(any(unlist(fit$EF) < 0))
  risk <- predict_multimodal_yfb(fit, data$validation$Y)$risk_scores
  assert_true(all(is.finite(risk)))
  bad <- tryCatch(fit_multimodal_yfb(data$training$Y, data$training$time,
                                     data$training$event, K = 2,
                                     control = list(prior_F = list(expression = "cauchy"))),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("prior_F", bad))
})

run_test("MMYFB-Prior-T7: signed priors accept centered data; point-exponential still rejects it", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 50,
                                        n_validation = 30, p_expression = 12,
                                        p_methylation = 12, seed = 32)
  # center each modality on the training means; apply the same means to validation
  mu <- lapply(data$training$Y, colMeans)
  center <- function(Y) Map(function(x, m) sweep(x, 2, m), Y, mu)
  Ytr <- center(data$training$Y); Yva <- center(data$validation$Y)
  fit <- fit_multimodal_yfb(Ytr, data$training$time, data$training$event, K = 2,
                            control = list(max_outer = 20L,
                              prior_F = list(expression = "point_laplace",
                                             methylation = "point_laplace")))
  assert_true(all(is.finite(predict_multimodal_yfb(fit, Yva)$risk_scores)))
  err <- tryCatch(fit_multimodal_yfb(Ytr, data$training$time, data$training$event, K = 2,
                                     control = list(max_outer = 2L)),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("nonnegative", err))
  # mixed: only the signed-prior modality may be signed
  mixed <- list(expression = data$training$Y$expression, methylation = Ytr$methylation)
  ok <- fit_multimodal_yfb(mixed, data$training$time, data$training$event, K = 2,
                           control = list(max_outer = 2L,
                             prior_F = list(expression = "point_exponential",
                                            methylation = "normal")))
  assert_true(inherits(ok, "multimodal_yfb_fit"))
})
