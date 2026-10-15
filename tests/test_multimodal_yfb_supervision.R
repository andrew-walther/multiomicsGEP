# ============================================================
# Script: test_multimodal_yfb_supervision.R
# Purpose: Tests for survival supervision in the multimodal YFB fitter
#          (Workstream C6): shared beta prior update, supervised start,
#          beta-prior warm-up.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: ebnm; code/multimodal_yfb_updates.R, code/fit_multimodal_yfb.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("multimodal_yfb_coordinate_kl")) source("code/multimodal_yfb_updates.R")
if (!exists("fit_multimodal_yfb")) source("code/fit_multimodal_yfb.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")

local({
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor", n_train = 80L,
                                    n_validation = 20L, p_expression = 20L,
                                    p_methylation = 20L, seed = 3L, truncate = FALSE,
                                    noise_scale = 0.5, beta = c(0, 0, 3))
  .mms <<- list(Y = d$training$Y, time = d$training$time, event = d$training$event)
})
.mms_ctrl <- list(intercept = TRUE, prior_update = "ebnm",
                  prior_F = list(expression = "point_laplace", methylation = "point_laplace"))

run_test("MMYFB-Sup-T1: shared-prior beta update is the posterior under the fixed prior", {
  set.seed(2); EZ <- rnorm(40); VZ <- runif(40, 0, 0.1); w <- runif(40, 0.2, 1)
  h <- 0.6 * w * EZ + rnorm(40, sd = 0.3)
  g <- ebnm::laplacemix(c(0.4, 0.6), c(0, 0), c(0, 0.8))
  res <- multimodal_yfb_update_beta_shared_k(EZ, VZ, w, h, g)
  A <- sum(w * (EZ^2 + VZ)); B <- sum(EZ * h)
  ref <- ebnm::ebnm(B / A, sqrt(1 / A), prior_family = "point_laplace", g_init = g, fix_g = TRUE)
  assert_near(res$mean, ref$posterior$mean, tol = 1e-10, msg = "posterior mean")
  assert_near(res$second, ref$posterior$sd^2 + ref$posterior$mean^2, tol = 1e-10,
              msg = "posterior second moment")
  assert_true(res$kl >= -1e-10, msg = "KL nonnegative")
  zero <- multimodal_yfb_update_beta_shared_k(rep(0, 5), rep(0, 5), rep(1, 5), rep(0, 5), g)
  assert_true(zero$mean == 0 && zero$kl == 0, msg = "no information gives beta = 0")
})

run_test("MMYFB-Sup-T2: supervised start needs signed loadings and aligns factor K with u0", {
  err <- tryCatch(fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                                     control = list(init_supervised = TRUE, max_outer = 1L)),
                  error = function(e) e)
  assert_true(inherits(err, "error"), msg = "point-exponential loadings rejected")
  fit <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                            control = c(.mms_ctrl, list(init_supervised = TRUE, max_outer = 1L,
                                                        max_inner = 1L)))
  u0 <- multimodal_yfb_cox_working(.mms$time, .mms$event, rep(0, length(.mms$time)))$u
  Yc <- lapply(.mms$Y, function(x) sweep(x, 2, colMeans(x)))
  zK <- Reduce(`+`, Map(function(y, f) y %*% f[, 3], Yc, fit$EF))
  assert_true(abs(cor(zK, u0)) > 0.4,
              msg = sprintf("factor 3 score correlates with u0 (|cor| = %.2f)", abs(cor(zK, u0))))
})

run_test("MMYFB-Sup-T3: shared beta prior is held during the warm-up and refit after it", {
  ctrl <- c(.mms_ctrl, list(beta_prior = "shared_point_laplace", beta_warmup = 3L,
                            max_inner = 1L))
  g0 <- ebnm::laplacemix(c(0.5, 0.5), c(0, 0), c(0, 1))
  held <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                             control = c(ctrl, list(max_outer = 3L)))
  assert_true(identical(held$g_beta, g0), msg = "prior unchanged within the warm-up")
  refit <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                              control = c(ctrl, list(max_outer = 4L)))
  assert_true(!identical(refit$g_beta, g0), msg = "prior refit after the warm-up")
  assert_true(all(is.finite(refit$EBeta)), msg = "finite coefficients")
})
