# ============================================================
# Script: test_multimodal_yfb_updates.R
# Purpose: Test exact modular updates for the isolated multimodal YFB model.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("multimodal_yfb_update_L_k")) source("code/multimodal_yfb_updates.R")

cat("\n========================================\n")
cat("  Tests: multimodal_yfb_updates.R\n")
cat("========================================\n\n")

.mm_prior <- list(pi = 0.6, rate = 1.2, point_mass = FALSE)

# T1: shared score update ----

run_test("MMYFB-Updates-T1: shared L A and B pool both modality blocks", {
  Y <- list(expression = matrix(c(3, 4), ncol = 1),
            methylation = matrix(c(5, 6), ncol = 1))
  R <- list(expression = matrix(c(2, 3), ncol = 1),
            methylation = matrix(c(4, 5), ncol = 1))
  update <- multimodal_yfb_update_L_k(Y, R,
    EF_k = list(expression = 2, methylation = 1),
    EF2_k = list(expression = 4.5, methylation = 1.25),
    Tau = list(expression = 3, methylation = 2), prior = .mm_prior)
  assert_near(update$A, 16, tol = 1e-12)
  assert_near(update$B, c(20, 28), tol = 1e-12)
  assert_near(update$x, c(1.25, 1.75), tol = 1e-12)
  assert_near(update$s2, rep(1 / 16, 2), tol = 1e-12)
  assert_true(all(update$second >= update$mean^2))
})

# T2: exact sequential loading update ----

run_test("MMYFB-Updates-T2: F update includes same-factor cross-feature term", {
  update <- multimodal_yfb_update_F_mk(
    Y_m = cbind(c(1, 2), c(3, 4)), Tau_m = c(2, 3),
    EL_k = c(1, 2), EL2_k = c(1.2, 4.5),
    R_mk = cbind(c(5, 6), c(7, 8)), w = c(0.5, 1),
    h_minus_k = c(2, 1), EZ_k = c(7, 10), VZ_k = c(0.3, 0.4),
    EF_k = c(1, 2), EF2_k = c(1.1, 4.2),
    EBeta_k = 0.5, EBeta2_k = 0.8, prior = .mm_prior)
  expected_A1 <- 2 * sum(c(1.2, 4.5)) + 0.8 * sum(c(0.5, 1) * c(1, 2)^2)
  expected_B1 <- 2 * sum(c(1, 2) * c(5, 6)) +
    sum(c(1, 2) * (0.5 * c(2, 1) - c(0.5, 1) * 0.8 * c(6, 8)))
  assert_near(update$details$A[1], expected_A1, tol = 1e-12)
  assert_near(update$details$B[1], expected_B1, tol = 1e-12)
  assert_near(update$details$x[1], expected_B1 / expected_A1, tol = 1e-12)
  assert_near(update$details$s2[1], 1 / expected_A1, tol = 1e-12)
  no_interaction <- 2 * sum(c(1, 2) * c(5, 6)) +
    sum(c(1, 2) * (0.5 * c(2, 1)))
  assert_true(abs(update$details$B[1] - no_interaction) > 1e-6,
              msg = "same-factor interaction must affect B")
})

run_test("MMYFB-Updates-T3: F update uses EBeta2 rather than squared EBeta", {
  args <- list(
    Y_m = cbind(c(1, 2), c(3, 4)), Tau_m = c(2, 3),
    EL_k = c(1, 2), EL2_k = c(1.2, 4.5),
    R_mk = cbind(c(5, 6), c(7, 8)), w = c(0.5, 1),
    h_minus_k = c(2, 1), EZ_k = c(7, 10), VZ_k = c(0.3, 0.4),
    EF_k = c(1, 2), EF2_k = c(1.1, 4.2), EBeta_k = 0.5,
    prior = .mm_prior
  )
  low <- do.call(multimodal_yfb_update_F_mk, c(args, list(EBeta2_k = 0.25)))
  high <- do.call(multimodal_yfb_update_F_mk, c(args, list(EBeta2_k = 0.8)))
  assert_true(high$details$A[1] > low$details$A[1])
})

# T3: prior, beta, and tau boundaries ----

run_test("MMYFB-Updates-T4: loading-prior update is modality-specific and handles zero slab mass", {
  prior_expression <- multimodal_yfb_update_loading_prior(c(1, 0.5), c(0.5, 1))
  prior_methylation <- multimodal_yfb_update_loading_prior(c(0, 0), c(0, 0))
  assert_near(prior_expression$pi, 0.75, tol = 1e-12)
  assert_near(prior_expression$rate, 1, tol = 1e-12)
  assert_true(prior_methylation$point_mass)
})

run_test("MMYFB-Updates-T4b: loading-prior update treats zero posterior magnitude as a point mass", {
  prior <- multimodal_yfb_update_loading_prior(c(0.5, 0.5), c(0, 0))
  assert_true(prior$point_mass)
})

run_test("MMYFB-Updates-T5: beta update includes projection uncertainty", {
  beta <- multimodal_yfb_update_beta_k(EZ_k = c(2, 3), VZ_k = c(1, 2),
                                        w = c(0.5, 1), h_minus_k = c(1, 2))
  assert_near(beta$A, 0.5 * 5 + 11, tol = 1e-12)
  assert_near(beta$B, 8, tol = 1e-12)
  assert_near(beta$x, 8 / 13.5, tol = 1e-12)
  expected_prior_var <- (8 / 13.5)^2 - 1 / 13.5
  expected_var <- 1 / (13.5 + 1 / expected_prior_var)
  assert_near(beta$variance, expected_var, tol = 1e-12)
  assert_near(beta$mean, expected_var * 8, tol = 1e-12)
  assert_true(beta$second >= beta$mean^2)
})

run_test("MMYFB-Updates-T6: tau is n divided by each expected residual sum", {
  tau <- multimodal_yfb_update_tau(
    Y = list(expression = matrix(c(2, 3), ncol = 1),
             methylation = matrix(c(4, 5), ncol = 1)),
    EL = matrix(c(1, 2), ncol = 1), EL2 = matrix(c(1.1, 4.2), ncol = 1),
    EF = list(expression = matrix(1, ncol = 1), methylation = matrix(2, ncol = 1)),
    EF2 = list(expression = matrix(1.2, ncol = 1), methylation = matrix(4.5, ncol = 1))
  )
  assert_near(tau$expected_residual2$expression[, 1], c(1.32, 2.04), tol = 1e-12)
  assert_near(tau$Tau$expression, 2 / 3.36, tol = 1e-12)
  assert_near(tau$Tau$methylation, 2 / 8.85, tol = 1e-12)
})

run_test("MMYFB-Updates-T7: zero-information updates retain their stated prior boundaries", {
  prior_mean <- multimodal_yfb_point_exponential_posterior(0, 0, .mm_prior)
  underflow_mean <- multimodal_yfb_point_exponential_posterior(0, -1e-238, .mm_prior)
  beta <- multimodal_yfb_update_beta_k(c(0, 0), c(0, 0), c(0, 0), c(0, 0),
                                        prior_variance = 2)
  assert_near(prior_mean$mean, 0.6 / 1.2, tol = 1e-12)
  assert_near(underflow_mean$mean, prior_mean$mean, tol = 1e-12)
  assert_near(prior_mean$second, 2 * 0.6 / 1.2^2, tol = 1e-12)
  assert_near(beta$second, 2, tol = 1e-12)
})

run_test("MMYFB-Updates-T8: tau reports zero-residual boundary", {
  error <- tryCatch({
    multimodal_yfb_update_tau(list(expression = matrix(1, 1, 1),
                                   methylation = matrix(1, 1, 1)),
                               matrix(1, 1, 1), matrix(1, 1, 1),
                               list(expression = matrix(1, 1, 1),
                                    methylation = matrix(1, 1, 1)),
                               list(expression = matrix(1, 1, 1),
                                    methylation = matrix(1, 1, 1)))
    NULL
  }, error = conditionMessage)
  assert_true(grepl("zero", error, ignore.case = TRUE))
})

run_test("MMYFB-Updates-T9: extreme negative pseudo-observations retain valid moments", {
  posterior <- multimodal_yfb_point_exponential_posterior(
    A = 1, B = -1e5, prior = .mm_prior
  )
  assert_true(is.finite(posterior$mean) && is.finite(posterior$second))
  assert_true(posterior$mean >= 0)
  assert_true(posterior$second >= posterior$mean^2)
})

run_test("MMYFB-Updates-T10: shared-score update refits its EB prior", {
  update <- multimodal_yfb_update_L_k(
    Y = list(expression = matrix(c(3, 4), ncol = 1),
             methylation = matrix(c(5, 6), ncol = 1)),
    R_minus_k = list(expression = matrix(c(2, 3), ncol = 1),
                       methylation = matrix(c(4, 5), ncol = 1)),
    EF_k = list(expression = 2, methylation = 1),
    EF2_k = list(expression = 4.5, methylation = 1.25),
    Tau = list(expression = 3, methylation = 2), prior = .mm_prior
  )
  expected <- multimodal_yfb_update_loading_prior(update$slab_prob, update$mean)
  assert_true(is.list(update$prior),
              msg = "shared-score update must return its refitted prior")
  assert_near(update$prior$pi, expected$pi, tol = 1e-12)
  assert_near(update$prior$rate, expected$rate, tol = 1e-12)
})


if (sys.nframe() == 0L) report_results("test_multimodal_yfb_updates.R")
