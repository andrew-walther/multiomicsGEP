# =============================================================================
# tests/test_empty_program_kl.R
# The beta KL of an empty program (all scores zero) must be zero.
#
# Found in attempt KL-1 (docs/plans/2026-11-06-kl-fix-pruning-rerun-plan.md).
# For A = B = 0 the β likelihood is flat, so q(β) = g and -KL = 0.
# Flooring A before solving introduced artificial information and exposed
# numerical instability in the point-Laplace solver.
# Run: Rscript tests/test_empty_program_kl.R
# =============================================================================

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(ebnm))
  source("tests/test_helpers.R")
  source("code/update_beta.R")
}

# Shared point-Laplace beta prior from a KL-1 fit (seed 20261016, noise 0.25):
# pi0 = 0.5, Laplace scale 0.363, so E_g[beta^2] = 0.5 * 2 * 0.363^2.
g_dead <- structure(data.frame(pi = c(.5, .5), mean = c(0, 0), scale = c(0, .363)),
                    class = c("laplacemix", "data.frame"))
prior_second <- 0.5 * 2 * 0.363^2

run_test("Empty program: q(beta_k) equals the shared prior", {
  set.seed(1)
  n <- 50
  res <- update_beta_k(w = runif(n), z_no_k = rnorm(n), EL_k = rep(0, n), EL2_k = rep(0, n),
                       prior_family = "point_laplace", alpha = 1, g_fixed = g_dead)
  assert_near(res$mean, 0, tol = 1e-12, msg = "posterior mean of an empty program")
  assert_near(res$second, prior_second, tol = 1e-6,
              msg = "posterior second moment must equal the prior's when the likelihood is flat")
})

run_test("Empty program: beta -KL is zero, never positive", {
  set.seed(1)
  n <- 50
  res <- update_beta_k(w = runif(n), z_no_k = rnorm(n), EL_k = rep(0, n), EL2_k = rep(0, n),
                       prior_family = "point_laplace", alpha = 1, g_fixed = g_dead)
  neg_kl <- res$neg_kl
  assert_true(res$A == 0 && is.infinite(res$s), "flat likelihood has zero precision")
  assert_true(neg_kl <= 1e-12, sprintf("-KL = %.3e is positive (KL must be >= 0)", neg_kl))
  assert_near(neg_kl, 0, tol = 1e-9, msg = "-KL of an empty program")
})

run_test("Weak beta information still goes through EBNM below the floor", {
  g <- ashr::normalmix(1, 0, 1)
  res <- update_beta_k(w = 1, z_no_k = 1e-5, EL_k = 1e-6, EL2_k = 1e-12,
                       prior_family = "normal", alpha = 1, g_fixed = g)
  assert_near(res$mean, 1e-11 / (1 + 1e-10), tol = 1e-14,
              msg = "positive unfloored information shifts the posterior mean")
  assert_true(is.finite(res$s), "weak likelihood contributes to the prior fit")
})

if (sys.nframe() == 0L) {
  r <- report_results("test_empty_program_kl.R")
  if (r$failed > 0) quit(status = 1)
}
