# =============================================================================
# tests/test_empty_program_kl.R
# The beta KL of an empty program (all scores zero) must be zero.
#
# Found in attempt KL-1 (docs/plans/2026-11-06-kl-fix-pruning-rerun-plan.md).
# When program k has zero scores, A_k is floored at 1e-10, so the EBNM noise
# SD is s_k = 1e5 and the likelihood for beta_k is flat. Then q(beta_k) = g and
#   -KL(q || g) = 0.
# ebnm 1.1.42's point-Laplace solver is inaccurate at s = 1e5 (posterior second
# moment 65535.87 against the prior's 0.1318; log-likelihood off by ~3e-6), so
# compute_ebnm_kl() returns -KL = +8.2e-6. prune_yf_nullcheck() then gives an
# empty program a deletion gain of about -1e-5, and its rule max(gain) >= 0
# keeps or removes empty programs according to the sign of that residue.
#
# Status: FAILS on the current code (2026-11-06). Not yet in run_tests.R, so
# the shared suite is not broken while the fix is decided.
# Run: Rscript tests/test_empty_program_kl.R
# =============================================================================

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(ebnm))
  source("tests/test_helpers.R")
  source("code/update_beta.R")
  source("code/compute_elbo.R")
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
  neg_kl <- compute_ebnm_kl(res$ebnm_result$log_likelihood, res$A, res$x,
                            res$mean, res$second)
  assert_true(neg_kl <= 1e-12, sprintf("-KL = %.3e is positive (KL must be >= 0)", neg_kl))
  assert_near(neg_kl, 0, tol = 1e-9, msg = "-KL of an empty program")
})

if (sys.nframe() == 0L) {
  r <- report_results("test_empty_program_kl.R")
  if (r$failed > 0) quit(status = 1)
}
