# ============================================================
# Script: test_multimodal_m1_equivalence.R
# Purpose: Check registered one-modality equivalence and KL edge cases.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: ebnm, survival, yaml and project update modules
# ============================================================
if (sys.nframe() == 0L) {
  source("tests/test_helpers.R")
  source("code/fit_ssbmf_yfb.R")
  source("code/fit_multimodal_yfb.R")
}

# Identical inputs, prior, weighting, projection centering and stopping options
# are required for the comparison. Unsupported arguments must fail loudly.
set.seed(20261112)
n <- 32L
Y <- matrix(rexp(n * 12L), n, 12L,
            dimnames = list(paste0("patient", seq_len(n)), paste0("feature", 1:12)))
time <- rexp(n) + seq_len(n) / 10000
status <- rep(c(1L, 1L, 0L, 1L), length.out = n)
names(time) <- names(status) <- rownames(Y)

#' Compare a single-modality fit with the M = 1 fit under matched settings.
#' @param strata NULL or a stratum label for each patient.
#' @return Invisible NULL; fails if a setting or result does not match.
#' @examples
#' # check_m1_equivalence(NULL)
check_m1_equivalence <- function(strata) {
  single <- fit_ssbmf_yfb(Y, time, status, K = 2, max_iter = 5, tol = 1e-5,
                           prior_LF = "point_exponential", prior_beta = "normal",
                           canonicalize = TRUE, w_S = 1, center_projection = TRUE,
                           strata_id = strata, verbose = FALSE)
  multi <- fit_multimodal_yfb(list(expression = Y), time, status, K = 2,
    control = list(max_outer = 5, max_inner = 1, tol = 1e-5,
                   prior_update = "ebnm", prior_F = list(expression = "point_exponential"),
                   survival_elbo = "corrected", w_S = 1, center_projection = TRUE),
    strata_id = strata)
  pairs <- list(EL = list(single$EL, multi$EL), EF = list(single$EF, multi$EF$expression),
                EBeta = list(single$EBeta, multi$EBeta),
                ELBO = list(tail(single$history$elbo_full, 1), multi$diagnostics$elbo))
  for (name in names(pairs)) {
    a <- pairs[[name]][[1]]; b <- pairs[[name]][[2]]
    assert_true(length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)), name)
    assert_true(max(abs(a - b)) <= 1e-6 * max(1, abs(a), abs(b)), name)
  }
  invisible(NULL)
}

run_test("M = 1 pooled, centered projection", { check_m1_equivalence(NULL) })
run_test("M = 1 stratified, centered projection", { check_m1_equivalence(rep(1:2, each = n / 2)) })

run_test("Scalar precision normalizer agrees with shared KL", {
  x <- c(-1, 0.2, 2); s2 <- 0.7
  mean <- c(-0.4, 0.1, 1); second <- mean^2 + 0.3
  log_ml <- dnorm(x, 0, sqrt(1 + s2), log = TRUE)
  actual <- multimodal_yfb_coordinate_kl(x, s2, mean, second, log_ml)
  expected <- -compute_ebnm_kl(sum(log_ml), 1 / s2, x, mean, second)
  assert_near(actual, expected, tol = 1e-10)
})
run_test("Flat beta likelihood retains shared prior second moment", {
  g <- ashr::normalmix(1, 0, 1)
  result <- multimodal_yfb_update_beta_shared_k(rep(0, 4), rep(0, 4), rep(1, 4), rep(0, 4), g, "normal")
  assert_near(result$mean, 0, tol = 1e-12)
  assert_near(result$second, 1, tol = 1e-12)
  assert_near(result$kl, 0, tol = 1e-12)
})
if (sys.nframe() == 0L) {
  r <- report_results("test_multimodal_m1_equivalence.R")
  if (r$failed > 0) quit(status = 1)
}
