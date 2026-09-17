# ============================================================
# Script: test_multimodal_yfb_fit.R
# Purpose: Integration tests for isolated multimodal YFB fitting and prediction.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("fit_multimodal_yfb")) source("code/fit_multimodal_yfb.R")
if (!exists("predict_multimodal_yfb")) source("code/predict_multimodal_yfb.R")

.mm_fit_fixture <- function() {
  set.seed(17)
  ids <- paste0("s", 1:12)
  L <- matrix(rexp(24), 12, 2)
  list(Y = list(expression = structure(L %*% t(matrix(c(1, .5, .2, 1), 2)) + matrix(rexp(24, 8), 12), dimnames = list(ids, paste0("g", 1:2))), methylation = structure(L %*% t(matrix(c(.4, .8, .6, .3), 2)) + matrix(rexp(24, 8), 12), dimnames = list(ids, paste0("cg", 1:2)))), time = setNames(rexp(12), ids), event = setNames(rep(c(1, 0), 6), ids))
}
run_test("MMYFB-Fit-T1: deterministic small fit has finite valid moments", {
  d <- .mm_fit_fixture()
  a <- fit_multimodal_yfb(d$Y, d$time, d$event, K = 2,
                           control = list(max_outer = 2, max_inner = 2))
  b <- fit_multimodal_yfb(d$Y, d$time, d$event, K = 2,
                           control = list(max_outer = 2, max_inner = 2))
  assert_near(a$EBeta, b$EBeta, tol = 1e-12)
  assert_true(all(is.finite(a$diagnostics$objective)))
  assert_true(all(is.finite(a$diagnostics$fixed_working_objective)))
  assert_true(all(a$EL2 >= a$EL^2))
  assert_true(all(vapply(a$EF2, function(x) all(is.finite(x)), logical(1))))
})
run_test("MMYFB-Fit-T2: frozen prediction ignores validation outcomes", {
  d <- .mm_fit_fixture()
  fit <- fit_multimodal_yfb(d$Y, d$time, d$event, K = 2,
                             control = list(max_outer = 1, max_inner = 1))
  pred <- predict_multimodal_yfb(fit, d$Y)
  expected <- as.vector((d$Y$expression %*% fit$EF$expression +
    d$Y$methylation %*% fit$EF$methylation) %*% fit$EBeta)
  assert_near(pred$risk_scores, expected, tol = 1e-12)
})
run_test("MMYFB-Fit-T3: diagnostics retain configured control and factor flags", {
  d <- .mm_fit_fixture()
  fit <- fit_multimodal_yfb(
    d$Y, d$time, d$event, K = 2,
    control = list(max_outer = 2, max_inner = 1, damping = 0.5,
                   tau_chunk_size = 1)
  )
  assert_true(identical(fit$diagnostics$controls$damping, 0.5))
  assert_true(length(fit$diagnostics$history) == fit$diagnostics$iterations)
  assert_true(length(fit$diagnostics$reconstruction_active) == 2)
  assert_true(length(fit$diagnostics$prognostic) == 2)
})
run_test("MMYFB-Fit-T4: chunked tau equals the exact unchunked update", {
  d <- .mm_fit_fixture()
  EL <- matrix(seq_len(nrow(d$Y$expression)), ncol = 1)
  EL2 <- EL^2 + 0.2
  EF <- lapply(d$Y, function(x) matrix(seq_len(ncol(x)) / 10, ncol = 1))
  EF2 <- lapply(EF, function(x) x^2 + 0.1)
  exact <- multimodal_yfb_update_tau(d$Y, EL, EL2, EF, EF2)
  chunked <- multimodal_yfb_update_tau_chunked(d$Y, EL, EL2, EF, EF2,
                                                chunk_size = 1)
  for (modality in names(d$Y)) {
    assert_near(chunked$Tau[[modality]], exact$Tau[[modality]], tol = 1e-12)
  }
})
if (sys.nframe() == 0L) report_results("test_multimodal_yfb_fit.R")
