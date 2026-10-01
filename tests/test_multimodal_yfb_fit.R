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
  assert_true(all(vapply(fit$diagnostics$history, function(x) {
    is.finite(x$max_relative_reconstruction_change) &&
      x$max_relative_reconstruction_change >= 0 &&
      is.finite(x$max_relative_eta_change) && x$max_relative_eta_change >= 0
  }, logical(1))))
  assert_true(length(fit$diagnostics$reconstruction_active) == 2)
  assert_true(length(fit$diagnostics$prognostic) == 2)
  assert_true(length(fit$diagnostics$survival_active) == 2)
  assert_true(fit$diagnostics$K_eff_reconstruction == fit$diagnostics$K_eff)
  assert_true(fit$diagnostics$K_eff_retained >= fit$diagnostics$K_eff_reconstruction)
  assert_true(fit$diagnostics$K_eff_retained >= fit$diagnostics$K_eff_survival)
})
run_test("MMYFB-Fit-T3b: identifiable changes are relative to the fitted scale", {
  change <- multimodal_yfb_relative_change(c(10, -20), c(10.01, -20.02))
  assert_near(change, 0.02 / 20.02, tol = 1e-12)
  assert_near(multimodal_yfb_relative_change(c(0, 0), c(0.002, -0.001)),
              0.002, tol = 1e-12)
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
run_test("MMYFB-Fit-T5: Cox warm start retains a finite prognostic coefficient", {
  set.seed(29)
  projection <- cbind(seq(-2, 2, length.out = 60), rnorm(60, sd = 0.1))
  time <- rexp(60, rate = exp(projection[, 1]))
  event <- rep(1L, 60)
  warm_start <- multimodal_yfb_cox_warm_start(projection, time, event)
  assert_true(all(is.finite(warm_start$mean)))
  assert_true(abs(warm_start$mean[1]) > 1e-6)
  assert_true(all(warm_start$second >= warm_start$mean^2))
})
run_test("MMYFB-Fit-T6: canonicalization preserves fitted products and prior scales", {
  EL <- matrix(c(2, 3, 1, 4), nrow = 2)
  EL2 <- EL^2 + 0.1
  EF <- list(
    expression = matrix(c(3, 4, 1, 2), nrow = 2),
    methylation = matrix(c(2, 1, 5, 3), nrow = 2)
  )
  EF2 <- lapply(EF, function(x) x^2 + 0.2)
  EBeta <- c(0.4, -0.3)
  EBeta2 <- EBeta^2 + 0.05
  prior_L <- rep(list(list(pi = 0.5, rate = 2, point_mass = FALSE)), 2)
  prior_F <- lapply(EF, function(x) {
    rep(list(list(pi = 0.5, rate = 3, point_mass = FALSE)), ncol(x))
  })
  Y <- list(expression = matrix(c(1, 3, 2, 4), nrow = 2),
            methylation = matrix(c(4, 2, 1, 3), nrow = 2))
  reconstruction_before <- lapply(names(EF), function(modality) {
    EL %*% t(EF[[modality]])
  })
  names(reconstruction_before) <- names(EF)
  predictor_before <- Reduce(`+`, Map(`%*%`, Y, EF)) %*% EBeta
  scale_before <- apply(Reduce(`+`, Map(`%*%`, Y, EF)), 2, stats::sd)
  canonical <- multimodal_yfb_canonicalize_factors(
    Y, EL, EL2, EF, EF2, EBeta, EBeta2, prior_L, prior_F
  )
  reconstruction_after <- lapply(names(EF), function(modality) {
    canonical$EL %*% t(canonical$EF[[modality]])
  })
  names(reconstruction_after) <- names(EF)
  predictor_after <- Reduce(`+`, Map(`%*%`, Y, canonical$EF)) %*% canonical$EBeta
  for (modality in names(EF)) {
    assert_near(reconstruction_after[[modality]], reconstruction_before[[modality]],
                tol = 1e-12)
  }
  assert_near(predictor_after, predictor_before, tol = 1e-12)
  assert_near(apply(Reduce(`+`, Map(`%*%`, Y, canonical$EF)), 2, stats::sd),
              c(1, 1), tol = 1e-12)
  assert_true(all(canonical$EL2 >= canonical$EL^2))
  for (modality in names(EF)) {
    assert_true(all(canonical$EF2[[modality]] >= canonical$EF[[modality]]^2))
  }
  assert_near(vapply(canonical$prior_L, `[[`, numeric(1), "rate"),
              2 / scale_before, tol = 1e-12)
  for (modality in names(EF)) {
    assert_near(vapply(canonical$prior_F[[modality]], `[[`, numeric(1), "rate"),
                3 * scale_before, tol = 1e-12)
  }
})
run_test("MMYFB-Fit-T7: an empty control list uses configured defaults", {
  d <- .mm_fit_fixture()
  fit <- fit_multimodal_yfb(d$Y, d$time, d$event, K = 2, control = list())
  assert_true(is.list(fit$diagnostics$controls))
  assert_true(fit$diagnostics$controls$max_outer == 500)
})
if (sys.nframe() == 0L) report_results("test_multimodal_yfb_fit.R")
