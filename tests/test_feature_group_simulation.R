# ============================================================
# Script: test_feature_group_simulation.R
# Purpose: Check training-only scaling, precision shares and recovery scoring.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: test helpers and feature-group simulation runner
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
if (file.exists("results/benchmark_sim/run_feature_group_simulation.R"))
  source("results/benchmark_sim/run_feature_group_simulation.R")

run_test("Feature-group simulation: validation scaling uses training means and SDs", {
  train <- list(rna = cbind(1:3, 2 * (1:3)), meth = matrix(c(4, 6, 8), 3))
  valid <- list(rna = matrix(c(4, 8), 1), meth = matrix(10, 1))
  out <- fgsim_scale(train, valid)
  assert_near(out$training$rna[, 1], c(-1, 0, 1), 0)
  assert_near(out$validation$rna, matrix(c(2, 2), 1), 0)
  assert_near(out$validation$meth, matrix(2, 1), 0)
  assert_near(out$sd$rna, c(1, 2), 0)
})

run_test("Feature-group simulation: modality precision shares match hand calculation", {
  fit <- list(Tau = c(2, 4, 3), EF2 = matrix(1:6, 3, 2))
  out <- fgsim_precision(fit, c("rna", "rna", "meth"))
  assert_near(out$precision[out$program == 1], c(10, 9), 0)
  assert_near(out$share[out$program == 1], c(10, 9) / 19, 1e-14)
  assert_near(out$precision[out$program == 2], c(28, 18), 0)
  fit$EF2[, 2] <- 0
  dead <- fgsim_precision(fit, c("rna", "rna", "meth"))
  assert_true(all(is.na(dead$share[dead$program == 2])))
})

# The serialized digest freezes every value and attribute of the pre-change
# generator output for this small synthetic example (serialization version 2).
run_test("Feature-group simulation: NULL methylation sparsity preserves generator output", {
  data <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 8L, n_validation = 6L, p_expression = 7L, p_methylation = 9L,
    seed = 42L, truncate = FALSE, methylation_active_fraction = NULL)
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  writeBin(serialize(data, NULL, version = 2), path)
  assert_equal(unname(tools::md5sum(path)), "9ab31fd9e2b70cee6240d9a40c22987a")
})

run_test("Feature-group simulation: methylation sparsity uses ceiling per program", {
  data <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 8L, n_validation = 6L, p_expression = 50L, p_methylation = 250L,
    seed = 11L, truncate = FALSE, methylation_active_fraction = 0.03)
  assert_equal(unname(colSums(data$truth$F$methylation != 0)), rep(8, 3))
  assert_equal(unname(colSums(data$truth$F$expression != 0)), rep(8, 3))
  smaller <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 8L, n_validation = 6L, p_methylation = 50L, seed = 11L,
    methylation_active_fraction = 0.03)
  assert_equal(unname(colSums(smaller$truth$F$methylation != 0)), rep(2, 3))
})

run_test("Feature-group simulation: shared truth and survival use training scaling", {
  data <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 30L, n_validation = 15L, p_expression = 7L, p_methylation = 9L,
    seed = 42L, truncate = FALSE)
  out <- fgsim_prepare(data, 42L)
  scaled <- fgsim_scale(data$training$Y, data$validation$Y)
  g <- unlist(Map(function(F, s) F[, 3] / s, data$truth$F, scaled$sd), use.names = FALSE)
  g <- g / sqrt(sum(g^2))
  assert_near(as.vector(out$g), g, 1e-12)
  assert_near(sqrt(sum(out$g^2)), 1, 1e-12)
  assert_near(out$Y, do.call(cbind, scaled$training), 1e-12)
  assert_near(out$Y_validation, do.call(cbind, scaled$validation), 1e-12)
  assert_near(as.vector(out$eta), as.vector(out$Y %*% g * out$a), 1e-12)
  assert_near(as.vector(out$eta_validation), as.vector(out$Y_validation %*% g * out$a), 1e-12)
  assert_near(stats::sd(out$eta), stats::sd(data$training$eta), 1e-12)
  assert_equal(out$groups, rep(names(scaled$training), vapply(scaled$training, ncol, integer(1))))
  assert_equal(out$scaled, scaled)
  assert_near(out$training_center, colMeans(out$Y), 1e-12)
  assert_equal(length(out$time), 30L)
  assert_equal(length(out$event), 30L)
  assert_equal(length(out$time_validation), 15L)
  assert_equal(length(out$event_validation), 15L)
  assert_true(all(out$time > 0) && all(out$time_validation > 0))
  assert_true(all(out$event %in% 0:1) && all(out$event_validation %in% 0:1))
  set.seed(42L + 12000L)
  for (cohort in c("training", "validation")) {
    eta <- if (cohort == "training") out$eta else out$eta_validation
    event_time <- (-log(runif(length(eta))) /
      (data$truth$weibull_scale * exp(eta)))^(1 / data$truth$weibull_shape)
    censor_base <- rexp(length(eta))
    censor_time <- censor_base * multimodal_yfb_calibrate_censoring(
      event_time, censor_base, data$truth$target_censoring)
    observed_time <- if (cohort == "training") out$time else out$time_validation
    observed_event <- if (cohort == "training") out$event else out$event_validation
    assert_near(as.vector(observed_time), as.vector(pmin(event_time, censor_time)), 1e-12)
    assert_equal(as.integer(observed_event), as.integer(event_time <= censor_time))
  }
  assert_equal(out, fgsim_prepare(data, 42L))
  other <- fgsim_prepare(data, 43L)
  assert_false(identical(out$time, other$time))
})

run_test("Feature-group simulation: arm inputs reject any changed data or initialization", {
  a <- list(Y = matrix(1:6, 3), Y_validation = matrix(1:4, 2),
            time = c(1, 2, 3), event = c(1L, 0L, 1L),
            time_validation = c(1, 2), event_validation = c(0L, 1L),
            EL_init = matrix(1:6, 3), EF_init = matrix(1:4, 2))
  fgsim_assert_inputs(a, a)
  for (field in names(a)) {
    b <- a
    b[[field]][1] <- 9
    err <- tryCatch({ fgsim_assert_inputs(a, b); NULL }, error = identity)
    assert_true(inherits(err, "error"), paste("changed", field, "must fail"))
  }
})

#' Build a small prepared data set for arm-input regression checks.
#' @return Prepared synthetic training and validation inputs.
#' @examples
#' \dontrun{fgsim_test_inputs()}
fgsim_test_inputs <- function() {
  data <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 20L, n_validation = 12L, p_expression = 8L,
    p_methylation = 10L, seed = 19L, truncate = FALSE)
  fgsim_prepare(data, 19L)
}

run_test("Feature-group simulation: arms share joint starts and true risk", {
  x <- fgsim_test_inputs()
  arms <- fgsim_arm_inputs(x, 2L)
  joint <- arms$grouped
  expression <- arms$expression_only
  idx <- which(x$groups == "expression")
  assert_equal(joint, arms$pooled)
  assert_equal(expression$Y, joint$Y[, idx, drop = FALSE])
  assert_equal(expression$Y_validation, joint$Y_validation[, idx, drop = FALSE])
  assert_equal(expression$EL_init, joint$EL_init)
  assert_equal(expression$EF_init, joint$EF_init[idx, , drop = FALSE])
  assert_equal(expression$training_center, joint$training_center[idx])
  assert_equal(expression$groups, joint$groups[idx])
  assert_near(as.vector(expression$g), as.vector(joint$g[idx] /
    sqrt(sum(joint$g[idx]^2))), 1e-12)
  assert_near(sqrt(sum(expression$g^2)), 1, 1e-12)
  assert_near(joint$EL_init %*% t(joint$EF_init),
    svd(x$Y, nu = 2L, nv = 2L)$u[, 1:2, drop = FALSE] %*%
      diag(svd(x$Y, nu = 2L, nv = 2L)$d[1:2]) %*%
      t(svd(x$Y, nu = 2L, nv = 2L)$v[, 1:2, drop = FALSE]), 1e-10)
  expected <- frozen_reverse_cindex(x$eta_validation, x$time_validation,
                                    x$event_validation)
  for (arm in arms) {
    assert_equal(arm$eta, x$eta)
    assert_equal(arm$eta_validation, x$eta_validation)
    assert_equal(fgsim_true_risk(arm), expected)
  }
  fgsim_assert_expression_inputs(joint, expression)
})

run_test("Feature-group simulation: true-risk concordance uses frozen validation risk", {
  x <- fgsim_test_inputs()
  expected <- frozen_reverse_cindex(x$eta_validation, x$time_validation,
                                    x$event_validation)
  assert_equal(fgsim_true_risk(x), expected)
})

run_test("Feature-group simulation: expression arm rejects mismatched inputs", {
  arms <- fgsim_arm_inputs(fgsim_test_inputs(), 2L)
  for (field in c("Y", "Y_validation", "time", "event", "time_validation",
                  "event_validation", "EL_init", "EF_init", "training_center",
                  "g", "groups", "eta", "eta_validation")) {
    bad <- arms$expression_only
    if (field == "groups") bad[[field]][1] <- "methylation" else
      bad[[field]][1] <- bad[[field]][1] + 1
    err <- tryCatch({ fgsim_assert_expression_inputs(arms$grouped, bad); NULL },
                    error = identity)
    assert_true(inherits(err, "error"), paste("changed", field, "must fail"))
  }
})

run_test("Feature-group simulation: recovery uses signed match and posterior variance", {
  g <- c(0, 1, 2, 3, 4, 6)
  groups <- rep(c("rna", "meth"), each = 3)
  fit <- list(EF = cbind(-g, c(0, 1, 0, 1, 0, 1), rep(0, 6)),
              EBeta = c(-1.96, 3, 0), VBeta = c(1, 1, 1),
              EBeta2 = c(100, 100, 100))
  result <- fgsim_recovery(fit, g, groups, 1.96)
  assert_true(result$recovered)
  assert_equal(result$matched_program, 1L)
  assert_equal(result$false_active, 1L)
  fit$EBeta[1] <- 1.96
  wrong_sign <- fgsim_recovery(fit, g, groups, 1.96)
  assert_false(wrong_sign$recovered)
  assert_equal(wrong_sign$false_active, 2L)
  fit$EBeta[1] <- -1.95
  inactive <- fgsim_recovery(fit, g, groups, 1.96)
  assert_false(inactive$recovered)
  assert_equal(inactive$false_active, 1L)
  fit$EF[, 1] <- rep(0, 6)
  fit$EBeta[1] <- -1.96
  no_match <- fgsim_recovery(fit, g, groups, 1.96)
  assert_false(no_match$recovered)
  assert_equal(no_match$false_active, 2L)
})

run_test("Feature-group simulation: failed fits remain in gate denominator", {
  metrics <- data.frame(feature_ratio = rep(1, 10), methylation_active_fraction = rep(0.03, 10),
    arm = rep("grouped", 10), recovered = c(rep(TRUE, 8), FALSE, FALSE),
    false_active = c(rep(1, 9), NA), converged = c(rep(TRUE, 9), FALSE),
    c_validation = c(rep(0.6, 8), NA, NA),
    c_true_risk = c(rep(0.7, 9), NA), empty_blocks = c(rep(1L, 9), NA),
    error = c(rep("", 9), "failed"))
  summary <- fgsim_gate(metrics)
  assert_equal(nrow(summary), 1L)
  assert_equal(summary$recovered, 8L)
  assert_equal(summary$converged, 9L)
  assert_equal(summary$c_defined, 8L)
  assert_near(summary$mean_c, 0.6, 1e-12)
  assert_equal(summary$true_c_defined, 9L)
  assert_near(summary$mean_true_c, 0.7, 1e-12)
  assert_equal(summary$total_empty_blocks, 9L)
  assert_false(summary$passed)
  metrics$error[10] <- ""
  metrics$false_active[10] <- 1
  assert_true(fgsim_gate(metrics)$passed)
})

run_test("Feature-group simulation: selection follows registered stop and comparator rule", {
  summary <- data.frame(feature_ratio = rep(c(1, 5), each = 2),
    methylation_active_fraction = rep(0.03, 4), arm = rep(c("grouped", "pooled"), 2),
    passed = c(TRUE, TRUE, TRUE, FALSE))
  selected <- fgsim_select(summary)
  assert_equal(selected$selected_prior, "grouped")
  assert_true(selected$step4_allowed)
  summary$passed <- c(TRUE, TRUE, FALSE, TRUE)
  selected <- fgsim_select(summary)
  assert_equal(selected$selected_prior, "pooled")
  assert_true(selected$step4_allowed)
  summary$passed <- c(FALSE, TRUE, FALSE, TRUE)
  selected <- fgsim_select(summary)
  assert_equal(selected$selected_prior, "pooled")
  assert_true(selected$step4_allowed)
  summary$passed <- rep(FALSE, 4)
  assert_false(fgsim_select(summary)$step4_allowed)
})

run_test("Feature-group simulation: grouped fits expose K priors per modality", {
  prior <- list(family = "point_laplace", pi = 0.5, rate = 1, point_mass = FALSE)
  fit <- list(prior_F = list(expression = rep(list(prior), 5),
                             methylation = rep(list(prior), 5)))
  fgsim_assert_priors(fit, 5L)
  for (bad in list(list(expression = vector("list", 5), methylation = vector("list", 5)),
                   list(expression = vector("list", 5)),
                   list(expression = vector("list", 5), methylation = vector("list", 4)))) {
    fit$prior_F <- bad
    err <- tryCatch({ fgsim_assert_priors(fit, 5L); NULL }, error = identity)
    assert_true(inherits(err, "error"))
  }
})

run_test("Feature-group simulation: pooled prior observer preserves fitter arithmetic", {
  data <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 20, n_validation = 10, p_expression = 8, p_methylation = 8, seed = 3)
  x <- fgsim_prepare(data, 3)
  s <- svd(x$Y, nu = 2, nv = 2)
  x$EL_init <- s$u %*% diag(sqrt(s$d[1:2]))
  x$EF_init <- s$v %*% diag(sqrt(s$d[1:2]))
  cfg <- yaml::read_yaml("config/globals.yml")
  cfg$benchmark$k_working_yfb <- 2L
  cfg$cavi$max_iter <- 2L
  before <- canonicalize_yf_factors
  observed <- fgsim_fit(x, NULL, cfg)
  direct <- fit_ssbmf_yfb(x$Y, x$time, x$event, K = 2L, max_iter = 2L,
    tol = cfg$cavi$tol, prior_LF = "point_laplace", prior_beta = "normal",
    canonicalize = TRUE, w_S = 1, alpha = cfg$benchmark$alpha,
    EL_init = x$EL_init, EF_init = x$EF_init, verbose = FALSE)
  assert_equal(observed$fit, direct)
  assert_equal(canonicalize_yf_factors, before)
  assert_equal(length(observed$priors), 2L)
  assert_true(all(vapply(observed$priors, function(p) !is.null(p$pi), logical(1))))
})

if (sys.nframe() == 0L) {
  results <- report_results("Feature-group simulation")
  if (results$failed > 0L) quit(status = 1L)
}
