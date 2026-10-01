# ============================================================
# Script: test_simulate_multimodal_yfb.R
# Purpose: Tests for the matched multimodal YFB simulation data contract.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R; survival; code/simulate_multimodal_yfb.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")
if (!exists("run_multimodal_yfb_replicate")) source("code/run_multimodal_yfb_simulation.R")

run_test("MMYFB-Sim-T1: scenarios generate matched nonnegative train/validation blocks", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 30,
                                        n_validation = 20, p_expression = 12,
                                        p_methylation = 10, seed = 23)
  assert_true(identical(names(data$training$Y), c("expression", "methylation")))
  assert_true(all(vapply(data$training$Y, function(x) all(is.finite(x) & x >= 0), logical(1))))
  assert_true(!any(rownames(data$training$Y$expression) %in%
                     rownames(data$validation$Y$expression)))
  assert_true(identical(colnames(data$training$Y$methylation),
                        colnames(data$validation$Y$methylation)))
})

run_test("MMYFB-Sim-T2: survival uses the stated raw joint predictor", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 20,
                                        n_validation = 20, p_expression = 10,
                                        p_methylation = 10, seed = 24)
  eta <- as.vector((data$training$Y$expression %*% data$truth$F$expression +
    data$training$Y$methylation %*% data$truth$F$methylation) %*% data$truth$beta)
  assert_near(eta, data$training$eta, tol = 1e-12)
})

run_test("MMYFB-Sim-T3: scenario-specific modality and feature designs apply", {
  expression_only <- simulate_multimodal_yfb_data("expression_only", n_train = 20,
    n_validation = 20, p_expression = 10, p_methylation = 10, seed = 25)
  unequal <- simulate_multimodal_yfb_data("unequal_methylation_features", n_train = 20,
    n_validation = 20, p_expression = 10, p_methylation = 10, seed = 26)
  assert_true(all(expression_only$truth$F$methylation[, expression_only$truth$beta != 0] == 0))
  assert_true(ncol(unequal$training$Y$methylation) == 75)
})

run_test("MMYFB-Sim-T4: factor matching recovers an exact permutation", {
  truth <- cbind(c(1, 0, 0), c(0, 1, 0))
  estimate <- truth[, 2:1]
  match <- multimodal_yfb_match_factors(truth, estimate)
  assert_true(all(abs(match$abs_correlation - 1) < 1e-12))
})

run_test("MMYFB-Sim-T5: all approved scenario names generate finite data", {
  for (scenario in names(multimodal_yfb_scenarios())) {
    data <- simulate_multimodal_yfb_data(scenario, n_train = 12, n_validation = 12,
                                          p_expression = 10, p_methylation = 10, seed = 30)
    assert_true(all(is.finite(data$training$time)))
    assert_true(all(data$training$event %in% 0:1))
  }
})

run_test("MMYFB-Sim-T5b: adverse-protective scenario contains both known directions", {
  scenario <- multimodal_yfb_scenarios()$adverse_protective
  assert_true(any(scenario$beta > 0))
  assert_true(any(scenario$beta < 0))
  assert_true(any(scenario$beta == 0))
})

run_test("MMYFB-Sim-T6: one replicate reports every required benchmark arm", {
  result <- run_multimodal_yfb_replicate(
    "both_informative", replicate_id = 1, seed = 31, n_train = 30,
    n_validation = 30, p_expression = 10, p_methylation = 10, K_init = 3,
    control = list(max_outer = 1L, max_inner = 1L)
  )
  expected <- c("joint_multimodal_yfb", "two_step_ebmf_cox", "expression_only_yfb",
                "methylation_only_yfb", "training_only_stacked_risk_cox")
  assert_true(identical(result$method, expected))
  assert_true(all(c("c_index", "heldout_partial_log_likelihood",
                    "true_risk_correlation", "runtime_seconds", "memory_mb") %in% names(result)))
  assert_true("prognostic_direction_recall" %in% names(result))
  joint <- result[result$method == "joint_multimodal_yfb", , drop = FALSE]
  assert_true(all(c("K_eff_reconstruction", "K_eff_survival",
                    "K_eff_retained") %in% names(joint)))
  assert_true(joint$K_eff_retained >= joint$K_eff_reconstruction)
  assert_true(joint$K_eff_retained >= joint$K_eff_survival)
})

run_test("MMYFB-Sim-T7: expression-only YFB baseline avoids a degenerate block", {
  data <- simulate_multimodal_yfb_data("expression_only", seed = 20260921L)
  fit <- fit_multimodal_yfb_single_modality(
    data$training$Y, data$training$time, data$training$event, "expression", K = 5,
    control = list(max_outer = 10L, max_inner = 2L)
  )
  prediction <- predict_multimodal_yfb_single_modality(fit, data$validation$Y)
  assert_true(all(is.finite(prediction$risk_scores)))
})

run_test("MMYFB-Sim-T8: null-risk correlation is explicitly undefined", {
  assert_true(is.na(multimodal_yfb_safe_correlation(c(1, 1), c(0, 0))))
})

run_test("MMYFB-Sim-T9: K sensitivity preserves every requested K", {
  output <- tempfile(fileext = ".csv")
  result <- run_multimodal_yfb_k_sensitivity(
    output_path = output, K_values = c(3L, 5L), replicates = 1L,
    scenarios = "both_informative", n_train = 20L, n_validation = 20L,
    p_expression = 10L, p_methylation = 10L,
    control = list(max_outer = 1L, max_inner = 1L)
  )
  assert_true(file.exists(output))
  assert_true(nrow(result) == 10L)
  assert_true(identical(sort(unique(result$K_init)), c(3L, 5L)))
})

run_test("MMYFB-Sim-T10: K selection chooses the smallest K on a one-SE C-index plateau", {
  results <- data.frame(
    method = rep(c("joint_multimodal_yfb", "expression_only_yfb"), each = 6),
    K_init = rep(rep(c(7L, 9L, 11L), each = 2), 2),
    c_index = c(0.72, 0.72, 0.73, 0.71, 0.74, 0.72,
                0.55, 0.56, 0.57, 0.58, 0.59, 0.60),
    heldout_partial_log_likelihood = seq_len(12),
    converged = TRUE
  )
  selected <- summarize_multimodal_yfb_k_selection(results)
  assert_true(selected$selected_K == 7L)
  assert_true(all(selected$summary$method == "joint_multimodal_yfb"))
})

run_test("MMYFB-Sim-T11: joint-only K-selection runner preserves paired seeds and K values", {
  path <- tempfile(fileext = ".csv")
  result <- run_multimodal_yfb_k_selection(
    output_path = path, K_values = c(3L, 4L), replicates = 1L,
    scenario = "adverse_protective", n_train = 30L, n_validation = 30L,
    p_expression = 10L, p_methylation = 10L,
    control = list(max_outer = 1L, max_inner = 1L)
  )
  assert_true(file.exists(path))
  assert_true(identical(sort(unique(result$K_init)), c(3L, 4L)))
  assert_true(length(unique(result$seed)) == 1L)
  assert_true(all(result$method == "joint_multimodal_yfb"))
  assert_true(all(c("prognostic_recall", "prognostic_direction_recall") %in%
    names(result)))
})

if (sys.nframe() == 0L) report_results("test_simulate_multimodal_yfb.R")
