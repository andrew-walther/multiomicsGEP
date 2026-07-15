# =============================================================================
# tests/run_tests.R
# Master test runner for multiomicsGEP
#
# Usage (from project root):
#   Rscript tests/run_tests.R
#
# Add new test files to the `test_files` vector below.
# =============================================================================

cat("============================================================\n")
cat(" multiomicsGEP Test Suite\n")
cat("============================================================\n\n")

# Source shared test infrastructure
source("tests/test_helpers.R")

# Source modules under test
source("code/update_beta.R")
source("code/update_L.R")
source("code/update_F.R")
source("code/update_tau.R")
source("code/compute_elbo.R")   # compute_ebnm_kl, compute_survival_elbo, compute_normal_kl
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R")
# fit_modular.R resets DATA_MODE <- "real" and has a runner block at the bottom
# that errors when real_Y is NULL.  Wrap in tryCatch: the function definition
# (lines 1-395) completes before the runner block fires, so the error is
# harmless — fit_supervised_mf_modular() is available after this call.
suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
# fit_cox_on_yf.R has the same runner-block pattern as fit_modular.R
suppressMessages(tryCatch(
  source("code/fit_cox_on_yf.R"),
  error = function(e) invisible(NULL)
))
source("code/predict_cox_on_yf.R")
source("code/select_alpha_cv.R")
source("code/preprocess_desurv.R")
source("results/multi_cohort_sim/fit_pca_cox.R")
source("code/pathway_enrichment.R")

# List all test files to run
test_files <- c(
  "tests/test_update_beta.R",
  "tests/test_update_L.R",
  "tests/test_update_F.R",
  "tests/test_update_tau.R",
  "tests/test_predict.R",
  "tests/test_elbo.R",
  "tests/test_select_alpha_cv.R",
  "tests/test_preprocess_desurv.R",
  "tests/test_multistart.R",
  "tests/test_yfb_multistart.R",
  "tests/test_warmstart_from_fit.R",
  "tests/test_deflation_init.R",
  "tests/test_select_K_cv.R",
  "tests/test_select_k_alpha_bo.R",
  "tests/test_update_F_cohort.R",
  "tests/test_fit_modular_cohort.R",
  "tests/test_fit_yf_cohort.R",
  "tests/test_fit_yf_frozen_f.R",
  "tests/test_normalization.R",
  "tests/test_fit_pca_cox.R",
  "tests/test_pathway_enrichment.R"
)

# Run each test file
total_passed <- 0L
total_failed <- 0L
all_errors <- character(0)

for (tf in test_files) {
  cat(sprintf("--- Running: %s ---\n", tf))
  reset_counts()
  source(tf, local = FALSE)
  file_results <- report_results(tf)
  total_passed <- total_passed + file_results$passed
  total_failed <- total_failed + file_results$failed
  all_errors   <- c(all_errors, file_results$errors)
}

# Final summary across all files
cat("============================================================\n")
cat(sprintf(" FINAL: %d passed, %d failed\n",
            total_passed, total_failed))
if (total_failed > 0) {
  cat("\nFailed tests:\n")
  for (e in all_errors) cat(" *", e, "\n")
  quit(status = 1)
} else {
  cat(" All tests PASSED.\n")
}
cat("============================================================\n")
