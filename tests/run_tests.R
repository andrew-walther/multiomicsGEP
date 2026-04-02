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
source("code/compute_elbo.R")   # compute_ebnm_kl, compute_survival_elbo
# fit_modular.R resets DATA_MODE <- "real" and has a runner block at the bottom
# that errors when real_Y is NULL.  Wrap in tryCatch: the function definition
# (lines 1-395) completes before the runner block fires, so the error is
# harmless — fit_supervised_mf_modular() is available after this call.
suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))

# List all test files to run
test_files <- c(
  "tests/test_update_beta.R",
  "tests/test_update_L.R",
  "tests/test_update_F.R",
  "tests/test_update_tau.R",
  "tests/test_predict.R",
  "tests/test_elbo.R"
)

# Run each test file
for (tf in test_files) {
  cat(sprintf("--- Running: %s ---\n", tf))
  reset_counts()
  source(tf, local = FALSE)
  report_results(tf)
}

# Final summary across all files
cat("============================================================\n")
cat(sprintf(" FINAL: %d passed, %d failed\n",
            .test_env$passed, .test_env$failed))
if (.test_env$failed > 0) {
  cat("\nFailed tests:\n")
  for (e in .test_env$errors) cat(" *", e, "\n")
  quit(status = 1)
} else {
  cat(" All tests PASSED.\n")
}
cat("============================================================\n")
