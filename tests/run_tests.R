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

# Source modules under test (add new modules here as they are created)
source("code/update_beta.R")

# List all test files to run
test_files <- c(
  "tests/test_update_beta.R"
  # future: "tests/test_update_L.R"
  # future: "tests/test_update_F.R"
  # future: "tests/test_update_tau.R"
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
