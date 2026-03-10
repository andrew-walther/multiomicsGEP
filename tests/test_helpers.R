# =============================================================================
# tests/test_helpers.R
# Lightweight test assertion framework for multiomicsGEP
#
# Usage:
#   source("tests/test_helpers.R")
#   run_test("My test", {
#     assert_near(actual, expected, tol = 1e-6)
#   })
#
# Pattern for future update modules: test_update_L.R, test_update_F.R, etc.
# will all source this file.
# =============================================================================

.test_env <- new.env(parent = emptyenv())
.test_env$passed <- 0L
.test_env$failed  <- 0L
.test_env$errors  <- character(0)

# -----------------------------------------------------------------------------
# Core assertion helpers
# -----------------------------------------------------------------------------

assert_true <- function(condition, msg = "") {
  if (!isTRUE(condition)) {
    stop(paste0("Assertion failed", if (nchar(msg) > 0) paste0(": ", msg) else ""))
  }
  invisible(NULL)
}

assert_false <- function(condition, msg = "") {
  assert_true(!condition, msg)
}

assert_near <- function(actual, expected, tol = 1e-6, msg = "") {
  diff <- max(abs(actual - expected))
  if (diff > tol) {
    stop(sprintf("Values differ by %.3e > tol %.3e%s",
                 diff, tol,
                 if (nchar(msg) > 0) paste0(": ", msg) else ""))
  }
  invisible(NULL)
}

assert_equal <- function(actual, expected, msg = "") {
  if (!identical(actual, expected)) {
    stop(paste0("Not identical",
                if (nchar(msg) > 0) paste0(": ", msg) else "",
                "\n  actual:   ", deparse(actual),
                "\n  expected: ", deparse(expected)))
  }
  invisible(NULL)
}

assert_length <- function(x, n, msg = "") {
  if (length(x) != n) {
    stop(sprintf("Expected length %d, got %d%s", n, length(x),
                 if (nchar(msg) > 0) paste0(": ", msg) else ""))
  }
  invisible(NULL)
}

assert_finite <- function(x, msg = "") {
  if (any(!is.finite(x))) {
    stop(paste0("Non-finite values found",
                if (nchar(msg) > 0) paste0(": ", msg) else ""))
  }
  invisible(NULL)
}

assert_positive <- function(x, msg = "") {
  if (any(x <= 0)) {
    stop(paste0("Non-positive values found",
                if (nchar(msg) > 0) paste0(": ", msg) else ""))
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Test runner
# -----------------------------------------------------------------------------

run_test <- function(name, expr) {
  result <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) e)

  if (isTRUE(result)) {
    .test_env$passed <- .test_env$passed + 1L
    cat(sprintf("  \033[32m[PASS]\033[0m %s\n", name))
  } else {
    .test_env$failed <- .test_env$failed + 1L
    msg <- conditionMessage(result)
    .test_env$errors <- c(.test_env$errors, sprintf("%s: %s", name, msg))
    cat(sprintf("  \033[31m[FAIL]\033[0m %s\n         %s\n", name, msg))
  }
  invisible(NULL)
}

# Call at end of a test file section
report_results <- function(section = "") {
  total <- .test_env$passed + .test_env$failed
  cat(sprintf("\n--- %s: %d/%d tests passed ---\n\n",
              if (nchar(section) > 0) section else "Results",
              .test_env$passed, total))
  invisible(list(passed = .test_env$passed, failed = .test_env$failed,
                 errors = .test_env$errors))
}

reset_counts <- function() {
  .test_env$passed <- 0L
  .test_env$failed  <- 0L
  .test_env$errors  <- character(0)
  invisible(NULL)
}
