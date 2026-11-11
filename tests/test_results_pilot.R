# ============================================================
# Script: test_results_pilot.R
# Purpose: Test deterministic pilot inputs and safe artifact/catalog handling.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: base R, yaml, digest; no targets installation or fit required
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
source("results/pilot/functions.R")

run_test("Pilot1: synthetic fixture is deterministic and uses training scaling", {
  spec <- pilot_settings("config/globals.yml")$fixture
  a <- pilot_simulate(spec); b <- pilot_simulate(spec)
  assert_equal(a, b)
  assert_near(colMeans(a$train), rep(0, ncol(a$train)), tol = 1e-12)
  assert_near(apply(a$train, 2, sd), rep(1, ncol(a$train)), tol = 1e-12)
  assert_length(a$time_test, spec$n_total - spec$n_train)
  assert_true(max(abs(colMeans(a$test))) > 1e-4,
              "Held-out rows must not have been centered on their own means")
})

run_test("Pilot2: the input family reports actual dimensions and event fraction", {
  x <- pilot_simulate(pilot_settings("config/globals.yml")$fixture)
  values <- with(pilot_input_summary(x), setNames(value, metric))
  assert_near(values[["training_event_fraction"]], mean(x$status_train))
  assert_equal(values[["n_train"]], as.numeric(nrow(x$train)))
  assert_equal(values[["n_test"]], as.numeric(nrow(x$test)))
})

run_test("Pilot3: package lock preserves raw hyphenated package versions", local({
  lock <- tempfile(fileext = ".csv"); on.exit(unlink(lock))
  version <- utils::packageDescription("survival", fields = "Version")
  write.csv(data.frame(Package = "survival", Version = version), lock, row.names = FALSE)
  assert_equal(pilot_runtime(lock)$packages$survival, version)
  write.csv(data.frame(Package = "survival", Version = "0.0-invalid"), lock, row.names = FALSE)
  err <- tryCatch(pilot_runtime(lock), error = identity)
  assert_true(inherits(err, "error"))
}))

run_test("Pilot4: immutable writes reject overwriting different content", local({
  folder <- tempfile(); dir.create(folder); on.exit(unlink(folder, recursive = TRUE))
  dest <- file.path(folder, "artifact")
  a <- tempfile(tmpdir = folder); writeLines("original", a); pilot_finalize(a, dest)
  b <- tempfile(tmpdir = folder); writeLines("replacement", b)
  err <- tryCatch(pilot_finalize(b, dest), error = identity)
  assert_true(inherits(err, "error"))
  assert_equal(readLines(dest), "original")
}))

run_test("Pilot5: artifacts resolve to a complete pilot-only manifest", local({
  output <- tempfile(); dir.create(output); on.exit(unlink(output, recursive = TRUE))
  family <- list(title = "Software test", figure_metric = "fraction")
  context <- list(run_id = "unit-example", source = list(commit = "unit", dirty = TRUE),
                  scientific_release = FALSE)
  files <- pilot_artifacts(data.frame(metric = "fraction", value = 0.5),
                            family, "example", context, output)
  assert_true(all(file.exists(files)))
  meta <- yaml::read_yaml(files[3])
  assert_equal(meta$run_id, context$run_id)
  release <- pilot_release(list(example = files), context, list(example = family), output)
  assert_true(all(file.exists(release)))
  manifest <- yaml::read_yaml(release[1])
  assert_false(manifest$scientific_release)
  catalog <- pilot_catalog(output, character())
  assert_equal(catalog$family, "example")
  assert_equal(catalog$run_id, context$run_id)
}))

run_test("Pilot6: stale and tampered artifacts cannot enter the fresh catalog", local({
  output <- tempfile(); dir.create(output); on.exit(unlink(output, recursive = TRUE))
  family <- list(title = "Software test", figure_metric = "fraction")
  context <- list(run_id = "unit-example", scientific_release = FALSE)
  files <- pilot_artifacts(data.frame(metric = "fraction", value = 0.5),
                            family, "example", context, output)
  pilot_release(list(example = files), context, list(example = family), output)
  err <- tryCatch(pilot_catalog(output, "model_fit"), error = identity)
  assert_true(inherits(err, "error"))
  writeLines("corrupt", files[1])
  err <- tryCatch(pilot_catalog(output, character()), error = identity)
  assert_true(inherits(err, "error"))
}))

run_test("Pilot7: incomplete release preserves the previous successful pointer", local({
  output <- tempfile(); dir.create(output); on.exit(unlink(output, recursive = TRUE))
  pointer <- file.path(output, "latest.yml"); writeLines("old-success", pointer)
  err <- tryCatch(pilot_release(list(example = "missing.csv"),
    list(run_id = "failed"), list(), output), error = identity)
  assert_true(inherits(err, "error"))
  assert_equal(readLines(pointer), "old-success")
}))

report_results("test_results_pilot.R")
if (sys.nframe() == 0L && .test_env$failed > 0L) quit(status = 1L)
