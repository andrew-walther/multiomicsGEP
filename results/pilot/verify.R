# ============================================================
# Script: verify.R
# Purpose: Check local pilot rebuild and release behavior in temporary copies.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: targets, yaml
# ============================================================

#' Exercise the pilot in isolated exported source trees
#' @return Invisibly returns the number of completed acceptance checks.
#' @examples
#' # Rscript results/pilot/verify.R
verify_pilot <- function() {
  root <- normalizePath(getwd())
  lib <- file.path(root, "results/pilot/library")
  .libPaths(c(normalizePath(lib, mustWork = FALSE), .libPaths()))
  old_lib <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  old_output <- Sys.getenv("SSBMF_PILOT_OUTPUT", unset = NA_character_)
  old_registry <- Sys.getenv("SSBMF_PILOT_REGISTRY", unset = NA_character_)
  scratch <- tempfile("pilot-acceptance-")
  dir.create(scratch)
  on.exit({
    setwd(root)
    if (is.na(old_lib)) Sys.unsetenv("R_LIBS_USER") else Sys.setenv(R_LIBS_USER = old_lib)
    if (is.na(old_output)) Sys.unsetenv("SSBMF_PILOT_OUTPUT") else Sys.setenv(SSBMF_PILOT_OUTPUT = old_output)
    if (is.na(old_registry)) Sys.unsetenv("SSBMF_PILOT_REGISTRY") else Sys.setenv(SSBMF_PILOT_REGISTRY = old_registry)
    unlink(scratch, recursive = TRUE)
  }, add = TRUE)
  Sys.unsetenv("SSBMF_PILOT_REGISTRY")
  Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
  Sys.setenv(SSBMF_PILOT_OUTPUT = file.path(scratch, "output"))
  dir.create(file.path(scratch, "config")); dir.create(file.path(scratch, "results"))
  stopifnot(file.copy(file.path(root, "config/globals.yml"), file.path(scratch, "config")),
            file.copy(file.path(root, "code"), scratch, recursive = TRUE),
            dir.create(file.path(scratch, "results/pilot")))
  pilot_files <- list.files(file.path(root, "results/pilot"), pattern = "\\.(R|csv|yml)$", full.names = TRUE)
  stopifnot(all(file.copy(pilot_files, file.path(scratch, "results/pilot"))))
  setwd(scratch)
  registry <- yaml::read_yaml("results/pilot/registry.yml")
  stopifnot(all(c("prediction", "input_audit") %in% names(registry$families)))
  first <- registry; first$families <- registry$families["prediction"]
  yaml::write_yaml(first, "results/pilot/registry.yml")
  output <- Sys.getenv("SSBMF_PILOT_OUTPUT")
  store <- file.path(output, "store")
  script <- "results/pilot/_targets.R"
  make <- function() targets::tar_make(script = script, store = store, envir = globalenv(), callr_function = NULL,
                                        reporter = "silent", use_crew = FALSE)
  stale <- function() targets::tar_outdated(script = script, store = store, envir = globalenv(),
                                            callr_function = NULL, reporter = "silent")
  fit_time <- function() targets::tar_meta(model_fit, time, store = store)$time
  catalog <- function() { source("results/pilot/functions.R"); pilot_catalog(output, stale()) }
  check <- function(ok, label) { if (!isTRUE(ok)) stop("Acceptance failed: ", label); cat("PASS ", label, "\n", sep = "") }
  append_line <- function(path, line) cat(line, "\n", file = path, append = TRUE)
  checks <- 0L
  make(); check(nrow(catalog()) == 1L && !length(stale()), "initial complete build"); checks <- checks + 1L
  initial <- read.csv(file.path(dirname(yaml::read_yaml(file.path(output, "latest.yml"))$manifest), "prediction.csv"))
  before <- fit_time(); make()
  check(identical(before, fit_time()), "unchanged build preserves fit timestamp"); checks <- checks + 1L
  # A distinct output store must reproduce the same aggregate values.
  second_output <- file.path(scratch, "independent-output")
  Sys.setenv(SSBMF_PILOT_OUTPUT = second_output); store <- file.path(second_output, "store")
  make()
  independent <- read.csv(file.path(dirname(yaml::read_yaml(file.path(second_output, "latest.yml"))$manifest), "prediction.csv"))
  check(identical(initial, independent), "independent stores reproduce aggregate"); checks <- checks + 1L
  Sys.setenv(SSBMF_PILOT_OUTPUT = output); store <- file.path(output, "store")
  functions <- readLines("results/pilot/functions.R")
  changed <- sub('main = family$title, ylab = "Fraction / concordance"', 'main = family$title, ylab = "Pilot fraction"', functions, fixed = TRUE)
  stopifnot(!identical(functions, changed)); writeLines(changed, "results/pilot/functions.R")
  check(!"model_fit" %in% stale(), "figure-only edit leaves fit fresh"); checks <- checks + 1L
  make(); check(identical(before, fit_time()), "figure-only build preserves fit timestamp"); checks <- checks + 1L
  append_line("code/fit_ssbmf_yfb.R", "# acceptance dependency change")
  check("model_fit" %in% stale(), "model dependency invalidates fit"); checks <- checks + 1L
  make()
  cfg <- yaml::read_yaml("config/globals.yml"); cfg$results_pilot$fixture$seed <- cfg$results_pilot$fixture$seed + 1L
  yaml::write_yaml(cfg, "config/globals.yml")
  check(all(c("simulation", "model_fit") %in% stale()), "fixture change invalidates simulation and fit"); checks <- checks + 1L
  make()
  writeLines("Temporary documentation edit", "README.md")
  check(!"model_fit" %in% stale(), "documentation edit leaves fit fresh"); checks <- checks + 1L
  first_functions <- readLines("results/pilot/functions.R")
  yaml::write_yaml(registry, "results/pilot/registry.yml")
  check(identical(first_functions, readLines("results/pilot/functions.R")), "second family requires no first-family change"); checks <- checks + 1L
  make(); check(setequal(catalog()$family, c("prediction", "input_audit")), "registry adds second family"); checks <- checks + 1L
  pointer <- yaml::read_yaml(file.path(output, "latest.yml"))
  artifact_names <- names(yaml::read_yaml(pointer$manifest)$artifact_md5)
  figure_name <- artifact_names[grepl("^prediction\\.(png|pdf)$", artifact_names)]
  stopifnot(length(figure_name) == 1L)
  artifact <- file.path(dirname(pointer$manifest), figure_name)
  stopifnot(file.remove(artifact))
  check(length(stale()) > 0L && inherits(try(catalog(), silent = TRUE), "try-error"), "missing artifact blocks catalog"); checks <- checks + 1L
  make(); check(file.exists(artifact) && nrow(catalog()) == 2L, "missing artifact repaired"); checks <- checks + 1L
  previous <- readLines(file.path(output, "latest.yml"))
  broken <- registry; broken$families$input_audit$command <- 'stop("injected family failure")'
  yaml::write_yaml(broken, "results/pilot/registry.yml")
  failed <- try(make(), silent = TRUE)
  check(inherits(failed, "try-error") && identical(previous, readLines(file.path(output, "latest.yml"))) &&
          inherits(try(catalog(), silent = TRUE), "try-error"), "family failure preserves pointer and blocks catalog"); checks <- checks + 1L
  yaml::write_yaml(registry, "results/pilot/registry.yml"); make()
  before_cli <- fit_time()
  cli <- system2(file.path(R.home("bin"), "Rscript"),
                 c("results/pilot/run.R", "catalog"), stdout = TRUE, stderr = TRUE)
  check(is.null(attr(cli, "status")) && any(grepl("prediction", cli)) &&
          any(grepl("input_audit", cli)) && identical(before_cli, fit_time()),
        "catalog CLI lists both families without fitting"); checks <- checks + 1L
  invisible(checks)
}

start <- proc.time()[["elapsed"]]
tryCatch({
  count <- verify_pilot()
  cat(sprintf("Completed %d checks in %.1f seconds.\n", count, proc.time()[["elapsed"]] - start))
}, error = function(e) {
  message("Pilot acceptance failed: ", conditionMessage(e))
  quit(status = 1L)
})
