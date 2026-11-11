# ============================================================
# Script: _targets.R
# Purpose: Declare the local synthetic pilot dependency graph.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: targets, yaml, digest, ebnm, survival, Rcpp
# ============================================================

source("results/pilot/functions.R")

targets::tar_option_set(seed = 1L, error = "stop", memory = "transient")
pilot_output <- Sys.getenv("SSBMF_PILOT_OUTPUT", "results/pilot/local")
pilot_registry_path <- Sys.getenv("SSBMF_PILOT_REGISTRY", "results/pilot/registry.yml")
pilot_registry <- yaml::read_yaml(pilot_registry_path)$families
if (!length(pilot_registry) || any(!grepl("^[a-z][a-z0-9_]*$", names(pilot_registry))))
  stop("Registry needs at least one family with a lower-case identifier.")

# Cheap global fingerprints are recalculated when the graph is read, even for
# status. They do not fit models or write artifacts, unlike an always-run target.
pilot_runtime_state <- pilot_runtime("results/pilot/package-versions.csv")
pilot_checkout_state <- pilot_source_state()
pilot_source_paths <- c(list.files("results/pilot", pattern = "\\.(R|csv)$", full.names = TRUE),
                        pilot_registry_path, "config/globals.yml")
pilot_source_hashes <- tools::md5sum(pilot_source_paths)

pipeline <- list(
  targets::tar_target(config_file, "config/globals.yml", format = "file"),
  targets::tar_target(settings, pilot_settings(config_file)),
  targets::tar_target(fixture_settings, settings$fixture),
  targets::tar_target(model_settings, settings$model),
  # Conservative pilot boundary: all model R/C++ source files. Later migration
  # may narrow this after auditing sourced helpers and compiled dependencies.
  targets::tar_target(model_files,
    list.files("code", pattern = "\\.(R|cpp)$", full.names = TRUE, recursive = TRUE), format = "file"),
  targets::tar_target(runtime, pilot_runtime_state),
  targets::tar_target(simulation, pilot_simulate(fixture_settings)),
  targets::tar_target(model_fit, pilot_fit(simulation, model_settings, model_files, runtime)),
  targets::tar_target(run_context, {
    pilot_source_hashes
    pilot_context(c(model_files, pilot_source_paths), pilot_checkout_state, runtime, settings)
  })
)

# Family commands are trusted version-controlled R expressions. Adding a producer
# plus a registry entry adds its targets without modifying earlier family code.
for (id in names(pilot_registry)) {
  family <- pilot_registry[[id]]
  if (length(family$command) != 1L || !is.character(family$command)) stop("Invalid family command.")
  family_command <- parse(text = family$command)
  if (length(family_command) != 1L) stop("Each family command must be one R expression.")
  result_name <- paste0(id, "_summary")
  pipeline <- c(pipeline, list(
    targets::tar_target_raw(result_name, family_command[[1]]),
    targets::tar_target_raw(paste0(id, "_artifacts"), substitute(
      pilot_artifacts(DATA, FAMILY, ID, run_context, pilot_output),
      list(DATA = as.name(result_name), FAMILY = family, ID = id)), format = "file")
  ))
}
artifact_names <- paste0(names(pilot_registry), "_artifacts")
artifact_call <- as.call(c(list(as.name("list")),
  setNames(lapply(artifact_names, as.name), names(pilot_registry))))
pipeline <- c(pipeline, list(targets::tar_target_raw("pilot_release_files", substitute(
  pilot_release(ARTIFACTS, run_context, REGISTRY, pilot_output),
  list(ARTIFACTS = artifact_call, REGISTRY = pilot_registry)), format = "file")))
pipeline
