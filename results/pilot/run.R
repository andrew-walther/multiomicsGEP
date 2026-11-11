# ============================================================
# Script: run.R
# Purpose: Inspect or explicitly build the small local results pilot.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: targets and the locked pilot library
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
command <- if (length(args)) args[[1]] else "status"
if (!command %in% c("status", "build", "catalog"))
  stop("Usage: Rscript results/pilot/run.R [status|build|catalog]")
if (dir.exists("results/pilot/library")) .libPaths(c("results/pilot/library", .libPaths()))
if (!requireNamespace("targets", quietly = TRUE))
  stop("Install dependencies first: Rscript results/pilot/install.R")
source("results/pilot/functions.R")
output <- Sys.getenv("SSBMF_PILOT_OUTPUT", "results/pilot/local")
store <- file.path(output, "store")
script <- "results/pilot/_targets.R"
if (command == "build") {
  # Serial local controller with an isolated fit process. No SLURM or real data.
  targets::tar_make(script = script, store = store, envir = globalenv(), callr_function = NULL,
                     reporter = "terse", use_crew = FALSE)
} else {
  outdated <- targets::tar_outdated(script = script, store = store, envir = globalenv(),
                                     callr_function = NULL, reporter = "silent")
  if (command == "catalog") {
    print(pilot_catalog(output, outdated), row.names = FALSE)
    cat("Artifacts under:", file.path(output, "runs"), "\n")
  } else if (length(outdated)) {
    cat("Stale or missing pilot targets (no fits run):\n", paste(outdated, collapse = "\n"), "\n")
  } else cat("Pilot is fresh. No fits run.\n")
}
