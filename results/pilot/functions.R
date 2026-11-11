# ============================================================
# Script: functions.R
# Purpose: Small synthetic result families and versioned pilot artifacts.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: targets, callr, yaml, digest, ebnm, survival, Rcpp
# ============================================================

# Configuration and provenance ----

#' Read smoke settings without changing the registered research settings
#' @param path Path to the global YAML configuration.
#' @return List with separate fixture and model settings.
#' @examples
#' pilot_settings("config/globals.yml")
#' @family results-pilot
pilot_settings <- function(path) {
  cfg <- yaml::read_yaml(path)
  x <- cfg$results_pilot
  x$model$K <- cfg$cavi$k_default
  x$model$max_iter <- cfg$cavi$max_iter
  x$model$tol <- cfg$cavi$tol
  x
}

#' Verify the exact tested R package versions and report the runtime
#' @param lock CSV listing Package and Version columns.
#' @return List containing R/platform identifiers and named package versions.
#' @examples
#' pilot_runtime("results/pilot/package-versions.csv")
#' @family results-pilot
pilot_runtime <- function(lock) {
  pins <- read.csv(lock, stringsAsFactors = FALSE)
  versions <- vapply(pins$Package, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
    utils::packageDescription(pkg, fields = "Version")
  }, character(1))
  if (any(versions != pins$Version))
    stop("Package versions differ from the pilot lock. Run Rscript results/pilot/install.R.")
  list(R = R.version.string, platform = R.version$platform,
       packages = as.list(setNames(versions, pins$Package)))
}

#' Identify the source checkout, including exploratory dirty state
#' @return List with the source commit and dirty flag, or an exported-source label.
#' @examples
#' pilot_source_state()
#' @family results-pilot
pilot_source_state <- function() {
  commit <- suppressWarnings(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(commit, "status")) || !length(commit))
    return(list(commit = "exported-source", dirty = TRUE))
  status <- system2("git", c("status", "--porcelain", "--untracked-files=normal"), stdout = TRUE)
  list(commit = commit[[1]], dirty = length(status) > 0L)
}

#' Hash declared source/configuration files for an immutable pilot run
#' @param files Existing source or configuration file paths.
#' @param source_state Commit and dirty-state list.
#' @param runtime Verified runtime list.
#' @param settings Complete smoke settings.
#' @return Provenance list with a content-derived run ID.
#' @examples
#' # context <- pilot_context(files, pilot_source_state(), runtime, settings)
#' @family results-pilot
pilot_context <- function(files, source_state, runtime, settings) {
  hashes <- as.list(setNames(unname(tools::md5sum(files)), files))
  x <- list(source = source_state, runtime = runtime, settings = settings,
            source_md5 = hashes, scope = "synthetic software pilot only",
            scientific_release = FALSE)
  x$run_id <- paste0("pilot-", substr(source_state$commit, 1, 8), "-",
                     substr(digest::digest(x, algo = "sha256"), 1, 16))
  x
}

# Synthetic families ----

#' Generate the README fixture with training-only standardization
#' @param spec Fixture settings from config/globals.yml.
#' @return List with standardized train/test matrices and survival vectors.
#' @examples
#' # x <- pilot_simulate(pilot_settings("config/globals.yml")$fixture)
#' @family results-pilot
pilot_simulate <- function(spec) {
  set.seed(spec$seed)
  n <- spec$n_total; p <- spec$p; K <- spec$k_true
  L <- matrix(rexp(n * K), n, K)
  F <- matrix(rexp(p * K) * rbinom(p * K, 1, spec$loading_probability), p, K)
  Y <- L %*% t(F) + matrix(rnorm(n * p, sd = spec$noise_sd), n, p)
  time <- rexp(n, exp(drop(L %*% unlist(spec$beta))))
  status <- rbinom(n, 1, spec$event_probability)
  tr <- seq_len(spec$n_train); te <- setdiff(seq_len(n), tr)
  train <- scale(Y[tr, , drop = FALSE])
  test <- scale(Y[te, , drop = FALSE], center = attr(train, "scaled:center"),
                scale = attr(train, "scaled:scale"))
  stopifnot(all(is.finite(train)), all(is.finite(test)))
  list(train = train, test = test, time_train = time[tr], status_train = status[tr],
       time_test = time[te], status_test = status[te])
}

#' Run the canonical fit in a clean local process
#' @param simulation Synthetic fixture list.
#' @param settings Explicit fitter settings.
#' @return Canonical fit with recorded warning messages.
#' @examples
#' # Called by pilot_fit(), not used as a new fitting algorithm.
#' @family results-pilot
#' @seealso pilot_fit
pilot_fit_process <- function(simulation, settings) {
  source("code/fit_ssbmf_yfb.R")
  messages <- character()
  fit <- withCallingHandlers(do.call(fit_ssbmf_yfb, c(list(Y = simulation$train,
    time = simulation$time_train, status = simulation$status_train, verbose = FALSE), settings)),
    warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
  fit$pilot_warnings <- messages
  fit
}

#' Fit only the canonical model on the software fixture, failing on a cap
#' @param simulation Synthetic fixture list.
#' @param settings Explicit fitter settings.
#' @param model_files Declared R/C++ model source dependencies.
#' @param runtime Verified package/runtime fingerprint, an explicit dependency.
#' @return Canonical fit object. No capped or nonfinite fit is returned.
#' @examples
#' # fit <- pilot_fit(simulation, settings$model, model_files, runtime)
#' @family results-pilot
#' @seealso fit_ssbmf_yfb
pilot_fit <- function(simulation, settings, model_files, runtime) {
  stopifnot(length(model_files) > 0L, length(runtime$packages) > 0L)
  # The canonical fitter compiles helpers and mutates its R session. Isolate those
  # effects from targets' graph environment while running only one job at a time.
  fit <- callr::r(pilot_fit_process, args = list(simulation, settings),
                  libpath = .libPaths(), system_profile = FALSE, user_profile = FALSE)
  for (message in fit$pilot_warnings) warning(message, call. = FALSE)
  if (!isTRUE(fit$history$converged)) stop("Pilot fit did not converge; no completed run produced.")
  if (any(!is.finite(c(fit$EF, fit$EBeta, fit$EL)))) stop("Nonfinite pilot fit.")
  fit
}

#' Produce aggregate checks from the fitted model and held-out synthetic rows
#' @param fit Canonical YFB fit.
#' @param simulation Fixture list with training and held-out data.
#' @return Data frame with metric and value columns.
#' @examples
#' # table <- pilot_prediction_summary(fit, simulation)
#' @family results-pilot
#' @seealso predict_ssbmf_yfb, frozen_reverse_cindex
pilot_prediction_summary <- function(fit, simulation) {
  scoring <- new.env(parent = globalenv())
  sys.source("code/predict_ssbmf_yfb.R", envir = scoring)
  sys.source("code/concordance_ci.R", envir = scoring)
  risk <- scoring$predict_ssbmf_yfb(simulation$test, fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
  cindex <- scoring$frozen_reverse_cindex(risk, simulation$time_test, simulation$status_test)
  rmse <- sqrt(mean((simulation$train - fit$EL %*% t(fit$EF))^2))
  if (any(!is.finite(c(risk, cindex, rmse)))) stop("Nonfinite pilot scoring output.")
  data.frame(metric = c("heldout_concordance", "reconstruction_rmse", "iterations"),
             value = c(cindex, rmse, fit$history$n_iter))
}

#' Produce a second result family that depends on inputs but not the fitter
#' @param simulation Fixture list.
#' @return Data frame with metric and value columns.
#' @examples
#' # table <- pilot_input_summary(simulation)
#' @family results-pilot
pilot_input_summary <- function(simulation) {
  data.frame(metric = c("training_event_fraction", "n_train", "n_test", "n_genes"),
    value = c(mean(simulation$status_train), nrow(simulation$train),
              nrow(simulation$test), ncol(simulation$train)))
}

# Immutable artifacts and catalog ----

#' Finalize an artifact without overwriting different historical content
#' @param temporary Completed temporary file.
#' @param destination Final artifact path.
#' @return Destination path, invisibly.
#' @examples
#' # pilot_finalize(temp_file, "results/pilot/local/runs/example/table.csv")
#' @family results-pilot
pilot_finalize <- function(temporary, destination) {
  on.exit(unlink(temporary), add = TRUE)
  if (file.exists(destination)) {
    if (!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(destination))))
      stop("Refusing to overwrite different content in immutable pilot artifact: ", destination)
  } else if (!file.rename(temporary, destination)) stop("Could not finalize ", destination)
  invisible(destination)
}

#' Write one family's aggregate table, figure and adjacent provenance
#' @param data Metric/value aggregate table.
#' @param family Registry entry with title and figure_metric.
#' @param id Family identifier.
#' @param context Run provenance.
#' @param output Local pilot output directory.
#' @return Character vector naming the CSV, PDF and YAML sidecar file targets.
#' @examples
#' # paths <- pilot_artifacts(table, family, "prediction", context, output)
#' @family results-pilot
pilot_artifacts <- function(data, family, id, context, output) {
  if (!all(c("metric", "value") %in% names(data)) || any(!is.finite(data$value)))
    stop("Invalid pilot aggregate schema or values.")
  if (!family$figure_metric %in% data$metric) stop("Figure metric missing from aggregate table.")
  folder <- file.path(output, "runs", context$run_id)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  table_path <- file.path(folder, paste0(id, ".csv"))
  temp <- tempfile(tmpdir = folder)
  write.csv(data, temp, row.names = FALSE)
  pilot_finalize(temp, table_path)
  figure_path <- file.path(folder, paste0(id, ".pdf"))
  temp <- tempfile(fileext = ".pdf", tmpdir = folder)
  grDevices::pdf(temp, width = 9, height = 5, timestamp = FALSE,
                  producer = FALSE, author = "Andrew Walther", title = family$title)
  tryCatch({
    value <- data$value[match(family$figure_metric, data$metric)]
    graphics::barplot(value, ylim = c(0, 1), names.arg = family$figure_metric,
                       main = family$title, ylab = "Fraction / concordance")
    graphics::mtext("Synthetic software pilot only", side = 3, line = 0.3)
  }, finally = grDevices::dev.off())
  pilot_finalize(temp, figure_path)
  sidecar <- file.path(folder, paste0(id, ".provenance.yml"))
  temp <- tempfile(tmpdir = folder)
  yaml::write_yaml(list(run_id = context$run_id, family = id,
    manifest = "manifest.yml", artifacts_md5 = as.list(setNames(
      unname(tools::md5sum(c(table_path, figure_path))), basename(c(table_path, figure_path))))), temp)
  pilot_finalize(temp, sidecar)
  c(table_path, figure_path, sidecar)
}

#' Write a completed software-pilot manifest and browsable aggregate catalog
#' @param artifacts Named list of family artifact path vectors.
#' @param context Complete provenance list.
#' @param registry Registry families.
#' @param output Pilot output root.
#' @return Paths to immutable manifest and catalog, then local latest-run pointer.
#' @examples
#' # pilot_release(artifacts, context, registry, output)
#' @family results-pilot
pilot_release <- function(artifacts, context, registry, output) {
  paths <- unlist(artifacts, use.names = FALSE)
  if (any(!file.exists(paths))) stop("Missing pilot artifact; release is incomplete.")
  folder <- file.path(output, "runs", context$run_id)
  catalog <- do.call(rbind, lapply(names(artifacts), function(id) {
    data.frame(family = id, title = registry[[id]]$title, run_id = context$run_id,
      table = basename(artifacts[[id]][1]), figure = basename(artifacts[[id]][2]),
      scope = "synthetic software pilot only", stringsAsFactors = FALSE)
  }))
  catalog_path <- file.path(folder, "catalog.csv")
  temp <- tempfile(tmpdir = folder); write.csv(catalog, temp, row.names = FALSE)
  pilot_finalize(temp, catalog_path)
  manifest <- c(context, list(status = "complete software pilot", families = names(artifacts),
    artifact_md5 = as.list(setNames(unname(tools::md5sum(c(paths, catalog_path))),
                                   basename(c(paths, catalog_path))))))
  manifest_path <- file.path(folder, "manifest.yml")
  temp <- tempfile(tmpdir = folder); yaml::write_yaml(manifest, temp)
  pilot_finalize(temp, manifest_path)
  # This is a pilot-only pointer, never results/current.yml or a scientific release.
  latest <- file.path(output, "latest.yml")
  temp <- tempfile(tmpdir = output)
  yaml::write_yaml(list(run_id = context$run_id, manifest = manifest_path,
                        scientific_release = FALSE), temp)
  if (!file.rename(temp, latest)) stop("Could not update the local pilot pointer.")
  c(manifest_path, catalog_path, latest)
}

#' Read only a fresh and checksum-verified pilot catalog
#' @param output Pilot output directory.
#' @param outdated Names of stale targets obtained without running fits.
#' @return Catalog data frame, or an explicit error if unavailable/stale/corrupt.
#' @examples
#' # pilot_catalog("results/pilot/local", targets::tar_outdated(...))
#' @family results-pilot
pilot_catalog <- function(output, outdated) {
  if (length(outdated)) stop("Pilot is stale or incomplete. Run the explicit build command.")
  latest <- file.path(output, "latest.yml")
  if (!file.exists(latest)) stop("No completed pilot. Run the explicit build command.")
  pointer <- yaml::read_yaml(latest)
  manifest <- yaml::read_yaml(pointer$manifest)
  paths <- file.path(dirname(pointer$manifest), names(manifest$artifact_md5))
  if (any(!file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), unname(unlist(manifest$artifact_md5))))
    stop("Pilot artifacts are missing or do not match the manifest.")
  read.csv(file.path(dirname(pointer$manifest), "catalog.csv"), stringsAsFactors = FALSE)
}
