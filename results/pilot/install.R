# ============================================================
# Script: install.R
# Purpose: Install the pilot's locked R packages in a local library.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: Base R
# ============================================================

# Run from the repository root: Rscript results/pilot/install.R

#' Install locked packages in dependency order
#'
#' @param lock Data frame with Package, Version, and Dependencies columns.
#' @param library_path Directory for locally installed packages.
#' @return Invisibly returns the names of packages processed.
#' @family pilot dependency setup
#' @examples
#' \dontrun{install_locked(lock, "results/pilot/library")}
install_locked <- function(lock, library_path) {
  dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
  library_path <- normalizePath(library_path, mustWork = TRUE)
  .libPaths(c(library_path, .libPaths()))
  if (anyDuplicated(lock$Package) || anyNA(lock$Package) ||
      anyNA(lock$Version) || any(!nzchar(lock$Package)) ||
      any(!nzchar(lock$Version))) stop("Invalid package lockfile")
  dependencies <- strsplit(lock$Dependencies, "|", fixed = TRUE)
  names(dependencies) <- lock$Package
  if (length(setdiff(unlist(dependencies), lock$Package))) {
    stop("Lockfile references an unlocked dependency")
  }
  remaining <- lock$Package
  done <- character()
  while (length(remaining)) {
    ready <- remaining[vapply(dependencies[remaining], function(x) {
      all(x %in% done)
    }, logical(1))]
    if (!length(ready)) stop("Dependency cycle in lockfile")
    for (pkg in ready) {
      version <- lock$Version[match(pkg, lock$Package)]
      visible <- packageDescription(pkg, fields = "Version")
      exact <- !is.na(visible) && identical(visible, version)
      if (!exact) {
        archive <- tempfile(fileext = ".tar.gz")
        on.exit(unlink(archive), add = TRUE)
        filename <- paste0(pkg, "_", version, ".tar.gz")
        urls <- c(paste0("https://cloud.r-project.org/src/contrib/Archive/", pkg, "/", filename),
                  paste0("https://cloud.r-project.org/src/contrib/", filename))
        downloaded <- FALSE
        for (url in urls) {
          status <- tryCatch(suppressWarnings(download.file(url, archive, mode = "wb", quiet = TRUE)),
                             error = function(e) 1L)
          if (identical(status, 0L)) { downloaded <- TRUE; break }
        }
        if (!downloaded) stop("Cannot download locked package: ", pkg, " ", version)
        install.packages(archive, lib = library_path, repos = NULL,
                         type = "source", dependencies = FALSE, quiet = TRUE)
        local <- installed.packages(lib.loc = library_path)
        if (!(pkg %in% rownames(local)) || local[pkg, "Version"] != version) {
          stop("Installation failed or version mismatch: ", pkg, " ", version)
        }
      }
      done <- c(done, pkg)
    }
    remaining <- setdiff(remaining, ready)
  }
  invisible(done)
}

lock <- read.csv("results/pilot/package-versions.csv", colClasses = "character",
                 na.strings = NULL, check.names = FALSE)
install_locked(lock, "results/pilot/library")
