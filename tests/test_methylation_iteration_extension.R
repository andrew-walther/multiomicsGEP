# ============================================================
# Script: test_methylation_iteration_extension.R
# Purpose: Verify isolated iteration-extension checks with synthetic fits.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R, yaml
# ============================================================

local({
  runner <- new.env(parent = globalenv())
  sys.source("results/multimodal_real/run_methylation_single_modality.R", runner)
  script <- "results/multimodal_real/run_methylation_iteration_extension.R"

  #' Create an isolated baseline and inject a synthetic fit boundary
  #' @param constant Whether risk and one program are constant.
  #' @param divergent Whether the synthetic extension changes the ELBO prefix.
  #' @return List containing temporary paths, evaluation environment, and call count.
  #' @examples
  #' # fixture <- make_fixture()
  make_fixture <- function(constant = FALSE, divergent = FALSE) {
    root <- tempfile("extension-test-")
    base <- file.path(root, "base")
    out <- file.path(root, "output")
    dir.create(file.path(base, "source", "code"), recursive = TRUE)
    dir.create(file.path(base, "source", "config"), recursive = TRUE)
    dir.create(out)
    writeLines("original source", file.path(base, "source", "code", "dummy.R"))
    yaml::write_yaml(list(cavi = list(max_iter = 300L, tol = 1e-4),
                          benchmark = list(k_merged_yfb_desurv = 2L, alpha = .5)),
                     file.path(base, "source", "config", "globals.yml"))
    saveRDS(list(toy = 1:6), file.path(base, "design.rds"))
    scores <- cbind(a = 1:6, b = c(1, -1, 0, 0, -1, 1))
    if (constant) scores[, 2] <- 1
    risk <- if (constant) rep(0, 6) else 1:6
    baseline <- list(row = data.frame(status = "capped", seconds = 3),
      fit = list(EF = diag(2), history = list(n_iter = 300L, elbo_full = -1000:-701)),
      scores = scores, risk = risk, activity = list(active = c(TRUE, FALSE)))
    saveRDS(baseline, file.path(base, "methylation_50000_PL_f0.rds"))
    state <- new.env(parent = emptyenv())
    state$calls <- 0L
    state$cap <- NA_integer_
    env <- new.env(parent = globalenv())
    env$source <- function(...) {
      env$C_LOCAL <- base
      env$C_OUT <- out
      env$c_match <- runner$c_match
      NULL
    }
    #' Save an extended synthetic checkpoint or reuse an existing checkpoint
    #' @param modality,count,arm,fold Registered fit identifiers.
    #' @return NULL invisibly after recording the effective cap.
    #' @examples
    #' # c_one("methylation", 50000, "PL", 0)
    env$c_one <- function(modality, count, arm, fold) {
      state$calls <- state$calls + 1L
      path <- file.path(env$C_LOCAL, "methylation_50000_PL_f0.rds")
      if (file.exists(path)) return(invisible(NULL))
      state$cap <- yaml::read_yaml(file.path(env$C_LOCAL, "source/config/globals.yml"))$cavi$max_iter
      extended <- baseline
      extended$row$status <- "converged"
      extended$fit$history$n_iter <- 350L
      extended$fit$history$elbo_full <- c(baseline$fit$history$elbo_full, -700:-651)
      if (divergent) extended$fit$history$elbo_full[1] <- 0
      extended$fit$history$delta_elbo_rel <- c(.01, .001)
      extended$scores <- cbind(-scores[, 2], scores[, 1])
      extended$risk <- 1.5 * risk
      extended$activity$active <- c(FALSE, TRUE)
      extended$warnings <- character()
      saveRDS(extended, path)
      invisible(NULL)
    }
    list(root = root, base = base, out = out, env = env, state = state)
  }

  #' Evaluate the extension script against an isolated synthetic fixture
  #' @param fixture Result of make_fixture.
  #' @return NULL invisibly, or an error from the extension script.
  #' @examples
  #' # run_fixture(fixture)
  run_fixture <- function(fixture) sys.source(script, envir = fixture$env)

  f <- make_fixture()
  paths <- file.path(f$base, c("design.rds", "source/config/globals.yml",
                              "methylation_50000_PL_f0.rds"))
  before <- tools::md5sum(paths)
  run_fixture(f)
  stopifnot(f$state$cap == 600L, identical(unname(before), unname(tools::md5sum(paths))))
  ext <- file.path(f$base, "iteration600")
  stopifnot(yaml::read_yaml(file.path(ext, "source/config/globals.yml"))$cavi$max_iter == 600L,
            identical(readRDS(paths[1]), readRDS(file.path(ext, "design.rds"))))
  checkpoint <- file.path(ext, "methylation_50000_PL_f0.rds")
  checkpoint_hash <- tools::md5sum(checkpoint)
  writeLines(c("  12.50 real 1.00 user 0.20 sys", "  10485760 maximum resident set size"),
             file.path(ext, "process.log"))
  run_fixture(f)
  metrics <- read.csv(file.path(f$out, "iteration_extension.csv"))
  programs <- read.csv(file.path(f$out, "iteration_extension_programs.csv"))
  expected_rms <- .5 * sqrt(mean((1:6)^2)) / sd(1:6)
  stopifnot(identical(unname(checkpoint_hash), unname(tools::md5sum(checkpoint))),
            identical(programs$extended_program, c(2L, 1L)),
            isTRUE(all.equal(programs$score_correlation, c(1, -1))),
            identical(programs$baseline_active, c(TRUE, FALSE)),
            identical(programs$extended_active, c(TRUE, FALSE)),
            abs(metrics$risk_RMS_change_over_baseline_SD - expected_rms) < 1e-12,
            f$state$calls == 2L, metrics$process_seconds == 12.5,
            metrics$peak_rss_MiB == 10, metrics$memory_measured,
            metrics$prefix_relative_elbo_error == 0,
            metrics$risk_correlation_defined, abs(metrics$risk_correlation_with_300 - 1) < 1e-12)
  unlink(f$root, recursive = TRUE)

  f <- make_fixture(constant = TRUE)
  run_fixture(f)
  metrics <- read.csv(file.path(f$out, "iteration_extension.csv"))
  programs <- read.csv(file.path(f$out, "iteration_extension_programs.csv"))
  stopifnot(!metrics$risk_correlation_defined, !metrics$risk_RMS_change_defined,
            is.na(metrics$risk_correlation_with_300), is.na(metrics$risk_RMS_change_over_baseline_SD),
            !programs$comparison_defined[2], is.na(programs$score_correlation[2]))
  unlink(f$root, recursive = TRUE)

  f <- make_fixture(divergent = TRUE)
  err <- tryCatch({ run_fixture(f); "" }, error = conditionMessage)
  stopifnot(grepl("does not reproduce the 300-iteration baseline", err, fixed = TRUE))
  unlink(f$root, recursive = TRUE)

  for (kind in c("config", "code")) {
    f <- make_fixture()
    run_fixture(f)
    ext <- file.path(f$base, "iteration600")
    if (kind == "config") {
      config_path <- file.path(ext, "source/config/globals.yml")
      cfg <- yaml::read_yaml(config_path)
      cfg$cavi$tol <- .1
      yaml::write_yaml(cfg, config_path)
      expected <- "Extension configuration differs from baseline."
    } else {
      writeLines("changed source", file.path(ext, "source/code/dummy.R"))
      expected <- "Extension model sources differ from baseline."
    }
    before_calls <- f$state$calls
    err <- tryCatch({ run_fixture(f); "" }, error = conditionMessage)
    stopifnot(grepl(expected, err, fixed = TRUE), f$state$calls == before_calls)
    unlink(f$root, recursive = TRUE)
  }
  cat("Methylation iteration extension tests passed.\n")
})
