# ============================================================
# Script: run_expression_iteration_extension_k5.R
# Purpose: Restart capped K = 5 expression fits with a verified 1,000-sweep cap.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, parallel; frozen K = 5 study dependencies
# ============================================================
# Run from the repository root: Rscript <this file> run
# One numerical thread per child, at most three children. For Longleaf request
# one CPU and several GB per child. Private matrices and fits stay under data/.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
E5_BASE <- "data/methylation_single_modality_2026-11-12/K05"
E5_OUT <- "results/multimodal_real/outputs/methylation_single_modality_2026-11-12/K05"
E5_SCRIPT <- "results/multimodal_real/run_expression_iteration_extension_k5.R"

#' Verify that an extended ELBO trajectory reproduces the saved prefix
#' @param old,new Numeric ELBO sequences, baseline and extended respectively.
#' @return Maximum relative prefix error; stops above 10^-10 or for invalid input.
#' @examples
#' e5_prefix(c(-10, -9), c(-10, -9, -8))
e5_prefix <- function(old, new) {
  if (!length(old) || length(new) < length(old) || any(!is.finite(c(old, new))))
    stop("K = 5 ELBO prefix is missing, short or nonfinite.")
  error <- max(abs(new[seq_along(old)] - old) / pmax(1, abs(old)))
  if (error > 1e-10) stop("K = 5 ELBO prefix failed to reproduce.")
  error
}

#' Change only the iteration cap in a copied configuration
#' @param cfg Frozen configuration list.
#' @return Configuration list with cavi$max_iter equal to 1000.
#' @examples
#' e5_config(list(cavi = list(max_iter = 300L)))
e5_config <- function(cfg) {
  cfg$cavi$max_iter <- 1000L
  cfg
}

#' Compute the mean and SE of five paired concordance differences
#' @param meth,expression Numeric fold-ordered vectors of length five.
#' @return Mean difference and sd(meth-expression)/sqrt(5), or NA if undefined.
#' @examples
#' e5_paired(rep(.6, 5), c(.5, .6, .7, .5, .6))
e5_paired <- function(meth, expression) {
  stopifnot(length(meth) == 5L, length(expression) == 5L)
  d <- meth - expression
  if (any(!is.finite(d))) return(c(mean_difference = NA_real_, paired_SE = NA_real_))
  c(mean_difference = mean(d), paired_SE = stats::sd(d) / sqrt(5))
}

#' Copy and verify the K = 5 source and design, allowing only the cap override
#' @return NULL invisibly; creates an isolated private iteration1000 directory.
#' @examples
#' # e5_prepare()
e5_prepare <- function() {
  ext <- file.path(E5_BASE, "iteration1000")
  src <- file.path(E5_BASE, "source")
  manifest <- read.csv(file.path(E5_OUT, "source_hashes.csv"))
  actual <- unname(tools::md5sum(file.path(src, manifest$path)))
  if (!identical(actual, manifest$md5)) stop("K = 5 baseline source hash mismatch.")
  dir.create(ext, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(file.path(ext, "source"))) {
    if (!file.copy(src, ext, recursive = TRUE) ||
        !file.copy(file.path(E5_BASE, "design.rds"), ext)) stop("K = 5 snapshot copy failed.")
    cfg <- yaml::read_yaml(file.path(src, "config/globals.yml"))
    yaml::write_yaml(e5_config(cfg), file.path(ext, "source/config/globals.yml"))
  }
  design <- c(file.path(E5_BASE, "design.rds"), file.path(ext, "design.rds"))
  if (anyNA(tools::md5sum(design)) || length(unique(unname(tools::md5sum(design)))) != 1L)
    stop("K = 5 copied design hash mismatch.")
  keep <- manifest$path != "config/globals.yml"
  if (!identical(unname(tools::md5sum(file.path(ext, "source", manifest$path[keep]))),
                 manifest$md5[keep])) stop("K = 5 copied source hash mismatch.")
  cfg <- yaml::read_yaml(file.path(src, "config/globals.yml"))
  changed <- yaml::read_yaml(file.path(ext, "source/config/globals.yml"))
  if (!identical(changed, e5_config(cfg))) stop("K = 5 cap is not the only config change.")
  paths <- c("design.rds", paste0("source/", manifest$path))
  write.csv(data.frame(K = 5L, path = paths,
    md5 = unname(tools::md5sum(file.path(ext, paths)))),
    file.path(E5_OUT, "expression_iteration1000_hashes.csv"), row.names = FALSE)
  invisible(NULL)
}

#' Load the frozen K = 5 runner and select the separate extension path
#' @return Environment containing unchanged fitting and matching functions.
#' @examples
#' # e5_runner()
e5_runner <- function() {
  env <- new.env(parent = globalenv())
  source(file.path(E5_BASE, "source/runner.R"), local = env)
  env$c_set_rank(5L)
  env$C_LOCAL <- file.path(E5_BASE, "iteration1000")
  env
}

#' Restart one capped K = 5 expression fit and enforce trajectory reproduction
#' @param arm PL or PE.
#' @param fold Zero for full fit, otherwise held-out fold 1 through 5.
#' @return NULL invisibly; saves fit privately and fails loudly on errors/prefix mismatch.
#' @examples
#' # e5_worker("PE", 0)
e5_worker <- function(arm, fold) {
  env <- e5_runner()
  name <- sprintf("expression_2218_%s_f%d.rds", arm, fold)
  old <- readRDS(file.path(E5_BASE, name))
  stopifnot(old$row$status == "capped", old$fit$history$n_iter == 300L,
            ncol(old$fit$EF) == 5L)
  env$c_one("expression", 2218L, arm, fold)
  new <- readRDS(file.path(env$C_LOCAL, name))
  if (is.null(new$fit) || new$row$status == "failed") stop("K = 5 expression refit failed: ", new$row$error)
  stopifnot(ncol(new$fit$EF) == 5L)
  e5_prefix(old$fit$history$elbo_full, new$fit$history$elbo_full)
  invisible(NULL)
}

#' Run the eight approved K = 5 restarts, then summarize without changing baseline tables
#' @return NULL invisibly; new aggregate iteration1000 tables under K05.
#' @examples
#' # e5_run()
e5_run <- function() {
  e5_prepare()
  jobs <- data.frame(arm = c(rep("PL", 3), rep("PE", 5)), fold = c(1L, 3L, 4L, 0L, 1L, 3L, 4L, 5L))
  # Batches stop before launching any later batch if a child fails.
  for (batch in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs))/3))) {
    result <- parallel::mclapply(batch, function(i) {
      log <- file.path(E5_BASE, "iteration1000", sprintf("expression_2218_%s_f%d.log", jobs$arm[i], jobs$fold[i]))
      system2(file.path(R.home("bin"), "Rscript"),
        c(shQuote(E5_SCRIPT), "one", jobs$arm[i], jobs$fold[i]), stdout = log, stderr = log)
    }, mc.cores = 3L)
    if (any(vapply(result, function(x) !is.numeric(x) || x != 0L, logical(1))))
      stop("K = 5 expression child failed. Inspect private iteration1000 logs.")
  }
  e5_summary()
  invisible(NULL)
}

#' Recompute K = 5 reference and paired differences using verified extended fits
#' @return NULL invisibly; writes fit diagnostics, matching, paired and program tables.
#' @examples
#' # e5_summary()
e5_summary <- function() {
  e5_prepare()
  env <- e5_runner()
  z <- readRDS(file.path(E5_BASE, "design.rds"))
  env$c_load()
  cfg <- yaml::read_yaml(file.path(env$C_LOCAL, "source/config/globals.yml"))
  diagnostics <- matching <- list()
  refs <- list()
  for (arm in c("PL", "PE")) {
    refs[[arm]] <- vector("list", 6)
    for (f in 0:5) {
      name <- sprintf("expression_2218_%s_f%d.rds", arm, f)
      old <- readRDS(file.path(E5_BASE, name))
      rerun <- old$row$status == "capped"
      new <- if (rerun) readRDS(file.path(env$C_LOCAL, name)) else old
      if (is.null(new$fit)) stop("K = 5 refit is missing its posterior.")
      prefix <- e5_prefix(old$fit$history$elbo_full, new$fit$history$elbo_full)
      map <- env$c_match(old$scores, new$scores)
      r <- stats::cor(old$scores, new$scores)
      refs[[arm]][[f+1L]] <- new
      C <- if (f > 0L) frozen_reverse_cindex(new$risk[new$te], z$time[new$te], z$event[new$te]) else NA_real_
      diagnostics[[length(diagnostics)+1L]] <- data.frame(K = 5L, arm, fold = f,
        rerun, baseline_status = old$row$status, baseline_iterations = old$row$iterations,
        status = new$row$status, iterations = new$row$iterations,
        iteration_cap = if (rerun) 1000L else 300L,
        final_relative_elbo_change = tail(new$fit$history$delta_elbo_rel, 1), tolerance = cfg$cavi$tol,
        prefix_relative_elbo_error = prefix, risk_correlation_with_300 = stats::cor(old$risk, new$risk),
        baseline_C = old$row$C, C, warning_count = length(new$warnings),
        seconds = new$row$seconds)
      matching[[length(matching)+1L]] <- data.frame(K = 5L, arm, fold = f,
        baseline_program = 1:5, extended_program = map,
        score_correlation = r[cbind(1:5, map)])
    }
  }
  diag <- do.call(rbind, diagnostics)
  write.csv(diag, file.path(E5_OUT, "expression_iteration1000.csv"), row.names = FALSE)
  write.csv(do.call(rbind, matching), file.path(E5_OUT, "expression_iteration1000_programs.csv"), row.names = FALSE)
  old_summary <- read.csv(file.path(E5_OUT, "cv_summary.csv"))
  cv <- read.csv(file.path(E5_OUT, "cv_fits.csv"))
  out <- lapply(seq_len(nrow(old_summary)), function(i) {
    s <- old_summary[i, ]; ex <- diag[diag$arm == s$arm & diag$fold > 0L, ]
    ex <- ex[match(1:5, ex$fold), ]
    m <- cv[cv$modality == "methylation" & cv$arm == s$arm & cv$count == s$count, ]
    paired <- e5_paired(m$C[match(1:5, m$fold)], ex$C)
    data.frame(K = 5L, count = s$count, arm = s$arm, cv_status = s$cv_status,
      methylation_mean_C = s$mean_C, expression_mean_C_300 = s$expression_mean_C,
      expression_mean_C_1000 = mean(ex$C), expression_converged_folds_1000 = sum(ex$status == "converged"),
      mean_difference_300 = s$mean_difference, paired_SE_300 = s$paired_SE,
      mean_difference_1000 = paired[1], paired_SE_1000 = paired[2])
  })
  write.csv(do.call(rbind, out), file.path(E5_OUT, "cv_summary_iteration1000.csv"), row.names = FALSE)
  programs <- read.csv(file.path(E5_OUT, "programs.csv"))
  for (i in seq_len(nrow(programs))) {
    p <- programs[i, ]; ex <- refs[[p$arm]][[1]]
    if (!ex$row$converged) next
    meth <- readRDS(file.path(E5_BASE, sprintf("methylation_%d_%s_f0.rds", p$count, p$arm)))
    programs$max_abs_expression_r[i] <- max(abs(stats::cor(meth$scores[, p$program], ex$scores)))
    folds <- cv[cv$modality == "methylation" & cv$count == p$count & cv$arm == p$arm, ]
    assess <- nrow(folds) == 5L && sum(folds$converged) >= 4L
    programs$methylation_only_assessable[i] <- assess
    programs$methylation_only[i] <- if (assess)
      p$full_survival_active && p$active_converged_folds >= 3 && programs$max_abs_expression_r[i] < .5 else NA
  }
  write.csv(programs, file.path(E5_OUT, "programs_iteration1000.csv"), row.names = FALSE)
  print(diag[, c("K", "arm", "fold", "rerun", "status", "iterations", "C")])
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[1] == "run") e5_run()
  else if (args[1] == "one") e5_worker(args[2], as.integer(args[3]))
  else if (args[1] == "summary") e5_summary()
  else stop("Use run, summary or one ARM FOLD. All fits are K = 5.")
}
