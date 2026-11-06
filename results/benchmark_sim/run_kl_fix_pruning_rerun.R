# ============================================================
# Script: run_kl_fix_pruning_rerun.R
# Purpose: Rerun the October 15 single-modality K_init = 7 pruning cells after
#          the compute_ebnm_kl() normalizer fix (issue #6, item 1)
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: survival, ebnm, yaml, parallel; run_k_selection_panel.R and
#               run_supervision_diagnostics.R helpers (make_testbed, score_fit)
# Run: Rscript results/benchmark_sim/run_kl_fix_pruning_rerun.R [cores] [output_dir] [n_seeds]
# ============================================================
# Registration: docs/plans/2026-11-06-kl-fix-pruning-rerun-plan.md and attempt
# KL-1 in docs/plans/2026-10-16-round2-attempts-log.md, committed before running.
#
# Test bed: unchanged October 15 calibrated bed (make_testbed(): n_train 400,
# n_validation 120, p 50 expression features, planted log HR beta_3 = 3.3,
# z-standardized cohorts), seeds 20261016:20261035, noise {0.25, 1}, K_init 7.
# Arms (fit settings copied from run_k_selection_panel.R, Task1-3):
#   reference  alpha_F = 0 (false-positive reference, as in Task1-3)
#   unpruned   selected canonical fit: canonicalize, w_S = 1, supervised start,
#              shared point-Laplace beta prior, warm-up 10; prune = FALSE
#   ordinary   unpruned + prune = TRUE (ELBO nullcheck pruning)
#   protected  ordinary + prune_protect_survival = TRUE
# All fits: point-Laplace L/F, max_iter and tol from config/globals.yml
# (300, 1e-5), converge_on = "elbo". No fitter default is changed.
# Criteria (Task1-1, unchanged): recovery >= 16/20; mean false-active <=
# reference + 0.5; 20/20 converged (ELBO flag of the returned fit); median
# retained K = 3 (pruning arms only; the unpruned arm keeps K = 7 by design).
# Parallelization: PSOCK workers, one serial fit per worker. 160 cells; the
# brief caps local workers at 4 while another job holds 9 cores. Each cell is
# checkpointed before aggregation; an existing output folder is never overwritten.

source("results/benchmark_sim/run_k_selection_panel.R")  # load_panel_code(), k_panel_summary()

#' Fit settings for each registered arm
#' @return Named list of argument lists passed to fit_ssbmf_yfb() on top of the shared ones.
#' @examples
#' names(kl_rerun_arms())
kl_rerun_arms <- function() {
  sel <- list(canonicalize = TRUE, w_S = 1, init_supervised = TRUE,
              beta_prior_mode = "shared", prior_beta = "point_laplace", beta_warmup = 10L)
  list(reference = list(alpha_F = 0),
       unpruned  = sel,
       ordinary  = c(sel, list(prune = TRUE)),
       protected = c(sel, list(prune = TRUE, prune_protect_survival = TRUE)))
}

#' Construct the 160-cell design
#' @param n_seeds Number of seeds from 20261016 (20 registered; fewer only for the smoke test).
#' @return Data frame with seed, noise, K, model, arm.
#' @examples
#' nrow(kl_rerun_grid())
kl_rerun_grid <- function(n_seeds = 20L) {
  expand.grid(seed = 20261015L + seq_len(n_seeds), noise = c(.25, 1), K = 7L,
              model = "single", arm = names(kl_rerun_arms()), stringsAsFactors = FALSE)
}

#' Fit and score one cell; errors and warnings are kept in the row
#' @param i Row of the design `grid` (exported to the worker).
#' @return One data-frame row; also written to <out_dir>/cells/cell_XXXX.csv.
#' @examples
#' # kl_rerun_cell(1L)
kl_rerun_cell <- function(i) {
  x <- grid[i, ]; t0 <- proc.time()[["elapsed"]]; fit_warnings <- character()
  answer <- withCallingHandlers(tryCatch({
    tb <- make_testbed(x$seed, x$noise)
    a <- c(list(K = x$K, max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = "point_laplace", verbose = FALSE, trace_supervision = TRUE,
                converge_on = "elbo"), kl_rerun_arms()[[x$arm]])
    fit <- do.call(fit_ssbmf_yfb, c(list(tb$Ytr, tb$ttr, tb$str), a))
    sc <- score_fit(fit, tb)
    # Diagnostic only: largest ELBO gain from deleting one program of the returned
    # fit (>= 0 means the nullcheck would prune again). Canonical arms only.
    gains <- if (isTRUE(a$canonicalize)) prune_yf_nullcheck(fit, tb$Ytr, tb$ttr, tb$str,
      protect_survival = isTRUE(a$prune_protect_survival)) else NA_real_
    gain <- max(gains)
    fmt <- function(v) paste(signif(v, 4), collapse = ";")
    cbind(x, error = NA_character_, runtime = proc.time()[["elapsed"]] - t0,
          data.frame(recovered = sc$recovered, false_active = sc$false_z,
            K_final = ncol(fit$EF), converged = fit$history$converged,
            elbo_converged = fit$history$elbo_converged,
            outputs_converged = fit$history$outputs_converged,
            c_validation = sc$c_validation, cor_planted = sc$cor_planted,
            n_active = sc$n_active_z, iterations = sc$n_iter,
            # Caps in intermediate pruning refits, reported apart from the returned flag
            pruning_caps = sum(!as.logical(fit$pruning$converged_before)),
            n_pruned = NROW(fit$pruning), n_programs = sc$n_programs,
            max_prune_gain = gain, prune_gains = fmt(gains),
            F_sq_norms = fmt(colSums(fit$EF^2)),
            max_abs_beta = sc$max_abs_beta))
  }, error = function(e) cbind(x, error = conditionMessage(e),
                               runtime = proc.time()[["elapsed"]] - t0)),
  warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w))
                          invokeRestart("muffleWarning") })
  answer$warnings <- paste(unique(fit_warnings), collapse = " | ")
  write.csv(answer, file.path(out_dir, "cells", sprintf("cell_%04d.csv", i)), row.names = FALSE)
  answer
}

#' Run the registered rerun
#' @param cores PSOCK worker count (at most 4 while S-15 runs).
#' @param output_dir Fresh output folder; refuses to overwrite.
#' @param n_seeds 20 for the registered run.
#' @return Invisibly the summary table; writes design, results, summary CSVs.
#' @examples
#' # run_kl_fix_pruning_rerun(4L)
run_kl_fix_pruning_rerun <- function(cores = 4L,
    output_dir = "results/benchmark_sim/outputs/kl_fix_pruning_rerun_2026-11-06",
    n_seeds = 20L) {
  out_dir <- output_dir
  if (dir.exists(file.path(out_dir, "cells"))) stop("Output already exists; do not overwrite")
  dir.create(file.path(out_dir, "cells"), recursive = TRUE)
  cfg <- yaml::read_yaml("config/globals.yml"); grid <- kl_rerun_grid(n_seeds)
  write.csv(grid, file.path(out_dir, "design.csv"), row.names = FALSE)
  # RUN_COMMIT is set when running from a `git archive` snapshot of the registered commit
  commit <- Sys.getenv("RUN_COMMIT", "")
  if (!nzchar(commit)) commit <- system("git rev-parse HEAD", intern = TRUE)
  writeLines(c(paste("commit:", commit),
               paste("started:", format(Sys.time())), paste("cores:", cores),
               paste("R:", R.version.string)), file.path(out_dir, "provenance.txt"))
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, c("load_panel_code", "kl_rerun_cell", "kl_rerun_arms",
                                "grid", "cfg", "out_dir"), envir = environment())
  invisible(parallel::clusterEvalQ(cl, load_panel_code()))
  rows <- tryCatch(parallel::parLapplyLB(cl, seq_len(nrow(grid)), kl_rerun_cell),
                   finally = parallel::stopCluster(cl))
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
  write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)
  summ <- k_panel_summary(res)
  summ$rank_pass[summ$arm %in% c("reference", "unpruned")] <- NA  # no rank selection without pruning
  write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
  print(summ, row.names = FALSE)
  if (any(!is.na(res$error))) warning("Cell errors retained in results.csv; those groups cannot pass")
  invisible(summ)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  run_kl_fix_pruning_rerun(if (length(a) >= 1) as.integer(a[1]) else 4L,
                           if (length(a) >= 2) a[2] else
                             "results/benchmark_sim/outputs/kl_fix_pruning_rerun_2026-11-06",
                           if (length(a) >= 3) as.integer(a[3]) else 20L)
}
