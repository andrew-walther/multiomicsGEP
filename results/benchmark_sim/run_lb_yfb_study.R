# ============================================================
# Script: run_lb_yfb_study.R
# Purpose: Run paired LB and YFB training and simulation comparisons.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm, yaml and the study fitting scripts
# Usage: Rscript results/benchmark_sim/run_lb_yfb_study.R --mode=training
#   --scaling=zscore --prior_LF=point_laplace --prior_beta=normal
#   --out=results/benchmark_sim/outputs/lb_yfb_run
#   --input=data/lf_prior_real_2026-11-12.rds
# ============================================================
# Serial execution, one BLAS thread. Patient-level objects stay under data/.

#' Parse and validate the explicit study settings
#' @param args Command-line arguments in --name=value form.
#' @return Named list of validated options.
#' @examples
#' # lb_yfb_runner_options(c('--mode=smoke', '--scaling=zscore',
#' #   '--prior_LF=point_laplace', '--prior_beta=normal', '--out=/tmp/run'))
lb_yfb_runner_options <- function(args) {
  if (any(!grepl('^--[^=]+=.+$', args))) stop('Use --name=value for every option')
  keys <- sub('^--([^=]+)=.*$', '\\1', args)
  vals <- sub('^--[^=]+=', '', args)
  allowed <- c('mode', 'scaling', 'prior_LF', 'prior_beta', 'out', 'input', 'K', 'w_S')
  if (anyDuplicated(keys) || any(!keys %in% allowed)) stop('Duplicate or unknown option')
  x <- as.list(stats::setNames(vals, keys))
  required <- c('mode', 'scaling', 'prior_LF', 'prior_beta', 'out')
  if (!all(required %in% keys)) stop('mode, scaling, prior_LF, prior_beta and out are required')
  if (!x$mode %in% c('smoke', 'simulation', 'training')) stop('Invalid mode')
  if (!x$scaling %in% c('zscore', 'sd')) stop('Invalid scaling')
  if (!x$prior_LF %in% c('point_laplace', 'point_exponential')) stop('Invalid prior_LF')
  if (!x$prior_beta %in% c('normal', 'point_normal', 'point_laplace')) stop('Invalid prior_beta')
  if (!is.null(x$K)) {
    k <- suppressWarnings(as.numeric(x$K))
    if (!is.finite(k) || k <= 0 || k != floor(k) || k > .Machine$integer.max) stop('K must be a positive integer')
    x$K <- as.integer(k)
  }
  if (is.null(x$w_S)) x$w_S <- 1 else {
    x$w_S <- suppressWarnings(as.numeric(x$w_S))
    if (!is.finite(x$w_S) || x$w_S <= 0) stop('w_S must be positive and finite')
  }
  if (is.null(x$input)) x$input <- 'data/lf_prior_real_2026-11-12.rds'
  x
}

#' Enumerate the registered datasets and survival weights
#' @return Data frame with 104 dataset-weight combinations.
#' @examples
#' # lb_yfb_runner_design()
lb_yfb_runner_design <- function() {
  a <- expand.grid(seed = 1:10, share = c(.0025, .01, .05, .0005), w_S = c(1, 10))
  a$design <- ifelse(a$share == .0005, 'S14_strong', 'S9')
  a$amp <- NA_real_
  b <- expand.grid(seed = 1:6, amp = c(.5, 1), w_S = c(1, 10))
  b$design <- 'matched'; b$share <- NA_real_
  rbind(a[, c('design', 'seed', 'share', 'amp', 'w_S')],
        b[, c('design', 'seed', 'share', 'amp', 'w_S')])
}

#' Generate the November 6 planted-factor dataset
#' @param seed Random seed. @param amp Signal amplitude.
#' @return List with patient matrix, outcomes, true loadings and split indices.
#' @examples
#' # lb_yfb_runner_planted(1, .5)
lb_yfb_runner_planted <- function(seed, amp) {
  set.seed(seed); n <- 300L; p <- 400L; K <- 4L
  L <- matrix(rexp(n * K), n, K); F <- matrix(0, p, K)
  for (k in seq_len(K)) F[sample.int(p, round(.1 * p)), k] <- amp
  tau <- rgamma(p, 2, 2)
  Y <- L %*% t(F) + sweep(matrix(rnorm(n * p), n, p), 2, sqrt(tau), '/') * .5
  beta <- c(1, -.8, .5, 0); eta <- drop(L %*% beta)
  ev <- (-log(runif(n)) / (.01 * exp(eta)))^(1 / 1.5)
  cen <- rexp(n) * median(ev) * 1.5
  list(X = Y, time = pmin(ev, cen), status = as.integer(ev <= cen),
       L = L, beta = beta, train = 1:200, validation = 201:300)
}

#' Write aggregate tables without patient-level objects
#' @param tables Named list of data frames. @param out Output directory.
#' @return Invisibly, output paths.
#' @examples
#' # lb_yfb_runner_write(list(metrics = data.frame(x = 1)), '/tmp/run')
lb_yfb_runner_write <- function(tables, out) {
  paths <- file.path(out, paste0(names(tables), '.csv'))
  for (i in seq_along(tables)) write.csv(tables[[i]], paths[i], row.names = FALSE)
  invisible(paths)
}

#' Express planted molecular loadings in the study-scaled feature units
#' @param truth Molecular loading on the generator's returned training matrix.
#' @param X Training matrix before study scaling.
#' @return Feature loading after division by training feature standard deviations.
#' @examples
#' # lb_yfb_runner_loading_truth(c(1, 2), matrix(c(1, 3, 2, 6), 2))
lb_yfb_runner_loading_truth <- function(truth, X) {
  truth / apply(X, 2, sd)
}

#' Preserve a failed simulation cell as two explicit model rows
#' @param cell One registered design row. @param expr Cell computation.
#' @return The computation result or failure rows and the error message.
#' @examples
#' # lb_yfb_runner_capture(cell, stop("example failure"))
lb_yfb_runner_capture <- function(cell, expr) {
  tryCatch(expr, error = function(e) {
    msg <- conditionMessage(e)
    rows <- do.call(rbind, lapply(c('LB', 'YFB'), function(model)
      data.frame(cell, model = model, converged = FALSE, n_iter = NA_integer_,
        foldin_converged = FALSE, recovery_success = FALSE,
        molecular_loading_correlation = NA_real_, matched_score_recovery = NA_real_,
        false_active = NA_integer_, validation_C = NA_real_,
        C_within = NA_real_, C_pooled = NA_real_, pairs_within = NA_real_,
        pairs_pooled = NA_real_, C_simulation = NA_real_, pairs_simulation = NA_real_, error = msg)))
    list(metrics = rows, error = msg)
  })
}

#' Run the registered simulation cells serially
#' @param opts Study options. @param K Fitted rank. @param max_iter Maximum sweeps.
#' @param tol Convergence tolerance. @param center Whether projection is centered.
#' @param checkpoint Ignored directory for patient-level fits.
#' @return Aggregate simulation, comparison and program tables.
#' @examples
#' # lb_yfb_runner_simulation(opts, 7, 300, 1e-5, FALSE, 'data/run')
lb_yfb_runner_simulation <- function(opts, K, max_iter, tol, center, checkpoint) {
  source('results/benchmark_sim/run_yfb_supervision_decision.R', local = .GlobalEnv)
  design <- lb_yfb_runner_design(); rows <- list(); comparisons <- list(); programs <- list()
  for (i in seq_len(nrow(design))) {
    cell <- design[i, ]; set.seed(cell$seed)
    result <- lb_yfb_runner_capture(cell, {
    if (cell$design == 'matched') {
      d <- lb_yfb_runner_planted(cell$seed, cell$amp)
      X <- d$X; tr <- d$train; va <- d$validation
      time <- d$time; status <- d$status
    } else {
      d <- yfb_decision_data(cell$seed, 2000L, cell$share, truth_mode = 'shared')
      # Historical generator SD-scales cohorts separately before returning them.
      # Further preprocessing below uses training statistics for both cohorts.
      X <- rbind(d$Y, d$Yv); tr <- seq_len(nrow(d$Y)); va <- nrow(d$Y) + seq_len(nrow(d$Yv))
      time <- c(d$time, d$tv); status <- c(d$event, d$ev)
    }
    cohort <- rep('simulation', nrow(X))
    Y <- lb_yfb_study_scale(X, cohort, train = tr, scaling = opts$scaling)
    pair <- lb_yfb_study_pair(Y[tr, , drop = FALSE], time[tr], status[tr],
      cohort[tr], K, opts$prior_LF, opts$prior_beta, cell$w_S, max_iter, tol, center)
    predictions <- list(); cell_rows <- list()
    for (model in c('LB', 'YFB')) {
      fit <- pair$fits[[model]]
      pred <- if (model == 'LB') predict_ssbmf_lb_foldin(Y[va, , drop = FALSE], fit) else
        predict_ssbmf_yfb(Y[va, , drop = FALSE], fit$EF, fit$EBeta,
                          EF_norms = fit$EF_norms, projection_center = fit$projection_center)
      predictions[[model]] <- pred
      sd_beta <- sqrt(if (!is.null(fit$VBeta)) fit$VBeta else pmax(fit$EBeta2 - fit$EBeta^2, 0))
      active <- abs(fit$EBeta) >= 1.96 * sd_beta & sd_beta > 0
      if (cell$design != 'matched') {
        # A molecular loading on Y/sd is the original loading divided by sd.
        truth <- lb_yfb_runner_loading_truth(d$truth, X[tr, , drop = FALSE])
        rr <- suppressWarnings(as.vector(cor(truth, fit$EF)))
        if (all(!is.finite(rr))) { recovery <- FALSE; false_active <- sum(active); corr <- 0 } else {
          k <- which.max(replace(abs(rr), !is.finite(rr), -Inf)); corr <- abs(rr[k])
          matched <- corr >= .7
          recovery <- matched && active[k] && fit$EBeta[k] * rr[k] > 0
          false_active <- sum(active[setdiff(seq_along(active), if (matched) k else integer())])
        }
      } else {
        scores <- if (model == 'LB') pred$L_test else pred$ZF_test
        r <- suppressWarnings(abs(cor(d$L[va, , drop = FALSE], scores)))
        # Dead estimated programs score zero recovery, as in the November 6 comparator.
        r[!is.finite(r)] <- 0
        # Greedy one-to-one matching follows the November 6 comparator.
        best <- rep(NA_real_, ncol(d$L)); match_k <- rep(NA_integer_, ncol(d$L))
        for (j in seq_len(min(dim(r)))) {
          valid <- is.finite(r); if (!any(valid)) break
          temp <- replace(r, !valid, -Inf); ij <- which(temp == max(temp), arr.ind = TRUE)[1, ]
          best[ij[1]] <- r[ij[1], ij[2]]; match_k[ij[1]] <- ij[2]
          r[ij[1], ] <- NA_real_; r[, ij[2]] <- NA_real_
        }
        corr <- if (all(is.finite(best[1:3]))) mean(best[1:3]) else NA_real_
        recovery <- NA
        # Only recovered true survival programs can exempt an estimated
        # program from false activity. Weak arbitrary matches do not qualify.
        recovered_k <- match_k[1:3][is.finite(best[1:3]) & best[1:3] >= .7 & !is.na(match_k[1:3])]
        false_active <- sum(active[setdiff(seq_along(active), recovered_k)])
      }
      concordance <- lb_yfb_study_concordance(pred$risk_scores, time[va], status[va], cohort[va])
      cell_rows[[model]] <- cbind(data.frame(cell, model = model,
        converged = isTRUE(fit$history$converged), n_iter = fit$history$n_iter,
        foldin_converged = if (model == 'LB') isTRUE(pred$converged) else TRUE,
        recovery_success = recovery,
        molecular_loading_correlation = if (cell$design == 'matched') NA_real_ else corr,
        matched_score_recovery = if (cell$design == 'matched') corr else NA_real_,
        false_active = false_active,
        validation_C = frozen_reverse_cindex(pred$risk_scores, time[va], status[va])),
        concordance, error = '')
    }
    list(metrics = do.call(rbind, cell_rows), comparison = cbind(cell, pair$comparison),
         program = cbind(cell, pair$program), pair = pair, validation_predictions = predictions,
         validation_time = time[va], validation_status = status[va], error = '')
    })
    rows[[i]] <- result$metrics
    comparisons[[i]] <- result$comparison
    programs[[i]] <- result$program
    saveRDS(list(cell = cell, result = result), file.path(checkpoint, sprintf('cell_%03d.rds', i)))
    if (any(!result$metrics$converged | !result$metrics$foldin_converged))
      warning(sprintf('Simulation cell %d/%d contains an unconverged fit or fold-in', i, nrow(design)))
    message(sprintf('Completed simulation cell %d/%d', i, nrow(design)))
  }
  comparison <- do.call(rbind, comparisons); program <- do.call(rbind, programs)
  if (is.null(comparison)) comparison <- data.frame()
  if (is.null(program)) program <- data.frame()
  list(simulation_metrics = do.call(rbind, rows), simulation_comparison = comparison,
       simulation_programs = program)
}

#' Execute the selected study using a fresh output directory
#' @param opts Validated options. @return Invisibly, aggregate tables.
#' @examples
#' # lb_yfb_runner_main(lb_yfb_runner_options(commandArgs(TRUE)))
lb_yfb_runner_main <- function(opts) {
  Sys.setenv(OPENBLAS_NUM_THREADS = '1', OMP_NUM_THREADS = '1', VECLIB_MAXIMUM_THREADS = '1')
  if (file.exists(opts$out)) stop('Output path already exists: ', opts$out)
  if (grepl('(^|/)data(/|$)', opts$out)) stop('Aggregate output must be outside data/')
  sources <- c('code/fit_ssbmf_lb.R', 'code/fit_ssbmf_yfb.R',
               'code/predict_ssbmf_lb.R', 'code/predict_ssbmf_yfb.R',
               'code/train_test_split.R', 'code/concordance_ci.R', 'code/lb_yfb_study.R')
  for (path in sources) source(path, local = .GlobalEnv)
  cfg <- yaml::read_yaml('config/globals.yml')
  K <- if (!is.null(opts$K)) opts$K else cfg$benchmark$k_merged_yfb_desurv
  max_iter <- cfg$cavi$max_iter; tol <- cfg$cavi$tol
  if (opts$mode == 'smoke') max_iter <- 2L
  center <- identical(opts$scaling, 'sd')
  run_id <- basename(normalizePath(opts$out, mustWork = FALSE))
  if (!grepl('^[A-Za-z0-9_.-]+$', run_id)) stop('Unsafe output directory name')
  checkpoint <- file.path('data', paste0('lb_yfb_study_', run_id))
  if (file.exists(checkpoint)) stop('Checkpoint path already exists: ', checkpoint)
  dir.create(opts$out, recursive = TRUE); dir.create(checkpoint, recursive = TRUE)
  if (opts$mode == 'simulation') sources <- c(sources, 'results/benchmark_sim/run_yfb_supervision_decision.R', 'code/simulate_multimodal_yfb.R', 'code/fit_yfb_supervised.R')
  provenance <- c(paste(names(opts)[!names(opts) %in% c('K', 'w_S')],
                        unlist(opts[!names(opts) %in% c('K', 'w_S')]), sep = '='),
                  paste0('w_S=', if (opts$mode == 'simulation') '1,10' else opts$w_S),
                  paste0('K=', K), paste0('max_iter=', max_iter), paste0('tol=', tol),
                  paste0('git_head=', system2('git', c('rev-parse', 'HEAD'), stdout = TRUE)),
                  paste0('git_status=', paste(system2('git', c('status', '--porcelain'), stdout = TRUE), collapse = '|')),
                  paste0('source_md5=', paste(names(tools::md5sum(sources)), tools::md5sum(sources), collapse = '|')))
  writeLines(provenance, file.path(opts$out, 'settings.txt'))
  if (opts$mode == 'simulation') {
    write.csv(lb_yfb_runner_design(), file.path(opts$out, 'design.csv'), row.names = FALSE)
    tables <- lb_yfb_runner_simulation(opts, K, max_iter, tol, center, checkpoint)
    lb_yfb_runner_write(tables, opts$out)
    if (any(nzchar(tables$simulation_metrics$error)))
      stop('One or more simulation cells failed; see simulation_metrics.csv and checkpoints')
    return(invisible(tables))
  }
  if (opts$mode == 'training') {
    d <- readRDS(opts$input)
  } else {
    d <- lb_yfb_runner_planted(1, .5)
    d$X <- d$X[1:40, 1:20, drop = FALSE]; d$time <- d$time[1:40]; d$status <- d$status[1:40]
    d$cohort <- rep(c('A', 'B'), each = 20); d$fold <- rep(1:5, 8)
  }
  needed <- c('X', 'time', 'status', 'cohort', 'fold')
  if (!all(needed %in% names(d))) stop('Input is missing: ', paste(setdiff(needed, names(d)), collapse = ', '))
  cv <- lb_yfb_study_cv(d$X, d$time, d$status, d$cohort, d$fold,
    opts$scaling, opts$prior_LF, opts$prior_beta, K, opts$w_S, max_iter, tol,
    checkpoint = checkpoint)
  saveRDS(cv, file.path(checkpoint, 'cv.rds'))
  lb_yfb_runner_write(list(cv_metrics = cv$metrics, cv_summary = cv$summary), opts$out)
  if (isTRUE(cv$failed)) stop('Training CV failed; see cv_metrics.csv and fold checkpoints')
  Y <- lb_yfb_study_scale(d$X, d$cohort, scaling = opts$scaling)
  pair <- tryCatch(lb_yfb_study_pair(Y, d$time, d$status, d$cohort, K, opts$prior_LF,
    opts$prior_beta, opts$w_S, max_iter, tol, center), error = identity)
  if (inherits(pair, 'error')) {
    saveRDS(list(error = conditionMessage(pair)), file.path(checkpoint, 'full_fit_error.rds'))
    lb_yfb_runner_write(list(full_metrics = data.frame(model = c('LB', 'YFB'),
      converged = FALSE, error = conditionMessage(pair))), opts$out)
    stop('Full-data fit failed after CV was saved: ', conditionMessage(pair))
  }
  saveRDS(list(cv = cv, pair = pair), file.path(checkpoint, 'patient_level.rds'))
  tables <- list(cv_metrics = cv$metrics, cv_summary = cv$summary,
                 full_metrics = pair$metrics, full_comparison = pair$comparison,
                 full_programs = pair$program)
  lb_yfb_runner_write(tables, opts$out)
  if (any(!cv$metrics$converged | !cv$metrics$foldin_converged) ||
      any(!pair$metrics$converged | !pair$metrics$foldin_converged))
    warning('At least one fit or fold-in did not converge; see the aggregate metrics')
  invisible(tables)
}

if (sys.nframe() == 0L) lb_yfb_runner_main(lb_yfb_runner_options(commandArgs(TRUE)))
