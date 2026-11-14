# ============================================================
# Script: run_share_targeted_weight.R
# Purpose: Fit and score share-targeted YFB simulation and training-data arms.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: survival, ebnm, yaml; S-9, S-13 and S-14 runners
# ============================================================
# At most two PSOCK workers and one BLAS thread per worker. Patient-level
# checkpoints belong in data/ and no external data are used here.

#' Allocate the fitted risk variance to false-active programs
#' @param components Patient by program matrix of fitted risk components.
#' @param false_active Logical vector identifying false-active programs.
#' @return Nonnegative variance ratio and signed covariance allocation.
#' @examples
#' # s16_false_variance(matrix(1:12, 4), c(TRUE, FALSE, FALSE))
s16_false_variance <- function(components, false_active) {
  stopifnot(is.matrix(components), length(false_active) == ncol(components),
            !anyNA(false_active), all(is.finite(components)))
  total <- rowSums(components)
  if (!any(false_active)) return(c(ratio = 0, signed = 0))
  denominator <- stats::var(total)
  if (!is.finite(denominator) || denominator <= 0)
    return(c(ratio = NA_real_, signed = NA_real_))
  false_score <- rowSums(components[, false_active, drop = FALSE])
  c(ratio = stats::var(false_score) / denominator,
    signed = stats::cov(false_score, total) / denominator)
}

#' Mark S-9 false-active programs, including all active programs under the null
#' @param active Logical activity vector.
#' @param correlations Signed correlations with the planted loading.
#' @param null Whether the training survival outcome is permuted.
#' @return Logical vector of false-active programs.
#' @examples
#' s16_false_active(c(TRUE, TRUE), c(.8, .2), FALSE)
s16_false_active <- function(active, correlations, null = FALSE) {
  stopifnot(length(active) == length(correlations), !anyNA(active))
  out <- active
  if (!null && any(is.finite(correlations))) {
    k <- which.max(replace(abs(correlations), !is.finite(correlations), -Inf))
    if (abs(correlations[k]) >= .7) out[k] <- FALSE
  }
  out
}

#' Fit a canonical S-16 posterior and retain its state for warm starts
#' @param Y,time,status Centered training matrix and outcomes.
#' @param weight Cox likelihood weight.
#' @param cfg Parsed project configuration with maximum sweeps set to 500.
#' @param initial Optional posterior state from a preceding fit.
#' @param fit_control Named extra fitter arguments for an approved opt-in control.
#' @return Fit, warnings and elapsed seconds.
#' @examples
#' # s16_fit(Y, time, status, 1, cfg)
s16_fit <- function(Y, time, status, weight, cfg, initial = NULL, fit_control = list()) {
  args <- list(Y = Y, time = time, status = status,
    K = cfg$benchmark$k_merged_yfb_desurv, max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
    prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
    w_S = weight, survival_block_weight = TRUE, beta_posterior = "meanfield",
    alpha = cfg$benchmark$alpha, verbose = FALSE, trace_supervision = TRUE,
    prune_preserve_posterior = TRUE, posterior_init = initial)
  ramp_omitted <- weight == 0 && !is.null(fit_control$survival_weight_ramp)
  if (ramp_omitted) fit_control$survival_weight_ramp <- NULL
  if (length(intersect(names(args), names(fit_control)))) stop("fit_control overrides fixed S-16 arguments")
  warnings <- if (ramp_omitted) "Survival weight ramp omitted for two-step weight zero" else character()
  start <- proc.time()[["elapsed"]]
  set.seed(42L)
  fit <- withCallingHandlers(do.call(fit_ssbmf_yfb, c(args, fit_control)), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(fit = fit, warnings = warnings, elapsed = proc.time()[["elapsed"]] - start)
}

#' Fit one target through at most three warm-started calibration passes
#' @param pilot Converged pilot fit.
#' @param Y,time,status Training input and outcomes.
#' @param target Target directional share.
#' @param cfg Project configuration.
#' @param fit_control Extra approved fitter controls.
#' @return Final fit and a record of all proposed and achieved shares.
#' @examples
#' # s16_calibrate(pilot, Y, time, status, .1, cfg)
s16_calibrate <- function(pilot, Y, time, status, target, cfg, fit_control = list(),
                          previous_cap = NULL) {
  proposed <- share_targeted_ws(pilot, Y, time, status, target, w_current = 1)
  if (proposed$no_signal_guard)
    return(list(fit = pilot, no_signal_guard = TRUE, capped = FALSE, passes = 0L,
                history = list(), achieved = NA_real_))
  fixed_loading <- proposed$pilot_loading
  records <- list(); current <- pilot
  for (pass in seq_len(3L)) {
    weight <- proposed$w_S
    if (proposed$capped && !is.null(previous_cap))
      return(list(capped_duplicate = TRUE, alias = previous_cap, history = records,
                  error = NULL))
    run <- s16_fit(Y, time, status, weight, cfg, current$posterior_state, fit_control)
    current <- run$fit
    row <- list(pass = pass, weight = weight, capped = proposed$capped,
                warnings = run$warnings, elapsed = run$elapsed,
                converged = isTRUE(current$history$converged))
    if (proposed$capped) row$cap_state <- current
    if (!row$converged) {
      records[[pass]] <- row
      return(list(fit = current, no_signal_guard = FALSE, capped = any(vapply(records, `[[`, logical(1), "capped")),
                  passes = pass, history = records, achieved = NA_real_, pilot_loading = fixed_loading))
    }
    next_proposal <- share_targeted_ws(current, Y, time, status, target,
      w_current = weight, pilot_loading = fixed_loading)
    row$ratio <- next_proposal$ratio
    row$achieved <- row$ratio / (1 + row$ratio)
    row$lost_activity <- next_proposal$lost_activity
    row$matching_correlation <- next_proposal$matching_correlation
    row$per_gene_share <- next_proposal$diagnostics$per_gene_median_share[next_proposal$selected_program]
    records[[pass]] <- row
    ratio_target <- target / (1 - target)
    if (row$ratio / ratio_target >= .5 && row$ratio / ratio_target <= 2) break
    proposed <- next_proposal
  }
  list(fit = current, pilot_loading = fixed_loading, no_signal_guard = FALSE,
       capped = any(vapply(records, `[[`, logical(1), "capped")), passes = length(records),
       history = records, achieved = tail(records, 1)[[1]]$achieved)
}

#' Score a simulated fit using registered recovery and false-active definitions
#' @param fit Fitted YFB state.
#' @param d Simulation data containing truth, Ytr, Yva and outcomes.
#' @param null Whether training survival is permuted.
#' @return One-row metrics including both risk-variance conventions.
#' @examples
#' # s16_sim_score(fit, d, FALSE)
s16_sim_score <- function(fit, d, null = FALSE) {
  F <- fit$EF; beta <- fit$EBeta
  sd_beta <- sqrt(pmax(fit$EBeta2 - beta^2, 0))
  active <- sd_beta > 0 & abs(beta) >= 1.96 * sd_beta
  correlations <- as.numeric(suppressWarnings(stats::cor(d$truth, F)))
  fa <- s16_false_active(active, correlations, null)
  component <- sweep(d$Ytr %*% F, 2, beta, "*")
  allocation <- s16_false_variance(component, fa)
  standardized <- abs(beta) * apply(d$Ytr %*% F, 2, stats::sd)
  risk <- as.numeric(d$Yva %*% F %*% beta)
  finite_corr <- is.finite(correlations)
  k <- if (any(finite_corr)) which.max(replace(abs(correlations), !finite_corr, -Inf)) else NA_integer_
  recovered <- !null && !is.na(k) && abs(correlations[k]) >= .7 && active[k] && beta[k] * correlations[k] > 0
  cv <- if (stats::sd(risk) > 0) frozen_reverse_cindex(risk, d$tv, d$ev) else NA_real_
  slope <- if (stats::sd(risk) > 0) unname(stats::coef(survival::coxph(
    survival::Surv(d$tv, d$ev) ~ risk, control = survival::coxph.control(timefix = FALSE)))) else NA_real_
  data.frame(recovered = recovered, cor_planted = if (is.na(k)) NA_real_ else abs(correlations[k]),
    false_active = sum(fa), false_risk_variance = unname(allocation["ratio"]),
    false_signed_allocation = unname(allocation["signed"]),
    max_false_standardized_beta = if (any(fa)) max(standardized[fa]) else 0,
    c_validation = cv, cal_slope = slope, converged = isTRUE(fit$history$converged))
}

#' Compute a two-step factorization and exact-times Cox reference
#' @param Y,time,status Centered training data and survival outcomes.
#' @param cfg Project configuration.
#' @param fit_control Approved fitter controls.
#' @return Fit, Cox coefficients, Wald activity and fit diagnostics.
#' @examples
#' # s16_twostep(Y, time, status, cfg)
s16_twostep <- function(Y, time, status, cfg, fit_control = list()) {
  run <- s16_fit(Y, time, status, 0, cfg, fit_control = fit_control)
  fit <- run$fit; live <- which(colSums(fit$EF^2) > 1e-8)
  if (!length(live)) stop("two-step fit has no live programs")
  Z <- Y %*% fit$EF[, live, drop = FALSE]
  warnings <- character()
  cf <- withCallingHandlers(survival::coxph(survival::Surv(time, status) ~ Z,
    control = survival::coxph.control(timefix = FALSE)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  coefficients <- as.numeric(stats::coef(cf)); se <- sqrt(diag(stats::vcov(cf)))
  if (any(!is.finite(coefficients)) || any(!is.finite(se)) || any(se <= 0))
    stop("two-step Cox coefficients or standard errors are not estimable")
  fit$EBeta[] <- 0; fit$EBeta2[] <- 0
  fit$w_S <- 0
  fit$EBeta[live] <- coefficients
  fit$EBeta2[live] <- coefficients^2 + se^2
  list(fit = fit, warnings = c(run$warnings, warnings), elapsed = run$elapsed,
       cox = cf)
}

#' Compute final-state directional and per-gene precision shares
#' @param fit Final posterior fit.
#' @param Y,time,status Exact training input and outcomes.
#' @param weight Effective Cox weight.
#' @return One row per program with activity and exact shares.
#' @examples
#' # s16_state_shares(fit, Y, time, status, 1)
s16_state_shares <- function(fit, Y, time, status, weight) {
  W <- calc_cox_taylor_yf(as.numeric(Y %*% fit$EF %*% fit$EBeta), time, status)$w
  scores <- Y %*% fit$EF
  beta_sd <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  active <- beta_sd > 0 & abs(fit$EBeta) >= 1.96 * beta_sd & colSums(fit$EF^2) > 1e-8
  tau <- rep(fit$Tau, length.out = ncol(Y))
  do.call(rbind, lapply(seq_len(ncol(fit$EF)), function(k) {
    f <- fit$EF[, k]; nf <- sqrt(sum(f^2)); l2 <- sum(fit$EL2[, k])
    if (nf == 0 || l2 == 0) return(data.frame(program = k, active = active[k], standardized_beta = fit$EBeta[k] * stats::sd(scores[, k]), directional_share = NA_real_, per_gene_share = NA_real_))
    f <- f / nf
    r <- weight * fit$EBeta2[k] * sum(W * as.numeric(Y %*% f)^2) / (l2 * sum(tau * f^2))
    rg <- weight * fit$EBeta2[k] * colSums(Y^2 * W) / (tau * l2)
    data.frame(program = k, active = active[k], standardized_beta = fit$EBeta[k] * stats::sd(scores[, k]), directional_share = r / (1 + r), per_gene_share = stats::median(rg / (1 + rg)))
  }))
}

#' Retain final-state quantities needed to audit precision shares
#' @param fit Posterior fit.
#' @param Y,time,status Exact centered input and outcomes.
#' @return Fit plus final Cox weights.
#' @examples
#' # s16_checkpoint(fit, Y, time, status)
s16_checkpoint <- function(fit, Y, time, status) {
  stopifnot(!is.null(fit$Tau), !is.null(fit$EL2), !is.null(fit$EF2),
            !is.null(fit$EBeta2))
  eta <- as.numeric(Y %*% fit$EF %*% fit$EBeta)
  list(fit = fit, cox_weights = calc_cox_taylor_yf(eta, time, status)$w,
       shares = s16_state_shares(fit, Y, time, status, if (is.null(fit$w_S)) 1 else fit$w_S))
}

#' Fit a single simulation dataset, preserving all target pass records
#' @param d Prepared S-13 or S-14 data list.
#' @param cfg Configuration.
#' @param targets Registered target shares, in ascending order.
#' @param null Permuted-survival null indicator.
#' @param fit_control Approved optional fitter controls.
#' @return Pilot, reference and target fits with metrics and warnings.
#' @examples
#' # s16_dataset(d, cfg, c(.01, .05, .1, .25, .5))
s16_dataset <- function(d, cfg, targets = c(.01, .05, .1, .25, .5),
                        null = FALSE, fit_control = list()) {
  stopifnot(identical(targets, sort(unique(targets))))
  safe <- function(expr) tryCatch(expr, error = function(e) list(error = conditionMessage(e)))
  pilot_run <- safe(s16_fit(d$Ytr, d$time, d$event, 1, cfg, fit_control = fit_control))
  pilot <- if (is.null(pilot_run$error)) pilot_run$fit else NULL
  pilot_entry <- if (is.null(pilot)) pilot_run else safe(list(
    checkpoint = s16_checkpoint(pilot, d$Ytr, d$time, d$event),
    metrics = s16_sim_score(pilot, d, null), run = pilot_run))
  reference <- safe(s16_twostep(d$Ytr, d$time, d$event, cfg, fit_control))
  reference_entry <- if (!is.null(reference$error)) reference else safe(list(
    run = reference, checkpoint = s16_checkpoint(reference$fit, d$Ytr, d$time, d$event),
    metrics = s16_sim_score(reference$fit, d, null)))
  result <- list(pilot = pilot_entry, twostep = reference_entry, targets = list())
  first_capped <- NULL
  for (target in targets) {
    key <- as.character(target)
    if (is.null(pilot) || !isTRUE(pilot$history$converged)) {
      result$targets[[key]] <- list(error = "pilot unavailable or nonconverged"); next
    }
    result$targets[[key]] <- safe({
      initial <- share_targeted_ws(pilot, d$Ytr, d$time, d$event, target, w_current = 1)
      if (initial$capped && !is.null(first_capped)) {
        list(capped_duplicate = TRUE, alias = first_capped)
      } else {
        calibrated <- s16_calibrate(pilot, d$Ytr, d$time, d$event, target, cfg, fit_control, first_capped)
        if (isTRUE(calibrated$capped_duplicate)) {
          calibrated
        } else {
        if (calibrated$capped) first_capped <- key
        calibrated$checkpoint <- s16_checkpoint(calibrated$fit, d$Ytr, d$time, d$event)
        calibrated$metrics <- s16_sim_score(calibrated$fit, d, null)
        rtruth <- suppressWarnings(as.numeric(stats::cor(d$truth, calibrated$fit$EF)))
        rk <- if (any(is.finite(rtruth))) which.max(replace(abs(rtruth), !is.finite(rtruth), -Inf)) else NA_integer_
        rt <- if (!is.na(rk) && is.null(reference$error)) suppressWarnings(abs(stats::cor(calibrated$fit$EF[, rk], reference$fit$EF))) else NA_real_
        calibrated$metrics$matched_twostep_correlation <- if (any(is.finite(rt))) max(rt, na.rm = TRUE) else NA_real_
        calibrated$capped_duplicate <- FALSE
        calibrated
        }
      }
    })
  }
  result
}

#' Fit secondary null arms at their paired signal dataset's final weights
#' @param signal Completed signal dataset result.
#' @param d Paired permuted-survival dataset.
#' @param cfg Configuration.
#' @param fit_control Approved optional fitter controls.
#' @return Target-indexed fit, checkpoint, diagnostics and metrics.
#' @examples
#' # s16_paired_null(signal, d, cfg)
s16_paired_null <- function(signal, d, cfg, fit_control = list()) {
  out <- list()
  for (key in names(signal$targets)) {
    arm <- signal$targets[[key]]
    if (!is.null(arm$error)) { out[[key]] <- arm; next }
    if (isTRUE(arm$capped_duplicate)) {
      out[[key]] <- list(capped_duplicate = TRUE, alias = arm$alias)
      next
    }
    if (!isTRUE(arm$fit$history$converged)) {
      out[[key]] <- list(error = "paired signal target did not converge")
      next
    }
    weight <- if (arm$no_signal_guard) 1 else tail(arm$history, 1)[[1]]$weight
    # Fixed paired weight: no recalibration on the null data.
    out[[key]] <- tryCatch({
      run <- s16_fit(d$Ytr, d$time, d$event, weight, cfg, fit_control = fit_control)
      list(weight = weight, run = run, pilot_loading = arm$pilot_loading,
        no_target = is.null(arm$pilot_loading), checkpoint = s16_checkpoint(run$fit, d$Ytr, d$time, d$event),
        metrics = s16_sim_score(run$fit, d, TRUE), capped_duplicate = FALSE)
    }, error = function(e) list(error = conditionMessage(e)))
  }
  out
}

#' Flatten a target result while retaining failures in the denominator
#' @param arm One calibrated or paired-null arm.
#' @param cell,seed,target,analysis Design identifiers.
#' @return One aggregate metrics row.
#' @examples
#' # s16_arm_row(arm, "weak", 1, .1, "primary")
s16_arm_row <- function(arm, cell, seed, target, analysis) {
  row <- data.frame(cell = cell, seed = seed, arm = if (is.na(target)) analysis else paste0("s", target), target = target, analysis = analysis,
    error = "", capped_duplicate = FALSE, alias = NA_character_, no_signal_guard = FALSE,
    capped = FALSE, passes = NA_integer_, weight = NA_real_, achieved = NA_real_,
    per_gene_share = NA_real_, share_program = NA_integer_, lost_activity = NA, elapsed = NA_real_, warnings = "",
    recovered = NA, cor_planted = NA_real_, false_active = NA_integer_,
    false_risk_variance = NA_real_, false_signed_allocation = NA_real_,
    max_false_standardized_beta = NA_real_, c_validation = NA_real_,
    cal_slope = NA_real_, matched_twostep_correlation = NA_real_, converged = FALSE, stringsAsFactors = FALSE)
  if (!is.null(arm$error)) { row$error <- arm$error; return(row) }
  if (isTRUE(arm$capped_duplicate)) {
    row$capped_duplicate <- TRUE; row$alias <- arm$alias
    row$error <- "capped duplicate; no additional fit"
    return(row)
  }
  row[names(arm$metrics)] <- arm$metrics
  if (!is.null(arm$checkpoint$shares)) {
    shares <- arm$checkpoint$shares
    if (!isTRUE(arm$no_target) && (!is.null(arm$pilot_loading) || any(shares$active))) {
      if (!is.null(arm$pilot_loading)) {
        cc <- suppressWarnings(abs(as.numeric(stats::cor(arm$pilot_loading, arm$checkpoint$fit$EF))))
        k <- which.max(replace(cc, !is.finite(cc), -Inf))
      } else k <- which.max(replace(abs(shares$standardized_beta), !shares$active, -Inf))
      row$share_program <- shares$program[k]; row$achieved <- shares$directional_share[k]
      row$per_gene_share <- shares$per_gene_share[k]
    }
  }
  if (!is.null(arm$run)) {
    row$weight <- arm$weight
    row$elapsed <- arm$run$elapsed
    row$warnings <- paste(arm$run$warnings, collapse = " | ")
  } else {
    row$no_signal_guard <- arm$no_signal_guard; row$capped <- arm$capped
    row$passes <- arm$passes
    if (arm$no_signal_guard) row$weight <- 1 else if (length(arm$history)) {
      last <- tail(arm$history, 1)[[1]]
      row$weight <- last$weight; if (length(arm$achieved) == 1L && is.finite(arm$achieved)) row$achieved <- arm$achieved
      if (!is.null(last$per_gene_share)) row$per_gene_share <- last$per_gene_share
      if (!is.null(last$lost_activity)) row$lost_activity <- last$lost_activity
      row$elapsed <- sum(vapply(arm$history, `[[`, numeric(1), "elapsed"))
      row$warnings <- paste(unlist(lapply(arm$history, `[[`, "warnings")), collapse = " | ")
    }
  }
  row
}

#' Evaluate the registered null gate with all ten replicates required
#' @param rows Aggregate rows for primary null arms.
#' @param targets Registered target shares.
#' @return One row per target, including counts and verdict.
#' @examples
#' # s16_null_gate(rows, c(.01, .05))
s16_null_gate <- function(rows, targets) {
  do.call(rbind, lapply(targets, function(target) {
    x <- rows[rows$cell == "weak_null" & rows$analysis == "primary" & !is.na(rows$target) & rows$target == target, , drop = FALSE]
    valid <- nrow(x) == 10L && length(unique(x$seed)) == 10L &&
      all(!x$capped_duplicate & x$error == "" &
            is.finite(x$false_active) & is.finite(x$false_risk_variance))
    mean_fa <- if (valid) mean(x$false_active) else NA_real_
    mean_var <- if (valid) mean(x$false_risk_variance) else NA_real_
    data.frame(target = target, n_rows = nrow(x), n_converged = sum(x$converged, na.rm = TRUE),
      mean_false_active = mean_fa, mean_false_risk_variance = mean_var,
      pass = valid && mean_fa <= 1 && mean_var < .01)
  }))
}

#' Atomically save a completed task checkpoint
#' @param object Task result with provenance.
#' @param path Destination in ignored data directory.
#' @return NULL, invisibly.
#' @examples
#' # s16_atomic_save(result, "data/result.rds")
s16_atomic_save <- function(object, path) {
  tmp <- tempfile(pattern = ".s16-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not atomically install checkpoint: ", path)
  invisible(NULL)
}

#' Execute a paired weak-signal/null task or low-variance task
#' @param task List with cell and seed.
#' @param cfg Configuration.
#' @param targets Allowed registered targets.
#' @param fit_control Opt-in fitter controls.
#' @param signature Checkpoint configuration and source signature.
#' @param local_dir Ignored checkpoint directory.
#' @return Checkpoint result, including explicit task errors.
#' @examples
#' # s16_screen_task(list(cell="weak", seed=1), cfg, .1, list(), "hash", "data/s16")
s16_screen_task <- function(task, cfg, targets, fit_control, signature, local_dir) {
  key <- sprintf("%s_s%02d", task$cell, task$seed)
  path <- file.path(local_dir, paste0(key, ".rds"))
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old$schema, 1L) || !identical(old$signature, signature) || !isTRUE(old$complete))
      stop("Stale or incomplete checkpoint: ", path)
    return(old)
  }
  start <- proc.time()[["elapsed"]]
  result <- tryCatch({
    if (task$cell == "weak") {
      signal_data <- s13_data(task$seed, "weak", FALSE)
      signal <- s16_dataset(signal_data, cfg, targets, FALSE, fit_control)
      null_data <- s13_data(task$seed, "weak", TRUE)
      null <- s16_dataset(null_data, cfg, targets, TRUE, fit_control)
      paired <- s16_paired_null(signal, null_data, cfg, fit_control)
      list(signal = signal, null = null, paired = paired)
    } else {
      d <- s14_data(task$seed, .0005, "strong")
      list(lowvariance = s16_dataset(d, cfg, targets, FALSE, fit_control))
    }
  }, error = function(e) list(error = conditionMessage(e)))
  out <- list(schema = 1L, signature = signature, complete = TRUE, task = task,
              elapsed = proc.time()[["elapsed"]] - start, result = result)
  s16_atomic_save(out, path)
  out
}

#' Report per-target rows from a completed task, including task-level failures
#' @param checkpoint Completed task checkpoint.
#' @param targets Target shares.
#' @return Aggregate data frame.
#' @examples
#' # s16_task_rows(checkpoint, c(.01, .05))
s16_task_rows <- function(checkpoint, targets) {
  task <- checkpoint$task; res <- checkpoint$result
  cells <- if (task$cell == "weak") c("weak", "weak_null", "weak_null") else "lowvariance"
  analyses <- if (task$cell == "weak") c("primary", "primary", "paired_signal_weight") else "primary"
  sets <- if (task$cell == "weak") c("signal", "null", "paired") else "lowvariance"
  do.call(rbind, lapply(seq_along(sets), function(i) do.call(rbind, lapply(targets, function(target) {
    arm <- if (!is.null(res$error)) list(error = res$error) else {
      z <- res[[sets[i]]]
      if (sets[i] != "paired") z <- z$targets
      z[[as.character(target)]]
    }
    if (is.null(arm)) arm <- list(error = "target result missing")
    s16_arm_row(arm, cells[i], task$seed, target, analyses[i])
  }))))
}

#' Execute an off-grid smoke dataset using the ordinary calibration path
#' @param d Prepared seed-99 synthetic dataset.
#' @param cfg Configuration.
#' @param targets Registered targets to exercise.
#' @param fit_control Approved opt-in controls.
#' @return The same dataset result as the screen.
#' @examples
#' # s16_smoke_run(d, cfg, .01, list())
s16_smoke_run <- function(d, cfg, targets, fit_control) {
  s16_dataset(d, cfg, targets, FALSE, fit_control)
}

#' Run the registered screen with reproducible cached dataset tasks
#' @param config_path Optional RDS containing fit_control and allowed_targets.
#' @param workers Integer number of PSOCK workers, at most two.
#' @return Invisible aggregate results.
#' @examples
#' # s16_screen(NULL, 2L)
s16_screen <- function(config_path = NULL, workers = 1L) {
  if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 2L || workers != as.integer(workers))
    stop("workers must be 1 or 2")
  settings <- if (is.null(config_path)) list() else readRDS(config_path)
  registered <- c(.01, .05, .1, .25, .5)
  targets <- if (is.null(settings$allowed_targets)) registered else settings$allowed_targets
  if (!is.numeric(targets) || anyNA(targets) || any(!targets %in% registered) ||
      !identical(targets, sort(unique(targets)))) stop("allowed_targets must be an ordered subset of registered targets")
  fit_control <- if (is.null(settings$fit_control)) list() else settings$fit_control
  if (!is.list(fit_control) || (length(fit_control) && (is.null(names(fit_control)) || any(!nzchar(names(fit_control))))))
    stop("fit_control must be a named list")
  source_files <- c(list.files("code", pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE), "results/benchmark_sim/run_share_targeted_weight.R", "code/share_targeted_ws.R",
    "code/fit_ssbmf_yfb.R", "results/benchmark_sim/run_yfb_survival_weight_p2000.R",
    "results/benchmark_sim/run_yfb_supervision_strength.R", "results/benchmark_sim/run_yfb_near_threshold.R",
    "results/benchmark_sim/run_yfb_supervision_decision.R",
    list.files("results/benchmark_sim", pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE),
    list.files("config", recursive = TRUE, full.names = TRUE))
  signature <- list(files = tools::md5sum(source_files), settings = settings, targets = targets)
  local_dir <- "data/share_targeted_ws_2026-11-13/screen"
  out_dir <- "results/benchmark_sim/outputs/share_targeted_ws_2026-11-13"
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- yaml::read_yaml("config/globals.yml"); cfg$cavi$max_iter <- 500L
  tasks <- c(lapply(1:10, function(seed) list(cell = "weak", seed = seed)),
             lapply(1:10, function(seed) list(cell = "lowvariance", seed = seed)))
  Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
  cl <- parallel::makePSOCKcluster(workers); on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
    source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
    source("code/concordance_ci.R"); source("code/share_targeted_ws.R")
    source("results/benchmark_sim/run_yfb_supervision_strength.R")
    source("results/benchmark_sim/run_yfb_near_threshold.R")
    source("results/benchmark_sim/run_share_targeted_weight.R"); NULL
  })
  parallel::clusterExport(cl, c("cfg", "targets", "fit_control", "signature", "local_dir"), envir = environment())
  checkpoints <- parallel::parLapplyLB(cl, tasks, function(task)
    s16_screen_task(task, cfg, targets, fit_control, signature, local_dir), chunk.size = 1L)
  rows <- do.call(rbind, lapply(checkpoints, s16_task_rows, targets = targets))
  references <- do.call(rbind, lapply(checkpoints, function(checkpoint) {
    task <- checkpoint$task; res <- checkpoint$result
    cells <- if (task$cell == "weak") c("weak", "weak_null") else "lowvariance"
    sets <- if (task$cell == "weak") c("signal", "null") else "lowvariance"
    do.call(rbind, lapply(seq_along(sets), function(i) do.call(rbind, lapply(c("pilot", "twostep"), function(which) {
      if (!is.null(res$error)) return(s16_arm_row(list(error = res$error), cells[i], task$seed, NA_real_, which))
      ref <- res[[sets[i]]][[which]]
      if (!is.null(ref$error)) return(s16_arm_row(ref, cells[i], task$seed, NA_real_, which))
      arm <- list(metrics = ref$metrics, run = ref$run, checkpoint = ref$checkpoint, weight = if (which == "pilot") 1 else 0)
      s16_arm_row(arm, cells[i], task$seed, NA_real_, which)
    }))))
  }))
  rows <- rbind(rows, references)
  warned <- rows[nzchar(rows$warnings) | nzchar(rows$error), , drop = FALSE]
  if (nrow(warned)) for (i in seq_len(nrow(warned)))
    message(sprintf("%s seed %s %s: %s %s", warned$cell[i], warned$seed[i], warned$analysis[i], warned$error[i], warned$warnings[i]))
  gate <- s16_null_gate(rows, targets)
  utils::write.csv(rows, file.path(out_dir, "screen_metrics.csv"), row.names = FALSE)
  utils::write.csv(gate, file.path(out_dir, "screen_null_gate.csv"), row.names = FALSE)
  utils::write.csv(data.frame(cell = vapply(tasks, `[[`, character(1), "cell"),
    seed = vapply(tasks, `[[`, integer(1), "seed"),
    elapsed = vapply(checkpoints, `[[`, numeric(1), "elapsed"),
    error = vapply(checkpoints, function(x) if (is.null(x$result$error)) "" else x$result$error, character(1))),
    file.path(out_dir, "screen_task_status.csv"), row.names = FALSE)
  invisible(list(metrics = rows, gate = gate))
}

#' Fit full or fold-specific real training arms using one training center
#' @param fold 0 for full training, 1 to 5 for held-out fold.
#' @param d Prepared S-15 training data.
#' @param cfg Configuration with working rank.
#' @param targets Eligible target shares.
#' @param fit_control Approved fitter controls.
#' @return Fits, full and held-out risks, and fold metrics.
#' @examples
#' # s16_real_task(1, d, cfg, .1, list())
s16_real_task <- function(fold, d, cfg, targets, fit_control) {
  tr <- if (fold == 0L) rep(TRUE, nrow(d$Y)) else d$outer != fold
  te <- !tr
  center <- colMeans(d$Y[tr, , drop = FALSE])
  Y <- sweep(d$Y, 2, center)
  ytr <- Y[tr, , drop = FALSE]
  safe <- function(expr) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) list(error = conditionMessage(e)))
    result$additional_warnings <- warnings
    result
  }
  pilot <- safe(s16_fit(ytr, d$time[tr], d$status[tr], 1, cfg, fit_control = fit_control))
  two <- safe(s16_twostep(ytr, d$time[tr], d$status[tr], cfg, fit_control))
  arms <- list(w1 = if (!is.null(pilot$error)) pilot else list(fit = pilot$fit, run = pilot, weight = 1),
               twostep = if (!is.null(two$error)) two else list(fit = two$fit, run = two, weight = 0))
  first_capped <- NULL
  for (target in targets) {
    key <- paste0("s", target)
    if (!is.null(pilot$error) || !isTRUE(pilot$fit$history$converged)) {
      arms[[key]] <- list(error = "pilot unavailable or nonconverged"); next
    }
    arms[[key]] <- safe({
      arm <- s16_calibrate(pilot$fit, ytr, d$time[tr], d$status[tr], target, cfg, fit_control, first_capped)
      if (isTRUE(arm$capped_duplicate)) {
        list(error = paste("capped duplicate of", arm$alias))
      } else {
        if (arm$capped) first_capped <- key
        arm$weight <- if (arm$no_signal_guard) 1 else tail(arm$history, 1)[[1]]$weight
        arm
      }
    })
  }
  for (key in names(arms)) {
    arm <- arms[[key]]
    if (!is.null(arm$error)) next
    checkpoint_result <- safe(s16_checkpoint(arm$fit, ytr, d$time[tr], d$status[tr]))
    if (!is.null(checkpoint_result$error)) {
      arm$error <- paste("checkpoint:", checkpoint_result$error)
      arm$scoring_warnings <- checkpoint_result$additional_warnings
      arms[[key]] <- arm; next
    }
    arm$checkpoint <- checkpoint_result
    evaluated <- safe({
      fit <- arm$fit
      risk <- as.numeric(Y %*% fit$EF %*% fit$EBeta)
      C <- if (fold > 0L && stats::sd(risk[te]) > 0)
        frozen_reverse_cindex(risk[te], d$time[te], d$status[te]) else NA_real_
      cvpl <- if (fold > 0L) calc_cox_taylor_yf(risk, d$time, d$status)$logPL -
        calc_cox_taylor_yf(risk[tr], d$time[tr], d$status[tr])$logPL else NA_real_
      list(risk = risk, C = C, cvpl = cvpl)
    })
    arm$scoring_warnings <- evaluated$additional_warnings
    if (!is.null(evaluated$error)) {
      arm$error <- paste("scoring:", evaluated$error)
    } else {
      arm$risk <- evaluated$risk
      arm$C <- evaluated$C; arm$cvpl <- evaluated$cvpl
    }
    arms[[key]] <- arm
  }
  list(fold = fold, center = center, training = tr, arms = arms)
}

#' Flatten one real-data fold result, preserving missing and failed arms
#' @param object Real task result or error.
#' @param fold Fold identifier.
#' @param arm_names Expected arm labels.
#' @return One row per arm.
#' @examples
#' # s16_real_rows(task, 1, c("w1", "twostep"))
s16_real_rows <- function(object, fold, arm_names) {
  do.call(rbind, lapply(arm_names, function(key) {
    arm <- if (!is.null(object$error)) object else object$arms[[key]]
    if (is.null(arm)) arm <- list(error = "missing arm")
    row <- data.frame(fold = fold, arm = key, weight = NA_real_, converged = FALSE,
      C = NA_real_, cvpl = NA_real_, achieved = NA_real_, per_gene_share = NA_real_,
      share_program = NA_integer_, n_iter = NA_integer_, passes = NA_integer_,
      capped = FALSE, no_signal_guard = FALSE, lost_activity = NA, warnings = "",
      elapsed = NA_real_, error = "", stringsAsFactors = FALSE)
    if (!is.null(arm$error)) {
      row$error <- arm$error
      if (!is.null(arm$fit)) {
        row$converged <- isTRUE(arm$fit$history$converged)
        row$n_iter <- arm$fit$history$n_iter
      }
      row$warnings <- paste(c(arm$additional_warnings, arm$scoring_warnings), collapse = " | ")
      return(row)
    }
    row$weight <- arm$weight; row$converged <- isTRUE(arm$fit$history$converged)
    row$C <- arm$C; row$cvpl <- arm$cvpl; row$n_iter <- arm$fit$history$n_iter
    sh <- arm$checkpoint$shares
    if (!is.null(arm$pilot_loading)) {
      corr <- suppressWarnings(abs(as.numeric(stats::cor(arm$pilot_loading, arm$fit$EF))))
      k <- if (any(is.finite(corr))) which.max(replace(corr, !is.finite(corr), -Inf)) else NA_integer_
    } else k <- if (any(sh$active)) which.max(replace(abs(sh$standardized_beta), !sh$active, -Inf)) else NA_integer_
    if (!is.na(k)) { row$share_program <- k; row$achieved <- sh$directional_share[k]; row$per_gene_share <- sh$per_gene_share[k] }
    if (!is.null(arm$run)) {
      row$warnings <- paste(arm$run$warnings, collapse = " | "); row$elapsed <- arm$run$elapsed
    } else {
      row$passes <- arm$passes; row$capped <- arm$capped
      row$no_signal_guard <- arm$no_signal_guard
      if (length(arm$history)) {
        row$warnings <- paste(unlist(lapply(arm$history, `[[`, "warnings")), collapse = " | ")
        row$elapsed <- sum(vapply(arm$history, `[[`, numeric(1), "elapsed"))
        last <- tail(arm$history, 1)[[1]]
        if (!is.null(last$lost_activity)) row$lost_activity <- last$lost_activity
      }
    }
    row$warnings <- paste(c(row$warnings[nzchar(row$warnings)], arm$additional_warnings, arm$run$additional_warnings, arm$scoring_warnings), collapse = " | ")
    row
  }))
}

#' Summarize registered real-data eligibility and select by mean fold CVPL
#' @param rows All full and five-fold arm metrics.
#' @return Summary and selected arm.
#' @examples
#' # s16_real_summary(rows)
s16_real_summary <- function(rows) {
  arms <- unique(rows$arm)
  summary <- do.call(rbind, lapply(arms, function(a) {
    x <- rows[rows$arm == a & rows$fold > 0, , drop = FALSE]
    full <- rows[rows$arm == a & rows$fold == 0, , drop = FALSE]
    complete <- nrow(x) == 5L && length(unique(x$fold)) == 5L &&
      all(is.finite(x$cvpl)) && all(is.finite(x$C)) && all(x$error == "")
    data.frame(arm = a, folds_converged = sum(x$converged),
      full_converged = nrow(full) == 1L && isTRUE(full$converged),
      eligible = complete && sum(x$converged) >= 4L && nrow(full) == 1L && isTRUE(full$converged) && full$error == "",
      mean_cvpl = if (complete) mean(x$cvpl) else NA_real_,
      mean_C = if (complete) mean(x$C) else NA_real_,
      se_C = if (complete) stats::sd(x$C) / sqrt(5) else NA_real_,
      median_weight = if (nrow(x) == 5 && all(is.finite(x$weight))) stats::median(x$weight) else NA_real_)
  }))
  candidate <- summary[summary$eligible & summary$arm != "twostep", , drop = FALSE]
  selected <- NA_character_
  if (nrow(candidate)) {
    near <- candidate[candidate$mean_cvpl >= max(candidate$mean_cvpl) - .5, , drop = FALSE]
    target <- suppressWarnings(as.numeric(sub("^s", "", near$arm)))
    target[near$arm == "w1"] <- 0
    selected <- near$arm[order(near$median_weight, target)][1]
  }
  list(summary = summary, selected = selected)
}

#' Compute paired held-out C differences from both fixed references
#' @param rows Real fold metrics.
#' @return Paired mean and standard error for each arm and reference.
#' @examples
#' # s16_paired_C(rows)
s16_paired_C <- function(rows) {
  do.call(rbind, lapply(unique(rows$arm), function(a) do.call(rbind, lapply(c("w1", "twostep"), function(ref) {
    x <- rows[rows$arm == a & rows$fold > 0, c("fold", "C")]
    y <- rows[rows$arm == ref & rows$fold > 0, c("fold", "C")]
    z <- merge(x, y, by = "fold", suffixes = c("", "_ref"))
    delta <- z$C - z$C_ref; ok <- nrow(z) == 5L && all(is.finite(delta))
    data.frame(arm = a, reference = ref, n = sum(is.finite(delta)),
      mean_difference = if (ok) mean(delta) else NA_real_,
      se_difference = if (ok) stats::sd(delta) / sqrt(5) else NA_real_)
  }))))
}

source("results/benchmark_sim/share_targeted_structure.R", local = TRUE)

#' Run real-data five-fold and full training fits, or describe eligible programs
#' @param mode Either real or structure.
#' @param config_path Settings RDS with fit_control and allowed_targets.
#' @param workers One or two PSOCK workers.
#' @return Invisible result tables.
#' @examples
#' # s16_real("real", "data/s16-settings.rds", 2L)
s16_real <- function(mode, config_path, workers = 1L) {
  if (is.null(config_path)) stop("real and structure modes require a reviewed settings RDS")
  if (length(workers) != 1L || is.na(workers) || !workers %in% 1:2) stop("workers must be 1 or 2")
  settings <- readRDS(config_path)
  registered <- c(.01, .05, .1, .25, .5)
  targets <- settings$allowed_targets
  if (!is.numeric(targets) || anyNA(targets) || any(!targets %in% registered) ||
      !identical(targets, sort(unique(targets)))) stop("settings must contain sorted allowed_targets")
  controls <- if (is.null(settings$fit_control)) list() else settings$fit_control
  if (length(targets)) {
    gate_path <- "results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/screen_null_gate.csv"
    if (!file.exists(gate_path)) stop("screen null gate is required before real targets")
    gate <- utils::read.csv(gate_path)
    if (!all(c("target", "pass") %in% names(gate)) || anyDuplicated(gate$target) ||
        anyNA(match(targets, gate$target)) || any(is.na(gate$pass[match(targets, gate$target)])) ||
        !all(gate$pass[match(targets, gate$target)])) stop("all real targets must pass the screen null gate")
    screen_paths <- file.path("data/share_targeted_ws_2026-11-13/screen", sprintf("weak_s%02d.rds", 1:10))
    if (!all(file.exists(screen_paths))) stop("all weak screen checkpoints are required to verify the gate")
    screen <- lapply(screen_paths, readRDS)
    if (any(!vapply(screen, function(x) {
      hashes <- x$signature$files
      length(hashes) > 0L && identical(hashes, tools::md5sum(names(hashes)))
    }, logical(1)))) stop("screen source changed before real-data gate verification")
    if (any(!vapply(screen, function(x) isTRUE(x$complete) && identical(x$schema, 1L) &&
      all(targets %in% x$signature$targets) &&
      identical(if (is.null(x$signature$settings$fit_control)) list() else x$signature$settings$fit_control, controls), logical(1))))
      stop("screen checkpoints do not match real fitter controls or targets")
    screen_rows <- do.call(rbind, lapply(screen, s16_task_rows, targets = targets))
    if (!all(s16_null_gate(screen_rows, targets)$pass)) stop("screen checkpoints do not reproduce the passing gate")
  }
  source("results/benchmark_sim/run_yfb_supervised_real_refit.R"); s11_load()
  source("code/share_targeted_ws.R")
  d <- s15_data()
  if (nrow(d$Y) != 273L || ncol(d$Y) != 2218L || sum(d$status) != 139L ||
      !identical(sort(unique(as.integer(d$outer))), 1:5))
    stop("prepared data do not match registered 273 patients, 139 events, 2218 genes")
  cfg <- yaml::read_yaml("config/globals.yml"); cfg$cavi$max_iter <- 500L
  if (cfg$benchmark$k_working_yfb != 5L) stop("working rank is not 5")
  cfg$benchmark$k_merged_yfb_desurv <- cfg$benchmark$k_working_yfb
  files <- c("results/benchmark_sim/run_share_targeted_weight.R", "code/share_targeted_ws.R",
             "code/fit_ssbmf_yfb.R", "results/benchmark_sim/run_yfb_supervised_real_refit.R", "config/globals.yml")
  files <- unique(c(files, list.files("code", pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE),
    list.files("results/benchmark_sim", pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE),
    list.files("config", recursive = TRUE, full.names = TRUE), LOCAL))
  signature <- list(files = tools::md5sum(files), settings = settings)
  local_dir <- "data/share_targeted_ws_2026-11-13/real"
  out_dir <- "results/benchmark_sim/outputs/share_targeted_ws_2026-11-13"
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(local_dir, sprintf("fold_%d.rds", 0:5))
  if (mode == "real") {
    Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
    cl <- parallel::makePSOCKcluster(workers); on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
      source("results/benchmark_sim/run_yfb_supervised_real_refit.R"); s11_load()
      source("code/share_targeted_ws.R"); source("results/benchmark_sim/run_share_targeted_weight.R"); NULL
    })
    parallel::clusterExport(cl, c("d", "cfg", "targets", "controls", "signature", "paths"), envir = environment())
    objects <- parallel::parLapplyLB(cl, 0:5, function(fold) {
      path <- paths[fold + 1L]
      if (file.exists(path)) {
        x <- readRDS(path)
        if (!identical(x$schema, 1L) || !identical(x$signature, signature) || !isTRUE(x$complete))
          stop("stale real checkpoint: ", path)
        return(x$result)
      }
      result <- tryCatch(s16_real_task(fold, d, cfg, targets, controls), error = function(e) list(error = conditionMessage(e)))
      s16_atomic_save(list(schema = 1L, signature = signature, complete = TRUE, result = result), path)
      result
    }, chunk.size = 1L)
    rows <- do.call(rbind, lapply(0:5, function(i) s16_real_rows(objects[[i + 1L]], i, c("w1", "twostep", paste0("s", targets)))))
    result <- s16_real_summary(rows)
    utils::write.csv(rows, file.path(out_dir, "real_fold_metrics.csv"), row.names = FALSE)
    utils::write.csv(result$summary, file.path(out_dir, "real_summary.csv"), row.names = FALSE)
    utils::write.csv(data.frame(selected = result$selected), file.path(out_dir, "real_selection.csv"), row.names = FALSE)
    utils::write.csv(s16_paired_C(rows), file.path(out_dir, "real_paired_c.csv"), row.names = FALSE)
    return(invisible(result))
  }
  if (!all(file.exists(paths))) stop("real full and fold checkpoints are required")
  ck <- lapply(paths, readRDS)
  if (any(!vapply(ck, function(x) identical(x$signature, signature) && isTRUE(x$complete), logical(1))))
    stop("real checkpoint signature mismatch")
  objects <- lapply(ck, `[[`, "result")
  rows <- do.call(rbind, lapply(0:5, function(i) s16_real_rows(objects[[i + 1L]], i, c("w1", "twostep", paste0("s", targets)))))
  summary <- s16_real_summary(rows)$summary
  eligible <- summary$arm[summary$eligible]
  if (!"twostep" %in% eligible) stop("two-step reference is not eligible for structure")
  structures <- lapply(eligible, function(a) tryCatch(s16_structure_arm(objects[[1]], objects[-1], d, a),
                                                       error = function(e) list(error = conditionMessage(e))))
  names(structures) <- eligible
  status <- data.frame(arm = eligible, error = vapply(structures, function(x) if (is.null(x$error)) "" else x$error, character(1)))
  utils::write.csv(status, file.path(out_dir, "structure_status.csv"), row.names = FALSE)
  valid <- structures[vapply(structures, function(x) is.null(x$error), logical(1))]
  if (length(valid)) {
    utils::write.csv(do.call(rbind, lapply(valid, `[[`, "program")), file.path(out_dir, "structure_programs.csv"), row.names = FALSE)
    for (a in names(valid)) {
      utils::write.csv(valid[[a]]$correlations, file.path(out_dir, paste0("structure_cor_", a, ".csv")))
      utils::write.csv(valid[[a]]$matches, file.path(out_dir, paste0("structure_match_", a, ".csv")), row.names = FALSE)
      utils::write.csv(valid[[a]]$fold_matches, file.path(out_dir, paste0("structure_folds_", a, ".csv")), row.names = FALSE)
    }
  }
  source("code/pathway_enrichment.R")
  sets <- s16_structure_genesets(
    "results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds",
    "data/PDAC_data/original/cmbSubtypes.RData")
  prepared <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  annotation <- do.call(rbind, lapply(names(valid), function(a) {
    x <- s16_structure_annotation(objects[[1]]$arms[[a]]$fit, prepared$sets$D4, sets)
    x$arm <- a
    x
  }))
  utils::write.csv(annotation, file.path(out_dir, "structure_annotation.csv"), row.names = FALSE)
  utils::write.csv(data.frame(signature = names(sets), n_signature = lengths(sets),
    n_measured = vapply(sets, function(x) length(intersect(x, prepared$sets$D4)), integer(1))),
    file.path(out_dir, "structure_signature_overlap.csv"), row.names = FALSE)
  if (any(nzchar(status$error))) stop("Program structure failed for an eligible arm. See structure_status.csv")
  invisible(structures)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || !args[1] %in% c("screen", "real", "structure", "smoke") || length(args) > 3L)
    stop("Usage: Rscript results/benchmark_sim/run_share_targeted_weight.R screen|real|structure [workers] [settings.rds]")
  workers <- if (length(args) >= 2L) as.integer(args[2]) else 1L
  settings <- if (length(args) >= 3L) args[3] else NULL
  if (args[1] == "screen") s16_screen(settings, workers) else if (args[1] == "smoke") {
    Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
    source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
    source("code/concordance_ci.R"); source("code/share_targeted_ws.R")
    source("results/benchmark_sim/run_yfb_supervision_strength.R")
    cfg <- yaml::read_yaml("config/globals.yml"); cfg$cavi$max_iter <- 500L
    options <- if (is.null(settings)) list() else readRDS(settings)
    controls <- if (is.null(options$fit_control)) list() else options$fit_control
    targets <- if (is.null(options$allowed_targets)) .01 else options$allowed_targets
    d <- s13_data(99L, "weak", FALSE)
    result <- s16_smoke_run(d, cfg, targets, controls)
    dir.create("data/share_targeted_ws_2026-11-13", recursive = TRUE, showWarnings = FALSE)
    s16_atomic_save(result, "data/share_targeted_ws_2026-11-13/smoke_seed99.rds")
  } else s16_real(args[1], settings, workers)
}
