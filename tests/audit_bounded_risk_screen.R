# ============================================================
# Script: audit_bounded_risk_screen.R
# Purpose: Audit completed S-17 screen checkpoints without fitting models.
# Author: Andrew Walther
# Created: 2026-11-14
# Dependencies: survival, ebnm, ashr, yaml; S-9, S-13, S-15 and S-17 runners
# ============================================================
# One CPU and one BLAS thread. The reconstructed inputs and checkpoints contain
# patient-level data and stay in memory or under the ignored data/ directory.
# Only case-level check outcomes and aggregate counts are written to outputs/.
#' Compare numeric objects to a strict relative tolerance
#' @param x,y Numeric objects to compare.
#' @return TRUE when dimensions and finite values agree within tolerance.
#' @examples
#' audit_equal(c(1, 2), c(1, 2))
audit_equal <- function(x, y) {
  is.numeric(x) && is.numeric(y) && identical(dim(x), dim(y)) && length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) && isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-10, check.attributes = FALSE))
}
#' Prepare the three registered training datasets without fitting a model
#' @return Named list of centered training matrices, event times and statuses.
#' @examples
#' # audit_inputs()
audit_inputs <- function() {
  source("results/benchmark_sim/run_yfb_survival_weight_p2000.R", local = .GlobalEnv)
  s9_load(); source("results/benchmark_sim/run_yfb_supervision_strength.R", local = .GlobalEnv)
  source("results/benchmark_sim/run_yfb_supervised_real_refit.R", local = .GlobalEnv)
  s11_load(); real <- s15_data()
  if (nrow(real$Y) != 273L || ncol(real$Y) != 2218L || sum(real$status) != 139L || any(!is.finite(real$Y))) stop("Reconstructed real training input does not match registration")
  s9 <- yfb_decision_data(1, 2000L, .01, truth_mode = "shared")
  weak <- s13_data(1, "weak", FALSE)
  list(real = list(Y = scale(real$Y, center = TRUE, scale = FALSE), time = real$time, status = real$status),
       s9 = list(Y = s9_zscore(s9$Y), time = s9$time, status = s9$event), weak = list(Y = weak$Ytr, time = weak$time, status = weak$event))
}
#' Track each fit against the selected pilot loading without refitting
#' @param object Completed S-17 case checkpoint.
#' @return Aggregate per-fit program matches and achieved shares.
#' @examples
#' # audit_track(checkpoint)
audit_track <- function(object) {
  dataset <- object$dataset; pilot <- object$arms$pilot
  base_loading <- NULL
  if (dataset != "s9" && length(pilot$fits) && length(pilot$checkpoints)) {
    p <- pilot$checkpoints[[1]]$shares
    eligible <- !is.na(p$active) & p$active & is.finite(p$standardized_beta)
    if (any(eligible)) {
      k <- p$program[which.max(replace(abs(p$standardized_beta), !eligible, -Inf))]
      base_loading <- pilot$fits[[1]]$EF[, k]
    }
  }
  rows <- list(); for (name in names(object$arms)) {
    arm <- object$arms[[name]]
    if (!length(arm$fits)) next
    for (i in seq_along(arm$fits)) {
      fit <- arm$fits[[i]]; shares <- arm$checkpoints[[i]]$shares
      reference <- if (!is.null(arm$calibration$pilot_loading)) arm$calibration$pilot_loading else base_loading
      reason <- if (dataset == "s9") "fixed-weight diagnostic has no pilot" else
        if (is.null(reference)) "pilot has no active selected program" else ""
      correlations <- rep(NA_real_, ncol(fit$EF))
      if (!nzchar(reason) && all(is.finite(reference)) && stats::sd(reference) > 0) {
        correlations <- vapply(seq_len(ncol(fit$EF)), function(k) {
          x <- fit$EF[, k]
          if (!all(is.finite(x)) || stats::sd(x) == 0) NA_real_ else abs(stats::cor(reference, x))
        }, numeric(1))
      }
      if (!nzchar(reason) && !any(is.finite(correlations))) reason <- "no finite loading correlation"
      selected <- if (nzchar(reason)) NA_integer_ else which.max(replace(correlations, !is.finite(correlations), -Inf))
      selected_row <- if (is.na(selected)) shares[FALSE, , drop = FALSE] else shares[shares$program == selected, , drop = FALSE]
      if (!is.na(selected) && nrow(selected_row) != 1L) { reason <- "selected program share missing or duplicated"; selected <- NA_integer_ }
      if (!is.na(selected) && !is.finite(selected_row$directional_share)) reason <- "directional share undefined"
      rows[[length(rows) + 1L]] <- data.frame(dataset = dataset, family = object$family, sigma = object$sigma,
        arm = name, pass_number = i, final = i == length(arm$fits), weight = arm$weights[i],
        guard = isTRUE(arm$guard), capped = arm$capped[i], selected_program = selected,
        correlation = if (is.na(selected)) NA_real_ else correlations[selected],
        active = if (is.na(selected)) NA else selected_row$active,
        directional_share = if (is.na(selected)) NA_real_ else selected_row$directional_share,
        per_gene_share = if (is.na(selected)) NA_real_ else selected_row$per_gene_share, reason = reason)
    }
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}
#' Audit all registered cases and save only non-patient-level results
#' @return Invisibly, audit rows; exits with an error after writing if checks fail.
#' @examples
#' # audit_screen()
audit_screen <- function() {
  source("results/benchmark_sim/run_bounded_risk_screen.R", local = .GlobalEnv)
  inputs <- audit_inputs()
  cfg <- yaml::read_yaml("config/globals.yml")
  output <- "results/benchmark_sim/outputs/bounded_risk_2026-11-14"
  private <- "data/bounded_risk_2026-11-14/screen"
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(dataset = c("real", "s9", "weak"), family = c("normal", "laplace"), sigma = c(.1, .25, .5, 1), stringsAsFactors = FALSE)
  rows <- list(); tracking <- list(); reference_signature <- NULL
  #' Append one aggregate case or fit check
  #' @param dataset,family,sigma Case identifiers.
  #' @param arm,pass_number Fit identifiers, or empty for a case check.
  #' @param check Check label.
  #' @param ok Logical outcome.
  #' @param detail Non-patient-level explanation or status.
  #' @return Nothing.
  #' @examples
  #' # record("real", "normal", .5, "pilot", 1, "rank", TRUE)
  record <- function(dataset, family, sigma, arm, pass_number, check, ok, detail = "") {
    rows[[length(rows) + 1L]] <<- data.frame(dataset = dataset, family = family, sigma = sigma, arm = arm, pass_number = pass_number, check = check, ok = isTRUE(ok), detail = detail, stringsAsFactors = FALSE)
  }
  for (j in seq_len(nrow(grid))) {
    dataset <- grid$dataset[j]; family <- grid$family[j]; sigma <- grid$sigma[j]
    #' Record a check for the current case.
    #' @param check,ok,detail Check name, result and explanation.
    #' @return Nothing.
    #' @examples
    #' # case_record("checkpoint_present", TRUE)
    case_record <- function(check, ok, detail = "")
      record(dataset, family, sigma, "", NA_integer_, check, ok, detail)
    key <- sprintf("%s_%s_%s", dataset, family, format(sigma, trim = TRUE, scientific = FALSE))
    path <- file.path(private, paste0(key, ".rds"))
    if (!file.exists(path)) { case_record("checkpoint_present", FALSE); next }
    object <- tryCatch(readRDS(path), error = function(e) e)
    if (inherits(object, "error")) { case_record("checkpoint_readable", FALSE, conditionMessage(object)); next }
    case_record("checkpoint_present", TRUE)
    tracking_result <- tryCatch(audit_track(object), error = function(e) e)
    if (inherits(tracking_result, "error")) case_record("tracking_error", FALSE, conditionMessage(tracking_result))
    else tracking[[length(tracking) + 1L]] <- tracking_result
    case_record("identity_and_complete", identical(object$schema, 1L) && isTRUE(object$complete) && identical(object$dataset, dataset) && identical(object$family, family) && identical(object$sigma, sigma) && is.null(object$error))
    sig <- object$signature; if (is.null(reference_signature)) reference_signature <- sig
    hashes <- sig$files; case_record("source_and_input_hashes", length(hashes) > 0L && !anyNA(hashes) && identical(hashes, tools::md5sum(names(hashes))) && identical(sig, reference_signature))
    versions <- tryCatch(vapply(names(sig$packages), function(p) as.character(utils::packageVersion(p)), character(1)), error = function(e) NULL)
    case_record("package_versions", !is.null(versions) && identical(versions, sig$packages))
    case_record("registered_settings", identical(sig$targets, c(.05, .25, .5)) && identical(sig$weights, c(40, 100, 1000)) && identical(sig$max_iter, 500L) &&
      identical(sig$seed, 1L) && identical(sig$families, c("normal", "laplace")) && identical(sig$sigma, c(.1, .25, .5, 1)))
    expected <- if (dataset == "s9") c("w40", "w100", "w1000") else c("pilot", "s0.05", "s0.25", "s0.5")
    case_record("expected_arms", setequal(names(object$arms), expected) && !anyDuplicated(names(object$arms)))
    pilot_converged <- dataset == "s9" || (length(object$arms$pilot$fits) > 0L && isTRUE(object$arms$pilot$fits[[1]]$history$converged))
    for (name in expected) {
      arm <- object$arms[[name]]; if (is.null(arm)) { record(dataset, family, sigma, name, NA_integer_, "arm_present", FALSE); next }
      if (!length(arm$fits)) {
        known <- dataset != "s9" && name != "pilot" && !pilot_converged && identical(arm$error, "pilot did not converge") && !isTRUE(arm$guard)
        record(dataset, family, sigma, name, 0L, "unavailable_target_recorded", known, if (known) "nonconverged pilot" else "missing fit without registered reason")
        next
      }
      nfits <- length(arm$fits); proper_count <- if (dataset == "s9" || name == "pilot" || isTRUE(arm$guard)) nfits == 1L else nfits <= 3L
      record(dataset, family, sigma, name, NA_integer_, "pass_count", proper_count)
      record(dataset, family, sigma, name, NA_integer_, "arm_error", is.null(arm$error))
      record(dataset, family, sigma, name, NA_integer_, "fit_checkpoint_counts", length(arm$checkpoints) == nfits && length(arm$weights) == nfits &&
               length(arm$warnings) == nfits && length(arm$elapsed) == nfits && length(arm$capped) == nfits)
      if (isTRUE(arm$guard))
        record(dataset, family, sigma, name, NA_integer_, "guard_reuses_pilot", dataset != "s9" && name != "pilot" && identical(arm$weights, 1) && identical(arm$fits[[1]], object$arms$pilot$fits[[1]]), "no-signal guard")
      if (dataset != "s9" && name != "pilot" && !isTRUE(arm$guard)) {
        history <- arm$calibration$history
        record(dataset, family, sigma, name, NA_integer_, "calibration_history", length(history) == nfits && nfits <= 3L && identical(vapply(history, `[[`, numeric(1), "weight"), arm$weights))
      }
      for (i in seq_along(arm$fits)) {
        #' Record a check for the current fit.
        #' @param check,ok,detail Check name, result and explanation.
        #' @return Nothing.
        #' @examples
        #' # fit_record("trace_shape", TRUE)
        fit_record <- function(check, ok, detail = "")
          record(dataset, family, sigma, name, i, check, ok, detail)
        tryCatch({
          fit <- arm$fits[[i]]; checkpoint <- arm$checkpoints[[i]]
          n <- fit$history$n_iter; K <- if (dataset == "real") 5L else 7L
          trace <- fit$history$s16; fit_record("sweeps_and_rank", length(n) == 1L && is.finite(n) && n >= 1L && n <= 500L && ncol(fit$EF) == K && ncol(fit$EL2) == K && length(fit$EBeta) == K)
          expected_trace <- is.data.frame(trace) && nrow(trace) == n * K && all(table(factor(trace$iteration, levels = seq_len(n))) == K) &&
            all(vapply(split(trace$program, trace$iteration), function(p) setequal(p, seq_len(K)), logical(1)))
          fit_record("trace_shape", expected_trace)
          fields <- c("Tau", "EL2", "EF2", "EBeta2", "CoxW")
          fit_record("finite_posterior_and_cox_weights", all(vapply(fields, function(f)
            length(fit[[f]]) > 0L && all(is.finite(fit[[f]])), logical(1))) && length(checkpoint$cox_weights) > 0L && all(is.finite(checkpoint$cox_weights)))
          prior <- s17_controls(family, sigma)$beta_prior_init
          prior_trace <- fit$history$prior_trace; prior_ok <- identical(fit$g_beta, prior) && identical(fit$posterior_state$g_beta, prior) && length(prior_trace) >= n && all(vapply(prior_trace[seq_len(n)], function(x)
              identical(x$prior_before, prior) && identical(x$prior_after, prior), logical(1)))
          fit_record("fixed_prior_every_sweep", prior_ok)
          fit_record("trace_canonical_beta", audit_equal(trace$beta_after, trace$canonical_scale * trace$beta_before))
          last <- trace[trace$iteration == n, , drop = FALSE]
          last <- last[order(last$program), , drop = FALSE]
          fit_record("final_trace_beta_matches_fit", audit_equal(last$beta_after, as.numeric(fit$EBeta)))
          fit_record("checkpoint_fit_matches", identical(checkpoint$fit, fit))
          dat <- inputs[[dataset]]; eta <- as.numeric(dat$Y %*% fit$EF %*% fit$EBeta)
          fresh_w <- calc_cox_taylor_yf(eta, dat$time, dat$status)$w
          fit_record("cox_weights_recomputed", audit_equal(as.numeric(fresh_w), as.numeric(fit$CoxW)) && audit_equal(as.numeric(fresh_w), as.numeric(checkpoint$cox_weights)))
          native <- is.finite(fit$history$delta_elbo_rel[n]) && fit$history$delta_elbo_rel[n] < cfg$cavi$tol
          fit_record("native_elbo_flag", identical(fit$history$native_elbo_converged, native))
          sweeps <- trace[!duplicated(trace$iteration), , drop = FALSE]
          sweeps <- sweeps[order(sweeps$iteration), , drop = FALSE]
          window <- if (n >= 50L) tail(sweeps, 50L) else NULL
          stable_window <- !is.null(window) && all(is.finite(window$eta_rel_change)) && is.finite(window$eta_start_max[1]) && is.finite(tail(window$eta_end_max, 1)) &&
            all(window$eta_rel_change < 1e-3) && abs(tail(window$eta_end_max, 1) - window$eta_start_max[1]) /
              max(window$eta_start_max[1], .Machine$double.eps) < .01
          posthoc <- s17_stability(fit, sigma)
          fit_record("posthoc_stability", identical(posthoc$assessable, n >= 50L) &&
            identical(posthoc$stable, isTRUE(fit$history$converged) && native && stable_window) && identical(isTRUE(fit$history$converged), native && stable_window) &&
            identical(posthoc$beta_pass, max(abs(fit$EBeta)) <= 4 * sigma))
          fit_record("weight_matches_fit", is.finite(arm$weights[i]) && identical(as.numeric(fit$w_S), as.numeric(arm$weights[i])))
        }, error = function(e) fit_record("unexpected_audit_error", FALSE, conditionMessage(e)))
      }
    }
  }
  checks <- do.call(rbind, rows)
  counts <- aggregate(rep(1L, nrow(checks)), checks[c("check", "ok")], sum)
  names(counts)[3] <- "count"
  tracked <- if (length(tracking)) do.call(rbind, tracking) else data.frame()
  utils::write.csv(tracked, file.path(output, "screen_loading_tracking.csv"), row.names = FALSE)
  utils::write.csv(checks, file.path(output, "screen_completeness_audit.csv"), row.names = FALSE)
  utils::write.csv(counts, file.path(output, "screen_completeness_audit_summary.csv"), row.names = FALSE)
  cat(sprintf("Audited %d expected cases; %d checks, %d failures.\n", nrow(grid), nrow(checks), sum(!checks$ok)))
  if (any(!checks$ok)) stop("Screen audit failed; see aggregate audit CSV")
  invisible(checks)
}
if (sys.nframe() == 0L) audit_screen()
