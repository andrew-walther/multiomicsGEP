# ============================================================
# Script: run_round2_prior_consistency_trace.R
# Purpose: Audit loading-prior KL consistency in an unchanged completed Single simulation.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing fixed-panel modules and prior-consistency helpers
# ============================================================
# Task2-5 observation only: follow its committed prospective attempt-log entry.
# This observer changes no shared fitting source and is not a remedy.
# Record analytic support conflicts separately from point-prior cases with
# positive native returned posterior second moments; do not infer native
# support solely from tiny reconstructed-probability discrepancies.
# Repeat ONLY completed Task2-4 cell21: Single seed20261026/noise .25,
# same bed400/120,p50,K7,framework Normal beta,no supervised start,
# profile_cox_intercept TRUE, original300 cap/ELBO criterion, no pruning.
# A process-local observer returns each original loading update unchanged.
# Hypothesis: molecular q and refitted g are out of sync; the stored KL
# uses the old prior rather than the returned one. Test actual loading
# updates, not only the synthetic support counterexample. Do NOT repair,
# damp, reclassify convergence or select a model using the adjusted monitor.
# Validity requires exact equality of the ENTIRE returned saved fit and
# native vs analytic posterior slab probabilities within1e-7. This is a
# software consistency tolerance, not an adoption/credible-program threshold.
# Report all mismatches, support conflicts and the F-prior-only correction
# to the original monitor. It is not claimed to correct every ELBO component.
# No external data or patient-level real-data input. Fresh output only.
# Runtime a few minutes; one serial local fit,1GB; no HPC submission.
# Run: Rscript results/benchmark_sim/run_round2_prior_consistency_trace.R
source("results/benchmark_sim/run_round2_fixed_prior.R")
source("results/benchmark_sim/round2_prior_consistency_helpers.R")
if(sys.nframe()==0L) {
  out <- "results/benchmark_sim/outputs/round2/prior_consistency_seed26"
  if(dir.exists(out))stop("Fresh prior-consistency diagnostic output required")
  dir.create(out,recursive=TRUE)
  load_panel_code()
  src <- "results/benchmark_sim/outputs/round2/profiled/cells/cell_021.rds"
  original_fit <- readRDS(src)
  a <- fixed_prior_controls("single","profiled")
  tb <- make_testbed(20261026L,.25)
  rows <- list()
  #' Save one factor-level diagnostic row in process-local memory
  #' @param row Aggregate prior/posterior diagnostic row.
  #' @return NULL, invisibly.
  #' @examples
  #' # collect(data.frame(iteration=1,factor=1))
  collect <- function(row) {
    rows[[length(rows)+1L]] <<- row
    invisible(NULL)
  }
  original_update <- multimodal_yfb_update_F_mk
  multimodal_yfb_update_F_mk <- round2_prior_audit_wrapper(original_update,collect)
  set.seed(20261026L);started <- proc.time()[["elapsed"]]
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    do.call(fit_cox_on_yf,c(list(tb$Ytr,tb$ttr,tb$str),a)),
    warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w))
    }),
    error=function(e) {
      writeLines(conditionMessage(e),file.path(out,"error.txt"));stop(e)
    },finally={
      multimodal_yfb_update_F_mk<-original_update
      write.csv(data.frame(message=warnings),file.path(out,"warnings.csv"),row.names=FALSE)
      if(length(rows))write.csv(do.call(rbind,rows),
        file.path(out,"loading_prior_trace.csv"),row.names=FALSE)
    })
  saveRDS(fit,file.path(out,"fit.rds"))
  d <- do.call(rbind,rows)
  write.csv(d,file.path(out,"loading_prior_trace.csv"),row.names=FALSE)
  identity <- data.frame(field=names(original_fit),identical=vapply(names(original_fit),
    function(n)identical(fit[[n]],original_fit[[n]]),logical(1)))
  write.csv(identity,file.path(out,"identity.csv"),row.names=FALSE)
  delta <- tapply(d$prior_step_delta,d$iteration,sum)
  monitor <- data.frame(iteration=seq_along(fit$history$elbo_full),
    original_monitor=fit$history$elbo_full,F_prior_log_density_change=NA_real_)
  monitor$F_prior_log_density_change[as.integer(names(delta))] <- as.numeric(delta)
  monitor$F_prior_KL_adjusted_monitor <- monitor$original_monitor+
    monitor$F_prior_log_density_change
  write.csv(monitor,file.path(out,"monitor_audit.csv"),row.names=FALSE)
  late <- d$iteration>=max(d$iteration)-49L
  summary <- data.frame(seed=20261026L,noise=.25,iterations=fit$history$n_iter,
    converged=fit$history$converged,loading_blocks=nrow(d),
    posterior_support_conflicts=sum(d$expected_log_prior_new==-Inf),
    last50_support_conflicts=sum(d$expected_log_prior_new[late]==-Inf),
    native_point_prior_conflicts=sum(d$new_pi==0 & d$posterior_second_moment_sum>0),
    last50_native_point_prior_conflicts=sum(d$new_pi[late]==0 & d$posterior_second_moment_sum[late]>0),
    max_probability_discrepancy=max(d$slab_probability_discrepancy),
    exact_fit_identity=identical(fit,original_fit),
    nonfinite_F_prior_adjusted_monitors=sum(!is.finite(monitor$F_prior_KL_adjusted_monitor)),
    runtime_seconds=proc.time()[["elapsed"]]-started)
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
  writeLines(paste("source_rds_md5",tools::md5sum(src)),file.path(out,"source.txt"))
  if(!all(identity$identical) || !identical(names(fit),names(original_fit)) ||
     !identical(fit,original_fit))stop("Observer changed the saved fit; interpretation blocked")
  if(anyNA(d$prior_step_delta) || summary$max_probability_discrepancy>1e-7)
    stop("Analytic/native posterior mismatch; interpretation blocked")
  print(summary)
}
