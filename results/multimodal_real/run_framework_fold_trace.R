# ============================================================
# Script: run_framework_fold_trace.R
# Purpose: Trace the unchanged TCGA framework fold1 that loses every survival coefficient.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing R1 runner, ebnm, yaml, survival
# ============================================================
# Task2-3 diagnostic, preregistered before fitting. No remedy or score change.
# Repeat EXACT R1 cell84: seed42, event-stratified fold1, unchanged cached
# TCGA training matrices, K7, per-factor Normal, no supervised initialization,
# original pruning, cap500 and output tolerance. Only trace=TRUE differs.
# Test hypothesis: per-factor EB variances reach zero when |x|sqrt(A)<=1,
# rather than the shared-prior absorption of the selected configuration.
# Report when loss/revival occurs and compare numerical fit with saved cell84.
# No adoption, criterion/default/preprocessing change, or external data.
# Read existing full-TCGA framework trace as a descriptive reference only.
# Fit must reproduce saved EL/EF/beta moments and history exactly; mismatch
# stops interpretation. Prior trace unit tests establish observational identity.
# One local serial fit, ~1hour, budget8GB. Never edit sources while running.
# Run: Rscript results/multimodal_real/run_framework_fold_trace.R
source("results/benchmark_sim/run_round2_fixed_prior.R")

#' Enable observation only for the exact historical framework control list
#' @return Multimodal control list identical to R1 except trace=TRUE.
#' @examples
#' framework_trace_controls()
framework_trace_controls <- function() {
  a <- fixed_prior_controls("multimodal","framework")
  a$trace <- TRUE
  a
}

#' Summarize zero posterior moments and pseudo-observation thresholds by sweep
#' @param x Factor-level trace from collapse_trace_rows(), in execution order.
#' @return Sweep data frame, including the first persistently all-zero sweep.
#' @examples
#' # framework_trace_boundaries(collapse_trace_rows(fit$diagnostics$trace))
framework_trace_boundaries <- function(x) {
  key <- paste(x$stage_outer,x$inner,sep=":")
  keys <- unique(key)
  z <- do.call(rbind,lapply(keys,function(k) {
    a <- x[key==k,,drop=FALSE]
    data.frame(stage_outer=a$stage_outer[1],inner=a$inner[1],
      all_zero=all(a$beta==0 & a$beta_sd==0),
      all_pseudo_subthreshold=all(is.finite(a$x) & a$A>0 & abs(a$x)*sqrt(a$A)<=1),
      max_abs_pseudo_z=max(abs(a$x)*sqrt(a$A)),
      nonzero_coefficients=sum(a$beta!=0),nondegenerate_coefficients=sum(a$beta_sd>0))
  }))
  live <- which(!z$all_zero)
  start <- if(all(!z$all_zero))NA_integer_ else {
    last <- if(length(live))max(live) else 0L
    if(last<nrow(z))last+1L else NA_integer_
  }
  z$persistent_zero_start <- start
  z
}

if(sys.nframe()==0L) {
  out <- "results/multimodal_real/outputs/round2/framework_fold_trace"
  if(dir.exists(out))stop("Fresh diagnostic output required")
  dir.create(out,recursive=TRUE)
  input <- "data/round2_tcga_training.rds"
  expected <- readLines("results/benchmark_sim/outputs/round2/fixed_prior/training_source.txt")
  stopifnot(identical(expected,paste("training_md5",tools::md5sum(input))))
  tr <- readRDS(input);va <- fixed_prior_folds(tr$event)==1L
  old <- readRDS("results/benchmark_sim/outputs/round2/fixed_prior/cells/cell_084.rds")
  train <- lapply(tr$Y,function(y)y[!va,,drop=FALSE])
  val <- lapply(tr$Y,function(y)y[va,,drop=FALSE])
  stopifnot(identical(rownames(train$expression),old$training_spec$subject_ids))
  load_panel_code();source("code/concordance_ci.R")
  warnings <- character();started <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(tryCatch({
    set.seed(42)
    fit_multimodal_yfb(train,tr$time[!va],tr$event[!va],7L,framework_trace_controls())
  },error=function(e) {
    writeLines(conditionMessage(e),file.path(out,"error.txt"));stop(e)
  }),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  saveRDS(fit,file.path(out,"fit.rds"))
  writeLines(unique(warnings),file.path(out,"warnings.txt"))
  writeLines(expected,file.path(out,"training_source.txt"))
  fields <- c("EL","EL2","EF","EF2","EBeta","EBeta2","Tau","mu","prior_L","prior_F","kl","g_beta")
  check <- data.frame(field=c(fields,"history"),
    identical=c(vapply(fields,function(n)identical(fit[[n]],old[[n]]),logical(1)),
      identical(fit$diagnostics$history,old$diagnostics$history)))
  write.csv(check,file.path(out,"identity.csv"),row.names=FALSE)
  x <- collapse_trace_rows(fit$diagnostics$trace)
  write.csv(x,file.path(out,"fold1_trace.csv"),row.names=FALSE)
  write.csv(framework_trace_boundaries(x),file.path(out,"fold1_boundaries.csv"),row.names=FALSE)
  full <- readRDS("results/multimodal_real/outputs/round2/diagnosis1/multimodal_framework.rds")
  write.csv(framework_trace_boundaries(collapse_trace_rows(full$diagnostics$trace)),
    file.path(out,"full_training_boundaries.csv"),row.names=FALSE)
  risk <- predict_multimodal_yfb(fit,val)$risk_scores
  sd <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
  row <- data.frame(source="round2_training_cv",model="multimodal",arm="framework_trace_fold1",
    n=sum(!va),K=ncol(fit$EL),active=sum(sd>0 & abs(fit$EBeta)>=1.96*sd),
    training_C=NA_real_,cv_C=frozen_reverse_cindex(risk,tr$time[va],tr$event[va]),
    risk_sd=stats::sd(risk),converged=fit$diagnostics$converged,
    iterations=fit$diagnostics$iterations,runtime_seconds=proc.time()[["elapsed"]]-started)
  bench <- read.csv("results/multimodal_real/outputs/round2/diagnosis1/comparison_all.csv")
  bench <- bench[bench$source=="benchmark_snapshot",]
  rows <- list(row,bench);cols <- unique(unlist(lapply(rows,names)))
  comparison <- do.call(rbind,lapply(rows,function(a){a[setdiff(cols,names(a))]<-NA;a[cols]}))
  write.csv(comparison,file.path(out,"comparison_with_benchmark.csv"),row.names=FALSE)
  if(!all(check$identical))stop("Trace did not reproduce saved fold; interpretation blocked")
}
