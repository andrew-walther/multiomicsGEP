# ============================================================
# Script: diagnose_methylation_moments.R
# Purpose: Reproduce the unchanged 25k CpG failure and observe its first invalid F-update inputs.
# Created: 2026-10-15
# Dependencies: yaml, existing multimodal fitter and aggregate diagnostic helper
# Run after Task5-M2 prospective commit: Rscript results/multimodal_real/diagnose_methylation_moments.R
# ============================================================
# One serial training-only cell, original seed42, K7, original500/5 caps and
# frozen M1 design. No altered control, tolerance, prior, preprocessing or
# variance clipping. Entry tracing observes inputs; it never replaces them.
# Original terminal rows stay intact; detailed patient state is ignored/local.
source("code/fit_multimodal_yfb.R")
source("code/predict_multimodal_yfb.R")
source("code/train_test_split.R")
source("code/feature_sensitivity.R")
source("code/cpg_feature_sensitivity.R")
source("code/multimodal_feature_sensitivity.R")
source("code/multimodal_moment_audit.R")
out<-"results/multimodal_real/outputs/feature_sensitivity/moment_diagnosis"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"status.csv")))stop("Diagnostic already has a terminal result; preserve it.")
base<-"results/multimodal_real/outputs/feature_sensitivity/methylation_count"
h<-read.csv(file.path(base,"input_hashes.csv"))
stopifnot(identical(unname(tools::md5sum(h$path)),h$md5))
terminal<-read.csv(file.path(base,"cell_07.csv"))
z<-readRDS("data/round2_tcga_feature_universe.rds")
design<-readRDS("data/round2_methylation_count_design.rds")
cell<-design$cells[design$cells$cell==7L,]
stopifnot(nrow(cell)==1L,cell$p==25000L,cell$fold==0L)
features<-z$manifest$features[[cell$set]]
Y<-z$methylation[,match(features,colnames(z$methylation)),drop=FALSE]
audit<-new.env(parent=emptyenv());audit$calls<-0L;audit$failed<-FALSE

#' Observe a traced call without changing its arguments or its return value
#' @param call_env Execution environment of the F update at entry.
#' @param frames Active execution frames used only to identify fitter indices.
#' @return Invisible NULL; aggregate audit records and ignored local failure state.
#' @examples
#' # Invoked by trace() below, never by the fitting equations.
record_moment_entry <- function(call_env,frames) {
  audit$calls<-audit$calls+1L
  a<-mget(c("w","h_minus_k","EZ_k","VZ_k"),call_env,inherits=FALSE)
  rows<-do.call(multimodal_moment_guard_summary,a)
  ids<-which(vapply(frames,function(e)exists("outer",e,inherits=FALSE) &&
    exists("inner",e,inherits=FALSE) && exists("modality",e,inherits=FALSE),logical(1)))
  stopifnot(length(ids)>0L)
  fit_env<-frames[[tail(ids,1L)]]
  rows$outer<-get("outer",fit_env);rows$inner<-get("inner",fit_env)
  rows$factor<-get("k",fit_env);rows$modality<-get("modality",fit_env)
  rows$call<-audit$calls
  if(any(rows$invalid)) {
    audit$failed<-TRUE;audit$failure<-rows
    state<-mget(c("EL","EL2","EF","EF2","Tau","EBeta","EBeta2",
      "working","projection","current_eta","predictor","prior_F","prior_L",
      "outer","inner","k","modality","mu"),fit_env,inherits=FALSE)
    saveRDS(list(arguments=a,state=state),"data/round2_methylation_moment_failure.rds")
    write.csv(rows,file.path(out,"failed_inputs.csv"),row.names=FALSE)
    if(exists("last",audit,inherits=FALSE))
      write.csv(audit$last,file.path(out,"previous_valid_inputs.csv"),row.names=FALSE)
  } else audit$last<-rows
  invisible(NULL)
}
trace("multimodal_yfb_update_F_mk",tracer=quote({
  record_moment_entry(environment(),sys.frames())
}),print=FALSE)
warning_messages<-character();started<-proc.time()[3];set.seed(42L)
message<-tryCatch(withCallingHandlers({
  fit<-fit_multimodal_yfb(list(methylation=Y),z$time,z$event,design$K,design$control)
  saveRDS(fit,"data/round2_methylation_moment_unexpected_fit.rds")
  "Unexpected successful fit: original failure not reproduced."
},warning=function(w)warning_messages<<-c(warning_messages,conditionMessage(w))),
error=conditionMessage,finally=untrace("multimodal_yfb_update_F_mk"))
status<-data.frame(cell=7L,p=25000L,seed=42L,calls=audit$calls,
  captured_invalid_inputs=audit$failed,original_error=terminal$error,
  diagnostic_error=message,same_error=identical(message,terminal$error),
  runtime_seconds=unname(proc.time()[3]-started),warnings=paste(unique(warning_messages),collapse=" | "))
write.csv(status,file.path(out,"status.csv"),row.names=FALSE)
stopifnot(identical(unname(tools::md5sum(h$path)),h$md5))
write.csv(h,file.path(out,"unchanged_input_hashes.csv"),row.names=FALSE)
file.copy(file.path(base,"results_with_benchmark.csv"),file.path(out,"original_results_with_benchmark.csv"),overwrite=FALSE)
print(status)
if(!status$same_error || !audit$failed)stop("Unchanged failure was not reproduced with invalid-input capture.")
