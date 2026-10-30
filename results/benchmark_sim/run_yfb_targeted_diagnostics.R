# ============================================================
# Script: run_yfb_targeted_diagnostics.R
# Purpose: Matched fixed/adaptive noise and guarded-update screen on saved synthetic cases.
# Created: 2026-10-30
# Dependencies: existing YFB runner; parallel PSOCK; no patient inputs
# ============================================================
# One case per socket worker; single BLAS thread and >=2 GB RAM per worker.
# Longleaf: allocate one CPU per worker, plus parent memory; no fork workers.
Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
source("results/benchmark_sim/run_yfb_supervision_decision.R")

#' Run four matched diagnostic variants for one registered synthetic case
#' @param variance Raw planted variance share; seed 1 and p=2000 are fixed.
#' @param out New checkpoint directory; existing results are never overwritten.
#' @return Status string; all arm errors are checkpointed and surfaced.
#' @examples
#' # Invoked by PSOCK workers after registration.
#' @family yfb_supervised
#' @seealso fit_yfb_supervised
 yfb_targeted_case <- function(variance,out) {
  cfg<-yaml::read_yaml("config/globals.yml")
  old<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30/attempt7"
  ref<-readRDS(file.path(old,sprintf("p2000_v%s_s01_unsupervised.rds",variance)))$fit
  saved<-readRDS(file.path(old,sprintf("p2000_v%s_s01_share0.05.rds",variance)))
  d<-yfb_decision_data(1,2000,variance,"shared")
  for(arm in c("adaptive","fixed","guard_adaptive","guard_fixed")) {
    path<-file.path(out,sprintf("v%s_%s.rds",variance,arm))
    if(file.exists(path))next
    warnings_seen<-character()
    ans<-tryCatch(withCallingHandlers({
      f<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,
        K=cfg$benchmark$k_merged_yfb_desurv,gamma=saved$metrics$gamma,
        rho=saved$metrics$rho,scaled=TRUE,tau_init=ref$Tau,
        tau_fixed=arm %in% c("fixed","guard_fixed"),guard=grepl("guard",arm))
      m<-yfb_decision_score(f,d)
      residual<-d$Y-f$EL%*%t(f$EF[[1]])
      extra<-f$VL%*%t(f$EF[[1]]^2)+f$EL^2%*%t(f$VF[[1]])+f$VL%*%t(f$VF[[1]])
      m$reconstruction_ratio<-sum(residual^2)/sum(d$Y^2)
      m$expected_reconstruction_ratio<-sum(residual^2+extra)/sum(d$Y^2)
      m$min_tau<-min(f$Tau);m$max_F<-max(abs(f$EF[[1]]))
      m$rejected_trials<-f$rejected_trials
      m$fixed_tau_difference<-if(f$tau_fixed)max(abs(f$Tau-ref$Tau)) else NA_real_
      m$success<-f$converged && !length(f$flags) && m$recovered &&
        m$false_active<=1 && m$expected_reconstruction_ratio<=1 && m$reconstruction_ratio<=1
      list(fit=f,metrics=m)
    },warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")}),
      error=function(e)list(fit=NULL,error=conditionMessage(e)))
    ans$variance<-variance;ans$arm<-arm;ans$gamma<-saved$metrics$gamma
    ans$warnings<-unique(warnings_seen);saveRDS(ans,path)
  }
  sprintf("variance=%s complete",variance)
}

if(sys.nframe()==0L) {
  out<-"results/benchmark_sim/outputs/yfb_targeted_2026-10-30"
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  design<-expand.grid(variance=c(.0025,.01,.05),arm=c("adaptive","fixed","guard_adaptive","guard_fixed"))
  write.csv(design,file.path(out,"design.csv"),row.names=FALSE)
  cl<-parallel::makePSOCKcluster(3,outfile=file.path(out,"workers.log"))
  parallel::clusterExport(cl,c("out","yfb_targeted_case"))
  parallel::clusterEvalQ(cl,{source("results/benchmark_sim/run_yfb_supervision_decision.R");NULL})
  status<-tryCatch(parallel::parLapplyLB(cl,c(.0025,.01,.05),yfb_targeted_case,out=out),
    finally=parallel::stopCluster(cl))
  writeLines(unlist(status),file.path(out,"task_status.txt"))
  files<-list.files(out,pattern="rds$",full.names=TRUE)
  if(length(files)!=nrow(design))stop("Incomplete targeted checkpoint inventory")
  rows<-lapply(files,function(path){
    a<-readRDS(path)
    if(!is.null(a$error))return(data.frame(variance=a$variance,arm=a$arm,error=a$error))
    cbind(variance=a$variance,arm=a$arm,a$metrics,warnings=paste(a$warnings,collapse=" | "))
  })
  errors<-vapply(rows,function(x)"error" %in% names(x) && nzchar(x$error[1]),TRUE)
  if(any(errors)) {
    write.csv(do.call(rbind,lapply(rows[errors],function(x)x[,c("variance","arm","error")])),file.path(out,"errors.csv"),row.names=FALSE)
    stop("Arm errors saved; inspect checkpoints before summarizing")
  }
  write.csv(do.call(rbind,rows),file.path(out,"summary.csv"),row.names=FALSE)
}
