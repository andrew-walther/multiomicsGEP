# ============================================================
# Script: summarize_yfb_supervision.R
# Purpose: Retain full registered denominators and evaluate Gate S from checkpoints.
# Created: 2026-10-30
# Dependencies: base R
# Run: Rscript results/benchmark_sim/summarize_yfb_supervision.R attempt1 [attempt2 ...]
# ============================================================

#' Enforce the registered conservative penalty without losing raw diagnostics
#' @param rows Per-fit records, including convergence, errors and flags.
#' @param K Number of fitted factors, used as the failed-fit false-active penalty.
#' @return Records with raw metrics preserved; failed factor fits are unrecovered.
#' @examples
#' # yfb_apply_failure_penalty(records,K)
yfb_apply_failure_penalty <- function(rows,K) {
 rows$raw_recovered<-rows$recovered;rows$raw_false_active<-rows$false_active
 rows$raw_converged<-rows$converged
 failed<-!rows$converged | nzchar(rows$error) |
   grepl("numerical_failure|nonfinite|all_zero|missing",rows$flags)
 rows$recovered[failed]<-FALSE;rows$converged[failed]<-FALSE
 rows$false_active[failed & rows$arm!="lasso"]<-K
 rows$false_active[rows$arm=="lasso"]<-NA_real_
 rows
}

#' Evaluate the pre-registered recovery, false-active and convergence gate
#' @param summary One row per feature-count, variance-cell and arm, with n,
#'   recovered, mean_false, converged, and errors.
#' @return List of cell verdicts, arm verdicts and the unchanged criteria.
#' @examples
#' # yfb_gate_cells(summary)
yfb_gate_cells <- function(summary) {
 criteria<-data.frame(genes=2000,seeds=10,minimum_gain=3,minimum_cells=2,
                      maximum_false_active=1,minimum_converged=9)
 cells<-summary[summary$p==criteria$genes & grepl("^share",summary$arm),]
 baseline<-summary[summary$p==criteria$genes & summary$arm=="unsupervised",]
 if(anyDuplicated(baseline$variance) || any(!cells$variance %in% baseline$variance))
   stop("Every gate cell requires one paired unsupervised baseline")
 cells$baseline_recovered<-baseline$recovered[match(cells$variance,baseline$variance)]
 cells$recovery_gain<-cells$recovered-cells$baseline_recovered
 cells$qualifies<-cells$n==criteria$seeds & cells$recovery_gain>=criteria$minimum_gain &
   cells$mean_false<=criteria$maximum_false_active &
   cells$converged>=criteria$minimum_converged
 gate<-do.call(rbind,lapply(split(cells,cells$arm),function(z)
   data.frame(arm=z$arm[1],qualifying_cells=sum(z$qualifies),
              passed=sum(z$qualifies)>=criteria$minimum_cells)))
 list(cells=cells,gate=gate,criteria=criteria)
}

#' Summarize one registered attempt, counting every missing or failed record
#' @param attempt Attempt directory name within the fixed output root.
#' @return Invisible summary table; writes metrics, gate cells and verdict CSVs.
#' @examples
#' # summarize_yfb_attempt("attempt2")
summarize_yfb_attempt <- function(attempt) {
 out<-file.path("results/benchmark_sim/outputs/yfb_supervision_2026-10-30",attempt)
 if(attempt %in% c("attempt4","attempt5","attempt6")) source("results/benchmark_sim/run_yfb_supervision_decision.R")
 design<-read.csv(file.path(out,"design.csv"))
 arms<-c("unsupervised","share0.05","share0.2","share0.5","canonical_w1","lasso")
 records<-list();j<-0L
 for(i in seq_len(nrow(design))) {
  if(attempt %in% c("attempt4","attempt5","attempt6")) dstate<-yfb_decision_data(design$seed[i],design$p[i],design$variance[i])
  for(arm in arms) {
  d<-design[i,];path<-file.path(out,sprintf("p%d_v%s_s%02d_%s.rds",d$p,d$variance,d$seed,arm))
  j<-j+1L
  saved<-if(file.exists(path))readRDS(path) else NULL
  records[[j]]<-if(!is.null(saved))saved$metrics else
    cbind(d,arm=arm,gamma=NA_real_,rho=if(attempt %in% c("attempt3","attempt4","attempt5","attempt6")).02 else .2,actual_share=NA_real_,elapsed=NA_real_,
      cor_planted=NA_real_,recovered=FALSE,false_active=7,converged=FALSE,
      c_validation=NA_real_,share=NA_real_,sweeps=NA_real_,decreases=NA_real_,
      max_beta=NA_real_,flags="missing",error="No checkpoint; inspect task_status.txt",warnings="")
  records[[j]]$raw_error<-records[[j]]$error
  records[[j]]$raw_flags<-records[[j]]$flags
  records[[j]]$lost_beta_variances<-0L
  records[[j]]$n_inactive<-NA_real_
  records[[j]]$all_zero_state<-NA
  records[[j]]$nonfinite_eta<-NA
  if(!is.null(saved$fit)) {
    ff<-saved$fit$EF
    if(is.list(ff) && !is.null(saved$fit$g_beta)) {
      gb<-saved$fit$g_beta
      nondegenerate<-any(gb$pi[gb$sd>0]>0)
      vv<-if(!is.null(saved$fit$VBeta))saved$fit$VBeta else
        saved$fit$EBeta2-saved$fit$EBeta^2
      lost<-if(nondegenerate)sum(!is.finite(vv)|vv<=0) else 0L
      records[[j]]$lost_beta_variances<-lost
      if(lost>0) {
        records[[j]]$error<-paste(records[[j]]$error,"Stored coefficient variance lost",sep="; ")
        records[[j]]$flags<-paste(records[[j]]$flags,"variance_loss",sep=";")
      }
    }
    if(is.list(ff))ff<-ff[[1]]
    records[[j]]$n_inactive<-sum(colSums(ff^2)==0 | colSums(saved$fit$EL^2)==0)
    records[[j]]$all_zero_state<-all(ff==0) || all(saved$fit$EL==0)
    if(!is.null(saved$fit$eta)) records[[j]]$nonfinite_eta<-any(!is.finite(saved$fit$eta))
  }
  records[[j]]$share_final_match<-NA_real_
  records[[j]]$share_final_live<-NA_real_
  if(attempt %in% c("attempt4","attempt5","attempt6") && !is.null(saved$fit) && is.list(saved$fit$EF)) {
    f<-saved$fit;sh<-yfb_precision_share_final(f,dstate$Y,dstate$time,dstate$event)
    rr<-suppressWarnings(cor(dstate$truth,f$EF[[1]]));rr[!is.finite(rr)]<-0
    records[[j]]$share_final_match<-sh[which.max(abs(rr))]
    records[[j]]$share_final_live<-if(all(is.na(sh)))NA_real_ else median(sh,na.rm=TRUE)
  }
 }
 }
 for(i in seq_along(records)) if(!"rho" %in% names(records[[i]]))
   records[[i]]$rho<-.2
 cols<-names(records[[1]])
 rows<-do.call(rbind,lapply(records,function(z)z[,cols]))
 rows<-yfb_apply_failure_penalty(rows,yaml::read_yaml("config/globals.yml")$benchmark$k_merged_yfb_desurv)
 write.csv(rows,file.path(out,"metrics.csv"),row.names=FALSE)
 groups<-split(rows,interaction(rows$p,rows$variance,rows$arm,drop=TRUE))
 mean_or_na<-function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE)
 summary<-do.call(rbind,lapply(groups,function(z)data.frame(p=z$p[1],variance=z$variance[1],
  arm=z$arm[1],n=nrow(z),recovered=sum(z$recovered),raw_recovered=sum(z$raw_recovered),
  raw_mean_false=mean_or_na(z$raw_false_active),raw_converged=sum(z$raw_converged),converged=sum(z$converged),
  mean_false=mean_or_na(z$false_active),mean_C=mean_or_na(z$c_validation),
  C_undefined=sum(!is.finite(z$c_validation)),errors=sum(nzchar(z$error)),
  missing=sum(z$flags=="missing"),warnings=sum(nzchar(z$warnings)),
  mean_share=mean_or_na(z$share),mean_share_matched=mean_or_na(z$share_final_match),
  mean_share_live=mean_or_na(z$share_final_live),mean_inactive=mean_or_na(z$n_inactive),decreases=sum(z$decreases,na.rm=TRUE),
  max_beta=if(all(is.na(z$max_beta)))NA_real_ else max(z$max_beta,na.rm=TRUE))))
 write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
 raw_summary<-summary
 raw_summary$recovered<-raw_summary$raw_recovered
 raw_summary$mean_false<-raw_summary$raw_mean_false
 raw_summary$converged<-raw_summary$raw_converged
 raw_evaluated<-yfb_gate_cells(raw_summary)
 write.csv(raw_evaluated$gate,file.path(out,"raw_gate_diagnostic.csv"),row.names=FALSE)
 evaluated<-yfb_gate_cells(summary)
 cells<-evaluated$cells;gate<-evaluated$gate
 write.csv(evaluated$criteria,file.path(out,"gate_criteria.csv"),row.names=FALSE)
 write.csv(cells,file.path(out,"gate_cells.csv"),row.names=FALSE)
 write.csv(gate,file.path(out,"gate.csv"),row.names=FALSE)
 counts<-data.frame(attempt=attempt,expected=nrow(rows),saved=nrow(rows)-sum(rows$flags=="missing"),
   missing=sum(grepl("missing",rows$flags)),raw_errors=sum(nzchar(rows$raw_error)),errors=sum(nzchar(rows$error)),
   nonconverged=sum(!rows$converged),undefined_C=sum(!is.finite(rows$c_validation)),
   warning_records=sum(nzchar(rows$warnings)),
   all_zero_returned_states=sum(rows$all_zero_state,na.rm=TRUE),
   inactive_factor_total=sum(rows$n_inactive,na.rm=TRUE),
   unavailable_factor_states=sum(is.na(rows$n_inactive)&rows$arm!="lasso"),
   nonfinite_returned_eta=sum(rows$nonfinite_eta,na.rm=TRUE),
   variance_loss_records=sum(rows$lost_beta_variances>0),gate_passed=any(gate$passed))
 write.csv(counts,file.path(out,"counts.csv"),row.names=FALSE)
 print(counts);print(gate)
 invisible(summary)
}
if(sys.nframe()==0L) for(a in commandArgs(TRUE))summarize_yfb_attempt(a)
