# ============================================================
# Script: summarize_yfb_supervision.R
# Purpose: Retain full registered denominators and evaluate Gate S from checkpoints.
# Created: 2026-10-30
# Dependencies: base R
# Run: Rscript results/benchmark_sim/summarize_yfb_supervision.R attempt1 [attempt2 ...]
# ============================================================

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
   cells$converged>=criteria$minimum_converged & cells$errors==0
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
 design<-read.csv(file.path(out,"design.csv"))
 arms<-c("unsupervised","share0.05","share0.2","share0.5","canonical_w1","lasso")
 records<-list();j<-0L
 for(i in seq_len(nrow(design))) for(arm in arms) {
  d<-design[i,];path<-file.path(out,sprintf("p%d_v%s_s%02d_%s.rds",d$p,d$variance,d$seed,arm))
  j<-j+1L
  records[[j]]<-if(file.exists(path))readRDS(path)$metrics else
    cbind(d,arm=arm,gamma=NA_real_,rho=if(attempt %in% c("attempt3","attempt4")).02 else .2,actual_share=NA_real_,elapsed=NA_real_,
      cor_planted=NA_real_,recovered=FALSE,false_active=7,converged=FALSE,
      c_validation=NA_real_,share=NA_real_,sweeps=NA_real_,decreases=NA_real_,
      max_beta=NA_real_,flags="missing",error="No checkpoint; inspect task_status.txt",warnings="")
 }
 for(i in seq_along(records)) if(!"rho" %in% names(records[[i]]))
   records[[i]]$rho<-.2
 cols<-names(records[[1]])
 rows<-do.call(rbind,lapply(records,function(z)z[,cols]))
 write.csv(rows,file.path(out,"metrics.csv"),row.names=FALSE)
 groups<-split(rows,interaction(rows$p,rows$variance,rows$arm,drop=TRUE))
 mean_or_na<-function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE)
 summary<-do.call(rbind,lapply(groups,function(z)data.frame(p=z$p[1],variance=z$variance[1],
  arm=z$arm[1],n=nrow(z),recovered=sum(z$recovered),converged=sum(z$converged),
  mean_false=mean_or_na(z$false_active),mean_C=mean_or_na(z$c_validation),
  C_undefined=sum(!is.finite(z$c_validation)),errors=sum(nzchar(z$error)),
  missing=sum(z$flags=="missing"),warnings=sum(nzchar(z$warnings)),
  mean_share=mean_or_na(z$share),decreases=sum(z$decreases,na.rm=TRUE),
  max_beta=if(all(is.na(z$max_beta)))NA_real_ else max(z$max_beta,na.rm=TRUE))))
 write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
 evaluated<-yfb_gate_cells(summary)
 cells<-evaluated$cells;gate<-evaluated$gate
 write.csv(evaluated$criteria,file.path(out,"gate_criteria.csv"),row.names=FALSE)
 write.csv(cells,file.path(out,"gate_cells.csv"),row.names=FALSE)
 write.csv(gate,file.path(out,"gate.csv"),row.names=FALSE)
 counts<-data.frame(attempt=attempt,expected=nrow(rows),saved=nrow(rows)-sum(rows$flags=="missing"),
   missing=sum(rows$flags=="missing"),errors=sum(nzchar(rows$error)),
   nonconverged=sum(!rows$converged),undefined_C=sum(!is.finite(rows$c_validation)),
   warning_records=sum(nzchar(rows$warnings)),gate_passed=any(gate$passed))
 write.csv(counts,file.path(out,"counts.csv"),row.names=FALSE)
 print(counts);print(gate)
 invisible(summary)
}
if(sys.nframe()==0L) for(a in commandArgs(TRUE))summarize_yfb_attempt(a)
