# ============================================================
# Script: run_supervised_multimodal_comparison.R
# Purpose: Audit the constant two-step rank and compare fixed TCGA-to-ICGC models
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, yaml, ebnm, flashier; existing model/loading helpers
# Run: Rscript results/multimodal_real/run_supervised_multimodal_comparison.R
# ============================================================
# Fixed before execution: K_init=7, fit seed=42, TCGA training and ICGC primary
# (primary validation, n=50) plus all donors (secondary, n=67).
# Arms: M4 selected multimodal; 10/1 unsupervised-start pruned framework;
# original stacked two-step EBMF -> Cox; selected single-modality expression.
# E did not adopt protection or another starting rank (ac572f6).
# Joint arms: cached raw expression log2 and methylation beta, signed F, intercept.
# Expression-only selected model: same 3000 expression genes, TCGA and ICGC
# each standardized once per platform; primary ICGC rows subset afterward.
# This corrects the initial subset-wise scaling before expression scores exist. No feature retuning.
# Metrics: frozen-orientation C, K_final, survival-active count, convergence.
# EBMF/Cox convergence not returned by the original helper: report NA explicitly.
# The selected joint arm is attempted first; each fit is checkpointed. Stop on
# a fit error without changing model controls, and retain the error record.
# Resume uses this experiment's own per-arm checkpoints, as the existing real runner
# does. Remove those checkpoints before rerunning with different inputs or code.
# C++ sweeps use existing defaults; no fork workers or new parallelism.

#' Return the actual retained rank instead of the requested greedy ceiling
#' @param fit Existing two-step EBMF-to-Cox object.
#' @return Integer retained rank, checked against the loading columns.
#' @examples
#' retained_two_step_rank(list(n_factors=2L,loadings=matrix(0,5,2)))
retained_two_step_rank <- function(fit) {
  if(!identical(as.integer(fit$n_factors),as.integer(ncol(fit$loadings))))
    stop("Two-step rank metadata disagrees with loading columns")
  as.integer(fit$n_factors)
}

#' Standardize ICGC expression once before selecting its primary-PDAC subset
#' @param all All-donor ICGC expression matrix with unique donor row names.
#' @param primary_ids Primary-PDAC donor IDs, in the required scoring order.
#' @return List of all-donor and primary matrices on the same platform scale.
#' @examples
#' # standardize_icgc_expression(Y_icgc, primary_ids)
standardize_icgc_expression <- function(all,primary_ids) {
  idx<-match(primary_ids,rownames(all))
  if(anyNA(idx))stop("Primary donor missing from the all-donor ICGC platform")
  z<-scale(all)
  list(validation_all=z,validation_primary=z[idx,,drop=FALSE])
}

#' Run the cached-rank audit and fixed four-arm comparison
#' @return Invisibly the score table; writes audit, fits, metrics and any failure.
#' @examples
#' # run_supervised_multimodal_comparison()
run_supervised_multimodal_comparison <- function() {
  suppressPackageStartupMessages(library(survival))
  for(f in c("code/multimodal_yfb_helpers.R","code/preprocess_multimodal_yfb.R",
    "code/multimodal_yfb_updates.R","code/fit_multimodal_yfb.R",
    "code/predict_multimodal_yfb.R","code/simulate_multimodal_yfb.R")) source(f,local=.GlobalEnv)
  tryCatch(source("code/fit_cox_on_yf.R",local=.GlobalEnv),error=function(e) {
    if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)
  })
  source("code/predict_cox_on_yf.R",local=.GlobalEnv)
  source("code/concordance_ci.R",local=.GlobalEnv)
  source("code/run_multimodal_yfb_simulation.R",local=.GlobalEnv)
  cfg<-yaml::read_yaml("config/globals.yml")
  out<-"results/multimodal_real/outputs/supervised";dir.create(out,recursive=TRUE,showWarnings=FALSE)
  d<-readRDS("data/processed_multiomics_tcga_icgc.rds");tr<-d$training
  ref<-NULL;audit<-list()
  for(k in 6:8) {
    p<-sprintf("results/multimodal_real/outputs/fits/fits_K%02d_point_laplace-point_laplace_pruned.rds",k)
    x<-readRDS(p)$two_step
    if(is.null(ref))ref<-x
    audit[[length(audit)+1L]]<-data.frame(K_requested=k,K_retained=retained_two_step_rank(x),
      max_loading_difference=max(abs(x$loadings-ref$loadings)),
      max_coefficient_difference=max(abs(x$coefficient-ref$coefficient)),
      max_training_risk_difference=max(abs(x$training_risk-ref$training_risk)))
  }
  write.csv(do.call(rbind,audit),file.path(out,"two_step_rank_audit.csv"),row.names=FALSE)
  ctrl<-list(intercept=TRUE,prior_update="ebnm",prune=TRUE,survival_elbo="plugin",
    prior_F=list(expression="point_laplace",methylation="point_laplace"))
  icgc_expression<-standardize_icgc_expression(d$validation_all$Y$expression,
    rownames(d$validation_primary$Y$expression))
  methods<-c("selected_multimodal","pruned_10_1","two_step_ebmf_cox","selected_expression")
  rows<-list()
  for(method in methods) {
    message("Fitting ",method);set.seed(42);t0<-proc.time()[["elapsed"]]
    fit_path<-file.path(out,paste0(method,".rds"))
    resumed<-file.exists(fit_path)
    previous<-if(file.exists(file.path(out,"comparison.csv")))
      read.csv(file.path(out,"comparison.csv")) else NULL
    fit<-tryCatch({
      if(resumed)readRDS(fit_path)
      else if(method=="selected_multimodal") fit_multimodal_yfb(tr$Y,tr$time,tr$event,7L,
        c(ctrl,list(init_supervised=TRUE,beta_prior="shared_point_laplace",beta_warmup=10L)))
      else if(method=="pruned_10_1")fit_multimodal_yfb(tr$Y,tr$time,tr$event,7L,ctrl)
      else if(method=="two_step_ebmf_cox")fit_multimodal_yfb_ebmf_cox(tr$Y,tr$time,tr$event,7L)
      else fit_cox_on_yf(scale(tr$Y$expression),tr$time,tr$event,K=7L,
        max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,canonicalize=TRUE,w_S=1,
        prior_LF="point_laplace",init_supervised=TRUE,beta_prior_mode="shared",
        prior_beta="point_laplace",beta_warmup=10L,verbose=FALSE)
    },error=function(e)e)
    if(inherits(fit,"error")) {
      write.csv(data.frame(method=method,error=conditionMessage(fit)),file.path(out,"failure.csv"),row.names=FALSE)
      stop("Stopped at ",method,": ",conditionMessage(fit))
    }
    secs<-if(resumed) attr(fit,"runtime_seconds") else proc.time()[["elapsed"]]-t0
    if(is.null(secs)) {
      old_secs<-if(!is.null(previous))previous$runtime[previous$method==method] else numeric()
      secs<-if(length(old_secs))old_secs[1] else NA_real_
    }
    if(!resumed) {
      attr(fit,"runtime_seconds")<-secs
      saveRDS(fit,fit_path)
    }
    if(method %in% c("selected_multimodal","pruned_10_1")) {
      K<-ncol(fit$EL);active<-fit$diagnostics$survival_active;conv<-fit$diagnostics$converged
    } else if(method=="selected_expression") {
      K<-ncol(fit$EF);sd<-sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
      active<-abs(fit$EBeta)>=1.96*sd & sd>0;conv<-fit$history$converged
    } else {
      K<-retained_two_step_rank(fit);X<-do.call(cbind,tr$Y)
      S<-X%*%fit$loadings%*%solve(crossprod(fit$loadings)+diag(1e-8,K))
      cx<-coxph(Surv(tr$time,tr$event)~S);z<-summary(cx)$coefficients[,"z"]
      if(any(!is.finite(z)))stop("Non-finite two-step Cox z statistics")
      active<-abs(z)>=1.96;conv<-NA
    }
    for(cohort in c("validation_primary","validation_all")) {
      va<-d[[cohort]]
      risk<-if(method %in% c("selected_multimodal","pruned_10_1")) predict_multimodal_yfb(fit,va$Y)$risk_scores
        else if(method=="selected_expression") predict_cox_on_yf(icgc_expression[[cohort]],fit$EF,fit$EBeta,
          EF_norms=fit$EF_norms)$risk_scores else predict_multimodal_yfb_ebmf_cox(fit,va$Y)
      C<-frozen_reverse_cindex(risk,va$time,va$event)
      rows[[length(rows)+1L]]<-data.frame(method=method,cohort=cohort,n=length(risk),K_init=7L,
        K_final=K,n_active=sum(active),active_programs=paste(which(active),collapse=";"),
        converged=conv,C=C,runtime=secs)
    }
    write.csv(do.call(rbind,rows),file.path(out,"comparison.csv"),row.names=FALSE)
  }
  ans<-do.call(rbind,rows);print(ans,row.names=FALSE);invisible(ans)
}
if(sys.nframe()==0L)run_supervised_multimodal_comparison()
