# ============================================================
# Script: run_yfb_supervision_decision.R
# Purpose: Registered paired simulation decision test of normalized YF supervision.
# Created: 2026-10-30
# Dependencies: survival, ebnm, glmnet, parallel, existing simulation generator
# Run: Rscript results/benchmark_sim/run_yfb_supervision_decision.R [cores] [seeds] [attempt]
# ============================================================
# Socket workers only; each holds n×p matrices and compiled sweep scratch.
# For Longleaf allocate one CPU per worker and at least 2 GB per worker;
# use single-threaded BLAS to avoid oversubscription. No patient data loaded.
Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
source("code/fit_yfb_supervised.R")
source("code/simulate_multimodal_yfb.R")
source("code/fit_cox_on_yf.R")
source("code/predict_cox_on_yf.R")
source("code/concordance_ci.R")

#' Make paired expression cohorts at a fixed planted centered variance share
#'
#' Extend the calibrated A1 generator (n=400, β=3.3, noise=.25) to larger p.
#' Rescale the third activity to the target variance; add feature offsets to
#' both cohorts so Gaussian noise is retained without negative reconstruction
#' inputs. Offsets cancel from the centered survival predictor. Outcomes are
#' regenerated from the expression-only projection, avoiding a methylation
#' contribution in the original two-modality generator.
#' @param seed,p RNG seed and gene count.
#' @param share Desired raw centered variance share of the planted rank-one signal.
#' @return List of cohorts, scaled truth, and independently regenerated outcomes.
#' @examples
#' yfb_decision_data(1,50,.01)
yfb_decision_data <- function(seed,p,share) {
  d<-simulate_multimodal_yfb_data("low_variance_prognostic_factor",n_train=400,
    n_validation=120,p_expression=p,p_methylation=50,seed=seed,
    truncate=FALSE,noise_scale=.25,beta=c(0,0,3.3))
  F<-d$truth$F$expression
  S<-outer(d$training$L[,3],F[,3]);Sc<-scale(S,scale=FALSE)
  base<-d$training$Y$expression-S;Bc<-scale(base,scale=FALSE)
  aa<-(1-share)*sum(Sc^2);bb<- -2*share*sum(Sc*Bc);cc<- -share*sum(Bc^2)
  a<-(-bb+sqrt(bb^2-4*aa*cc))/(2*aa)
  tr<-base+a*S
  va<-d$validation$Y$expression+(a-1)*outer(d$validation$L[,3],F[,3])
  # Offsets define the simulated measurement scale and use no outcomes.
  offset<-pmax(-pmin(apply(tr,2,min),apply(va,2,min)),0)+.1
  tr<-sweep(tr,2,offset,"+");va<-sweep(va,2,offset,"+")
  sdtr<-apply(tr,2,sd);sdva<-apply(va,2,sd)
  truth<-F[,3]*sdtr
  # Keep the calibrated coefficient per unit raw projected expression. A
  # centered projection changes only the baseline hazard, not relative risk.
  outcomes<-function(X) {
    eta<-drop(scale(X,scale=FALSE)%*%F[,3])*3.3
    et<-(-log(runif(nrow(X)))/(.01*exp(eta)))^(1/1.5)
    cb<-rexp(nrow(X));cs<-multimodal_yfb_calibrate_censoring(et,cb,.3)
    list(time=pmin(et,cb*cs),event=as.numeric(et<=cb*cs),eta=eta)
  }
  ot<-outcomes(tr);ov<-outcomes(va)
  list(Y=yfb_scale(tr,rep(1,nrow(tr))),Yv=yfb_scale(va,rep(1,nrow(va))),
    time=ot$time,event=ot$event,tv=ov$time,ev=ov$event,
    truth=truth,oracle=ov$eta,actual_share=sum((a*Sc)^2)/sum(scale(tr,scale=FALSE)^2))
}

#' Recompute per-factor loading precision shares at a saved final state
#' @param fit Experimental fitter result with named-list EF and its moments.
#' @param Y Scaled, nonnegative patient-by-feature training matrix.
#' @param time,event Training survival outcomes.
#' @return Factor vector of median gene shares; NA for dead columns.
#' @examples
#' # yfb_precision_share_final(fit,Y,time,event)
yfb_precision_share_final <- function(fit,Y,time,event) {
  F<-fit$EF[[1]];live<-colSums(F^2)>0 & colSums(fit$EL2)>0
  yc<-scale(Y,scale=FALSE)
  pr<-yfb_projection(yc,F,fit$EF2[[1]],VF=if(!is.null(fit$VF))fit$VF[[1]] else NULL)
  tay<-yfb_cox(time,event,drop(pr$Z%*%fit$EBeta))
  rec<-outer(fit$Tau,colSums(fit$EL2))
  sur<-outer(colSums(yc^2*tay$w),fit$EBeta2*pr$s^2)*fit$gamma*fit$rho
  out<-rep(NA_real_,ncol(F))
  out[live]<-apply((sur/(rec+sur))[,live,drop=FALSE],2,median)
  out
}

#' Score loading recovery and posterior survival activity without reorientation
#' @param fit Fitter result, with matrix EF for the canonical comparator.
#' @param d Generated paired cohorts.
#' @param canonical Whether fit uses the historical matrix interface.
#' @return One metric row, with failures retained in the denominator.
#' @examples
#' # Called for every registered factor-model arm.
yfb_decision_score <- function(fit,d,canonical=FALSE) {
  F<-if(canonical)fit$EF else fit$EF[[1]]
  r<-suppressWarnings(cor(d$truth,F));r[!is.finite(r)]<-0;k<-which.max(abs(r))
  sd<-sqrt(if(!is.null(fit$VBeta))fit$VBeta else pmax(fit$EBeta2-fit$EBeta^2,0))
  active<-abs(fit$EBeta)>=1.96*sd & sd>0
  matched<-abs(r[k])>=.7
  risk<-if(canonical)predict_cox_on_yf(d$Yv,F,fit$EBeta,EF_norms=fit$EF_norms)$risk_scores else {
    Z<-sweep(scale(d$Yv,scale=FALSE)%*%F,2,fit$s,"*")
    drop(Z%*%fit$EBeta)
  }
  data.frame(cor_planted=abs(r[k]),recovered=matched && active[k] && fit$EBeta[k]*r[k]>0,
    false_active=sum(active[setdiff(seq_along(active),if(matched)k else integer())]),
    converged=if(canonical)isTRUE(fit$history$converged) else fit$converged,
    c_validation=frozen_reverse_cindex(risk,d$tv,d$ev),
    share=if(canonical)tail(fit$history$trace$shareA_median[,k],1) else tail(fit$history$share_median,1),
    sweeps=if(canonical)fit$history$n_iter else nrow(fit$history),
    decreases=if(canonical)sum(diff(fit$history$elbo_full)< -1e-8*abs(head(fit$history$elbo_full,-1))) else sum(fit$history$objective_decreased),
    max_beta=max(abs(fit$EBeta)),flags=if(canonical)"" else paste(fit$flags,collapse=";"),
    error=if(!canonical && !is.null(fit$failure))fit$failure$message else "",stringsAsFactors=FALSE)
}

#' Calibrate a fixed multiplier on identifiable, live loading columns
#' @param rec,surv Gene-by-factor nonnegative reconstruction and survival precision.
#' @param live Logical factor mask; dead columns carry no defined relative precision.
#' @param target Desired median survival share in (0,1).
#' @param rho Fixed damping applied to the survival loading message.
#' @return Positive gamma solving median(ργA_S/(A_R+ργA_S)) = target.
#' @examples
#' yfb_calibrate_gamma(matrix(2,3,2),matrix(.1,3,2),c(TRUE,FALSE),.2,.2)
yfb_calibrate_gamma <- function(rec,surv,live,target,rho) {
  ratio<-as.vector(surv[,live,drop=FALSE]/rec[,live,drop=FALSE])
  if(!length(ratio) || any(!is.finite(ratio)) || any(ratio<=0))
    stop("Cannot calibrate positive precision on live factors")
  # Solve directly: even-length sample medians need not commute exactly
  # with the nonlinear odds transformation.
  root<-uniroot(function(logg)median(plogis(log(ratio)+log(rho)+logg))-target,
                c(-700,700),tol=1e-10)$root
  exp(root)
}

#' Fit all paired arms for one generated dataset and checkpoint each result
#' @param row Dataset design row.
#' @param out Output directory; never overwrites a successful checkpoint.
#' @return Small status string; errors are saved as rows, not omitted.
#' @examples
#' # Invoked by socket workers.
yfb_decision_task <- function(row,out,fit_rho=.2) {
  cfg<-yaml::read_yaml("config/globals.yml")
  d<-yfb_decision_data(row$seed,row$p,row$variance)
  prefix<-sprintf("p%d_v%s_s%02d",row$p,row$variance,row$seed)
  warnings_seen<-character()
  record<-function(arm,expr,gamma=NA_real_) {
    path<-file.path(out,paste0(prefix,"_",arm,".rds"))
    if(file.exists(path)) return(readRDS(path)$fit)
    started<-proc.time()[3]
    ans<-tryCatch(withCallingHandlers(force(expr),warning=function(w) {
      warnings_seen<<-c(warnings_seen,conditionMessage(w))
      invokeRestart("muffleWarning")
    }),error=function(e)e)
    if(inherits(ans,"error")) {
      metrics<-data.frame(cor_planted=NA_real_,recovered=FALSE,false_active=cfg$benchmark$k_merged_yfb_desurv,
        converged=FALSE,c_validation=NA_real_,share=NA_real_,sweeps=NA_integer_,
        decreases=NA_integer_,max_beta=NA_real_,flags="error",error=conditionMessage(ans))
      fit<-NULL
    } else {fit<-ans$fit;metrics<-ans$metrics}
    metrics<-cbind(row,arm=arm,gamma=gamma,rho=fit_rho,actual_share=d$actual_share,
      elapsed=proc.time()[3]-started,metrics,warnings=paste(unique(warnings_seen),collapse=" | "))
    saveRDS(list(metrics=metrics,fit=fit),path)
    fit
  }
  u<-record("unsupervised",{
    f<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=cfg$benchmark$k_merged_yfb_desurv,gamma=0,scaled=TRUE)
    list(fit=f,metrics=yfb_decision_score(f,d))
  },0)
  # γ chosen once per dataset from the unsupervised fit's precision audit,
  # with the damping accounted for; γ remains fixed through every fitted arm.
  if(is.null(u)) stop("Unsupervised reference failed; cannot calibrate gamma")
  pr<-yfb_projection(scale(d$Y,scale=FALSE),u$EF[[1]],u$EF2[[1]],VF=if(!is.null(u$VF))u$VF[[1]] else NULL)
  tay<-yfb_cox(d$time,d$event,u$eta)
  rec<-outer(u$Tau,colSums(u$EL2))
  surv<-outer(colSums(scale(d$Y,scale=FALSE)^2*tay$w),u$EBeta2*pr$s^2)
  live<-colSums(u$EF[[1]]^2)>0 & colSums(u$EL2)>0
  for(target in c(.05,.2,.5)) {
    gamma<-yfb_calibrate_gamma(rec,surv,live,target,.2)
    record(paste0("share",target),{
      f<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=cfg$benchmark$k_merged_yfb_desurv,gamma=gamma,rho=fit_rho,scaled=TRUE)
      list(fit=f,metrics=yfb_decision_score(f,d))
    },gamma)
  }
  record("canonical_w1",{
    f<-fit_cox_on_yf(d$Y,d$time,d$event,K=cfg$benchmark$k_merged_yfb_desurv,max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,
      prior_LF="point_laplace",prior_beta="normal",canonicalize=TRUE,w_S=1,
      alpha=cfg$benchmark$alpha,verbose=FALSE,trace_supervision=TRUE,sign_correction=FALSE)
    list(fit=f,metrics=yfb_decision_score(f,d,TRUE))
  },1)
  record("lasso",{
    set.seed(row$seed+9000)
    f<-glmnet::cv.glmnet(d$Y,survival::Surv(d$time,d$event),family="cox",alpha=1,nfolds=cfg$cavi$n_cv_folds)
    risk<-drop(predict(f,d$Yv,s="lambda.min",type="link"))
    cf<-drop(as.matrix(coef(f,s="lambda.min")))
    rr<-suppressWarnings(cor(cf,d$truth));if(!is.finite(rr))rr<-0
    list(fit=NULL,metrics=data.frame(cor_planted=abs(rr),recovered=rr>=.7,
      false_active=NA_real_,converged=TRUE,c_validation=frozen_reverse_cindex(risk,d$tv,d$ev),
      share=NA_real_,sweeps=NA_integer_,decreases=NA_integer_,max_beta=max(abs(cf)),flags="",error=""))
  })
  prefix
}

# Main ----
if(sys.nframe()==0L) {
 args<-commandArgs(TRUE);cores<-if(length(args))as.integer(args[1]) else 4L
 seeds<-if(length(args)>1)as.integer(args[2]) else 10L
 attempt<-if(length(args)>2)args[3] else "attempt1"
 out<-file.path("results/benchmark_sim/outputs/yfb_supervision_2026-10-30",attempt)
 fit_rho<-if(attempt %in% c("attempt3","attempt4","attempt5","attempt6")) .02 else .2
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 design<-expand.grid(seed=seq_len(seeds),p=c(50L,2000L),variance=c(.0025,.01,.05))
 write.csv(design,file.path(out,"design.csv"),row.names=FALSE)
 cl<-parallel::makePSOCKcluster(cores,outfile=file.path(out,"workers.log"))
 on.exit(parallel::stopCluster(cl))
 parallel::clusterExport(cl,c("out","fit_rho","yfb_decision_data","yfb_decision_score","yfb_decision_task","yfb_calibrate_gamma"))
 parallel::clusterEvalQ(cl,{
   Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
   source("code/fit_yfb_supervised.R");source("code/simulate_multimodal_yfb.R")
   source("code/fit_cox_on_yf.R");source("code/predict_cox_on_yf.R");source("code/concordance_ci.R")
   multimodal_yfb_use_cpp();NULL
 })
 status<-parallel::parLapplyLB(cl,split(design,seq_len(nrow(design))),function(row)
   tryCatch(yfb_decision_task(row,out,fit_rho),error=function(e)paste("TASK ERROR",conditionMessage(e))),chunk.size=1L)
 writeLines(unlist(status),file.path(out,"task_status.txt"))
}
