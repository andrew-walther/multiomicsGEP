# ============================================================
# Script: run_supervised_robustness.R
# Purpose: Evaluate fixed-start, prior and training-only feature robustness
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, yaml, ebnm, parallel; existing test-bed and CV helpers
# Run: Rscript results/benchmark_sim/run_supervised_robustness.R [cores]
# ============================================================
# Fixed before execution: 20 seeds (20261015+1:20), noise {.25,1}, calibrated
# make_testbed(n_train=400,beta3=3.3). K=7; E did not adopt protection or a
# different K (ac572f6); no pruning; canonical w_S=1, supervised start,
# shared point-Laplace beta prior, warm-up10. Defaults/stopping unchanged.
# Each dataset: default point-Laplace L/F; 20 custom starts from the signed
# SVD plus supervised program, each jittered in BOTH EL and EF by Gaussian
# SD=.1*its column SD; normal L/F default start. Jitter RNG=seed+100000+j.
# Choose highest final monitored history$elbo_full, not recovery; compare
# its recovery to default. Retain every start, convergence and final ELBO.
# Then run run_supervised_train_cv.R unchanged (original five folds/seed2026),
# copying original CSVs as *_prefix.csv first. Finally feature top-n per cohort
# {1000,2000,3000,5000}, same training preprocessing and folds, no external data.
# Each count also gets one full-training fit (seed42) for its survival-active count.
# Training CV preprocessing is done once before folds, as in the old runner;
# absolute CV may be optimistic. Do not select a method from external labels.
# Stop the prompt on a fit/loading error, retaining explicit failures and
# completed starts; do not tune. All socket workers write distinct seed files.

#' Construct signed SVD plus supervised initialization, with reproducible jitter
#' @param tb make_testbed() list with Ytr,ttr,str and seed.
#' @param start Integer 1:20; zero returns the unjittered constructed start.
#' @return List EL_init (n by 7) and EF_init (p by 7), before canonicalization.
#' @examples
#' # make_jitter_start(tb,1L)
make_jitter_start <- function(tb,start) {
  Y<-tb$Ytr;K<-7L;s<-svd(Y,nu=K,nv=K);d<-sqrt(pmax(s$d[1:K],0))
  L<-s$u%*%diag(d,K,K);F<-s$v%*%diag(d,K,K)
  u0<-calc_cox_taylor_yf(rep(0,nrow(Y)),tb$ttr,tb$str)$u
  # F_K = Yᵀu₀; L_K = R_{−K} F_K / ‖F_K‖², before canonicalization.
  direction<-as.vector(crossprod(Y,u0));F[,K]<-direction
  R<-Y-L[,-K,drop=FALSE]%*%t(F[,-K,drop=FALSE])
  L[,K]<-as.vector(R%*%direction)/sum(direction^2)
  if(start>0L) {
    set.seed(tb$seed+100000L+start)
    L<-L+sweep(matrix(rnorm(length(L)),nrow(L),K),2,.1*apply(L,2,sd),"*")
    F<-F+sweep(matrix(rnorm(length(F)),nrow(F),K),2,.1*apply(F,2,sd),"*")
  }
  list(EL_init=L,EF_init=F)
}

#' Select a start solely by its final monitored ELBO
#' @param rows Data frame with finite final_elbo and recovery fields.
#' @return Selected one-row data frame; first start wins an exact tie.
#' @examples
#' best_elbo_start(data.frame(start=1:2,final_elbo=c(-10,-5),recovered=c(TRUE,FALSE)))
best_elbo_start <- function(rows) {
  if(any(!is.finite(rows$final_elbo)))stop("Non-finite final ELBO among starts")
  rows[which.max(rows$final_elbo),,drop=FALSE]
}

#' Load original single-modality test-bed and scoring functions
#' @return Nothing; original definitions loaded in the global environment.
#' @examples
#' load_robustness_code()
load_robustness_code <- function() {
  suppressPackageStartupMessages(library(survival))
  for(e in parse("results/benchmark_sim/run_supervision_diagnostics.R"))
    if(is.call(e)&&identical(e[[1]],as.name("<-"))&&
       as.character(e[[2]])%in%c("load_code","make_testbed","score_fit"))eval(e,.GlobalEnv)
  load_code();if(!exists("fit_cox_on_yf"))stop("Fitter unavailable")
  assign("N_TRAIN",400L,.GlobalEnv);assign("BETA3",3.3,.GlobalEnv)
}

#' Fit all pre-specified starts/prior arms for one simulated dataset
#' @param seed,noise Fixed dataset seed and noise multiplier.
#' @return Status row; complete starts saved in the dataset's own CSV.
#' @examples
#' # robustness_dataset(20261016L,.25)
robustness_dataset <- function(seed,noise) {
  tb<-make_testbed(seed,noise);tb$seed<-seed;rows<-list()
  path<-file.path(OUT,sprintf("starts_seed%d_noise%s.csv",seed,noise))
  cfg<-yaml::read_yaml("config/globals.yml")
  tryCatch({
    for(j in 0:21) {
      arm<-if(j==0L)"default"else if(j==21L)"normal_prior"else "jitter"
      a<-list(K=7L,max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,prior_LF=if(j==21L)"normal"else"point_laplace",
        canonicalize=TRUE,w_S=1,init_supervised=TRUE,beta_prior_mode="shared",
        prior_beta="point_laplace",beta_warmup=10L,verbose=FALSE,trace_supervision=TRUE)
      if(j%in%1:20)a<-c(a,list(init_method="custom"),make_jitter_start(tb,j))
      set.seed(seed)
      fit_warnings<-character()
      fit<-withCallingHandlers(do.call(fit_cox_on_yf,c(list(tb$Ytr,tb$ttr,tb$str),a)),
        warning=function(w)fit_warnings<<-c(fit_warnings,conditionMessage(w)))
      sc<-score_fit(fit,tb)
      rows[[length(rows)+1L]]<-cbind(seed=seed,noise=noise,arm=arm,start=j,
        final_elbo=tail(fit$history$elbo_full,1),warnings=paste(unique(fit_warnings),collapse=" | "),sc)
      write.csv(do.call(rbind,rows),path,row.names=FALSE)
    }
    data.frame(seed=seed,noise=noise,error=NA_character_)
  },error=function(e)data.frame(seed=seed,noise=noise,error=conditionMessage(e)))
}

#' Execute robustness in the required order without loading external cohorts
#' @param cores Number of local socket workers for simulation only.
#' @return Nothing; CSVs and unchanged-runner outputs are saved under OUT.
#' @examples
#' # run_supervised_robustness(6L)
run_supervised_robustness <- function(cores=6L) {
  OUT<<-"results/benchmark_sim/outputs/supervised_robustness"
  dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
  grid<-expand.grid(seed=20261015L+1:20,noise=c(.25,1))
  cl<-parallel::makeCluster(cores)
  parallel::clusterExport(cl,c("load_robustness_code","make_jitter_start","robustness_dataset","OUT"))
  parallel::clusterEvalQ(cl,load_robustness_code())
  status<-tryCatch(parallel::parLapplyLB(cl,seq_len(nrow(grid)),function(i,g)
    robustness_dataset(g$seed[i],g$noise[i]),g=grid),finally=parallel::stopCluster(cl))
  status<-do.call(rbind,status);write.csv(status,file.path(OUT,"dataset_status.csv"),row.names=FALSE)
  if(any(!is.na(status$error)))stop("Robustness simulation failed; see dataset_status.csv and retained start rows. No tuning.")
  files<-file.path(OUT,sprintf("starts_seed%d_noise%s.csv",grid$seed,grid$noise))
  all<-do.call(rbind,lapply(files,read.csv));write.csv(all,file.path(OUT,"all_starts.csv"),row.names=FALSE)
  if(any(!is.na(all$warnings) & nzchar(all$warnings)))warning("Simulation warnings are retained in all_starts.csv")
  selected<-do.call(rbind,lapply(split(all,list(all$seed,all$noise),drop=TRUE),function(x) {
    best<-best_elbo_start(x[x$arm=="jitter",]);def<-x[x$arm=="default",]
    best$recovery_differs<-best$recovered!=def$recovered
    best$default_recovered<-def$recovered
    best
  }))
  write.csv(selected,file.path(OUT,"best_starts.csv"),row.names=FALSE)
  panel<-rbind(all[all$arm%in%c("default","normal_prior"),],selected[names(all)])
  panel$arm[panel$arm=="jitter"]<-"best_of_20"
  sm<-do.call(rbind,lapply(split(panel,list(panel$arm,panel$noise),drop=TRUE),function(x)
    data.frame(arm=x$arm[1],noise=x$noise[1],n=nrow(x),recovered=sum(x$recovered),
      false_active=mean(x$false_z),nonconverged=sum(!x$converged),C=mean(x$c_validation))))
  sm$best_recovery_differs<-vapply(seq_len(nrow(sm)),function(i)
    if(sm$arm[i]=="best_of_20")sum(selected$recovery_differs[selected$noise==sm$noise[i]])else NA_integer_,integer(1))
  write.csv(sm,file.path(OUT,"simulation_summary.csv"),row.names=FALSE)
  # Preserve old files exactly before the explicitly unchanged CV runner writes.
  old<-"results/benchmark_sim/outputs/supervised_train_cv"
  for(nm in c("folds","summary")) {
    from<-file.path(old,paste0(nm,".csv"));to<-file.path(OUT,paste0(nm,"_prefix.csv"))
    if(!file.exists(from))stop("Missing pre-fix CV output: ",from)
    if(!file.copy(from,to,overwrite=FALSE))stop("Could not preserve pre-fix output: ",to)
  }
  result<-system2(file.path(R.home("bin"),"Rscript"),"results/benchmark_sim/run_supervised_train_cv.R",
                  stdout=file.path(OUT,"train_cv.log"),stderr=file.path(OUT,"train_cv_errors.log"))
  if(result!=0L)stop("Unchanged training CV runner failed; see log")
  for(nm in c("folds","summary")) {
    if(!file.copy(file.path(old,paste0(nm,".csv")),
      file.path(OUT,paste0(nm,"_postfix.csv")),overwrite=FALSE))
      stop("Could not preserve post-fix CV output: ",nm)
  }
  cv<-read.csv(file.path(OUT,"folds_postfix.csv"))
  if(any(!is.finite(cv$C)))stop("Unchanged CV runner has failed/undefined folds; see saved outputs")
  before<-read.csv(file.path(OUT,"summary_prefix.csv"))
  after<-read.csv(file.path(OUT,"summary_postfix.csv"))
  comparison<-merge(before,after,by="config",suffixes=c("_prefix","_postfix"))
  comparison$change_in_mean_C<-comparison$mean_fold_C_postfix-comparison$mean_fold_C_prefix
  write.csv(comparison,file.path(OUT,"cv_before_after.csv"),row.names=FALSE)
  # Reuse its loading and fold prefix, ending before common fit configurations.
  env<-globalenv()
  for(e in parse("results/benchmark_sim/run_supervised_train_cv.R")) {
    if(is.call(e)&&identical(e[[1]],as.name("<-"))&&identical(e[[2]],as.name("common")))break
    eval(e,env)
  }
  feature_rows<-list();full_rows<-list()
  for(top in c(1000L,2000L,3000L,5000L)) {
    pp<-preprocess_merged_cohorts(env$train_raw,log_transform_flags=env$PLATFORM_LOG_TRANSFORM[env$TRAIN],
      top_n=top,rank_transform=FALSE,per_platform_standardize=TRUE,normalize_method="none",
      selection_per_cohort=TRUE,selection_method="combined_rank")
    Y<-pp$Y
    for(f in 1:5) {
      tr<-which(env$fold!=f);te<-which(env$fold==f);set.seed(42)
      fit<-fit_cox_on_yf(Y[tr,,drop=FALSE],env$time_train[tr],env$status_train[tr],K=7L,
        max_iter=env$cfg$cavi$max_iter,tol=env$cfg$cavi$tol,prior_LF="point_laplace",canonicalize=TRUE,
        w_S=1,init_supervised=TRUE,beta_prior_mode="shared",prior_beta="point_laplace",beta_warmup=10L,
        verbose=FALSE)
      risk<-predict_cox_on_yf(Y[te,,drop=FALSE],fit$EF,fit$EBeta,EF_norms=fit$EF_norms)$risk_scores
      sd<-sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
      feature_rows[[length(feature_rows)+1L]]<-data.frame(top_n=top,p=ncol(Y),fold=f,
        C=frozen_reverse_cindex(risk,env$time_train[te],env$status_train[te]),
        n_active=sum(abs(fit$EBeta)>=1.96*sd & sd>0),converged=fit$history$converged)
      write.csv(do.call(rbind,feature_rows),file.path(OUT,"feature_folds.csv"),row.names=FALSE)
    }
    set.seed(42)
    full<-fit_cox_on_yf(Y,env$time_train,env$status_train,K=7L,
      max_iter=env$cfg$cavi$max_iter,tol=env$cfg$cavi$tol,prior_LF="point_laplace",canonicalize=TRUE,
      w_S=1,init_supervised=TRUE,beta_prior_mode="shared",prior_beta="point_laplace",beta_warmup=10L,
      verbose=FALSE)
    full_sd<-sqrt(pmax(full$EBeta2-full$EBeta^2,0))
    full_active<-abs(full$EBeta)>=1.96*full_sd & full_sd>0
    full_rows[[length(full_rows)+1L]]<-data.frame(top_n=top,p=ncol(Y),
      n_active_full=sum(full_active),active_programs=paste(which(full_active),collapse=";"),
      full_converged=full$history$converged)
    write.csv(do.call(rbind,full_rows),file.path(OUT,"feature_full_fits.csv"),row.names=FALSE)
  }
  ft<-do.call(rbind,feature_rows)
  if(any(!is.finite(ft$C)))stop("Undefined feature-CV C; partial fold outputs retained")
  fs<-do.call(rbind,lapply(split(ft,ft$top_n),function(x)data.frame(top_n=x$top_n[1],p=x$p[1],
    mean_C=mean(x$C),SE=sd(x$C)/sqrt(5),mean_active=mean(x$n_active),nonconverged=sum(!x$converged))))
  fs<-merge(fs,do.call(rbind,full_rows),by=c("top_n","p"))
  write.csv(fs,file.path(OUT,"feature_summary.csv"),row.names=FALSE)
  print(sm,row.names=FALSE);print(fs,row.names=FALSE)
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE);run_supervised_robustness(if(length(a))as.integer(a[1])else 6L)
}
