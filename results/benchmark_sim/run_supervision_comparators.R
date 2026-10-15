# ============================================================
# Script: run_supervision_comparators.R
# Purpose: Compare unsupervised and two-step methods on the calibrated planted bed
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, ebnm, flashier, glmnet, yaml, parallel
# Run: Rscript results/benchmark_sim/run_supervision_comparators.R [cores]
# ============================================================
# Fixed design: 20 seeds (20261015+1:20), noise {.25,1}, n_train=400,
# beta3=3.3, n_validation=120, p_expression=50 (make_testbed unchanged).
# Arms: EBMF (point-Laplace both sides, greedy/backfit, K<=7) alone;
# EBMF -> Cox; EBMF -> lasso-Cox (10-fold training CV, lambda.min);
# PCA -> Cox (K=7); selected canonical shared point-Laplace model, K=7,
# beta_warmup=10, no pruning. E did not adopt a change (ac572f6).
# DeSurv excluded: Prompt A found no installed/vendored implementation.
# Two-step projections YF use training F, matching the existing single-modality
# comparator contract. The candidate program is the largest |unpenalized Cox z|.
# Lasso has no coefficient z: use that same stage-2 Cox ranking, then test
# the lasso coefficient's sign. Lasso nonzero coefficients are not credible
# intervals: report their count separately; false survival-active is NA.
# Recovery: selected |cor(F_true3,F_k)|>=.7 with correct beta sign. EBMF
# alone reports best absolute loading correlation, without a survival sign test.
# C uses frozen fitted Cox orientation, never reoriented on validation.
# Folds for lasso fixed by dataset seed, stratified by event; no parallel BLAS
# or fork workers requested. Socket cluster, one dataset per task.

#' Load the original test-bed functions and required packages
#' @return Nothing; installs original definitions in the global environment.
#' @examples
#' load_comparator_code()
load_comparator_code <- function() {
  suppressPackageStartupMessages({library(survival);library(ebnm);library(flashier);library(glmnet)})
  wanted <- c("load_code","make_testbed","score_fit")
  for (e in parse("results/benchmark_sim/run_supervision_diagnostics.R"))
    if (is.call(e) && identical(e[[1]],as.name("<-")) && as.character(e[[2]]) %in% wanted)
      eval(e,.GlobalEnv)
  load_code()
  if (!exists("fit_cox_on_yf")) stop("Fitter definitions unavailable")
  assign("N_TRAIN",400L,.GlobalEnv);assign("BETA3",3.3,.GlobalEnv)
}

#' Evaluate loading recovery of a two-step method using its Cox-z-selected program
#' @param F p by K loading matrix; rows in original test-bed gene order.
#' @param beta,z K-vectors of coefficients and unpenalized Cox z statistics.
#' @param tb Original make_testbed() list, including F_true and beta_true.
#' @param unsupervised If TRUE, select by correlation and omit sign/active metrics.
#' @return One data-frame row: recovery, selected correlation/program and false-active count.
#' @examples
#' # comparator_metrics(F, beta, z, tb)
comparator_metrics <- function(F,beta,z,tb,unsupervised=FALSE) {
  r <- suppressWarnings(cor(tb$F_true[,3],F));r[!is.finite(r)]<-0
  if (unsupervised) {
    k<-which.max(abs(r))
    return(data.frame(program=k,cor_planted=abs(r[k]),recovered=abs(r[k])>=.7,
                      false_active=NA_integer_))
  }
  if(any(!is.finite(c(beta,z)))) stop("Non-finite Cox coefficients or z statistics")
  k<-which.max(abs(z));matched<-which(abs(r)>=.7)
  data.frame(program=k,cor_planted=abs(r[k]),
    recovered=abs(r[k])>=.7 && sign(beta[k]*r[k])==sign(tb$beta_true[3]),
    false_active=sum(abs(z[setdiff(seq_along(z),matched)])>=1.96))
}

#' Fit every fixed comparator for one dataset, sharing its unsupervised fit
#' @param seed,noise Fixed seed and noise multiplier.
#' @return Data frame of all five methods, or an explicit dataset error row.
#' @examples
#' # run_comparator_dataset(20261016L,.25)
run_comparator_dataset <- function(seed,noise) {
  tb<-make_testbed(seed,noise);t0<-proc.time()[["elapsed"]]
  tryCatch({
    set.seed(seed)
    ff<-flashier::flash(tb$Ytr,var_type=2,greedy_Kmax=7L,backfit=TRUE,verbose=0,
                        ebnm_fn=ebnm::ebnm_point_laplace)
    F<-as.matrix(flashier::ldf(ff,type="2")$F)
    if(ncol(F)<1L) stop("EBMF retained no programs")
    S<-tb$Ytr%*%F;Sv<-tb$Yva%*%F
    cx<-survival::coxph(survival::Surv(tb$ttr,tb$str)~S)
    b<-as.numeric(coef(cx));z<-summary(cx)$coefficients[,"z"]
    if(any(!is.finite(c(b,z)))) stop("EBMF Cox fit is singular/non-finite")
    uns<-comparator_metrics(F,b,z,tb,TRUE)
    rows<-list(cbind(method="ebmf_alone",uns,C=NA_real_,K_final=ncol(F),
                     n_active=NA_integer_,converged=NA))
    sc<-comparator_metrics(F,b,z,tb)
    rows[[2]]<-cbind(method="ebmf_cox",sc,
      C=frozen_reverse_cindex(as.vector(Sv%*%b),tb$tva,tb$sva),K_final=ncol(F),
      n_active=sum(abs(z)>=1.96),converged=NA)
    set.seed(seed);fold<-integer(nrow(S))
    for(st in unique(tb$str)) {
      idx<-which(tb$str==st);fold[idx]<-sample(rep(1:10,length.out=length(idx)))
    }
    # glmnet requires at least two predictor columns; a one-program fit cannot
    # satisfy this arm as specified and is reported rather than padded silently.
    cv<-glmnet::cv.glmnet(S,survival::Surv(tb$ttr,tb$str),family="cox",alpha=1,foldid=fold)
    bl<-as.numeric(coef(cv,s="lambda.min"))
    sc<-comparator_metrics(F,bl,z,tb);sc$false_active<-NA_integer_
    rows[[3]]<-cbind(method="ebmf_lasso_cox",sc,
      C=frozen_reverse_cindex(as.vector(Sv%*%bl),tb$tva,tb$sva),K_final=ncol(F),
      n_active=sum(bl!=0),converged=NA)
    pc<-svd(tb$Ytr,nu=0,nv=7L)$v;Sp<-tb$Ytr%*%pc
    cx<-survival::coxph(survival::Surv(tb$ttr,tb$str)~Sp)
    bp<-as.numeric(coef(cx));zp<-summary(cx)$coefficients[,"z"]
    sc<-comparator_metrics(pc,bp,zp,tb)
    rows[[4]]<-cbind(method="pca_cox",sc,
      C=frozen_reverse_cindex(as.vector(tb$Yva%*%pc%*%bp),tb$tva,tb$sva),K_final=7L,
      n_active=sum(abs(zp)>=1.96),converged=NA)
    set.seed(seed)
    fit<-fit_cox_on_yf(tb$Ytr,tb$ttr,tb$str,K=7L,max_iter=cfg$cavi$max_iter,
      tol=cfg$cavi$tol,prior_LF="point_laplace",canonicalize=TRUE,w_S=1,
      init_supervised=TRUE,beta_prior_mode="shared",prior_beta="point_laplace",
      beta_warmup=10L,verbose=FALSE,trace_supervision=TRUE)
    sf<-score_fit(fit,tb)
    rows[[5]]<-data.frame(method="supervised",program=sf$match_k,cor_planted=sf$cor_planted,
      recovered=sf$recovered,false_active=sf$false_z,C=sf$c_validation,K_final=ncol(fit$EF),
      n_active=sf$n_active_z,converged=sf$converged)
    cbind(seed=seed,noise=noise,error=NA_character_,runtime=proc.time()[["elapsed"]]-t0,
          do.call(rbind,rows))
  },error=function(e) data.frame(seed=seed,noise=noise,error=conditionMessage(e),
                                 runtime=proc.time()[["elapsed"]]-t0))
}

#' Run the fixed 40-dataset socket-cluster experiment and report all failures
#' @param cores Number of local socket workers.
#' @return Invisibly, raw result rows; CSV results and summaries are written.
#' @examples
#' # run_comparator_panel(6L)
run_comparator_panel <- function(cores=6L) {
  out<-"results/benchmark_sim/outputs/supervision_comparators"
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  cfg<<-yaml::read_yaml("config/globals.yml")
  grid<-expand.grid(seed=20261015L+1:20,noise=c(.25,1))
  cl<-parallel::makeCluster(cores);on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl,c("load_comparator_code","comparator_metrics",
                               "run_comparator_dataset","cfg"))
  parallel::clusterEvalQ(cl,load_comparator_code())
  rows<-parallel::parLapplyLB(cl,seq_len(nrow(grid)),function(i,g) {
    run_comparator_dataset(g$seed[i],g$noise[i])
  },g=grid)
  cols<-unique(unlist(lapply(rows,names)))
  res<-do.call(rbind,lapply(rows,function(x){x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(res,file.path(out,"results.csv"),row.names=FALSE)
  if(any(!is.na(res$error))) stop("Comparator dataset failure(s); inspect results.csv. No tuning or retry.")
  sm<-do.call(rbind,lapply(split(res,list(res$method,res$noise),drop=TRUE),function(x)
    data.frame(method=x$method[1],noise=x$noise[1],n=nrow(x),recovered=sum(x$recovered),
      mean_cor=mean(x$cor_planted),mean_C=if(all(is.na(x$C))) NA_real_ else mean(x$C),
      false_active=if(all(is.na(x$false_active))) NA_real_ else mean(x$false_active),
      median_K=median(x$K_final),nonconverged=if(all(is.na(x$converged))) NA_integer_ else sum(!x$converged))))
  write.csv(sm,file.path(out,"summary.csv"),row.names=FALSE);print(sm,row.names=FALSE)
  invisible(res)
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE);run_comparator_panel(if(length(a))as.integer(a[1])else 6L)
}
