# ============================================================
# Script: run_supervision_comparators.R
# Purpose: Compare unsupervised and two-step methods on the calibrated planted bed
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, ebnm, flashier, glmnet, yaml, parallel
# Run: Rscript results/benchmark_sim/run_supervision_comparators.R [cores]
# Round2: per-arm errors, zero/one-factor handling, full20-seed denominators.
# Outputs are fresh round2/ checkpoints; original results are never overwritten.
# ============================================================
# Fixed design: 20 seeds (20261015+1:20), noise {.25,1}, n_train=400,
# beta3=3.3, n_validation=120, p_expression=50 (make_testbed unchanged).
# Arms: EBMF (point-Laplace both sides, greedy/backfit, K<=7) alone;
# EBMF -> Cox; EBMF -> lasso-Cox (10-fold training CV, lambda.min);
# PCA -> Cox (K=7); selected canonical shared point-Laplace model, K=7,
# beta_warmup=10, no pruning. E did not adopt a change (ac572f6).
# Andrew explicitly approved omitting DeSurv from this signed-Gaussian panel.
# Rank preprocessing would change the loading interpretation; reported top270
# clamps to all50 features. Task3 real-data head-to-head is deferred for review.
# Final five-arm design is committed before scientific execution.
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

# Small-rank and per-arm handling ----

#' Construct an explicit unavailable comparator row without dropping a dataset
#' @param method Comparator name.
#' @param K Actual retained rank, when available.
#' @param applicable Whether the method applies to the retained rank.
#' @param reason,error Structural inapplicability or exact fitting error.
#' @return One row with recovery FALSE and undefined metrics retained as NA.
#' @examples
#' comparator_empty_row("ebmf_lasso_cox",1L,FALSE,"one factor")
comparator_empty_row <- function(method,K=NA_integer_,applicable=TRUE,
  reason=NA_character_,error=NA_character_) {
  data.frame(method=method,program=NA_integer_,cor_planted=NA_real_,
    recovered=FALSE,false_active=NA_real_,C=NA_real_,K_final=K,
    n_active=NA_integer_,converged=NA,iterations=NA_integer_,
    applicable=applicable,reason=reason,error=error,warnings="")
}

#' Capture one comparator's error and warnings while preserving all other arms
#' @param method Comparator name.
#' @param expr Expression producing one scored data-frame row.
#' @param K Actual rank, when known before the arm.
#' @return Standard scored row, or a failure row with the exact message.
#' @examples
#' comparator_try_arm("bad",stop("planned"))
comparator_try_arm <- function(method,expr,K=NA_integer_) {
  warnings <- character();base <- comparator_empty_row(method,K)
  x <- withCallingHandlers(tryCatch({
    ans <- force(expr)
    stopifnot(is.data.frame(ans),nrow(ans)==1L)
    ans$method <- method
    for(n in setdiff(names(base),names(ans)))ans[[n]] <- base[[n]]
    ans[names(base)]
  },error=function(e)comparator_empty_row(method,K,error=conditionMessage(e))),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  x$warnings <- paste(unique(warnings),collapse=" | ");x
}

#' Evaluate the three EBMF arms with fixed zero/one-factor rules
#' @param F Actual p-by-K EBMF loadings, including p-by-zero.
#' @param tb Unchanged simulation bed with train/validation data and truth.
#' @param seed Fixed dataset seed for the event-stratified lasso folds.
#' @return Three rows; no arm or dataset is removed from the denominator.
#' @examples
#' # comparator_ebmf_rows(F,tb,20261016L)
comparator_ebmf_rows <- function(F,tb,seed) {
  K <- ncol(F);methods <- c("ebmf_alone","ebmf_cox","ebmf_lasso_cox")
  if(K==0L) {
    rows <- lapply(methods,function(m) {
      x <- comparator_empty_row(m,0L,applicable=m!="ebmf_lasso_cox",
        reason="zero factors: constant risk, C undefined")
      x$cor_planted <- 0;x$n_active <- 0L
      if(m=="ebmf_cox")x$false_active <- 0
      x
    })
    return(do.call(rbind,rows))
  }
  S <- tb$Ytr%*%F;Sv <- tb$Yva%*%F
  rows <- list(comparator_try_arm("ebmf_alone",{
    cbind(comparator_metrics(F,NULL,NULL,tb,TRUE),C=NA_real_,K_final=K,
      n_active=NA_integer_,converged=NA)
  },K))
  cox_warnings <- character()
  cx <- withCallingHandlers(tryCatch(
    survival::coxph(survival::Surv(tb$ttr,tb$str)~S),error=function(e)e),
    warning=function(w) {cox_warnings <<- c(cox_warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  rows[[2L]] <- comparator_try_arm("ebmf_cox",{
    if(inherits(cx,"error"))stop(conditionMessage(cx))
    b <- as.numeric(coef(cx));z <- as.numeric(summary(cx)$coefficients[,"z"])
    cbind(comparator_metrics(F,b,z,tb),
      C=frozen_reverse_cindex(as.vector(Sv%*%b),tb$tva,tb$sva),K_final=K,
      n_active=sum(abs(z)>=1.96),converged=NA,iterations=cx$iter)
  },K)
  if(K==1L) {
    rows[[3L]] <- comparator_empty_row("ebmf_lasso_cox",K,applicable=FALSE,
      reason="one factor: ordinary Cox only; lasso not applicable")
  } else rows[[3L]] <- comparator_try_arm("ebmf_lasso_cox",{
    if(inherits(cx,"error"))stop("Cox-z ranking unavailable: ",conditionMessage(cx))
    z <- as.numeric(summary(cx)$coefficients[,"z"])
    set.seed(seed);fold <- integer(nrow(S))
    for(st in unique(tb$str)) {
      idx <- which(tb$str==st);fold[idx] <- sample(rep(1:10,length.out=length(idx)))
    }
    cv <- glmnet::cv.glmnet(S,survival::Surv(tb$ttr,tb$str),family="cox",alpha=1,foldid=fold)
    b <- as.numeric(coef(cv,s="lambda.min"));sc <- comparator_metrics(F,b,z,tb)
    sc$false_active <- NA_real_
    cbind(sc,C=frozen_reverse_cindex(as.vector(Sv%*%b),tb$tva,tb$sva),
      K_final=K,n_active=sum(b!=0),converged=NA)
  },K)
  for(j in 2:3)rows[[j]]$warnings <- paste(unique(c(cox_warnings,
    rows[[j]]$warnings[nzchar(rows[[j]]$warnings)])),collapse=" | ")
  do.call(rbind,rows)
}

#' Fit all fixed methods with errors isolated to their dependent arms
#' @param seed,noise Fixed simulation seed and noise multiplier.
#' @return Five rows, including all unavailable/error cases.
#' @examples
#' # run_comparator_dataset(20261016L,.25)
run_comparator_dataset <- function(seed,noise) {
  tb <- make_testbed(seed,noise);t0 <- proc.time()[["elapsed"]]
  set.seed(seed)
  ebmf_warnings <- character()
  ff <- withCallingHandlers(tryCatch({
    f <- flashier::flash(tb$Ytr,var_type=2,greedy_Kmax=7L,backfit=TRUE,
      verbose=0,ebnm_fn=ebnm::ebnm_point_laplace)
    # ldf() throws for a valid zero-factor fit: use its actual retained count.
    if(f$n_factors==0L)matrix(numeric(),ncol(tb$Ytr),0L) else
      as.matrix(flashier::ldf(f,type="2")$F)
  },error=function(e)e),warning=function(w) {
    ebmf_warnings <<- c(ebmf_warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  if(inherits(ff,"error")) {
    rows <- lapply(c("ebmf_alone","ebmf_cox","ebmf_lasso_cox"),function(m)
      comparator_empty_row(m,error=paste("EBMF:",conditionMessage(ff))))
  } else {
    rows <- list(comparator_ebmf_rows(ff,tb,seed))
  }
  rows <- list(do.call(rbind,rows))
  rows[[1L]]$warnings <- vapply(rows[[1L]]$warnings,function(w)
    paste(unique(c(ebmf_warnings,w[nzchar(w)])),collapse=" | "),character(1))
  rows[[length(rows)+1L]] <- comparator_try_arm("pca_cox",{
    pc <- svd(tb$Ytr,nu=0,nv=7L)$v;Sp <- tb$Ytr%*%pc
    cx <- survival::coxph(survival::Surv(tb$ttr,tb$str)~Sp)
    b <- as.numeric(coef(cx));z <- as.numeric(summary(cx)$coefficients[,"z"])
    cbind(comparator_metrics(pc,b,z,tb),
      C=frozen_reverse_cindex(as.vector(tb$Yva%*%pc%*%b),tb$tva,tb$sva),
      K_final=ncol(pc),n_active=sum(abs(z)>=1.96),converged=NA,iterations=cx$iter)
  },7L)
  rows[[length(rows)+1L]] <- comparator_try_arm("supervised",{
    set.seed(seed)
    fit <- fit_cox_on_yf(tb$Ytr,tb$ttr,tb$str,K=7L,max_iter=cfg$cavi$max_iter,
      tol=cfg$cavi$tol,prior_LF="point_laplace",canonicalize=TRUE,w_S=1,
      init_supervised=TRUE,beta_prior_mode="shared",prior_beta="point_laplace",
      beta_warmup=10L,verbose=FALSE,trace_supervision=TRUE)
    sf <- score_fit(fit,tb)
    data.frame(program=sf$match_k,cor_planted=sf$cor_planted,
      recovered=sf$recovered,false_active=sf$false_z,C=sf$c_validation,
      K_final=ncol(fit$EF),n_active=sf$n_active_z,converged=sf$converged,
      iterations=fit$history$n_iter)
  },7L)
  cbind(seed=seed,noise=noise,runtime=proc.time()[["elapsed"]]-t0,do.call(rbind,rows))
}

# Fixed-denominator reporting ----

#' Summarize every dataset, distinguishing undefined C from removed observations
#' @param res All method-by-dataset rows, including errors/inapplicability.
#' @return Per-method/noise summary with the full n and metric-specific counts.
#' @examples
#' # comparator_panel_summary(results)
comparator_panel_summary <- function(res) {
  do.call(rbind,lapply(split(res,list(res$method,res$noise),drop=TRUE),function(x) {
    defined_C <- is.finite(x$C);defined_cor <- is.finite(x$cor_planted)
    data.frame(method=x$method[1],noise=x$noise[1],n=nrow(x),
      recovered=sum(x$recovered),errors=sum(!is.na(x$error)),
      n_applicable=sum(x$applicable),n_C_defined=sum(defined_C),
      mean_C_defined=if(any(defined_C))mean(x$C[defined_C]) else NA_real_,
      n_cor_defined=sum(defined_cor),mean_cor_defined=if(any(defined_cor))mean(x$cor_planted[defined_cor]) else NA_real_,
      false_active=if(anyNA(x$false_active))NA_real_ else mean(x$false_active),
      median_K=if(anyNA(x$K_final))NA_real_ else median(x$K_final),
      n_convergence_reported=sum(!is.na(x$converged)),
      nonconverged=if(all(is.na(x$converged)))NA_integer_ else sum(!x$converged,na.rm=TRUE))
  }))
}

#' Count actual zero/one/multiple-factor EBMF fits without treating failures as zero
#' @param res All comparator rows; EBMF-alone identifies each dataset once.
#' @return Counts per noise level, including method failures.
#' @examples
#' # comparator_ebmf_rank_counts(results)
comparator_ebmf_rank_counts <- function(res) {
  res <- res[res$method=="ebmf_alone",]
  do.call(rbind,lapply(split(res,res$noise),function(x)
    data.frame(noise=x$noise[1],n=nrow(x),n_zero=sum(x$K_final==0,na.rm=TRUE),
      n_one=sum(x$K_final==1,na.rm=TRUE),n_two_or_more=sum(x$K_final>=2,na.rm=TRUE),
      n_error=sum(!is.na(x$error)))))
}

#' Run the unchanged 40-dataset panel into fresh per-dataset checkpoints
#' @param cores Number of local serial socket workers.
#' @param out_dir Fresh Round2 output directory; original results stay intact.
#' @return Invisibly all results; writes full-denominator summary and rank counts.
#' @examples
#' # run_comparator_panel(6L)
run_comparator_panel <- function(cores=6L,
  out_dir="results/benchmark_sim/outputs/supervision_comparators/round2") {
  if(dir.exists(out_dir))stop("Fresh output folder required: ",out_dir)
  dir.create(file.path(out_dir,"cells"),recursive=TRUE)
  cfg <- yaml::read_yaml("config/globals.yml")
  grid <- expand.grid(seed=20261015L+1:20,noise=c(.25,1))
  write.csv(grid,file.path(out_dir,"design.csv"),row.names=FALSE)
  cl <- parallel::makeCluster(cores);on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl,c("load_comparator_code","comparator_metrics",
    "comparator_empty_row","comparator_try_arm","comparator_ebmf_rows",
    "run_comparator_dataset","cfg"),envir=environment())
  invisible(parallel::clusterEvalQ(cl,load_comparator_code()))
  rows <- parallel::parLapplyLB(cl,seq_len(nrow(grid)),function(i,g,out) {
    x <- run_comparator_dataset(g$seed[i],g$noise[i])
    write.csv(x,file.path(out,"cells",sprintf("cell_%03d.csv",i)),row.names=FALSE);x
  },g=grid,out=out_dir,chunk.size=1L)
  res <- do.call(rbind,rows)
  write.csv(res,file.path(out_dir,"results.csv"),row.names=FALSE)
  sm <- comparator_panel_summary(res)
  write.csv(sm,file.path(out_dir,"summary.csv"),row.names=FALSE)
  write.csv(comparator_ebmf_rank_counts(res),file.path(out_dir,"ebmf_rank_counts.csv"),row.names=FALSE)
  print(sm,row.names=FALSE)
  if(any(!is.na(res$error)))stop("Per-arm failure(s) retained with all other results; inspect results.csv")
  invisible(res)
}
if(sys.nframe()==0L) {
  a <- commandArgs(TRUE);run_comparator_panel(if(length(a))as.integer(a[1]) else 6L)
}
