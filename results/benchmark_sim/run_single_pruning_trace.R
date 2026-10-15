# ============================================================
# Script: run_single_pruning_trace.R
# Purpose: Observe single-modality pruning transitions without changing fitting.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing single-modality fitter and calibrated panel definitions
# ============================================================
# Task1-4 diagnosis, committed before execution. Same selected protected model,
# first prescribed seed20261016, K7, both noise levels; no new statistical arm.
# Record the unmodified fit at each pruning nullcheck and its actual gains.
# Determine the first transition at which loading recovery or credible survival
# protection is lost. No threshold, data, update, initialization or criterion is
# changed. This observes the failure; it cannot adopt a remedy. No external data.
# Local runtime minutes, memory under 1GB; two serial fits, no nested workers.
# Run: Rscript results/benchmark_sim/run_single_pruning_trace.R

#' Observe the fit states presented to the unmodified pruning nullcheck
#' @param Y Numeric subject-by-feature matrix.
#' @param time,status Training survival vectors.
#' @param fit_args Named arguments for the unmodified single-modality fitter.
#' @return List of the unchanged fit and observed pruning-stage states.
#' @examples
#' # trace_single_pruning(Y,time,event,list(K=7L,prune=TRUE,canonicalize=TRUE))
trace_single_pruning <- function(Y,time,status,fit_args) {
  original <- get("prune_yf_nullcheck",envir=.GlobalEnv)
  stages <- list()
  # Only this R process's function binding is wrapped; no source file is edited.
  # The original nullcheck receives every argument and its return is unchanged.
  observer <- function(fit,...) {
    delta <- original(fit,...)
    stages[[length(stages)+1L]] <<- list(fit=fit,delta=delta)
    delta
  }
  assign("prune_yf_nullcheck",observer,envir=.GlobalEnv)
  on.exit(assign("prune_yf_nullcheck",original,envir=.GlobalEnv))
  fit <- do.call(fit_cox_on_yf,c(list(Y,time,status),fit_args))
  list(fit=fit,stages=stages)
}

#' Trace the first fixed panel seed at both noise levels and retain every failure
#' @param out_dir Fresh diagnosis folder, separate from the ongoing rank panel.
#' @return Stage-level data frame; raw posterior states saved outside git.
#' @examples
#' # run_single_pruning_trace()
run_single_pruning_trace <- function(out_dir="results/benchmark_sim/outputs/k_selection/round2/attempt4_trace") {
  if(dir.exists(out_dir))stop("Diagnosis folder exists; do not overwrite")
  dir.create(out_dir,recursive=TRUE)
  source("results/benchmark_sim/run_k_selection_panel.R",local=.GlobalEnv)
  load_panel_code()
  cfg <- yaml::read_yaml("config/globals.yml")
  a <- list(K=7L,max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,
    prior_LF="point_laplace",canonicalize=TRUE,w_S=1,init_supervised=TRUE,
    beta_prior_mode="shared",prior_beta="point_laplace",beta_warmup=10L,
    prune=TRUE,prune_protect_survival=TRUE,converge_on="elbo",
    trace_supervision=TRUE,verbose=FALSE)
  rows <- list()
  for(noise in c(.25,1)) {
    warnings <- character()
    tb <- make_testbed(20261016L,noise)
    ans <- withCallingHandlers(tryCatch({
      tr <- trace_single_pruning(tb$Ytr,tb$ttr,tb$str,a)
      saveRDS(tr,file.path(out_dir,paste0("states_noise",noise,".rds")))
      states <- tr$stages
      # The terminal K=1 fit has no nullcheck; include that returned state too.
      states[[length(states)+1L]] <- list(fit=tr$fit,delta=NULL)
      do.call(rbind,lapply(seq_along(states),function(i) {
        f <- states[[i]]$fit;delta <- states[[i]]$delta
        rho <- suppressWarnings(cor(tb$F_true[,3],f$EF));rho[!is.finite(rho)] <- 0
        k <- which.max(abs(rho));sd <- sqrt(pmax(f$EBeta2[k]-f$EBeta[k]^2,0))
        g <- f$g_beta
        data.frame(seed=20261016L,noise=noise,stage=i,K=ncol(f$EF),
          terminal=i==length(states),cor_planted=abs(rho[k]),
          recovered=abs(rho[k])>=.7 && sign(rho[k]*f$EBeta[k])==1,
          matched_factor=k,beta=f$EBeta[k],beta_sd=sd,
          protected=sd>0 && abs(f$EBeta[k])>=1.96*sd,
          matched_delta=if(length(delta))delta[k] else NA_real_,
          dropped=if(length(delta) && max(delta)>=0)which.max(delta) else NA_integer_,
          prior_pi=sum(g$pi[g$scale>0]),prior_scale=max(g$scale),
          converged=f$history$converged,iterations=f$history$n_iter,error=NA_character_)
      }))
    },error=function(e)data.frame(seed=20261016L,noise=noise,error=conditionMessage(e))),
      warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
    ans$warnings <- paste(unique(warnings),collapse=" | ")
    rows[[length(rows)+1L]] <- ans
  }
  cols <- unique(unlist(lapply(rows,names)))
  out <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(out,file.path(out_dir,"stages.csv"),row.names=FALSE)
  print(out,row.names=FALSE);invisible(out)
}
if(sys.nframe()==0L)run_single_pruning_trace()
