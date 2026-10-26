# ============================================================
# Script: run_round2_collapse_diagnosis.R
# Purpose: Diagnose shared-prior collapse on unchanged TCGA training data only.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing single/multimodal fitters, survival, ebnm, yaml, parallel
# Run: Rscript results/multimodal_real/run_round2_collapse_diagnosis.R [cell|all]
# ============================================================
# Task2-1 design, to be committed BEFORE fitting. No remedy is selected here.
# K7, seed42, original model-specific iteration caps/tolerances. Five arms per
# model: selected; starting shared Laplace prior held fixed; shared Normal with
# the same ten-iteration warmup; historical per-factor Normal without supervised
# start; selected without supervised start. Each contrast changes one option;
# the historical framework is a separate reference. Single remains unpruned;
# multimodal retains its specified ordinary pruning. No protection was adopted.
# Compare prior before/after each sweep with coefficient extinction and the actual
# reconstruction/eta changes; a fitted point mass followed by persistent zero
# coefficients supports the absorbing-prior hypothesis. Do not alter convergence.
# This is diagnosis, not a simulation/CV adoption attempt. No external predictions.
# Input: unchanged $training extracted from the existing cache with Andrew's
# explicit approval; extraction avoids re-imputation of the screened CpGs.
# Runtime: 500 outer x 5 inner sweeps on 144 x (3000+10000) can exceed one hour
# per multimodal arm. Memory: budget 8GB per serial arm (dense matrices and
# temporary projections); six local PSOCK workers, no nested fit parallelism.
# Longleaf array: results/multimodal_real/round2_collapse_diagnosis.slurm.
# Every arm saves its fit, trace and status independently; errors remain rows.

#' Return the five fixed diagnostic specifications for either model
#' @param model Single or multimodal model name.
#' @param arm One of selected, fixed, normal, framework, no_start.
#' @return Arguments for fit_cox_on_yf or the multimodal control list.
#' @examples
#' collapse_arm_controls("multimodal", "fixed")
collapse_arm_controls <- function(model, arm) {
  stopifnot(model %in% c("single","multimodal"),
    arm %in% c("selected","fixed","normal","framework","no_start"))
  cfg <- yaml::read_yaml("config/globals.yml")
  if(model=="single") {
    a <- list(K=7L,max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,
      prior_LF="point_laplace",canonicalize=TRUE,w_S=1,init_supervised=TRUE,
      beta_prior_mode="shared",prior_beta="point_laplace",beta_warmup=10L,
      trace_supervision=TRUE,converge_on="elbo",verbose=FALSE)
    if(arm=="fixed")a$beta_prior_fix <- TRUE
    if(arm=="normal")a$prior_beta <- "normal"
    if(arm=="framework") {
      a$init_supervised <- FALSE;a$beta_prior_mode <- "per_program";a$prior_beta <- "normal"
    }
  } else {
    a <- list(intercept=TRUE,prior_update="ebnm",prune=TRUE,survival_elbo="plugin",
      prior_F=list(expression="point_laplace",methylation="point_laplace"),
      init_supervised=TRUE,beta_prior="shared_point_laplace",beta_warmup=10L,trace=TRUE)
    if(arm=="fixed")a$beta_prior_fix <- TRUE
    if(arm=="normal")a$beta_prior <- "shared_normal"
    if(arm=="framework") {a$init_supervised <- FALSE;a$beta_prior <- "eb_normal"}
  }
  if(arm=="no_start")a$init_supervised <- FALSE
  a
}

#' Flatten every recorded multimodal sweep without changing its numerical scale
#' @param trace List of outer records, including all recursive pruning stages.
#' @return Data frame with one row per stage, outer iteration, inner sweep and factor.
#' @examples
#' # collapse_trace_rows(fit$diagnostics$trace)
collapse_trace_rows <- function(trace) {
  rows <- list()
  for(j in seq_along(trace)) {
    x <- trace[[j]]
    for(i in seq_along(x$prior_sweeps)) {
      s <- x$prior_sweeps[[i]]
      # ebnm mixture weights multiply zero-centered components: only positive
      # scale components contribute to slab probability. A zero-scale component
      # is a point mass even when the mixture has a single weight equal to one.
      gb <- s$prior_before;ga <- s$prior_after
      sb <- if(is.null(gb)) numeric() else if(!is.null(gb$scale))gb$scale else gb$sd
      sa <- if(is.null(ga)) numeric() else if(!is.null(ga$scale))ga$scale else ga$sd
      rows[[length(rows)+1L]] <- data.frame(stage_outer=j,K=x$K,outer=x$iteration,
        inner=i,factor=seq_len(x$K),
        prior_pi_before=if(is.null(gb))NA_real_ else sum(gb$pi[sb>0]),
        prior_pi_after=if(is.null(ga))NA_real_ else sum(ga$pi[sa>0]),
        prior_scale_before=if(length(sb))max(sb) else NA_real_,
        prior_scale_after=if(length(sa))max(sa) else NA_real_,
        beta=s$beta_mean,beta_sd=s$beta_sd,x=s$x,A=s$A,
        canonical_scale=s$canonical_scale,live_programs=x$live_programs,
        reconstruction=x$elbo_components["reconstruction"],
        survival=x$elbo_components["survival"],KL=x$elbo_components["kl"],
        reconstruction_change=max(x$reconstruction_change),
        full_reconstruction_change=max(x$full_reconstruction_change),
        centered_reconstruction_change=max(x$centered_reconstruction_change),eta_change=x$eta_change,
        centered_eta_change=x$centered_eta_change,eta_mean=x$eta_mean,eta_sd=x$eta_sd,
        row.names=NULL)
    }
  }
  do.call(rbind,rows)
}

#' Fit one specified training-only diagnosis arm and retain errors and warnings
#' @param i Integer cell in 1:10 (five single, then five multimodal arms).
#' @param out_dir Fresh output folder for this committed attempt.
#' @param training_path TCGA-only training cache; defaults to the approved extraction.
#' @return Status row; writes arm-specific fit, trace and status checkpoints.
#' @examples
#' # run_collapse_cell(6L, "results/multimodal_real/outputs/round2/diagnosis1")
run_collapse_cell <- function(i,out_dir,training_path="data/round2_tcga_training.rds") {
  stopifnot(length(i)==1L,!is.na(i),i %in% 1:10)
  arms <- c("selected","fixed","normal","framework","no_start")
  model <- if(i<=5L)"single" else "multimodal";arm <- arms[(i-1L)%%5L+1L]
  stem <- file.path(out_dir,paste(model,arm,sep="_"))
  if(file.exists(paste0(stem,"_status.csv")))stop("Cell already exists; do not overwrite: ",stem)
  warnings <- character();t0 <- proc.time()[["elapsed"]];set.seed(42)
  ans <- withCallingHandlers(tryCatch({
  source("code/fit_multimodal_yfb.R",local=.GlobalEnv)
  source("code/concordance_ci.R",local=.GlobalEnv)
  tryCatch(source("code/fit_cox_on_yf.R",local=.GlobalEnv),error=function(e) {
    if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)
  })
  tr <- readRDS(training_path)
  spec <- readRDS("results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned.rds")$joint$training_spec
  stopifnot(identical(rownames(tr$Y$expression),spec$subject_ids),
    identical(lapply(tr$Y,colnames),spec$feature_names),
    identical(names(tr$time),spec$subject_ids),identical(names(tr$event),spec$subject_ids))
    a <- collapse_arm_controls(model,arm)
    fit <- if(model=="single")do.call(fit_cox_on_yf,c(list(scale(tr$Y$expression),tr$time,tr$event),a)) else
      fit_multimodal_yfb(tr$Y,tr$time,tr$event,7L,a)
    saveRDS(fit,paste0(stem,".rds"))
    if(model=="multimodal") {
      z <- collapse_trace_rows(fit$diagnostics$trace)
      write.csv(z,paste0(stem,"_trace.csv"),row.names=FALSE)
      h <- tail(fit$diagnostics$history,1)[[1]]
      eta <- multimodal_yfb_projection_moments(
        Map(function(y,mu)sweep(y,2,mu),tr$Y,fit$mu),fit$EF,fit$EF2)$EZ %*% fit$EBeta
      converged <- fit$diagnostics$converged;iterations <- fit$diagnostics$iterations
      rc <- h$max_relative_reconstruction_change;ec <- h$max_relative_eta_change
      caps <- sum(!as.logical(fit$diagnostics$pruning$converged_before))
    } else {
      eta <- scale(tr$Y$expression) %*% fit$EF %*% fit$EBeta
      converged <- fit$history$converged;iterations <- fit$history$n_iter
      rc <- NA_real_;ec <- tail(fit$history$output_change,1)
      caps <- sum(!as.logical(fit$pruning$converged_before))
      # Same raw-prior diagnostic schema as multimodal, one sweep per iteration.
      pr <- fit$history$prior_trace
      prior_rows <- lapply(seq_along(pr), function(j) {
        x <- pr[[j]]
        #' Effective probability assigned to positive-scale prior components.
        #' @param g EBNM mixture or NULL for the per-factor reference.
        #' @return Numeric probability, or NA when no shared prior exists.
        #' @examples
        #' # gpi(ebnm::laplacemix(c(.5,.5),c(0,0),c(0,1)))
        gpi <- function(g) {
          if(is.null(g)) return(NA_real_)
          sc <- if(!is.null(g$scale)) g$scale else g$sd
          sum(g$pi[sc>0])
        }
        data.frame(outer=j,factor=seq_along(x$x),x=x$x,A=x$A,
          beta=x$beta_mean,beta_sd=x$beta_sd,
          prior_pi_before=gpi(x$prior_before),prior_pi_after=gpi(x$prior_after))
      })
      if(length(prior_rows))
        write.csv(do.call(rbind,prior_rows),paste0(stem,"_prior_trace.csv"),row.names=FALSE)
      bt <- fit$history$trace
      z <- data.frame(iteration=rep(seq_len(iterations),each=ncol(fit$EF)),
        factor=rep(seq_len(ncol(fit$EF)),times=iterations),
        beta=as.vector(t(bt$beta_mean)),beta_second=as.vector(t(bt$beta_second)))
      write.csv(z,paste0(stem,"_trace.csv"),row.names=FALSE)
    }
    sd <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0));active <- sd>0 & abs(fit$EBeta)>=1.96*sd
    data.frame(source="round2",model=model,arm=arm,n=length(tr$time),K=ncol(fit$EL),
      active=sum(active),risk_sd=stats::sd(as.vector(eta)),
      training_C=frozen_reverse_cindex(as.vector(eta),tr$time,tr$event),
      converged=converged,iterations=iterations,pruning_caps=caps,
      reconstruction_change=rc,eta_change=ec,error=NA_character_)
  },error=function(e)data.frame(source="round2",model=model,arm=arm,
      error=conditionMessage(e))),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  ans$warnings <- paste(unique(warnings),collapse=" | ")
  ans$runtime_seconds <- proc.time()[["elapsed"]]-t0
  write.csv(ans,paste0(stem,"_status.csv"),row.names=FALSE)
  ans
}

#' Run the fixed diagnosis locally or as one Longleaf array cell
#' @param cell "all" or an integer in 1:10.
#' @param out_dir New attempt folder; cells cannot overwrite existing checkpoints.
#' @return Invisibly status rows with immutable benchmark training rows appended.
#' @examples
#' # run_round2_collapse_diagnosis("all")
run_round2_collapse_diagnosis <- function(cell="all",
  out_dir="results/multimodal_real/outputs/round2/diagnosis1") {
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  ids <- if(identical(cell,"all"))1:10 else as.integer(cell)
  stopifnot(all(ids %in% 1:10),!anyNA(ids))
  if(length(ids)>1L) {
    cl <- parallel::makeCluster(6L)
    on.exit(parallel::stopCluster(cl))
    invisible(parallel::clusterEvalQ(cl,source("results/multimodal_real/run_round2_collapse_diagnosis.R")))
    rows <- parallel::parLapplyLB(cl,ids,run_collapse_cell,out_dir=out_dir)
  } else rows <- list(run_collapse_cell(ids,out_dir))
  # Every reported training result carries the frozen benchmark on the same
  # TCGA training donors; no external predictions are computed in this task.
  b <- read.csv("results/benchmark_snapshot_2026-10-14/multimodal/real_sweep_summary.csv")
  b <- b[b$K_init==7L & b$cohort=="tcga_train" &
    ((b$config=="pruned" & b$method=="joint_multimodal_yfb") |
     (b$config=="pruned_expression_only" & b$method %in%
       c("expression_only_yfb","joint_multimodal_yfb"))),]
  benchmark_model <- ifelse(b$config=="pruned","multimodal_framework",
    ifelse(b$method=="joint_multimodal_yfb","expression_only_framework","single_alpha_F0"))
  ref <- data.frame(source="benchmark_snapshot",model=benchmark_model,arm=b$config,n=144L,
    K=b$K_final,active=b$n_survival_active,training_C=b$c_index,converged=b$converged)
  rows <- c(rows,list(ref));cols <- unique(unlist(lapply(rows,names)))
  ans <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(ans,file.path(out_dir,paste0("comparison_",cell,".csv")),row.names=FALSE)
  print(ans,row.names=FALSE);invisible(ans)
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_round2_collapse_diagnosis(if(length(args))args[1] else "all")
}
