# ============================================================
# Script: run_round2_fixed_prior.R
# Purpose: Evaluate fixed/preregistered candidate panels on fixed simulation and TCGA CV beds.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing fitters and benchmark helpers, yaml, survival, parallel
# Run: Rscript results/benchmark_sim/run_round2_fixed_prior.R [cores|prepare|finalize]
# ============================================================
# R1 remains the default100-cell fixed-prior design below. The approved opt-in
# profiled route is entered by run_round2_profiled_cox.R, with its separately
# committed one-change design,92 new fits and10 unchanged CV references.
# Task2-2, committed BEFORE execution. One change from selected: hold the
# initial shared point-Laplace beta prior (pi=.5,scale1) fixed beyond warmup.
# R1 full-TCGA diagnostic already fails credible-program/output-convergence
# requirements. This completes its specified simulation/CV evaluation; no
# candidate can be adopted by ignoring those failed necessary requirements.
# Simulation: unchanged calibrated beds, n400/validation120, p50 per modality,
# beta3=3.3, seeds20261016:20261035, noise .25/1, K7. Single unpruned selected;
# multimodal ordinary unprotected pruning. Original caps/tolerances/ELBO and
# output rules. Recovery>=16/20 both noises, null-active<=reference+.5, every
# returned fit converged. K-final gate awaits K-selection adoption (not adopted).
# Simulation references: saved Task1-3 K7 alpha_F0 Single and10/1 multimodal.
# TCGA: same144 patients and imputed matrices as diagnosis1, event-stratified
# fivefold CV (seed42), training-only centering/scaling for Single, same selected
# and framework controls as the committed diagnosis. Both models use the same
# fixed folds. Fit candidate and its model's diagnostic-framework reference
# in each fold, seed42. Compare paired fold C: candidate not lower by >oneSE.
# Cross-validation fit convergence/caps/errors are all reported, never hidden.
# No external data, no default/criterion changes, no inference from validation
# C selects a simulation arm. Full-TCGA diagnostic rows include frozen benchmark
# rows; CV reports also include these labelled non-CV benchmark rows, not refits.
# Fresh outputs only; every fit checkpoints errors/warnings independently.
# Longleaf: round2_fixed_prior.slurm,100 cells,8GB/cell,12h; never submit here.
# Dense TCGA fold fits can exceed1hour locally. One serial fit/socket worker;
# budget8GB/worker, no nested BLAS/fit parallelism, six local workers by default.
source("results/benchmark_sim/run_k_selection_panel.R")
source("results/multimodal_real/run_round2_collapse_diagnosis.R")

#' Assign every patient once to fixed event-stratified training-only CV folds
#' @param event Binary event vector in frozen TCGA patient order.
#' @return Integer fold vector of the same length, with values1:5.
#' @examples
#' fixed_prior_folds(rep(c(0,1),20))
fixed_prior_folds <- function(event) {
  stopifnot(length(event)>=5L,all(event %in% c(0,1)))
  set.seed(42);fold <- integer(length(event))
  for(st in unique(event)) {
    idx <- which(event==st);fold[idx] <- sample(rep(1:5,length.out=length(idx)))
  }
  fold
}

#' Create the fixed80 simulation and20 TCGA fold fit specifications
#' @param candidate Fixed or profiled; the latter adds two full-training fits.
#' @return100 unique rows for fixed,102 for profiled; each row identifies exactly one independent fit.
#' @examples
#' nrow(fixed_prior_design())
fixed_prior_design <- function(candidate="fixed") {
  stopifnot(candidate %in% c("fixed","profiled"))
  sim <- expand.grid(model=c("single","multimodal"),seed=20261015L+1:20,
    noise=c(.25,1),stringsAsFactors=FALSE)
  sim$kind <- "simulation";sim$arm <- candidate;sim$fold <- NA_integer_
  cv <- expand.grid(model=c("single","multimodal"),arm=c(candidate,"framework"),
    fold=1:5,stringsAsFactors=FALSE)
  cv$kind <- "training_cv";cv$seed <- 42L;cv$noise <- NA_real_
  d <- rbind(sim,cv[names(sim)])
  if(candidate=="profiled") {
    full <- data.frame(model=c("single","multimodal"),seed=42L,noise=NA_real_,
      kind="training_full",arm=candidate,fold=NA_integer_)
    d <- rbind(d,full[names(d)])
  }
  d
}

#' Select fixed-prior and framework settings without changing other fit options
#' @param model Single or multimodal model name.
#' @param arm Fixed, framework reference, or profiled remedy.
#' @return Existing diagnostic argument list with trace disabled only.
#' @examples
#' fixed_prior_controls("multimodal","fixed")
fixed_prior_controls <- function(model,arm) {
  a <- collapse_arm_controls(model,if(arm=="profiled")"framework" else arm)
  if(arm=="profiled")a$profile_cox_intercept <- TRUE
  if(model=="single")a$trace_supervision <- FALSE else a$trace <- FALSE
  a
}

#' Scale validation expression using training-fold statistics only
#' @param train,validation Patient-by-feature expression matrices.
#' @return List of training/validation matrices in the same training-defined units.
#' @examples
#' fixed_prior_scale(matrix(1:12,4,3),matrix(3:8,2,3))
fixed_prior_scale <- function(train,validation) {
  tr <- scale(train);s <- attr(tr,"scaled:scale");mu <- attr(tr,"scaled:center")
  if(any(!is.finite(s)) || any(s<=0))stop("Zero/nonfinite training-fold feature SD")
  list(training=tr,validation=sweep(sweep(validation,2,mu),2,s,"/"))
}

#' Fit one pre-specified simulation or TCGA fold with isolated failure records
#' @param i Integer row of fixed_prior_design().
#' @param candidate Fixed or profiled curvature remedy.
#' @param out_dir Fresh attempt output folder.
#' @return One metrics/error row; saves a checkpoint and ignored fit RDS.
#' @examples
#' # fixed_prior_cell(1L,"results/benchmark_sim/outputs/round2/fixed_prior")
fixed_prior_cell <- function(i,out_dir,candidate="fixed") {
  d <- fixed_prior_design(candidate);stopifnot(i %in% seq_len(nrow(d)))
  x <- d[i,];path <- file.path(out_dir,"cells",sprintf("cell_%03d",i))
  if(file.exists(paste0(path,".csv")))stop("Cell exists; do not overwrite")
  warnings <- character();t0 <- proc.time()[["elapsed"]]
  row <- withCallingHandlers(tryCatch({
    set.seed(x$seed)
    if(x$kind=="simulation") {
      if(x$model=="single") {
        tb <- make_testbed(x$seed,x$noise)
        train <- tb$Ytr;val <- tb$Yva;time <- tb$ttr;event <- tb$str
        vt <- tb$tva;ve <- tb$sva;truth <- tb$F_true
      } else {
        tb <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",n_train=400L,
          n_validation=120L,p_expression=50L,p_methylation=50L,seed=x$seed,
          truncate=FALSE,noise_scale=x$noise,beta=c(0,0,3.3))
        train <- tb$training$Y;val <- tb$validation$Y;time <- tb$training$time
        event <- tb$training$event;vt <- tb$validation$time;ve <- tb$validation$event
        truth <- do.call(rbind,tb$truth$F)
      }
    } else {
      tr <- readRDS("data/round2_tcga_training.rds")
      spec <- readRDS("results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned.rds")$joint$training_spec
      stopifnot(identical(rownames(tr$Y$expression),spec$subject_ids),
        identical(lapply(tr$Y,colnames),spec$feature_names),
        identical(names(tr$time),spec$subject_ids),identical(names(tr$event),spec$subject_ids))
      fold <- fixed_prior_folds(tr$event);va <- if(x$kind=="training_full")rep(FALSE,length(fold)) else fold==x$fold
      time <- tr$time[!va];event <- tr$event[!va];vt <- tr$time[va];ve <- tr$event[va]
      validation_rows <- if(x$kind=="training_full")rep(TRUE,length(va)) else va
      vt <- tr$time[validation_rows];ve <- tr$event[validation_rows]
      if(x$model=="single") {
        z <- fixed_prior_scale(tr$Y$expression[!va,,drop=FALSE],tr$Y$expression[validation_rows,,drop=FALSE])
        train <- z$training;val <- z$validation
      } else {
        train <- lapply(tr$Y,function(y)y[!va,,drop=FALSE])
        val <- lapply(tr$Y,function(y)y[validation_rows,,drop=FALSE])
      }
    }
    a <- fixed_prior_controls(x$model,x$arm);set.seed(x$seed)
    if(x$model=="single") {
      fit <- do.call(fit_cox_on_yf,c(list(train,time,event),a))
      risk <- predict_cox_on_yf(val,fit$EF,fit$EBeta,fit$EF_norms)$risk_scores
      conv <- fit$history$converged;iters <- fit$history$n_iter
      caps <- sum(!as.logical(fit$pruning$converged_before));F <- fit$EF
    } else {
      fit <- fit_multimodal_yfb(train,time,event,7L,a)
      risk <- predict_multimodal_yfb(fit,val)$risk_scores
      conv <- fit$diagnostics$converged;iters <- fit$diagnostics$iterations
      caps <- sum(!as.logical(fit$diagnostics$pruning$converged_before));F <- do.call(rbind,fit$EF)
    }
    prior_metrics <- NULL
    if(candidate=="profiled" && x$arm=="profiled") {
      used <- if(x$model=="single")isTRUE(fit$profile_cox_intercept) else
        isTRUE(fit$diagnostics$controls$profile_cox_intercept)
      if(!used)stop("Fit did not use the requested profiled Cox curvature")
      prior_metrics <- data.frame(profiled_curvature=TRUE)
    }
    saveRDS(fit,paste0(path,".rds"))
    sd <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
    active <- sd>0 & abs(fit$EBeta)>=1.96*sd
    rho <- NA_real_;recovered <- NA;false <- NA_real_
    if(x$kind=="simulation") {
      rr <- suppressWarnings(cor(truth[,3],F));rr[!is.finite(rr)] <- 0
      k <- which.max(abs(rr));rho <- rr[k]
      recovered <- abs(rho)>=.7 && sign(rho*fit$EBeta[k])==1
      false <- sum(active[setdiff(seq_along(active),if(abs(rho)>=.7)k else integer())])
    }
    metrics <- data.frame(error=NA_character_,n_train=length(time),n_validation=length(vt),
      C=frozen_reverse_cindex(risk,vt,ve),risk_sd=stats::sd(risk),
      K_final=ncol(fit$EL),n_active=sum(active),recovered=recovered,
      cor_planted=abs(rho),false_active=false,converged=conv,iterations=iters,pruning_caps=caps)
    if(!is.null(prior_metrics))metrics <- cbind(metrics,prior_metrics)
    metrics
  },error=function(e)data.frame(error=conditionMessage(e))),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  # Error-only batches still have the complete metric schema for aggregation.
  metric_names <- c("n_train","n_validation","C","risk_sd","K_final","n_active",
    "recovered","cor_planted","false_active","converged","iterations","pruning_caps")
  for(n in setdiff(metric_names,names(row)))row[[n]] <- NA
  row <- cbind(x,runtime_seconds=proc.time()[["elapsed"]]-t0,warnings=paste(unique(warnings),collapse=" | "),row)
  write.csv(row,paste0(path,".csv"),row.names=FALSE);row
}

#' Summarize fixed criteria without removing failed seeds or undefined folds
#' @param res All100 checkpoint rows, with missing metric columns filled NA.
#' @param reference Saved original K7 simulation reference rows.
#' @param candidate Fixed or profiled curvature remedy.
#' @return Simulation/CV tables and optional full-training checks; failures cannot pass.
#' @examples
#' # fixed_prior_summary(results,reference)
fixed_prior_summary <- function(res,reference,candidate="fixed") {
  sim <- res[res$kind=="simulation",]
  sm <- do.call(rbind,lapply(split(sim,list(sim$model,sim$noise),drop=TRUE),function(x) {
    ref <- reference[reference$model==x$model[1] & reference$noise==x$noise[1],]
    fp <- if(anyNA(x$false_active))NA_real_ else mean(x$false_active)
    rc <- sum(x$recovered %in% TRUE);nc <- sum(x$converged %in% TRUE)
    valid <- nrow(x)==20L && all(is.na(x$error))
    data.frame(model=x$model[1],noise=x$noise[1],n=nrow(x),recovered=rc,
      false_active=fp,reference_false_active=mean(ref$false_active),converged=nc,
      errors=sum(!is.na(x$error)),recovery_pass=valid && rc>=16L,
      false_positive_pass=valid && is.finite(fp) && fp<=mean(ref$false_active)+.5,
      convergence_pass=valid && nc==20L)
  }))
  cv <- res[res$kind=="training_cv",]
  cm <- do.call(rbind,lapply(split(cv,cv$model),function(x) {
    a <- x[x$arm==candidate,];b <- x[x$arm=="framework",]
    z <- merge(a,b,by="fold",suffixes=c("_fixed","_reference"))
    valid <- nrow(z)==5L && all(is.na(z$error_fixed)) && all(is.na(z$error_reference)) &&
      all(is.finite(z$C_fixed)) && all(is.finite(z$C_reference))
    delta <- if(valid)z$C_fixed-z$C_reference else NA_real_
    m <- if(valid)mean(delta) else NA_real_;se <- if(valid)stats::sd(delta)/sqrt(5) else NA_real_
    data.frame(model=x$model[1],n_folds=nrow(z),candidate_C=if(valid)mean(z$C_fixed) else NA_real_,
      reference_C=if(valid)mean(z$C_reference) else NA_real_,delta=m,paired_SE=se,
      pass=valid && m>= -se,converged_candidate=sum(a$converged %in% TRUE),
      converged_reference=sum(b$converged %in% TRUE),errors=sum(!is.na(x$error)))
  }))
  full <- res[res$kind=="training_full",]
  if(nrow(full)) {
    full$no_collapse <- is.na(full$error) & !is.na(full$n_active) & full$n_active>=1L &
      is.finite(full$risk_sd) & full$risk_sd>0
  }
  if(candidate=="profiled") {
    cm$candidate_finite_folds <- cm$reference_finite_folds <- 0L
    cm$candidate_mean_if_all5 <- cm$reference_mean_if_all5 <- NA_real_
    for(i in seq_len(nrow(cm))) {
      z <- cv[cv$model==cm$model[i],]
      a <- z[z$arm==candidate,];b <- z[z$arm=="framework",]
      cm$candidate_finite_folds[i] <- sum(is.finite(a$C) & is.na(a$error))
      cm$reference_finite_folds[i] <- sum(is.finite(b$C) & is.na(b$error))
      if(nrow(a)==5L && cm$candidate_finite_folds[i]==5L)
        cm$candidate_mean_if_all5[i] <- mean(a$C)
      if(nrow(b)==5L && cm$reference_finite_folds[i]==5L)
        cm$reference_mean_if_all5[i] <- mean(b$C)
    }
  }
  list(simulation=sm,training_cv=cm,training_full=full)
}

#' Compare checkpoint keys without CSV's all-NA integer/logical type ambiguity
#' @param row One loaded checkpoint row.
#' @param design_row The corresponding row of the fixed design.
#' @return TRUE only when every required key value, including NA, is unchanged.
#' @examples
#' fixed_prior_checkpoint_matches(fixed_prior_design()[1,],fixed_prior_design()[1,])
fixed_prior_checkpoint_matches <- function(row,design_row) {
  nrow(row)==1L && all(names(design_row) %in% names(row)) &&
    all(vapply(names(design_row),function(n)
      identical(as.character(row[[n]]),as.character(design_row[[n]])),logical(1)))
}

#' Reuse the ten unchanged TCGA CV framework references from completed R1
#' @param out_dir Prepared fresh profiled-curvature attempt folder.
#' @param source_dir Completed fixed-prior attempt folder.
#' @return Invisibly the ten checkpoint source/MD5 records, copied without refitting.
#' @examples
#' # fixed_prior_reuse_cv("results/benchmark_sim/outputs/round2/profiled")
fixed_prior_reuse_cv <- function(out_dir,
  source_dir="results/benchmark_sim/outputs/round2/fixed_prior") {
  if(!file.exists(file.path(source_dir,"results.csv")))stop("R1 must complete before reference reuse")
  if(!identical(readLines(file.path(out_dir,"training_source.txt")),
    readLines(file.path(source_dir,"training_source.txt"))) ||
    !identical(unname(tools::md5sum(file.path(out_dir,"folds.csv"))),
      unname(tools::md5sum(file.path(source_dir,"folds.csv")))))
    stop("CV reference training cache or folds differ")
  d <- fixed_prior_design("profiled");old <- fixed_prior_design()
  complete <- read.csv(file.path(source_dir,"results.csv"))
  if(nrow(complete)!=nrow(old) ||
    !all(vapply(seq_len(nrow(old)),function(i)
      fixed_prior_checkpoint_matches(complete[i,],old[i,]),logical(1))))
    stop("R1 aggregate is incomplete or differs from its fixed design")
  ids <- which(d$kind=="training_cv" & d$arm=="framework")
  records <- list()
  for(i in ids) {
    src <- file.path(source_dir,"cells",sprintf("cell_%03d.csv",i))
    dst <- file.path(out_dir,"cells",sprintf("cell_%03d.csv",i))
    row <- read.csv(src)
    if(!fixed_prior_checkpoint_matches(row,old[i,]) ||
      !fixed_prior_checkpoint_matches(row,d[i,]) || !is.na(row$error))
      stop("Invalid cached framework reference at cell",i)
    if(file.exists(dst))stop("Reference checkpoint exists; never overwrite")
    row$reference_source <- src;row$reference_md5 <- unname(tools::md5sum(src))
    write.csv(row,dst,row.names=FALSE)
    records[[length(records)+1L]] <- data.frame(cell=i,source=src,md5=row$reference_md5)
  }
  records <- do.call(rbind,records)
  write.csv(records,file.path(out_dir,"reused_cv_references.csv"),row.names=FALSE)
  invisible(records)
}

#' Execute all100 fixed-prior evaluation cells after preregistration
#' @param cores Positive number of serial PSOCK workers.
#' @param out_dir Fresh attempt folder, never overwrites earlier fits.
#' @param mode All locally, prepare Longleaf inputs only, or finalize all checkpoints.
#' @param candidate Fixed or profiled curvature remedy; default retains R1 design.
#' @return Invisibly aggregate results; writes tables before reporting any errors.
#' @examples
#' # run_fixed_prior_panel(6L)
run_fixed_prior_panel <- function(cores=6L,
  out_dir="results/benchmark_sim/outputs/round2/fixed_prior",mode="all",candidate="fixed") {
  stopifnot(candidate %in% c("fixed","profiled"))
  if(candidate=="profiled" && out_dir=="results/benchmark_sim/outputs/round2/fixed_prior")
    stop("A fresh profiled output folder is required")
  if(!mode %in% c("all","prepare","finalize"))stop("mode must be all, prepare or finalize")
  if(mode!="finalize" && dir.exists(out_dir))stop("Fresh output folder required")
  if(mode=="finalize" && (!dir.exists(out_dir) || file.exists(file.path(out_dir,"results.csv"))))
    stop("Finalization requires existing cells and no previous aggregate")
  stopifnot(length(cores)==1L,is.finite(cores),cores>=1L,cores==as.integer(cores))
  reference <- read.csv("results/benchmark_sim/outputs/k_selection/round2/attempt3/results.csv")
  reference <- reference[reference$K==7L & reference$arm=="reference",]
  stopifnot(nrow(reference)==80L,all(is.na(reference$error)),!anyDuplicated(reference[c("model","seed","noise")]))
  d <- fixed_prior_design(candidate)
  if(mode!="finalize") {
  dir.create(file.path(out_dir,"cells"),recursive=TRUE)
  write.csv(d,file.path(out_dir,"design.csv"),row.names=FALSE)
  write.csv(reference,file.path(out_dir,"simulation_references.csv"),row.names=FALSE)
  tr <- readRDS("data/round2_tcga_training.rds")
  write.csv(data.frame(patient_id=rownames(tr$Y$expression),event=tr$event,
    fold=fixed_prior_folds(tr$event)),file.path(out_dir,"folds.csv"),row.names=FALSE)
  writeLines(paste("training_md5",tools::md5sum("data/round2_tcga_training.rds")),
    file.path(out_dir,"training_source.txt"))
  }
  if(candidate=="profiled" && mode!="finalize")fixed_prior_reuse_cv(out_dir)
  if(mode=="prepare")return(invisible(d))
  expected_md5 <- paste("training_md5",tools::md5sum("data/round2_tcga_training.rds"))
  if(!identical(readLines(file.path(out_dir,"training_source.txt")),expected_md5))
    stop("TCGA training cache changed after preparation")
  if(mode=="all") {
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl,c("load_panel_code","collapse_arm_controls","fixed_prior_design",
    "fixed_prior_controls","fixed_prior_folds","fixed_prior_scale","fixed_prior_cell"),envir=environment())
  invisible(parallel::clusterEvalQ(cl,{
    load_panel_code();source("code/concordance_ci.R")
  }))
  ids <- if(candidate=="profiled")which(!(d$kind=="training_cv" & d$arm=="framework")) else seq_len(nrow(d))
  fresh <- tryCatch(parallel::parLapplyLB(cl,ids,fixed_prior_cell,out_dir=out_dir,candidate=candidate,chunk.size=1L),
    finally=parallel::stopCluster(cl))
  rows <- lapply(file.path(out_dir,"cells",sprintf("cell_%03d.csv",seq_len(nrow(d)))),read.csv)
  } else {
    paths <- file.path(out_dir,"cells",sprintf("cell_%03d.csv",seq_len(nrow(d))))
    if(!all(file.exists(paths)))stop("All design cell checkpoints required before finalizing")
    rows <- lapply(paths,read.csv)
    for(i in seq_along(rows)) {
      if(!fixed_prior_checkpoint_matches(rows[[i]],d[i,]))
        stop("Checkpoint design mismatch at cell",i)
    }
  }
  cols <- unique(unlist(lapply(rows,names)))
  res <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(res,file.path(out_dir,"results.csv"),row.names=FALSE)
  s <- fixed_prior_summary(res,reference,candidate)
  write.csv(s$simulation,file.path(out_dir,"simulation_summary.csv"),row.names=FALSE)
  old <- read.csv("results/multimodal_real/outputs/round2/diagnosis1/comparison_all.csv")
  # A frozen training C is a labelled reference, not a fivefold CV estimate.
  cm <- s$training_cv;cm$source <- "round2_training_cv";cm$arm <- paste0(candidate,"_vs_framework")
  benchmark <- old[old$source=="benchmark_snapshot",c("source","model","arm","n","K","active","training_C")]
  rows <- list(cm,benchmark);cols <- unique(unlist(lapply(rows,names)))
  cm <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(cm,file.path(out_dir,"training_cv_summary.csv"),row.names=FALSE)
  cv <- res[res$kind=="training_cv",];cv$source <- "round2_training_cv";cv$cv_C <- cv$C
  rows <- list(cv,old[old$source=="benchmark_snapshot",]);cols <- unique(unlist(lapply(rows,names)))
  comparison <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(comparison,file.path(out_dir,"training_cv_with_benchmark.csv"),row.names=FALSE)
  if(nrow(s$training_full)) {
    full <- s$training_full;full$source <- "round2_training_full";full$training_C <- full$C
    rows <- list(full,old[old$source=="benchmark_snapshot" | old$arm %in% c("selected","framework"),]);cols <- unique(unlist(lapply(rows,names)))
    full <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
    write.csv(full,file.path(out_dir,"training_full_with_benchmark.csv"),row.names=FALSE)
  }
  print(s)
  if(any(!is.na(res$error)))stop("Fit errors retained in checkpoints and summaries; inspect results.csv")
  invisible(res)
}
if(sys.nframe()==0L) {
  a <- commandArgs(TRUE)
  if(length(a) && a[1]=="cell") {
    load_panel_code();source("code/concordance_ci.R")
    fixed_prior_cell(as.integer(a[2]),"results/benchmark_sim/outputs/round2/fixed_prior")
  } else if(length(a) && a[1] %in% c("prepare","finalize"))run_fixed_prior_panel(mode=a[1])
  else run_fixed_prior_panel(if(length(a))as.integer(a[1]) else 6L)
}
