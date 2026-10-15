# ============================================================
# REVIEW STATUS: Andrew deferred the real-data head-to-head on 2026-10-15.
# This is prepared software only. Do not execute real-data fits or scoring
# until Andrew authorizes resuming Task3 after review.
# Script: run_desurv_h2h.R
# Purpose: Compare authors' DeSurv with frozen D4 and selected supervised Single.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: DeSurv, cvwrapr, survival, yaml, ebnm, parallel; project helpers
# Run: Rscript results/benchmark_sim/run_desurv_h2h.R [all|prepare|finalize|bo] [workers]
# ============================================================
# Package SHA: cc9d4c2486b3125704ab71eea61f35e040a40080 (DeSurv1.0.1).
# Paper SHA: 0a870a78673296a3802aeddfee0514e4a8b0f9ae (rashidlab/DeSurv-paper).
# cvwrapr SHA: 9188a9e0b5fb94f227bdd19a8b9e4f46b08a8723.
# PRE-RUN design: reported TCGA+CPTAC parameters, unchanged author preprocessing,
# same eligible patients as D4. Full BO costs4-8h on30 cores (paper03); use
# Andrew's authorized reported-parameter fallback, not a reduced local search.
# Reported parameter RDS SHA256:
# 23ec3fe34cfd306cabfef64ee64edd25d60eddfd0adf756be5d6d3e977d842ec.
# Fit100 independent full seed fits then the paper04 consensus/final fit.
# Record every seed error/warning/convergence flag; never select on external C.
# Fixed positive-hazard orientation for DeSurv/selected Single; D4 uses the frozen
# documented -eta orientation and cached risks. Single is refit with seed42.
# Per-cohort C, mean C and paired cohort-stratified bootstrap B2000/seed1.
# Every new real-data report includes frozen benchmark rows, never a D4 refit.
# Runtime: author production100-start fits can exceed1h locally; seed fits are
# checkpointed independently. Budget16GB/6 serial socket workers on Longleaf;
# no nested fit parallelism. Full BO can take several hours on30 cores.
# No running analysis script is edited and no Longleaf jobs are submitted here.

# Reported settings ----

#' Return the paper's reported parameters and full independent-seed controls
#' @param seed Integer independent author seed.
#' @return Named argument list for DeSurv::desurv_fit.
#' @examples
#' desurv_reported_controls(1L)
#' @family desurv_author_pipeline
#' @seealso DeSurv::desurv_fit
desurv_reported_controls <- function(seed) {
  list(k=3L,alpha=.333935556177224,lambda=.349448737941329,
    nu=.0556512168045738,lambdaW=0,lambdaH=0,seed=seed,
    tol=1e-7,tol_init=1e-5,maxit=4000L,imaxit=4000L,ninit=1L,
    parallel_init=FALSE,verbose=FALSE)
}

# Independent seed fits ----

#' Fit one full author seed and retain its exact failure or convergence status
#' @param i Integer seed from the fixed100-start protocol.
#' @param training Author-preprocessed list with ex and sampInfo time/event.
#' @param out_dir Fresh directory for independent seed checkpoints.
#' @return Status row; writes an ignored fitted RDS and a CSV status.
#' @examples
#' # desurv_seed_fit(1L,training,"results/benchmark_sim/outputs/desurv_h2h/round2/seeds")
#' @family desurv_author_pipeline
#' @seealso desurv_reported_controls
desurv_seed_fit <- function(i,training,out_dir) {
  stopifnot(length(i)==1L,!is.na(i),i %in% 1:100)
  stem <- file.path(out_dir,sprintf("seed_%03d",i))
  if(file.exists(paste0(stem,"_status.csv")))stop("Seed checkpoint exists: ",stem)
  warnings <- character();t0 <- proc.time()[["elapsed"]]
  row <- withCallingHandlers(tryCatch({
    f <- do.call(DeSurv::desurv_fit,c(list(X=training$ex,
      y=training$sampInfo$time,d=training$sampInfo$event),desurv_reported_controls(i)))
    saveRDS(f,paste0(stem,".rds"))
    if(isTRUE(f$nan_flag) || !all(is.finite(unlist(f$loss))))
      stop("Author optimizer returned a non-finite fit or nan_flag=TRUE")
    data.frame(seed=i,converged=f$convergence,nan_flag=f$nan_flag,
      iterations=length(f$lossit),loss=f$loss$loss,error=NA_character_)
  },error=function(e)data.frame(seed=i,converged=NA,nan_flag=NA,
    iterations=NA,loss=NA_real_,error=conditionMessage(e))),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  row$warnings <- paste(unique(warnings),collapse=" | ")
  row$runtime_seconds <- proc.time()[["elapsed"]]-t0
  write.csv(row,paste0(stem,"_status.csv"),row.names=FALSE)
  row
}

# Full training-only tuning ----

#' Return the paper's full training-only Bayesian optimization specification
#' @param training Raw unscreened author-merged ex and aligned sampInfo.
#' @param workers Number of parallel grid workers.
#' @return Arguments for DeSurv::desurv_cv_bayesopt_refine; preprocessing is in-fold.
#' @examples
#' # desurv_bo_controls(raw_training,30L)
#' @family desurv_author_pipeline
#' @seealso DeSurv::desurv_cv_bayesopt_refine
desurv_bo_controls <- function(training,workers) {
  cfg <- yaml::read_yaml("config/globals.yml")
  list(X=training$ex,y=training$sampInfo$time,d=training$sampInfo$event,
    dataset=training$sampInfo$dataset,coarse_bounds=list(
      k_grid=list(lower=2L,upper=12L,type="integer"),
      alpha_grid=list(lower=0,upper=1,type="continuous"),
      lambda_grid=list(lower=1e-3,upper=1e3,scale="log10"),
      nu_grid=list(lower=0,upper=1,type="continuous"),
      ntop=list(lower=50L,upper=300L,type="integer")),
    bo_fixed=list(n_starts=30L,ngene=cfg$preprocessing$top_n_genes_desurv,
      lambdaW_grid=0,lambdaH_grid=0),preprocess=TRUE,method_trans_train="rank",
    engine="warmstart",nfolds=5L,tol=1e-5,maxit=4000L,
    max_refinements=0L,tol_gain=.002,plateau=1L,top_k=10L,
    shrink_base=.3,importance_gain=.1,coarse_control=list(n_init=50L,
      n_iter=100L,candidate_pool=4000L,exploration_weight=.01,seed=123),
    refine_control=NULL,verbose=TRUE,parallel_grid=workers>1L,ncores_grid=workers)
}

#' Apply the unchanged paper03 full-BO postprocessing
#' @param result Completed full DeSurv BO result; training data only.
#' @return List with the authors' rank assessment and final reported parameters.
#' @examples
#' # desurv_bo_summary(full_bo_result)
#' @family desurv_author_pipeline
#' @seealso desurv_bo_controls
desurv_bo_summary <- function(result) {
  path <- "data/round2_desurv_paper_bo_helpers.R"
  if(!file.exists(path) || !identical(unname(tools::md5sum(path)),
    "6f02c072e5b64c5d7582ff76773bee8b"))
    stop("Provide the unchanged bo_helpers.R from the pinned paper at ",path)
  env <- new.env(parent=globalenv());source(path,local=env)
  rank <- env$select_bo_k_by_cv_se(result)
  params <- env$standardize_bo_params(result$overall_best$params)
  if(is.null(rank$k_selected))stop("Author rank selection failed: ",rank$reason)
  # Match paper03 exactly: overwrite k in overall_best, retain other parameters.
  params$k <- rank$k_selected
  list(k_selection=rank,params=params)
}

# Head-to-head execution ----

#' Execute the fixed head-to-head or the authors' full training-only BO search
#' @param mode all for the comparison, prepare for a seed array, finalize after
#'   all100 seed checkpoints exist, or bo for full training-only Longleaf tuning.
#' @param workers Number of serial PSOCK workers; no nested seed parallelism.
#' @return Invisibly the real-data report, or the full training-only BO object.
#' @examples
#' # run_desurv_h2h("all",2L)
#' @family desurv_author_pipeline
#' @seealso desurv_author_prepare, desurv_author_risk
run_desurv_h2h <- function(mode="all",workers=2L) {
  if(!mode %in% c("all","prepare","finalize","bo"))
    stop("mode must be all, prepare, finalize or bo")
  if(dir.exists("data/round2_R_library"))
    .libPaths(c(normalizePath("data/round2_R_library"),.libPaths()))
  if(!requireNamespace("DeSurv",quietly=TRUE) ||
    !identical(packageDescription("DeSurv")$RemoteSha,
      "cc9d4c2486b3125704ab71eea61f35e040a40080"))
    stop("Install the approved pinned DeSurv package before running")
  source("results/benchmark_sim/desurv_author_adapter.R")
  source("code/concordance_ci.R")
  cfg <- yaml::read_yaml("config/globals.yml")
  source("results/benchmark_sim/benchmark_helpers.R",local=environment())
  root <- getwd();tmp <- tempfile("desurv-author-");dir.create(tmp)
  if(!file.symlink(PDAC_DATA_ROOT,file.path(tmp,"data")))stop("Cannot link author data root")
  on.exit({setwd(root);unlink(file.path(tmp,"data"));unlink(tmp)},add=TRUE)
  # This local file was verified byte-identical to the pinned paper loader.
  env <- new.env(parent=globalenv())
  loader <- file.path(PDAC_DATA_ROOT,"load_data_internal.R")
  if(!identical(unname(tools::md5sum(loader)),"3792de03e1649c4fbc1987770905aab9"))
    stop("Local author loader differs from the pinned paper pipeline")
  source(loader,local=env)
  setwd(tmp)
  raw <- lapply(setNames(cfg$pdac$training_cohorts,cfg$pdac$training_cohorts),env$load_data_internal)
  setwd(root)
  train <- desurv_author_prepare(raw)
  project_raw <- lapply(setNames(cfg$pdac$training_cohorts,cfg$pdac$training_cohorts),
    function(ds)load_pdac_raw(ds,PDAC_DATA_ROOT))
  time <- unlist(lapply(project_raw,`[[`,"time"),use.names=FALSE)
  event <- unlist(lapply(project_raw,`[[`,"status"),use.names=FALSE)
  stopifnot(identical(as.numeric(train$sampInfo$time),as.numeric(time)),
    identical(as.numeric(train$sampInfo$event),as.numeric(event)))
  benchmark_fit <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_fits.rds")$D4
  stopifnot(ncol(train$ex)==nrow(benchmark_fit$EL))
  out <- file.path("results/benchmark_sim/outputs/desurv_h2h",if(mode=="bo")"round2_full_bo" else "round2")
  if(mode=="finalize") {
    if(!file.exists(file.path(out,"training.rds")))stop("Prepare the seed array first")
    stopifnot(identical(readRDS(file.path(out,"training.rds")),train))
  } else {
    if(dir.exists(out))stop("Fresh output folder required: ",out)
    dir.create(out,recursive=TRUE)
    saveRDS(train,file.path(out,"training.rds"))
    write.csv(data.frame(cohort=names(raw),n=vapply(raw,function(x)sum(x$sampInfo$keep==1),integer(1))),
      file.path(out,"training_patient_counts.csv"),row.names=FALSE)
  }
  if(mode=="bo") {
    # Exact paper03 full DeSurv search, excluding unrelated NMF baseline searches.
    # Supply unscreened raw values so preprocessing remains inside each fold.
    raw_training <- desurv_author_merge(raw)
    ans <- do.call(DeSurv::desurv_cv_bayesopt_refine,
      desurv_bo_controls(raw_training,workers))
    saveRDS(ans,file.path(out,"full_bo.rds"))
    selected <- desurv_bo_summary(ans)
    saveRDS(selected,file.path(out,"full_bo_selection.rds"))
    write.csv(as.data.frame(selected$params),file.path(out,"full_bo_params.csv"),row.names=FALSE)
    return(invisible(ans))
  }
  seed_dir <- file.path(out,"seeds")
  if(mode!="finalize")dir.create(seed_dir)
  if(mode=="prepare")return(invisible(train))
  if(mode=="finalize") {
    files <- file.path(seed_dir,sprintf("seed_%03d_status.csv",1:100))
    if(!all(file.exists(files)))stop("All100 seed statuses, including errors, must exist before finalizing")
    rows <- lapply(files,read.csv)
    if(!identical(vapply(rows,function(x)as.integer(x$seed),integer(1)),1:100))
      stop("Seed checkpoint identities disagree with the fixed protocol")
  } else if(workers>1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(if(!is.null(cl))parallel::stopCluster(cl),add=TRUE)
    invisible(parallel::clusterEvalQ(cl,{
      if(dir.exists("data/round2_R_library")) .libPaths(c(normalizePath("data/round2_R_library"),.libPaths()))
      source("results/benchmark_sim/run_desurv_h2h.R")
    }))
    rows <- parallel::parLapplyLB(cl,1:100,desurv_seed_fit,training=train,out_dir=seed_dir,chunk.size=1L)
    parallel::stopCluster(cl);cl <- NULL
  } else rows <- lapply(1:100,desurv_seed_fit,training=train,out_dir=seed_dir)
  statuses <- do.call(rbind,rows)
  write.csv(statuses,file.path(out,"seed_status.csv"),row.names=FALSE)
  keep <- statuses$seed[is.na(statuses$error)]
  if(!length(keep))stop("All author seed fits failed; inspect every seed status")
  fits <- lapply(keep,function(i)readRDS(file.path(seed_dir,sprintf("seed_%03d.rds",i))))
  init <- DeSurv::desurv_consensus_seed(fits,train$ex,ntop=270L,k=3L,
    min_frequency=max(1L,as.integer(ceiling(.3*length(fits)))))
  p <- desurv_reported_controls(NULL)
  final <- DeSurv::desurv_fit(train$ex,time,event,k=p$k,alpha=p$alpha,
    lambda=p$lambda,nu=p$nu,lambdaW=0,lambdaH=0,W0=init$W0,H0=init$H0,
    beta0=init$beta0,seed=NULL,tol=1e-7,maxit=4000L,verbose=FALSE)
  saveRDS(final,file.path(out,"desurv_final.rds"))
  write.csv(data.frame(model="DeSurv",n_success=length(keep),n_error=sum(!is.na(statuses$error)),
    seed_caps=sum(!statuses$converged,na.rm=TRUE),converged=final$convergence,
    iterations=length(final$lossit),nan_flag=final$nan_flag),file.path(out,"fit_status.csv"),row.names=FALSE)
  if(isTRUE(final$nan_flag) || !all(is.finite(unlist(final$loss))))stop("Final DeSurv fit is non-finite")
  source("code/preprocess_desurv.R")
  tryCatch(source("code/fit_cox_on_yf.R"),error=function(e) {
    if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)
  })
  source("code/predict_cox_on_yf.R")
  pp <- preprocess_merged_cohorts(project_raw,PLATFORM_LOG_TRANSFORM[names(raw)],
    top_n=cfg$preprocessing$top_n_genes_desurv,rank_transform=FALSE,
    per_platform_standardize=TRUE,normalize_method="none",selection_per_cohort=TRUE,
    selection_method="combined_rank")
  genes <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/d4_gene_names.rds")
  stopifnot(identical(pp$gene_names,genes))
  set.seed(42)
  single <- fit_cox_on_yf(pp$Y,time,event,K=cfg$benchmark$k_merged_yfb_desurv,
    max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,canonicalize=TRUE,w_S=1,
    prior_LF="point_laplace",init_supervised=TRUE,beta_prior_mode="shared",
    prior_beta="point_laplace",beta_warmup=10L,verbose=FALSE)
  saveRDS(single,file.path(out,"selected_single.rds"))
  frozen <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_riskscores.rds")$D4
  risks <- list();results <- list()
  for(co in cfg$pdac$external_cohorts) {
    setwd(tmp);author_raw <- env$load_data_internal(co);setwd(root)
    val <- desurv_author_prepare(setNames(list(author_raw),co),genes=rownames(final$W))
    raw_val <- load_pdac_raw(co,PDAC_DATA_ROOT)
    stopifnot(identical(as.numeric(val$sampInfo$time),as.numeric(raw_val$time)),
      identical(as.numeric(val$sampInfo$event),as.numeric(raw_val$status)),
      identical(raw_val$time,frozen[[co]]$time),identical(raw_val$status,frozen[[co]]$status))
    ve <- preprocess_desurv_cohort(raw_val$Y,raw_val$gene_names,top_n=NULL,
      log_transform=PLATFORM_LOG_TRANSFORM[[co]],cohort_name=co,
      rank_transform=FALSE,per_platform_standardize=TRUE)
    common <- genes[genes %in% ve$gene_names];idx <- match(common,genes)
    sv <- predict_cox_on_yf(ve$Y[,match(common,ve$gene_names),drop=FALSE],
      single$EF[idx,,drop=FALSE],single$EBeta,EF_norms=single$EF_norms)$risk_scores
    rr <- list(DeSurv=desurv_author_risk(final,val$ex,270L),selected_single=sv,
      D4=-frozen[[co]]$risk,time=raw_val$time,event=raw_val$status)
    risks[[co]] <- rr
    for(m in c("DeSurv","selected_single"))results[[length(results)+1L]] <- data.frame(
      source="round2",model=m,cohort=co,n=length(rr$time),
      C=frozen_reverse_cindex(rr[[m]],rr$time,rr$event),
      converged=if(m=="DeSurv")final$convergence else single$history$converged)
  }
  saveRDS(risks,file.path(out,"riskscores.rds"))
  b <- read.csv("results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_results.csv")
  b <- b[b$model=="D4",]
  ref <- data.frame(source="benchmark_snapshot",model="D4",cohort=b$cohort,
    n=vapply(b$cohort,function(co)length(risks[[co]]$time),integer(1)),C=b$c_index,
    converged=benchmark_fit$history$converged)
  table <- rbind(do.call(rbind,results),ref)
  means <- do.call(rbind,lapply(split(table,table$model),function(x)
    data.frame(source=x$source[1],model=x$model[1],cohort="MEAN_OF_COHORTS",
      n=sum(x$n),C=mean(x$C),converged=x$converged[1])))
  table <- rbind(table,means)
  write.csv(table,file.path(out,"cindex_comparison.csv"),row.names=FALSE)
  diffs <- list()
  for(ref_name in c("D4","selected_single")) {
    for(co in names(risks)) {
      r <- risks[[co]]
      ci <- bootstrap_concordance_diff_ci(r$DeSurv,r[[ref_name]],r$time,r$event,
        B=2000L,seed=1L,flip_a=FALSE,flip_b=FALSE)
      diffs[[length(diffs)+1L]] <- data.frame(cohort=co,reference=ref_name,
        delta=ci$estimate,lower=ci$lower,upper=ci$upper)
    }
    ci <- bootstrap_concordance_diff_ci_stratified(lapply(risks,`[[`,"DeSurv"),
      lapply(risks,`[[`,ref_name),lapply(risks,`[[`,"time"),lapply(risks,`[[`,"event"),B=2000L,seed=1L)
    diffs[[length(diffs)+1L]] <- data.frame(cohort="MEAN_OF_COHORTS",reference=ref_name,
      delta=ci$estimate,lower=ci$lower,upper=ci$upper)
  }
  write.csv(do.call(rbind,diffs),file.path(out,"paired_differences.csv"),row.names=FALSE)
  print(table,row.names=FALSE);invisible(table)
}
# Command line ----

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_desurv_h2h(if(length(args))args[1] else "all",if(length(args)>1L)as.integer(args[2]) else 2L)
}
