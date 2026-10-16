# ============================================================
# Script: run_feature_expansion_multimodal.R
# Purpose: Evaluate the full-range methylation-only curve at fixed K under the historical framework.
# Created: 2026-10-15
# Dependencies: yaml, survival, existing multimodal fitter and feature helpers
# Run from project root after prospective commit: Rscript .../run_feature_expansion_multimodal.R design|all|CELL
# ============================================================
# First stage is methylation alone; joint comparisons are a separately logged
# next stage, not silently substituted. No filter/weight/K/criterion selection.
# The whole-universe imputation and screening are fixed before patient CV.
# Historical framework: raw beta input, fitted intercept, point-Laplace F,
# point-exponential L, EB Normal beta, original working curvature and caps.
# Fixed K7/no deletion is the approved sensitivity arm, not adopted parsimony.
# Serial local all; independent HPC cells: one CPU,8GB, no nested parallelism.
# Large cells may need Longleaf; do not submit jobs or silently downsample.
source("code/fit_multimodal_yfb.R")
source("code/predict_multimodal_yfb.R")
source("code/train_test_split.R")
source("code/feature_sensitivity.R")
source("code/cpg_feature_sensitivity.R")
source("code/multimodal_feature_sensitivity.R")
cfg <- yaml::read_yaml("config/globals.yml")
args <- commandArgs(TRUE);stopifnot(length(args)==1L)
mode <- args[1]
out <- "results/multimodal_real/outputs/feature_sensitivity/methylation_count"
local <- "data/round2_methylation_count_design.rds"
input <- "data/round2_tcga_feature_universe.rds"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
z <- readRDS(input)

# Freeze the design without fitting ----
if(mode=="design") {
  if(file.exists(local))stop("Design already exists; do not overwrite.")
  validate_cpg_reconstruction(z$methylation,z$baseline$Y$methylation)
  stopifnot(identical(z$manifest$design$actual,c(10000L,25000L,50000L,100000L,305352L)))
  folds <- create_stratified_folds(z$event,cfg$cavi$n_cv_folds,seed=42L)
  cells <- expand.grid(fold=0:cfg$cavi$n_cv_folds,set=seq_along(z$manifest$features))
  cells$cell <- seq_len(nrow(cells));cells$p <- z$manifest$design$actual[cells$set]
  control <- modifyList(multimodal_yfb_default_control(),list(intercept=TRUE,
    prior_update="ebnm",prune=FALSE,survival_elbo="plugin",tau_model="feature",
    prior_F=list(methylation="point_laplace"),init_supervised=FALSE,
    beta_prior="eb_normal",profile_cox_intercept=FALSE))
  design <- list(cells=cells,folds=folds,control=control,K=cfg$benchmark$k_merged_yfb_desurv)
  stopifnot(design$K==7L,nrow(cells)==30L,control$max_outer==cfg$multimodal_yfb$max_outer)
  saveRDS(design,local)
  write.csv(cells,file.path(out,"cells.csv"),row.names=FALSE)
  # Historical references share144 TCGA patients but use different molecular
  # inputs; these training rows provide context, not a paired methylation test.
  paths <- c(joint="results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned.rds",
    expression_only="results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned_expression-only.rds")
  references <- list()
  for(nm in names(paths)) {
    fit <- readRDS(paths[[nm]])$joint
    Y <- z$baseline$Y[fit$training_spec$modalities]
    pred <- predict_multimodal_yfb(fit,Y)
    score <- multimodal_feature_metrics(z$time,z$event,pred$risk_scores)
    references[[nm]] <- data.frame(source="benchmark_snapshot",model=nm,
      evaluation="training",p=if(nm=="joint")10000L else 0L,p_expression=3000L,
      C=score$C,logPL=score$logPL,risk_sd=score$risk_sd,
      credible=sum(fit$diagnostics$prognostic),retained=fit$diagnostics$K_eff_retained,
      converged=fit$diagnostics$converged,iterations=fit$diagnostics$iterations)
  }
  write.csv(do.call(rbind,references),file.path(out,"benchmark_rows.csv"),row.names=FALSE)
  paths <- c(input,local,"config/globals.yml","code/fit_multimodal_yfb.R",
    "code/multimodal_yfb_helpers.R","code/multimodal_yfb_updates.R",
    "code/multimodal_yfb_sweep.cpp","code/preprocess_multimodal_yfb.R",
    "code/predict_multimodal_yfb.R","code/profile_cox_intercept.R",
    "code/train_test_split.R","code/feature_sensitivity.R","code/cpg_feature_sensitivity.R",
    "code/multimodal_feature_sensitivity.R","results/multimodal_real/run_feature_expansion_multimodal.R")
  write.csv(data.frame(path=paths,md5=unname(tools::md5sum(paths))),file.path(out,"input_hashes.csv"),row.names=FALSE)
  writeLines(c("Count-only methylation sensitivity. Conditional CV after fixed imputation/screening.",
    "No mean/proportion/SD filter; those require a separate prospective comparison.",
    "No external score, model adoption, altered cap or constant-risk substitution.",
    capture.output(str(control)),capture.output(sessionInfo())),file.path(out,"design.txt"))
  quit(status=0)
}

# Evaluate immutable cells ----
design <- readRDS(local)
h <- read.csv(file.path(out,"input_hashes.csv"))
stopifnot(identical(unname(tools::md5sum(h$path)),h$md5))
if(mode=="all")ids <- design$cells$cell else {
  ids <- suppressWarnings(as.integer(mode))
  if(length(ids)!=1L || is.na(ids) || !ids%in%design$cells$cell)stop("Use design, all or cell1-30.")
}
for(id in ids) {
  cell <- design$cells[design$cells$cell==id,]
  prefix <- file.path(out,sprintf("cell_%02d",id))
  if(file.exists(paste0(prefix,".csv"))) {
    cat("Preserved terminal cell",id,"\n");next
  }
  features <- z$manifest$features[[cell$set]]
  Y <- z$methylation[,match(features,colnames(z$methylation)),drop=FALSE]
  test <- if(cell$fold==0L)seq_len(nrow(Y)) else design$folds$folds[[cell$fold]]
  train <- if(cell$fold==0L)test else setdiff(seq_len(nrow(Y)),test)
  row <- data.frame(source="methylation_count",model="methylation_only",cell=id,
    p=length(features),p_expression=0L,fold=cell$fold,
    evaluation=if(cell$fold==0L)"training" else "conditional_CV",
    C=NA_real_,logPL=NA_real_,risk_sd=NA_real_,credible=NA_integer_,retained=NA_integer_,
    converged=FALSE,iterations=NA_integer_,reconstruction_rmse=NA_real_,
    runtime_seconds=NA_real_,gc_peak_mb=NA_real_,warnings="",error="")
  started <- proc.time()[3];warnings <- character();gc(reset=TRUE)
  cat("Starting methylation cell",id,"p",length(features),"fold",cell$fold,"\n")
  tryCatch(withCallingHandlers({
    set.seed(42L)
    fit <- fit_multimodal_yfb(list(methylation=Y[train,,drop=FALSE]),
      z$time[train],z$event[train],design$K,design$control)
    pred <- predict_multimodal_yfb(fit,list(methylation=Y[test,,drop=FALSE]))
    score <- multimodal_feature_metrics(z$time[test],z$event[test],pred$risk_scores)
    row$C <- score$C;row$logPL <- score$logPL;row$risk_sd <- score$risk_sd
    sd_beta <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
    row$credible <- sum(sd_beta>0 & abs(fit$EBeta)>=cfg$multimodal_yfb$prognostic_z_threshold*sd_beta)
    row$retained <- fit$diagnostics$K_eff_retained
    row$converged <- fit$diagnostics$converged;row$iterations <- fit$diagnostics$iterations
    residual <- Y[train,,drop=FALSE]-sweep(fit$EL%*%t(fit$EF$methylation),2,fit$mu$methylation,"+")
    row$reconstruction_rmse <- sqrt(mean(residual^2))
    saveRDS(fit,paste0(prefix,".rds"))
    write.csv(data.frame(factor=seq_along(fit$EBeta),beta=fit$EBeta,sd_beta=sd_beta,
      beta_per_projection_sd=fit$EBeta*apply(Y[train,,drop=FALSE]%*%fit$EF$methylation,2,stats::sd),
      class=fit$diagnostics$factor_class,pve=fit$diagnostics$factor_pve),
      paste0(prefix,"_programs.csv"),row.names=FALSE)
  },warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e){row$error <<- conditionMessage(e)})
  row$warnings <- paste(warnings,collapse=" | ")
  row$runtime_seconds <- proc.time()[3]-started;row$gc_peak_mb <- feature_gc_peak_mb(gc())
  write.csv(row,paste0(prefix,".csv"),row.names=FALSE);print(row)
}
files <- list.files(out,"^cell_[0-9]+[.]csv$",full.names=TRUE)
rows <- c(lapply(files,feature_read_terminal),list(read.csv(file.path(out,"benchmark_rows.csv"))))
columns <- unique(unlist(lapply(rows,names)))
for(i in seq_along(rows))for(nm in setdiff(columns,names(rows[[i]])))rows[[i]][[nm]] <- NA
comparison <- do.call(rbind,lapply(rows,function(x)x[,columns,drop=FALSE]))
write.csv(comparison,file.path(out,"results_with_benchmark.csv"),row.names=FALSE)
if(any(nzchar(comparison$error[comparison$source=="methylation_count"])))
  stop("At least one fit failed; terminal diagnostics retained.")
