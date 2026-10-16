# ============================================================
# Script: run_feature_expansion.R
# Purpose: Evaluate nested expression feature counts at fixed K using benchmark settings.
# Created: 2026-10-15
# Dependencies: yaml, survival, existing YFB fitter and preprocessing
# Run from repository root: Rscript results/benchmark_sim/run_feature_expansion.R prepare|all|CELL
# ============================================================
# Training only; no external loading or selection. Screening/standardization
# precede patient folds, so CV is explicitly conditional on those fixed inputs.
# Cold starts at every count are authorized after the frozen warm inverse
# proved singular. No pruning, supervised initialization or weight selection.
# Longleaf: independent CELL processes need one CPU and approximately 4GB;
# prepare once, then array1-24. Do not share output cells between live workers.
# Local all is serial; no forked workers or threaded parallel fitting assumed.

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
mode <- args[1]
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")
source("code/train_test_split.R")
# The legacy source defines functions before its demonstration runner stops.
# Accept only that known NULL-input assertion, never suppress another error.
tryCatch(source("code/fit_cox_on_yf.R"),error=function(e) {
  if(conditionMessage(e)!="!is.null(real_Y) is not TRUE")stop(e)
})
source("code/predict_cox_on_yf.R")
source("code/select_K.R")
source("code/feature_sensitivity.R")
stopifnot(exists("fit_cox_on_yf"),exists("calc_cox_taylor_yf"))
out <- "results/benchmark_sim/outputs/feature_expansion/cold_primary"
local <- "data/round2_feature_cold_primary.rds"
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# Freeze inputs and cell design without fitting ----
if(mode=="prepare") {
  if(file.exists(local))stop("Prepared input already exists; do not overwrite a frozen design.")
  train <- cfg$pdac$training_cohorts
  raw <- lapply(setNames(train,train),load_pdac_raw,pdac_root=PDAC_DATA_ROOT)
  common <- Reduce(intersect,lapply(raw,function(z)z$gene_names))
  ranked <- lapply(train,function(ds) {
    x <- raw[[ds]]$Y[,match(common,raw[[ds]]$gene_names),drop=FALSE]
    if(PLATFORM_LOG_TRANSFORM[[ds]])x <- log2_plus1_transform(x)
    colnames(x) <- common;x
  });names(ranked) <- train
  frozen_genes <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/d4_gene_names.rds")
  manifest <- feature_count_manifest(ranked,common,
    c(length(frozen_genes),5000L,10000L,length(common)),cfg$preprocessing$top_n_genes_desurv)
  stopifnot(identical(manifest$genes[[1]],frozen_genes),
    identical(manifest$design$actual_shared,c(2064L,5001L,10000L,17900L)))
  all <- preprocess_merged_cohorts(raw,PLATFORM_LOG_TRANSFORM[train],top_n=length(common),
    rank_transform=FALSE,per_platform_standardize=TRUE,normalize_method="none",
    selection_per_cohort=TRUE,selection_method="combined_rank")
  baseline <- preprocess_merged_cohorts(raw,PLATFORM_LOG_TRANSFORM[train],
    top_n=cfg$preprocessing$top_n_genes_desurv,rank_transform=FALSE,
    per_platform_standardize=TRUE,normalize_method="none",selection_per_cohort=TRUE,
    selection_method="combined_rank")
  stopifnot(identical(baseline$gene_names,frozen_genes),
    identical(unname(baseline$Y),unname(all$Y[,match(frozen_genes,all$gene_names)])),
    all(is.finite(all$Y)))
  time <- unlist(lapply(train,function(ds)raw[[ds]]$time))
  status <- unlist(lapply(train,function(ds)raw[[ds]]$status))
  folds <- create_stratified_folds(status,cfg$cavi$n_cv_folds,seed=42L)
  cells <- expand.grid(fold=0:cfg$cavi$n_cv_folds,set=seq_along(manifest$genes))
  cells$cell <- seq_len(nrow(cells));cells$p <- manifest$design$actual_shared[cells$set]
  controls <- list(K=cfg$benchmark$k_merged_yfb_desurv,max_iter=cfg$cavi$max_iter,
    tol=cfg$cavi$tol,alpha=cfg$benchmark$alpha,alpha_F=0,
    prior_LF="point_exponential",prior_beta="normal",prune=FALSE,
    canonicalize=FALSE,init_supervised=FALSE,verbose=FALSE)
  stopifnot(controls$K==7L,nrow(all$Y)==273L,nrow(cells)==24L)
  saveRDS(list(Y=all$Y,genes=all$gene_names,manifest=manifest,time=time,status=status,
    folds=folds,cells=cells,controls=controls),local)
  write.csv(manifest$design,file.path(out,"feature_counts.csv"),row.names=FALSE)
  write.csv(cells,file.path(out,"cells.csv"),row.names=FALSE)
  tracked <- c("config/globals.yml","code/fit_cox_on_yf.R","code/feature_sensitivity.R",
    "code/preprocess_desurv.R","code/train_test_split.R",
    "results/benchmark_sim/run_feature_expansion.R",local)
  write.csv(data.frame(path=tracked,md5=unname(tools::md5sum(tracked))),
    file.path(out,"input_hashes.csv"),row.names=FALSE)
  frozen <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_fits.rds")$D4
  fp <- predict_cox_on_yf(baseline$Y,frozen$EF,frozen$EBeta,frozen$EF_norms)
  # Historical D4 predates the sign fix: frozen risk=-eta (DECISIONS 2026-10-02).
  fm <- feature_survival_metrics(time,status,-fp$risk_scores)
  frozen_rows <- data.frame(source="benchmark_snapshot",p=length(frozen_genes),K=7,
    fold=0,evaluation="training",C=fm$C,logPL=fm$logPL,risk_sd=fm$risk_sd,
    converged=frozen$history$converged,iterations=frozen$history$n_iter)
  write.csv(frozen_rows,file.path(out,"benchmark_rows.csv"),row.names=FALSE)
  writeLines(c("Conditional training CV: screening/scaling fixed before folds.",
    "Local RDS contains ordered patients/genes/folds; never commit it.",
    "Existing zero-variance columns retained; no feature filter changed.",
    capture.output(str(controls)),capture.output(sessionInfo())),file.path(out,"design.txt"))
  quit(status=0)
}

# Fit each immutable cell and preserve every error/warning ----
z <- readRDS(local)
if(mode=="all")ids <- z$cells$cell else {
  ids <- suppressWarnings(as.integer(mode))
  if(length(ids)!=1L || is.na(ids) || !ids%in%z$cells$cell)stop("Use prepare, all or cell1-24.")
}
hashes <- read.csv(file.path(out,"input_hashes.csv"))
stopifnot(identical(unname(tools::md5sum(hashes$path)),hashes$md5))
for(id in ids) {
  cell <- z$cells[z$cells$cell==id,]
  prefix <- file.path(out,sprintf("cell_%02d",id))
  if(file.exists(paste0(prefix,".csv"))) {
    cat("Preserved existing terminal cell",id,"\n");next
  }
  genes <- z$manifest$genes[[cell$set]]
  Y <- z$Y[,match(genes,z$genes),drop=FALSE]
  test <- if(cell$fold==0L)seq_len(nrow(Y)) else z$folds$folds[[cell$fold]]
  train <- if(cell$fold==0L)test else setdiff(seq_len(nrow(Y)),test)
  warnings <- character();started <- proc.time()[3];gc(reset=TRUE)
  row <- data.frame(source="feature_sensitivity",cell=id,p=length(genes),K=z$controls$K,
    fold=cell$fold,evaluation=if(cell$fold==0L)"training" else "conditional_CV",
    C=NA_real_,logPL=NA_real_,risk_sd=NA_real_,active=NA_integer_,retained=NA_integer_,
    converged=FALSE,iterations=NA_integer_,constant_columns=sum(apply(Y,2,stats::sd)==0),
    runtime_seconds=NA_real_,gc_peak_mb=NA_real_,warnings="",error="")
  cat("Starting cell",id,"p",length(genes),"fold",cell$fold,"\n")
  tryCatch(withCallingHandlers({
    set.seed(42L)
    controls <- z$controls;controls$sign_correction <- cell$fold==0L
    fit <- do.call(fit_cox_on_yf,c(list(Y=Y[train,,drop=FALSE],
      time=z$time[train],status=z$status[train]),controls))
    training <- predict_cox_on_yf(Y[train,,drop=FALSE],fit$EF,fit$EBeta,fit$EF_norms)
    # Match existing CV orientation using training outcomes only.
    orientation <- feature_survival_metrics(z$time[train],z$status[train],training$risk_scores)$C
    if(cell$fold!=0L && is.finite(orientation) && orientation<.5)fit$EBeta <- -fit$EBeta
    pred <- predict_cox_on_yf(Y[test,,drop=FALSE],fit$EF,fit$EBeta,fit$EF_norms)
    metrics <- feature_survival_metrics(z$time[test],z$status[test],pred$risk_scores)
    classes <- classify_factors(fit,Y[train,,drop=FALSE],
      beta_thresh=cfg$k_selection$beta_threshold,pve_thresh=cfg$k_selection$pve_threshold)
    row$C <- metrics$C;row$logPL <- metrics$logPL;row$risk_sd <- metrics$risk_sd
    row$active <- sum(classes$surv_active);row$retained <- sum(classes$category!="dead")
    row$converged <- fit$history$converged;row$iterations <- fit$history$n_iter
    classes$EBeta_by_cohort <- NULL
    write.csv(classes,paste0(prefix,"_programs.csv"),row.names=FALSE)
    saveRDS(list(fit=fit,genes=genes,train=train,test=test),paste0(prefix,".rds"))
  },warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e){row$error <<- conditionMessage(e)})
  row$warnings <- paste(unique(warnings),collapse=" | ")
  row$runtime_seconds <- proc.time()[3]-started
  # R allocation high-water mark, not operating-system resident memory.
  memory <- gc();row$gc_peak_mb <- sum(memory[,6])
  write.csv(row,paste0(prefix,".csv"),row.names=FALSE)
  print(row)
}

# Publish all terminal cells alongside the labelled frozen training reference.
# Training reference C is not a CV estimate; no external cohort is scored.
files <- list.files(out,pattern="^cell_[0-9]+[.]csv$",full.names=TRUE)
rows <- c(lapply(files,read.csv,stringsAsFactors=FALSE),
  list(read.csv(file.path(out,"benchmark_rows.csv"),stringsAsFactors=FALSE)))
columns <- unique(unlist(lapply(rows,names)))
for(i in seq_along(rows))for(nm in setdiff(columns,names(rows[[i]])))rows[[i]][[nm]] <- NA
comparison <- do.call(rbind,lapply(rows,function(x)x[,columns,drop=FALSE]))
write.csv(comparison,file.path(out,"results_with_benchmark.csv"),row.names=FALSE)
if(any(nzchar(comparison$error[comparison$source=="feature_sensitivity"])))
  stop("At least one cell failed; all terminal rows and diagnostics were preserved.")
