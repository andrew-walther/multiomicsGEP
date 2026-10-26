# ============================================================
# Script: prepare_feature_universe.R
# Purpose: Reconstruct the original training CpG universe and verify the frozen screen before expansion.
# Created: 2026-10-15
# Dependencies: yaml, impute, existing multiomics loader and feature helpers
# Run only after prospective commit: Rscript results/multimodal_real/prepare_feature_universe.R
# ============================================================
# Data preparation only: no model, CV scoring, feature/default adoption or
# external outcomes. The original shared-probe rule uses validation probe IDs
# only. No validation expression, survival, clinical data or scores are read.
# Imputation follows the existing helper's exact submatrix/default arguments,
# while exposing warnings instead of suppressing them. Fit inputs must match
# frozen values bit-for-bit. No fallback/tolerance or re-imputation of 10k alone.
# Serial; raw input loading may peak at several GB. Longleaf: one CPU,8GB.
# Local full matrices/patient manifests stay ignored; only aggregate audit is safe.
source("code/load_multiomics_data.R")
source("code/cpg_feature_sensitivity.R")
source("results/multimodal_real/summarize_cpg_variability.R")
cfg <- yaml::read_yaml("config/globals.yml")
out <- "results/multimodal_real/outputs/feature_sensitivity/universe_preflight"
local <- "data/round2_tcga_feature_universe.rds"
if(file.exists(local))stop("Prepared universe already exists; do not overwrite.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
warnings <- character();started <- proc.time()[3]
tryCatch(withCallingHandlers({
  root <- cfg$multiomics_data$root
  meth_path <- file.path(root,"tcga","tcga_meth_filtered_sex.rds")
  probe_path <- file.path(root,"icgc","icgc_meth_filtered_sex.rds")
  expr_path <- file.path(root,"tcga","expr_log_normalized.rds")
  surv_path <- file.path(cfg$multiomics_data$pdac_root,"original","TCGA_PAAD.survival_data.rds")
  meth <- readRDS(meth_path)$meth
  validation_probe_ids <- rownames(readRDS(probe_path)$meth)
  raw_probe_count <- nrow(meth)
  cpgs <- intersect(rownames(meth),validation_probe_ids)
  rm(validation_probe_ids);gc()
  expr <- readRDS(expr_path)
  surv <- readRDS(surv_path)
  cache <- readRDS("data/round2_tcga_training.rds")
  stopifnot(identical(sort(colnames(expr)),sort(colnames(meth))),
    all(colnames(expr)%in%surv$sampID))
  s <- surv[match(colnames(expr),surv$sampID),]
  usable <- filter_usable_survival(colnames(expr),s$time,s$event,"TCGA")
  ids <- colnames(expr)[usable$keep]
  stopifnot(identical(ids,rownames(cache$Y$methylation)),
    identical(usable$time,cache$time),identical(usable$event,cache$event))
  meth <- meth[cpgs,ids,drop=FALSE]
  stopifnot(!any(is.infinite(meth)),all(meth[!is.na(meth)]>=0),all(meth[!is.na(meth)]<=1))
  observed_summary <- cpg_variability_table(meth)
  filter_counts <- summarize_cpg_filters(observed_summary,colnames(cache$Y$methylation))
  missing <- sum(is.na(meth));idx <- which(rowSums(is.na(meth))>0)
  if(length(idx)) {
    block <- meth[idx,,drop=FALSE]
    # Same input and default arguments as impute_methylation_knn(). The outer
    # warning observer retains diagnostics; no warning suppression is used.
    imputed <- impute::impute.knn(block)$data
    observed <- !is.na(block)
    stopifnot(identical(as.numeric(imputed[observed]),as.numeric(block[observed])))
    meth[idx,] <- imputed
    rm(block,imputed,observed);gc()
  }
  Y <- t(meth);rm(meth);gc()
  validate_cpg_reconstruction(Y,cache$Y$methylation)
  stopifnot(identical(unname(t(expr[colnames(cache$Y$expression),ids,drop=FALSE])),
    unname(cache$Y$expression)))
  eligible <- sum(apply(Y,2,stats::var)>0)
  manifest <- cpg_count_manifest(Y,c(ncol(cache$Y$methylation),25000L,50000L,100000L,eligible),
    colnames(cache$Y$methylation))
  stopifnot(nrow(Y)==144L,ncol(Y)==305352L,all(is.finite(Y)),
    all(is.finite(expr[,ids,drop=FALSE])))
  saveRDS(list(methylation=Y,expression_tcga_raw=t(expr[,ids,drop=FALSE]),
    baseline=cache,manifest=manifest,observed_summary=observed_summary,
    time=usable$time,event=usable$event),local)
  write.csv(manifest$design,file.path(out,"counts.csv"),row.names=FALSE)
  write.csv(filter_counts,file.path(out,"filter_counts.csv"),row.names=FALSE)
  write.csv(data.frame(raw_TCGA_probes=raw_probe_count,eligible_shared_probes=ncol(Y),
    excluded_by_original_shared_rule=raw_probe_count-ncol(Y),n=nrow(Y),
    missing_before=missing,imputed_probes=length(idx),missing_after=sum(is.na(Y)),
    frozen_values_identical=TRUE,frozen_rank_identical=TRUE,expression_identical=TRUE,
    runtime_seconds=proc.time()[3]-started),file.path(out,"validation.csv"),row.names=FALSE)
  paths <- c(meth_path,probe_path,expr_path,surv_path,"data/round2_tcga_training.rds",
    "config/globals.yml","code/load_multiomics_data.R","code/cpg_feature_sensitivity.R",
    "results/multimodal_real/summarize_cpg_variability.R",
    "results/multimodal_real/prepare_feature_universe.R",local)
  write.csv(data.frame(path=paths,md5=unname(tools::md5sum(paths))),file.path(out,"input_hashes.csv"),row.names=FALSE)
  writeLines(c("PASS: no fit or external outcome/score read; original shared eligibility retained.",
    "Original KNN missing-row submatrix/default arguments; warnings retained.",
    "Frozen methylation values, order, variance rank and expression subset match exactly.",
    capture.output(sessionInfo())),file.path(out,"status.txt"))
},warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
error=function(e){
  write.csv(data.frame(warning=warnings),file.path(out,"warnings.csv"),row.names=FALSE)
  writeLines(conditionMessage(e),file.path(out,"error.txt"));stop(e)
})
write.csv(data.frame(warning=warnings),file.path(out,"warnings.csv"),row.names=FALSE)
cat("Universe verified without fitting; warnings:",length(warnings),"\n")
