# ============================================================
# Script: run_posterior_k_panel.R
# Purpose: Test Andrew's approved posterior-state pruning refit remedy.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: unchanged calibrated panel definitions, parallel, yaml
# Run: Rscript results/benchmark_sim/run_posterior_k_panel.R [cores]
# ============================================================
# Task1-5 design, committed before execution. One change: the selected protected
# single-modality model carries its full retained posterior into pruning refits.
# All updates, priors, thresholds, calibrated data, seeds, K_init and stopping
# rules are unchanged. This extra remedy was explicitly approved by Andrew after
# the case trace showed loss during a refit following deletion of another program.
# Fit all160 single protected cells: 20seeds x2noises x4K. Reuse, never refit or
# overwrite, the160 reference and160 mean-only protected rows from the complete
# original800-cell Task1-3 panel. Require that source panel to be complete/error-free.
# Adoption: every K/noise recovery>=16/20, false<=reference+.5,20/20 ELBO-converged,
# medianK=3. Missing/error cells cannot pass. Report intermediate caps separately.
# No external data, criterion changes, tuning or benchmark replacement here.
# Six PSOCK workers, one serial fit each, chunk.size1 to avoid coarse scheduling
# tails. Budget4GB/worker; local runtime expected minutes, possibly>1hour if many
# refits cap. Longleaf array: round2_posterior_k_selection.slurm; never submit here.
source("results/benchmark_sim/run_k_selection_panel.R")

#' Run every fixed single protected cell with posterior-state refitting
#' @param cores Positive socket-worker count.
#' @param out_dir New output folder; refuses any existing attempt folder.
#' @return Invisibly summary table, with new cells and frozen source rows retained.
#' @examples
#' # run_posterior_k_panel(6L)
run_posterior_k_panel <- function(cores=6L,
  out_dir="results/benchmark_sim/outputs/k_selection/round2/attempt5_posterior") {
  if(dir.exists(out_dir))stop("Attempt folder exists; do not overwrite")
  source_path <- "results/benchmark_sim/outputs/k_selection/round2/attempt3/results.csv"
  old <- read.csv(source_path)
  required <- k_panel_grid()
  keys <- c("seed","noise","K","model","arm")
  actual <- merge(required,old,by=keys)
  if(nrow(old)!=800L || nrow(actual)!=800L || any(!is.na(old$error)) ||
    anyDuplicated(old[keys]))stop("Task1-3 must have all800 unique error-free cells first")
  if(length(cores)!=1L || !is.finite(cores) || cores<1L || cores!=as.integer(cores))
    stop("cores must be a positive integer")
  dir.create(file.path(out_dir,"cells"),recursive=TRUE)
  grid <- required
  ids <- which(grid$model=="single" & grid$arm=="protected")
  design <- grid[ids,];design$posterior_start <- TRUE
  write.csv(design,file.path(out_dir,"design.csv"),row.names=FALSE)
  writeLines(c(paste("source",source_path),paste("source_md5",tools::md5sum(source_path))),
    file.path(out_dir,"reference_source.txt"))
  cfg <- yaml::read_yaml("config/globals.yml")
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl,c("load_panel_code","panel_cell","grid","cfg","out_dir"),
    envir=environment())
  invisible(parallel::clusterEvalQ(cl,load_panel_code()))
  rows <- tryCatch(parallel::parLapplyLB(cl,ids,panel_cell,posterior_start=TRUE,chunk.size=1L),
    finally=parallel::stopCluster(cl))
  cols <- unique(unlist(lapply(rows,names)))
  fresh <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(fresh,file.path(out_dir,"results_new.csv"),row.names=FALSE)
  ref <- old[old$model=="single",]
  ref$arm[ref$arm=="protected"] <- "mean_only_protected"
  rows <- list(ref,fresh);cols <- unique(unlist(lapply(rows,names)))
  res <- do.call(rbind,lapply(rows,function(x) {x[setdiff(cols,names(x))]<-NA;x[cols]}))
  write.csv(res,file.path(out_dir,"results.csv"),row.names=FALSE)
  summ <- k_panel_summary(res)
  write.csv(summ,file.path(out_dir,"summary.csv"),row.names=FALSE)
  flags <- c("recovery_pass","false_positive_pass","convergence_pass","rank_pass")
  x <- summ[summ$arm=="protected",]
  pass <- nrow(x)==8L && all(!is.na(as.matrix(x[flags])) & as.matrix(x[flags]))
  failed <- vapply(flags,function(f)any(is.na(x[[f]]) | !x[[f]]),logical(1))
  adoption <- data.frame(model="single",adopted=pass,
    initialization="full_posterior",reason=if(pass)"all_simulation_criteria" else
      paste(flags[failed],collapse=";"))
  write.csv(adoption,file.path(out_dir,"adoption.csv"),row.names=FALSE)
  print(summ,row.names=FALSE);print(adoption,row.names=FALSE);invisible(summ)
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_posterior_k_panel(if(length(args))as.integer(args[1]) else 6L)
}
