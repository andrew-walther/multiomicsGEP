# ============================================================
# Script: run_working_cycle_audit.R
# Purpose: Observe one unchanged registered Single simulation working cycle trajectory.
# Created: 2026-10-16
# Dependencies: existing panel modules and working_cycle_audit_helpers.R
# ============================================================
# Task2-6/Task2-7: run ONLY after the corresponding prospective code/log commit. Andrew approved
# the complete-cycle observer on2026-10-16. This is one unchanged cell021
# replay, NOT a remedy or adoption panel. No external/real data are loaded.
# See the attempt log for fixed arms, original criteria and stopping rule.
# Serial local; ~minutes,1GB. No nested workers or HPC submission.
# Detailed telemetry remains ignored in data/; commit compact aggregates only.
# Corrected Task2-7 run: Rscript results/benchmark_sim/run_working_cycle_audit.R corrected
source("results/benchmark_sim/run_round2_fixed_prior.R")
source("results/benchmark_sim/working_cycle_audit_helpers.R")

#' Select the fixed paths for one of the two registered diagnostic attempts
#' @param mode Empty character vector for historical Task2-6, or "corrected" for Task2-7.
#' @return Named out/detailed paths; any other mode fails explicitly.
#' @examples
#' cycle_audit_paths("corrected")
cycle_audit_paths <- function(mode=character()) {
  if(!identical(mode,character()) && !identical(mode,"corrected"))
    stop("Only the original or registered corrected attempt is permitted")
  suffix<-if(length(mode))"_corrected" else ""
  list(out=paste0("results/benchmark_sim/outputs/round2/working_cycle_seed26",suffix),
    detailed=paste0("data/round2_working_cycle_seed26",suffix,".rds"))
}

if(sys.nframe()==0L) {
  paths_for_attempt<-cycle_audit_paths(commandArgs(trailingOnly=TRUE))
  out<-paths_for_attempt$out;detailed<-paths_for_attempt$detailed
  if(dir.exists(out) || file.exists(detailed))stop("Fresh audit output required")
  source_fit<-"results/benchmark_sim/outputs/round2/profiled/cells/cell_021.rds"
  paths<-c(source_fit,"code/fit_cox_on_yf.R","code/multimodal_yfb_updates.R",
    "code/update_L_surv_YFB.R","code/update_beta.R","code/update_tau.R",
    "code/profile_cox_intercept.R","config/globals.yml",
    "results/benchmark_sim/working_cycle_audit_helpers.R",
    "results/benchmark_sim/run_working_cycle_audit.R")
  hashes<-tools::md5sum(paths)
  if(anyNA(hashes))stop("Required audit source missing")
  reference<-readRDS(source_fit)
  dir.create(out,recursive=TRUE)
  dir.create(dirname(detailed),showWarnings=FALSE)
  write.csv(data.frame(path=paths,md5=unname(hashes)),file.path(out,"input_hashes.csv"),row.names=FALSE)
  load_panel_code()
  tb<-make_testbed(20261026L,.25)
  controls<-fixed_prior_controls("single","profiled")
  audit<-cycle_observer();warnings<-character();started<-proc.time()[["elapsed"]]
  set.seed(20261026L)
  fit<-tryCatch(withCallingHandlers(audit$run(function()
    do.call(fit_cox_on_yf,c(list(tb$Ytr,tb$ttr,tb$str),controls))),
    warning=function(w)warnings<<-c(warnings,conditionMessage(w))),
    error=function(e){writeLines(conditionMessage(e),file.path(out,"error.txt"));stop(e)},
    finally={
      saveRDS(audit$rows,detailed)
      write.csv(data.frame(message=warnings),file.path(out,"warnings.csv"),row.names=FALSE)
    })
  saveRDS(fit,file.path(out,"fit.rds"))
  identity<-data.frame(field=names(reference),identical=vapply(names(reference),
    function(nm)identical(fit[[nm]],reference[[nm]]),logical(1)))
  write.csv(identity,file.path(out,"identity.csv"),row.names=FALSE)
  if(!identical(fit,reference))stop("Entire-fit identity failed; interpretation blocked")
  if(!identical(hashes,tools::md5sum(paths)))stop("Audit inputs changed during execution")
  d<-do.call(rbind,audit$rows)
  # Adjacent finite differences are diagnostic; infinite support transitions
  # remain explicit. They cannot be converted into a finite ascent statement.
  previous<-c(NA_real_,head(d$reference_objective,-1L))
  d$finite_pair<-is.finite(previous)&is.finite(d$reference_objective)
  d$delta<-NA_real_;d$delta[d$finite_pair]<-d$reference_objective[d$finite_pair]-previous[d$finite_pair]
  d$decrease<-d$finite_pair & !is.na(d$delta) & d$delta < -1e-8
  saveRDS(d,detailed)
  # Keep one row per outer iteration and a compact transition summary in git.
  terminal<-d[d$event=="canonical",]
  write.csv(terminal,file.path(out,"iteration_summary.csv"),row.names=FALSE)
  groups<-split(d,interaction(d$event,d$iteration>=251,drop=TRUE))
  summary<-do.call(rbind,lapply(groups,function(x)data.frame(event=x$event[1],
    last50=x$iteration[1]>=251,rows=nrow(x),complete=sum(x$complete),
    nonfinite=sum(!is.finite(x$reference_objective)),finite_pairs=sum(x$finite_pair),
    decreases=sum(x$decrease),min_delta=if(any(x$finite_pair))min(x$delta[x$finite_pair]) else NA_real_,
    max_delta=if(any(x$finite_pair))max(x$delta[x$finite_pair]) else NA_real_)))
  write.csv(summary,file.path(out,"transition_summary.csv"),row.names=FALSE)
  result<-data.frame(seed=20261026L,noise=.25,iterations=fit$history$n_iter,
    converged=fit$history$converged,exact_fit_identity=identical(fit,reference),
    rows=nrow(d),last50_complete=all(d$complete[d$iteration>=251]),
    max_moment_gap=max(d$moment_gap,na.rm=TRUE),
    implemented_common_prior=FALSE,runtime_seconds=proc.time()[["elapsed"]]-started)
  write.csv(result,file.path(out,"summary.csv"),row.names=FALSE)
  writeLines(c("One unchanged simulation replay; no remedy or scientific adoption.",
    "L/beta reference priors are ephemeral, transformed for density bookkeeping only.",
    "The fitter does not carry these priors through canonicalization; no common implemented objective is certified.",
    "Infinite KL transitions remain nonfinite; adjacent finite differences are descriptive.",
    "The anchored expected quadratic is not exact E_q logPL or a global lower bound.",
    "Original300 cap and convergence failure remain unchanged. Stop and regroup."),
    file.path(out,"interpretation.txt"))
  if(!result$last50_complete)stop("Incomplete late-cycle density coverage; interpretation blocked")
  print(result)
}
