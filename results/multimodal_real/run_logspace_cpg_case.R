# ============================================================
# Script: run_logspace_cpg_case.R
# Purpose: Test stable Cox arithmetic on the unchanged 25k CpG case without tuning.
# Created: 2026-10-15
# Dependencies: survival, existing multimodal fitter/predictor
# Run only after Task5-M3 prospective commit: Rscript results/multimodal_real/run_logspace_cpg_case.R
# ============================================================
# One analytical change: control$cox_arithmetic="logspace". Same cell7,
# original seed42,144 patients,K7,25k raw-beta CpGs, imputation/ranking,
# priors, intercept, precision,500/5 caps and convergence/classification rules.
# No damping, clipping, survival weight, altered initialization or external
# score. This case tests numeric usability, not complete scientific adoption.
source("code/fit_multimodal_yfb.R")
source("code/predict_multimodal_yfb.R")
source("code/multimodal_moment_audit.R")
out<-"results/multimodal_real/outputs/feature_sensitivity/logspace_case"
if(file.exists(file.path(out,"case.csv")))stop("Preserve the terminal arithmetic case.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base<-"results/multimodal_real/outputs/feature_sensitivity/methylation_count"
h<-read.csv(file.path(base,"input_hashes.csv"))
# Only these two sources have an opt-in dispatch added; their original branch
# is independently checked against the protected base by the test suite.
changed_sources<-c("code/fit_multimodal_yfb.R","code/multimodal_yfb_helpers.R")
unchanged<-h[!h$path%in%changed_sources,]
stopifnot(identical(unname(tools::md5sum(unchanged$path)),unchanged$md5))
paths<-c(h$path,"code/cox_working_logspace.R","code/multimodal_moment_audit.R",
 "data/round2_methylation_moment_failure.rds","results/multimodal_real/run_logspace_cpg_case.R")
hashes<-tools::md5sum(paths)
z<-readRDS("data/round2_tcga_feature_universe.rds")
design<-readRDS("data/round2_methylation_count_design.rds")
cell<-design$cells[design$cells$cell==7L,];stopifnot(cell$p==25000L,cell$fold==0L)
features<-z$manifest$features[[cell$set]]
Y<-z$methylation[,match(features,colnames(z$methylation)),drop=FALSE]
control<-design$control;control$cox_arithmetic<-"logspace"
cfg<-yaml::read_yaml("config/globals.yml")

# Replay the previously captured working point, without fitting ----
saved<-readRDS("data/round2_methylation_moment_failure.rds")
eta<-saved$state$predictor$mean
old<-multimodal_yfb_cox_working(z$time,z$event,eta)
new<-multimodal_yfb_cox_working(z$time,z$event,eta,arithmetic="logspace")
replay<-data.frame(arithmetic=c("original","logspace"),
 working_nonfinite=c(sum(!is.finite(c(old$w,old$h))),sum(!is.finite(c(new$w,new$h)))),
 logPL=c(old$log_partial_likelihood,new$log_partial_likelihood),
 total_working_weight=c(sum(old$w),sum(new$w)),events=sum(z$event))
write.csv(replay,file.path(out,"captured_point_replay.csv"),row.names=FALSE)
stopifnot(replay$working_nonfinite[1]>0,replay$working_nonfinite[2]==0,
 abs(sum(new$w)-sum(z$event))<1e-8)

# Observe progress without changing any fit argument ----
progress<-new.env(parent=emptyenv());progress$last_outer<-0L
#' Record aggregate predictors at each outer sweep's first working evaluation
#' @param call_env Environment of the working-function entry.
#' @param frames Active frames used only to locate the current fit indices.
#' @return Invisible NULL; appends aggregate progress, no argument modifications.
#' @examples
#' # Invoked by the process-local trace below.
observe_logspace_case <- function(call_env,frames) {
 ids<-which(vapply(frames,function(e)exists("outer",e,inherits=FALSE) &&
   exists("settings",e,inherits=FALSE),logical(1)))
 if(length(ids)) {
  f<-frames[[tail(ids,1L)]];outer<-get("outer",f)
  if(outer>progress$last_outer) {
   x<-get("eta",call_env);row<-data.frame(outer=outer,eta_min=min(x),eta_max=max(x),
     nonfinite=sum(!is.finite(x)))
   file<-file.path(out,"progress.csv")
   write.table(row,file,sep=",",row.names=FALSE,col.names=!file.exists(file),append=file.exists(file))
   cat("Outer",outer,"eta range",range(x),"\n");flush.console()
   progress$last_outer<-outer
  }
 }
 invisible(NULL)
}
trace("multimodal_yfb_cox_working",tracer=quote({
 observe_logspace_case(environment(),sys.frames())
}),print=FALSE)
row<-data.frame(source="logspace_25k",model="methylation_only",p=25000L,
 evaluation="training",C=NA_real_,logPL=NA_real_,risk_sd=NA_real_,credible=NA_integer_,
 retained=NA_integer_,converged=FALSE,iterations=NA_integer_,runtime_seconds=NA_real_,
 warnings="",error="")
warning_messages<-character();started<-proc.time()[3];set.seed(42L)
tryCatch(withCallingHandlers({
 fit<-fit_multimodal_yfb(list(methylation=Y),z$time,z$event,design$K,control)
 saveRDS(fit,file.path(out,"fit.rds"))
 risk<-predict_multimodal_yfb(fit,list(methylation=Y))$risk_scores
 stopifnot(all(is.finite(risk)))
 row$risk_sd<-sd(risk)
 row$C<-if(sd(risk)==0)NA_real_ else survival::concordance(survival::Surv(z$time,z$event)~risk,reverse=TRUE)$concordance
 row$logPL<-multimodal_yfb_cox_working(z$time,z$event,risk,arithmetic="logspace")$log_partial_likelihood
 row$converged<-fit$diagnostics$converged;row$iterations<-fit$diagnostics$iterations
 bsd<-sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
 row$credible<-sum(bsd>0 & abs(fit$EBeta)>=cfg$multimodal_yfb$prognostic_z_threshold*bsd)
 row$retained<-fit$diagnostics$K_eff_retained
},warning=function(w)warning_messages<<-c(warning_messages,conditionMessage(w))),
error=function(e)row$error<<-conditionMessage(e),finally=untrace("multimodal_yfb_cox_working"))
row$runtime_seconds<-unname(proc.time()[3]-started);row$warnings<-paste(unique(warning_messages),collapse=" | ")
write.csv(row,file.path(out,"case.csv"),row.names=FALSE)
reference<-read.csv(file.path(base,"results_with_benchmark.csv"))
cols<-union(names(reference),names(row))
for(nm in setdiff(cols,names(reference)))reference[[nm]]<-NA
for(nm in setdiff(cols,names(row)))row[[nm]]<-NA
write.csv(rbind(reference[,cols],row[,cols]),file.path(out,"results_with_benchmark.csv"),row.names=FALSE)
stopifnot(identical(hashes,tools::md5sum(paths)))
write.csv(data.frame(path=paths,md5=unname(hashes)),file.path(out,"input_hashes.csv"),row.names=FALSE)
print(row)
if(nzchar(row$error))stop("Arithmetic case failed; terminal diagnostics preserved.")
