# ============================================================
# Script: audit_loading_balance.R
# Purpose: Describe loading-update likelihood balance in saved TCGA fits without refitting.
# Created: 2026-10-15
# Dependencies: base R, existing multimodal helpers, yaml via existing fitter
# Run only after Task2-B1 prospective commit: Rscript results/multimodal_real/audit_loading_balance.R
# ============================================================
# Four saved terminal states; all capped. These are isolated coordinate
# conditionals at one frozen state, not sequential sweep telemetry, causal
# evidence for collapse, new prediction scores or a selected survival weight.
# No global objective ratio substitutes for per-coordinate information shares.
source("code/fit_multimodal_yfb.R")
source("code/multimodal_program_agreement.R")
source("code/loading_balance_audit.R")
out<-"results/multimodal_real/outputs/round2/loading_balance"
if(file.exists(file.path(out,"likelihood_shares.csv")))stop("Preserve the completed audit.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths<-c(benchmark_joint="results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned.rds",
 benchmark_expression="results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned_expression-only.rds",
 profiled_joint="results/benchmark_sim/outputs/round2/profiled/cells/cell_102.rds",
 methylation_10k="results/multimodal_real/outputs/feature_sensitivity/methylation_count/cell_01.rds")
inputs<-c(unname(paths),"data/round2_tcga_training.rds","config/globals.yml",
 "code/fit_multimodal_yfb.R","code/multimodal_yfb_helpers.R","code/multimodal_yfb_updates.R",
 "code/profile_cox_intercept.R","code/multimodal_program_agreement.R","code/loading_balance_audit.R")
hashes<-tools::md5sum(inputs)
training<-readRDS("data/round2_tcga_training.rds")
rows<-components<-list()
for(nm in names(paths)) {
 fit<-readRDS(paths[[nm]])
 if(grepl("^benchmark",nm))fit<-fit$joint
 a<-multimodal_feature_adapter(fit,training$Y)
 stopifnot(!fit$diagnostics$converged,fit$diagnostics$iterations==500L)
 Y<-Map(function(x,m)sweep(x,2,m),a$blocks,fit$mu[names(a$blocks)])
 raw<-multimodal_yfb_projection_moments(Y,fit$EF,fit$EF2)
 predictor<-multimodal_yfb_predictor_moments(raw$EZ,raw$VZ,fit$EBeta,fit$EBeta2)
 working<-multimodal_yfb_cox_working(training$time,training$event,predictor$mean)
 profile<-isTRUE(fit$diagnostics$controls$profile_cox_intercept)
 survival_Y<-if(profile)lapply(Y,profile_cox_center,w=working$w) else Y
 if(profile)working$h<-working$u+working$w*profile_cox_center(predictor$mean,working$w)
 projection<-multimodal_yfb_projection_moments(survival_Y,fit$EF,fit$EF2)
 current_eta<-multimodal_yfb_predictor_moments(projection$EZ,projection$VZ,fit$EBeta,fit$EBeta2)$mean
 for(k in seq_along(fit$EBeta)) {
  residual<-multimodal_yfb_residual_minus_k(Y,fit$EL,fit$EF,k)
  h_minus<-working$h-working$w*(current_eta-fit$EBeta[k]*projection$EZ[,k])
  for(m in names(Y)) {
   terms<-loading_balance_terms(Y[[m]],survival_Y[[m]],fit$Tau[[m]],
     fit$EL[,k],fit$EL2[,k],residual[[m]],projection$EZ[,k],fit$EF[[m]][,k],
     fit$EBeta[k],fit$EBeta2[k],working$w,h_minus)
   row<-data.frame(source=nm,converged=FALSE,iterations=500L,profiled=profile,
     factor=k,modality=m,features=ncol(Y[[m]]),beta=fit$EBeta[k],beta_second=fit$EBeta2[k],
     beta_point_mass_zero=fit$EBeta[k]==0 && fit$EBeta2[k]==0,
     undefined_precision=sum(is.na(terms$precision_share)),undefined_signal=sum(is.na(terms$signal_share)))
   for(field in names(terms)) {
    x<-terms[[field]];x<-x[is.finite(x)]
    quantiles<-if(length(x))stats::quantile(x,c(.5,.9,1),names=FALSE) else rep(NA_real_,3)
    for(j in 1:3)row[[paste(field,c("median","p90","max")[j],sep="_")]]<-quantiles[j]
   }
   rows[[length(rows)+1L]]<-row
  }
 }
 saved<-attr(fit$diagnostics$elbo,"components")
 stopifnot(all(c("reconstruction","survival","kl")%in%names(saved)))
 components[[nm]]<-data.frame(source=nm,converged=FALSE,iterations=500L,
   saved_reconstruction=unname(saved["reconstruction"]),saved_survival=unname(saved["survival"]),
   saved_kl=unname(saved["kl"]),saved_removal_monitor=as.numeric(fit$diagnostics$elbo))
}
stopifnot(identical(hashes,tools::md5sum(inputs)))
write.csv(do.call(rbind,rows),file.path(out,"likelihood_shares.csv"),row.names=FALSE)
write.csv(do.call(rbind,components),file.path(out,"saved_objective_components.csv"),row.names=FALSE)
write.csv(data.frame(path=inputs,md5=unname(hashes)),file.path(out,"input_hashes.csv"),row.names=FALSE)
writeLines(c("Saved capped states only: no fitting, new C, external outcome or weight selection.",
 "Likelihood shares exclude prior penalties; point-prior gating needs separate consideration.",
 "Each coordinate is evaluated independently with other moments fixed, not as a sequential sweep.",
 "Signal share uses absolute numerator terms before cancellation, not variance explained or causal influence.",
 "Cached objective components can use earlier priors; they are not a verified common current-state ELBO.",
 capture.output(sessionInfo())),file.path(out,"interpretation.txt"))
cat("Saved-state loading-balance audit completed without refitting.\n")
