# ============================================================
# Script: run_yfb_saved_fit_diagnostics.R
# Purpose: Replay fixed registered fits with factor-wise live precision traces.
# Created: 2026-10-30
# Dependencies: registered runner and saved synthetic checkpoints
# ============================================================
source("results/benchmark_sim/run_yfb_supervision_decision.R")
out<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
rows<-list()
# Inspect the smallest share-target arm, seed 1, in every large-p variance
# cell. These are exact diagnostic replays, not additional gate datasets.
for(variance in c(.0025,.01,.05)) {
  path<-file.path(out,"attempt4",sprintf("p2000_v%s_s01_share0.05.rds",variance))
  saved<-readRDS(path);d<-yfb_decision_data(1,2000,variance)
  f<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,
    K=yaml::read_yaml("config/globals.yml")$benchmark$k_merged_yfb_desurv,
    gamma=saved$metrics$gamma,rho=saved$metrics$rho,scaled=TRUE)
  saveRDS(f,file.path(out,sprintf("diagnostic_v%s.rds",variance)))
  write.csv(f$history,file.path(out,sprintf("diagnostic_trace_v%s.csv",variance)),row.names=FALSE)
  m<-yfb_decision_score(f,d)
  shares<-yfb_precision_share_final(f,d$Y,d$time,d$event)
  m$variance<-variance;m$share_final_live<-median(shares,na.rm=TRUE)
  m$largest_mean_difference<-if(is.null(saved$fit))NA_real_ else
    max(abs(f$EF[[1]]-saved$fit$EF[[1]]))
  m$first_decrease<-if(any(f$history$objective_decreased))
    min(f$history$iter[f$history$objective_decreased]) else NA_integer_
  rows[[length(rows)+1]]<-m
}
write.csv(do.call(rbind,rows),file.path(out,"diagnostic_replays.csv"),row.names=FALSE)
