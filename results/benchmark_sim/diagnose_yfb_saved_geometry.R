# ============================================================
# Script: diagnose_yfb_saved_geometry.R
# Purpose: Measure reconstruction and scale changes in saved synthetic replays.
# Created: 2026-10-30
# Dependencies: existing simulation generator; saved diagnostic RDS/CSV files
# ============================================================
source("results/benchmark_sim/run_yfb_supervision_decision.R")
out<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
reps<-read.csv(file.path(out,"diagnostic_replays.csv"))
rows<-list()
for(v in reps$variance) {
  f<-readRDS(file.path(out,sprintf("diagnostic_v%s.rds",v)))
  d<-yfb_decision_data(1,2000,v)
  residual<-d$Y-f$EL%*%t(f$EF[[1]])
  uncertainty<-f$EL2%*%t(f$EF2[[1]])-f$EL^2%*%t(f$EF[[1]]^2)
  row<-data.frame(variance=v,mean_reconstruction_loss=sum(residual^2)/sum(d$Y^2),
    expected_reconstruction_loss=sum(residual^2+uncertainty)/sum(d$Y^2),
    F_norm_ratio=tail(f$history$norm_F,1)/f$history$norm_F[1],
    L_norm_ratio=tail(f$history$norm_L,1)/f$history$norm_L[1],
    median_tau=median(f$Tau),max_beta=max(abs(f$EBeta)))
  stopifnot(all(is.finite(as.matrix(row))),row$expected_reconstruction_loss>=row$mean_reconstruction_loss)
  rows[[length(rows)+1]]<-row
}
write.csv(do.call(rbind,rows),file.path(out,"diagnostic_geometry.csv"),row.names=FALSE)
print(do.call(rbind,rows))
