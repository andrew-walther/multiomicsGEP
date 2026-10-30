# ============================================================
# Script: diagnose_yfb_saved_geometry.R
# Purpose: Measure reconstruction and scale changes in saved synthetic replays.
# Created: 2026-10-30
# Dependencies: existing simulation generator; saved diagnostic RDS/CSV files
# ============================================================
source("results/benchmark_sim/run_yfb_supervision_decision.R")
base<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
args<-commandArgs(TRUE);attempt<-if(length(args))args[1] else "attempt7"
out<-file.path(base,paste0(attempt,"_diagnostics"))
reps<-read.csv(file.path(out,"diagnostic_replays.csv"))
rows<-list();factors<-list();inference<-list();instability<-list()
for(v in reps$variance) {
  f<-readRDS(file.path(out,sprintf("diagnostic_v%s.rds",v)))
  d<-yfb_decision_data(1,2000,v,truth_mode=if(attempt=="attempt7")"shared" else "legacy")
  residual<-d$Y-f$EL%*%t(f$EF[[1]])
  uncertainty<-if(!is.null(f$VL))
    f$VL%*%t(f$EF[[1]]^2)+f$EL^2%*%t(f$VF[[1]])+f$VL%*%t(f$VF[[1]]) else
    f$EL2%*%t(f$EF2[[1]])-f$EL^2%*%t(f$EF[[1]]^2)
  row<-data.frame(variance=v,mean_reconstruction_loss=sum(residual^2)/sum(d$Y^2),
    expected_reconstruction_loss=sum(residual^2+uncertainty)/sum(d$Y^2),
    F_norm_ratio=tail(f$history$norm_F,1)/f$history$norm_F[1],
    L_norm_ratio=tail(f$history$norm_L,1)/f$history$norm_L[1],
    median_tau=median(f$Tau),min_tau=min(f$Tau),
    max_gene_relative_loss=max(colSums(residual^2+uncertainty)/colSums(d$Y^2)),
    max_loading_energy_fraction=max(rowSums(f$EF[[1]]^2))/sum(f$EF[[1]]^2),
    max_beta=max(abs(f$EBeta)))
  stopifnot(all(is.finite(as.matrix(row))),row$expected_reconstruction_loss>=row$mean_reconstruction_loss)
  rows[[length(rows)+1]]<-row
  # These descriptors use the saved posterior only; no coefficient refit or
  # alternative uncertainty definition enters recovery or the gate.
  h<-f$history
  stopifnot(nrow(h)>1,all(is.finite(h$norm_F)),all(h$norm_F>0))
  jumps<-diff(log(h$norm_F));jf<-which.max(jumps)+1L
  jb<-which.max(diff(h$max_beta))+1L;jo<-which.min(diff(h$objective))+1L
  instability[[length(instability)+1]]<-data.frame(variance=v,
    largest_F_jump_sweep=h$iter[jf],F_jump_ratio=exp(jumps[jf-1L]),
    largest_beta_increase_sweep=h$iter[jb],
    largest_objective_drop_sweep=h$iter[jo],
    largest_objective_drop=diff(h$objective)[jo-1L])
  F<-f$EF[[1]];rr<-drop(suppressWarnings(cor(d$truth,F)))
  vz<-if(!is.null(f$VBeta))f$VBeta else f$EBeta2-f$EBeta^2
  z<-f$EBeta/sqrt(vz)
  factors[[length(factors)+1]]<-data.frame(variance=v,factor=seq_along(z),
    loading_correlation=rr,beta=f$EBeta,beta_variance=vz,z=z,
    survival_active=is.finite(z)&abs(z)>=1.96,
    activity_second=colMeans(f$EL2),loading_norm=sqrt(colSums(F^2)),
    orphan=colSums(f$EL2)==0 & colSums(F^2)>0)
  Z<-f$Z;live<-apply(Z,2,sd)>0
  correlations<-if(sum(live)>1)cor(Z[,live,drop=FALSE]) else matrix(NA_real_,1,1)
  off<-abs(correlations[upper.tri(correlations)])
  inference[[length(inference)+1]]<-data.frame(variance=v,
    gamma=f$gamma,rho=f$rho,max_projection_correlation=if(length(off))max(off) else NA_real_,
    min_beta_sd=min(sqrt(vz)),max_abs_z=max(abs(z)),
    live_projection_columns=sum(live))
}
write.csv(do.call(rbind,rows),file.path(out,"diagnostic_geometry.csv"),row.names=FALSE)
print(do.call(rbind,rows))

write.csv(do.call(rbind,factors),file.path(out,"diagnostic_factors.csv"),row.names=FALSE)
write.csv(do.call(rbind,inference),file.path(out,"diagnostic_inference.csv"),row.names=FALSE)

write.csv(do.call(rbind,instability),file.path(out,"diagnostic_instability.csv"),row.names=FALSE)
