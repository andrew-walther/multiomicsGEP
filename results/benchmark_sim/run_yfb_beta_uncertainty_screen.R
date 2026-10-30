# ============================================================
# Script: run_yfb_beta_uncertainty_screen.R
# Purpose: Separate tempering and mean-field effects in a fixed Gaussian working problem.
# Created: 2026-10-30
# Dependencies: existing YFB posterior helpers; base R; synthetic inputs only
# ============================================================
Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
source("code/fit_yfb_supervised.R")
source("code/yfb_joint_beta.R")
out<-"results/benchmark_sim/outputs/yfb_targeted_2026-10-30"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"beta_metrics.csv")))stop("Existing coefficient screen: do not overwrite")
cfg<-yaml::read_yaml("config/globals.yml");K<-cfg$benchmark$k_merged_yfb_desurv;n<-400
saved<-readRDS("results/benchmark_sim/outputs/yfb_supervision_2026-10-30/attempt7/p2000_v0.0025_s01_share0.05.rds")
set.seed(491);base<-matrix(rnorm(n*K),n,K)
Xorth<-qr.Q(qr(scale(base,scale=FALSE)))*sqrt(n)
set.seed(492);common<-rnorm(n)
Xcorr<-scale(sqrt(.9)*common+sqrt(.1)*base)
Xcorr<-sweep(Xcorr,2,sqrt(colSums(Xcorr^2)/n),"/")
rows<-list();j<-0L
for(design in c("orthogonal","correlated")) {
 X<-if(design=="orthogonal")Xorth else Xcorr
 for(signal in c("null","planted")) for(seed in 1:50) {
  beta<-rep(0,K);if(signal=="planted")beta[1]<-1
  set.seed(8000+seed);y<-drop(X%*%beta)+rnorm(n)
  for(weight in c(1,saved$metrics$gamma)) {
   joint<-yfb_joint_beta(X,y,weight,1/K,1)
   mf<-yfb_meanfield_beta(X,y,weight,1/K,1)
   for(method in c("joint","mean_field")) {
    q<-if(method=="joint")joint else mf
    active<-abs(q$mean)>=1.96*sqrt(q$variance) & q$variance>0
    nulls<-if(signal=="null")seq_len(K) else 2:K
    j<-j+1L;rows[[j]]<-data.frame(design=design,signal=signal,seed=seed,
      weight=weight,method=method,false_active=sum(active[nulls]),
      recovered=if(signal=="planted")active[1] && q$mean[1]>0 else NA,
      mean_sd=mean(sqrt(q$variance)),max_abs_z=max(ifelse(q$variance>0,abs(q$mean)/sqrt(q$variance),0)),
      converged=if(method=="joint")TRUE else q$converged,
      iterations=if(method=="joint")NA_integer_ else q$iterations)
   }
  }
 }
}
metrics<-do.call(rbind,rows);write.csv(metrics,file.path(out,"beta_metrics.csv"),row.names=FALSE)
groups<-split(metrics,interaction(metrics$design,metrics$signal,metrics$weight,metrics$method,drop=TRUE))
summary<-do.call(rbind,lapply(groups,function(d)data.frame(design=d$design[1],signal=d$signal[1],weight=d$weight[1],method=d$method[1],
 n=nrow(d),mean_false=mean(d$false_active),se_false=sd(d$false_active)/sqrt(nrow(d)),
 recovered=if(all(is.na(d$recovered)))NA_integer_ else sum(d$recovered),
 mean_sd=mean(d$mean_sd),nonconverged=sum(!d$converged))))
write.csv(summary,file.path(out,"beta_summary.csv"),row.names=FALSE)
