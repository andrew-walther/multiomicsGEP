# ============================================================
# Script: plot_yfb_targeted_diagnostics.R
# Purpose: Compare saved reconstruction and precision trajectories without refitting.
# Created: 2026-10-30
# Dependencies: base R graphics, yaml; synthetic targeted checkpoints
# ============================================================
b<-"results/benchmark_sim/outputs/yfb_targeted_2026-10-30"
design<-read.csv(file.path(b,"design.csv"));vs<-sort(unique(design$variance))
max_iter<-yaml::read_yaml("config/globals.yml")$cavi$max_iter
arms<-c("adaptive","fixed","guard_adaptive","guard_fixed")
labels<-c("Adaptive noise","Fixed noise","Guarded adaptive","Guarded fixed")
cols<-c("#B3482A","#24577A","#418361","#8A5AA0")
traces<-lapply(vs,function(v)lapply(arms,function(a)
  readRDS(file.path(b,sprintf("v%s_%s.rds",v,a)))$fit$history))
for(z in traces)for(h in z)stopifnot(all(is.finite(h$expected_reconstruction_ratio)),
 all(h$expected_reconstruction_ratio>0),all(h$min_tau>0))
vals<-c("expected_reconstruction_ratio","min_tau")
ylabels<-c("log10 expected error ratio","log10 minimum precision")
png(file.path(b,"trajectories.png"),width=1800,height=1050,res=160,
 type=if(Sys.info()[["sysname"]]=="Darwin")"quartz" else getOption("bitmapType"))
par(mfrow=c(2,3),mar=c(3.2,4.5,2.4,.8),oma=c(3.6,0,1,0),bty="l",las=1)
for(metric in seq_along(vals)) {
 limits<-range(unlist(lapply(traces,function(z)lapply(z,function(h)log10(h[[vals[metric]]])))))
 for(j in seq_along(vs)) {
  plot(NA,xlim=c(1,max_iter),
   ylim=limits,xlab="",ylab=if(j==1)ylabels[metric] else "",
   main=if(metric==1)sprintf("Raw variance share: %g%%",100*vs[j]) else "")
  for(k in seq_along(arms))lines(traces[[j]][[k]]$iter,
    log10(traces[[j]][[k]][[vals[metric]]]),col=cols[k],lty=k,lwd=1.7)
 }
}
par(fig=c(0,1,0,.10),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0))
plot.new();text(.5,.78,"Accepted sweep",cex=1)
legend(.5,.2,xjust=.5,yjust=.5,legend=labels,col=cols,lty=seq_along(arms),lwd=1.7,
 horiz=TRUE,bty="n",cex=.95)
dev.off()
