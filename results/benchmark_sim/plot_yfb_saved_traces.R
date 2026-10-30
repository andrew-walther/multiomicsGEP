# ============================================================
# Script: plot_yfb_saved_traces.R
# Purpose: Plot saved representative supervision trajectories without refitting.
# Created: 2026-10-30
# Dependencies: base R graphics
# ============================================================

#' Draw comparable panels from fixed diagnostic replay CSVs
#' @param out Directory containing diagnostic_replays.csv and its trace CSVs.
#' @return Invisibly NULL; draws on the current device without modifying data.
#' @examples
#' # plot_yfb_saved_traces("results/benchmark_sim/outputs/yfb_supervision_2026-10-30")
#' @family yfb_supervised
#' @seealso run_yfb_saved_fit_diagnostics.R
plot_yfb_saved_traces <- function(out) {
  reps<-read.csv(file.path(out,"diagnostic_replays.csv"))
  reps<-reps[order(reps$variance),]
  traces<-lapply(reps$variance,function(v)
    read.csv(file.path(out,sprintf("diagnostic_trace_v%s.csv",v))))
  values<-lapply(traces,function(h) {
    vals<-list((h$objective-h$objective[1])/max(1,abs(h$objective[1])),
               h$share_median,h$max_beta,log10(h$norm_F/h$norm_F[1]))
    if(any(!is.finite(unlist(vals)))) stop("Nonfinite diagnostic plot value")
    vals
  })
  labels<-c("Objective","Live share","max |beta|","log10 F ratio")
  old<-par(mfrow=c(4,length(traces)),mar=c(2.2,4.2,2,0.7),oma=c(2,1,3.5,0),
           cex=.9,bty="l",las=1)
  on.exit(par(old))
  xmax<-max(vapply(traces,function(h)max(h$iter),0.))
  for(metric in seq_len(4)) {
    limits<-range(unlist(lapply(values,function(v)v[[metric]])))
    if(metric==2)limits<-c(0,1)
    for(j in seq_along(traces)) {
      h<-traces[[j]];y<-values[[j]][[metric]]
      plot(h$iter,y,type="l",lwd=1.5,col="#24577A",xlim=c(1,xmax),ylim=limits,
           xlab="",ylab=if(j==1)labels[metric] else "",
           main=if(metric==1)sprintf("Variance share: %g%%",100*reps$variance[j]) else "")
      if(metric==1) {
        ii<-which(h$objective_decreased)
        points(h$iter[ii],y[ii],pch=16,cex=.4,col="#C45A00")
      }
    }
  }
  mtext("Sweep",side=1,outer=TRUE,line=.5)
  mtext("Stable posterior evaluation does not guarantee stable joint fitting",side=3,
        outer=TRUE,line=1.3,cex=1.1)
  mtext("Dots mark recorded objective decreases",side=3,outer=TRUE,line=.2,cex=.9)
  invisible(NULL)
}

if(sys.nframe()==0L) {
  base<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
  args<-commandArgs(TRUE);attempt<-if(length(args))args[1] else "attempt7"
  out<-file.path(base,paste0(attempt,"_diagnostics"))
  png(file.path(out,"diagnostic_traces.png"),width=1760,height=1600,res=160,
      type=if(Sys.info()[["sysname"]]=="Darwin")"quartz" else getOption("bitmapType"))
  plot_yfb_saved_traces(out)
  dev.off()
}
