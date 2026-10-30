# ============================================================
# Script: audit_yfb_variance_storage.R
# Purpose: Audit stored coefficient uncertainty without fitting or omitting failures.
# Created: 2026-10-30
# Dependencies: base R; synthetic attempt checkpoints
# ============================================================
out<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
rows<-list()
for(attempt in paste0("attempt",seq_len(6))) {
  files<-list.files(file.path(out,attempt),pattern="rds$",full.names=TRUE)
  for(path in files) {
    saved<-readRDS(path);f<-saved$fit
    if(is.null(f) || !is.list(f$EF) || is.null(f$g_beta))next
    g<-f$g_beta;nondegenerate<-any(g$pi[g$sd>0]>0)
    diff<-f$EBeta2-f$EBeta^2
    vv<-if(!is.null(f$VBeta))f$VBeta else diff
    rows[[length(rows)+1]]<-data.frame(attempt=attempt,p=saved$metrics$p,
      variance=saved$metrics$variance,seed=saved$metrics$seed,arm=saved$metrics$arm,
      separate_variance=!is.null(f$VBeta),nondegenerate_prior=nondegenerate,
      min_subtracted_variance=min(diff),min_carried_variance=min(vv),
      loss_count=if(nondegenerate)sum(!is.finite(vv)|vv<=0) else 0L,
      max_beta=max(abs(f$EBeta)))
  }
}
write.csv(do.call(rbind,rows),file.path(out,"variance_storage_audit.csv"),row.names=FALSE)
