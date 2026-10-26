# ============================================================
# Script: summarize_feature_expansion.R
# Purpose: Report fixed-K feature performance and program agreement with frozen references.
# Created: 2026-10-15
# Dependencies: yaml, survival, existing feature sensitivity helpers; base graphics
# Run: Rscript results/benchmark_sim/summarize_feature_expansion.R
# ============================================================
# No fitting, selection or external scoring. Undefined fold/constant-factor
# quantities remain NA. Patient matrices are read only from ignored local RDS.
source("code/feature_sensitivity.R")
out <- "results/benchmark_sim/outputs/feature_expansion/cold_primary"
z <- readRDS("data/round2_feature_cold_primary.rds")
files <- file.path(out,sprintf("cell_%02d.csv",z$cells$cell))
stopifnot(all(file.exists(files)))
rows <- do.call(rbind,lapply(files,feature_read_terminal))
stopifnot(nrow(rows)==24L,!any(nzchar(rows$error)))
reference <- read.csv(file.path(out,"benchmark_rows.csv"))
counts <- z$manifest$design$actual_shared
summary <- vector("list",length(counts))
for(i in seq_along(counts)) {
  a <- rows[rows$p==counts[i]&rows$fold>0,]
  b <- rows[rows$p==counts[i]&rows$fold==0,]
  stopifnot(nrow(a)==5L,nrow(b)==1L)
  complete <- all(is.finite(a$C))
  summary[[i]] <- data.frame(source="feature_sensitivity",p=counts[i],K=7,
    training_C=b$C,conditional_CV_C=if(complete)mean(a$C) else NA_real_,
    conditional_CV_SE=if(complete)stats::sd(a$C)/sqrt(nrow(a)) else NA_real_,
    finite_CV_folds=sum(is.finite(a$C)),converged_CV_folds=sum(a$converged),
    training_converged=b$converged,retained=b$retained,active=b$active,
    total_seconds=sum(c(a$runtime_seconds,b$runtime_seconds)),
    constant_columns=b$constant_columns)
}
summary <- do.call(rbind,summary)
benchmark <- summary[1,,drop=FALSE];benchmark[] <- NA
benchmark$source <- "benchmark_snapshot";benchmark$p <- reference$p
benchmark$K <- reference$K;benchmark$training_C <- reference$C
benchmark$training_converged <- reference$converged
write.csv(rbind(summary,benchmark),file.path(out,"summary_with_benchmark.csv"),row.names=FALSE)

# Full-training agreement ----
fits <- lapply(seq_along(counts),function(i) {
  id <- z$cells$cell[z$cells$set==i&z$cells$fold==0]
  readRDS(file.path(out,sprintf("cell_%02d.rds",id)))
})
pairs <- list();risk <- list();matrices <- list();k <- 0L
for(i in seq_along(counts)) {
  refs <- if(i==1L)"frozen_D4" else c("fresh_baseline","preceding_count")
  for(ref in refs) {
    if(ref=="frozen_D4") {
      old <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_fits.rds")$D4
      # Previously frozen historical risk orientation, not data-driven here.
      old$EBeta <- -old$EBeta
      old_genes <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/d4_gene_names.rds")
    } else {
      j <- if(ref=="fresh_baseline")1L else i-1L
      old <- fits[[j]]$fit;old_genes <- fits[[j]]$genes
      stopifnot(identical(fits[[j]]$train,fits[[i]]$train))
    }
    new <- fits[[i]]
    a <- feature_program_agreement(old,old_genes,new$fit,new$genes,
      z$Y[,match(old_genes,z$genes),drop=FALSE],
      z$Y[,match(new$genes,z$genes),drop=FALSE])
    k <- k+1L
    pairs[[k]] <- cbind(reference=ref,p=counts[i],a$pairs)
    risk[[k]] <- data.frame(reference=ref,p=counts[i],risk_correlation=a$risk_correlation)
    matrices[[k]] <- a
  }
}
write.csv(do.call(rbind,pairs),file.path(out,"program_agreement.csv"),row.names=FALSE)
write.csv(do.call(rbind,risk),file.path(out,"risk_agreement.csv"),row.names=FALSE)
saveRDS(matrices,file.path(out,"agreement_matrices.rds"))

# Performance plot ----
png(file.path(out,"feature_count_performance.png"),width=1800,height=850,res=160)
par(mfrow=c(1,2),mar=c(5,5,3,1))
limits <- range(c(summary$training_C,reference$C,
  summary$conditional_CV_C-summary$conditional_CV_SE,
  summary$conditional_CV_C+summary$conditional_CV_SE),finite=TRUE)
plot(counts,summary$conditional_CV_C,type="b",pch=19,col="#176b91",ylim=limits,
  xlab="Shared expression genes",ylab="Concordance",main="K=7, alpha_F=0: conditional training CV")
arrows(counts,summary$conditional_CV_C-summary$conditional_CV_SE,
  counts,summary$conditional_CV_C+summary$conditional_CV_SE,
  angle=90,code=3,length=.05,col="#176b91")
lines(counts,summary$training_C,type="b",pch=17,col="#ad5a10")
abline(h=reference$C,lty=2,col="#555555")
legend("bottomleft",c("Five-fold CV mean +/- SE","Full-training C","Frozen D4 training C"),
  col=c("#176b91","#ad5a10","#555555"),pch=c(19,17,NA),lty=c(1,1,2),bty="n",cex=.8)
plot(counts,summary$retained,type="b",pch=19,ylim=c(0,7),
  xlab="Shared expression genes",ylab="Operational program count",main="Shrinkage at fixed seven available slots")
lines(counts,summary$active,type="b",pch=17,col="#ad5a10")
legend("topleft",c("Retained","Survival-active (legacy threshold)"),
  col=c("black","#ad5a10"),pch=c(19,17),lty=1,bty="n",cex=.8)
dev.off()
writeLines(c("CV conditions on screening/scaling before patient folds; not nested preprocessing CV.",
  "Error bars are fold-level SE, not an external performance confidence interval.",
  "No feature count/default adopted; no external scoring.",
  "Operational activity is not Bayesian credibility; constant/dead correlations remain NA.",
  "Frozen external mean C approximately .627 is historical, not the training/CV comparison."),
  file.path(out,"interpretation.txt"))
print(summary)

# Agreement curve on full-training fits ----
# Show the matching denominator: a high correlation among surviving matches
# does not establish preservation of every original program.
a <- do.call(rbind,pairs)
r <- do.call(rbind,risk)
comparison <- data.frame(p=counts[-1],matched=0L,activity=NA_real_,loading=NA_real_,
  projection=NA_real_,risk=NA_real_)
for(i in seq_len(nrow(comparison))) {
  take <- a[a$reference=="fresh_baseline"&a$p==comparison$p[i],]
  matched <- !is.na(take$old_factor)&!is.na(take$new_factor)
  comparison$matched[i] <- sum(matched)
  for(field in c("activity","loading","projection")) {
    values <- take[[paste0(field,"_correlation")]][matched]
    values <- values[is.finite(values)]
    comparison[i,field] <- if(length(values))stats::median(values) else NA_real_
  }
  comparison$risk[i] <- r$risk_correlation[r$reference=="fresh_baseline"&r$p==comparison$p[i]]
}
write.csv(comparison,file.path(out,"agreement_curve.csv"),row.names=FALSE)
png(file.path(out,"feature_count_agreement.png"),width=1500,height=850,res=160)
par(mar=c(5,5,3,1))
matplot(comparison$p,comparison[,c("activity","loading","projection","risk")],
  type="b",pch=1:4,lty=1:4,col=c("#176b91","#ad5a10","#3d7d41","#733b88"),
  ylim=c(-1,1),xlab="Shared expression genes",ylab="Correlation with fresh 2,064-gene fit",
  main="Full-training agreement: median over matched programs")
text(comparison$p,rep(-.50,nrow(comparison)),labels=paste("Matched:",comparison$matched),pos=c(4,3,2),cex=.8)
legend("bottomleft",c("Activity median","Shared-loading median","Projection median","Risk score"),
  col=c("#176b91","#ad5a10","#3d7d41","#733b88"),pch=1:4,lty=1:4,bty="n",cex=.8)
dev.off()

# Loading shrinkage, without inventing a gene-retention threshold ----
loading_rows <- list()
for(i in 0:nrow(rows)) {
  if(i==0L) {
    fit <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_fits.rds")$D4
    source_label <- "benchmark_snapshot";cell <- NA_integer_;fold <- 0L
  } else {
    cell <- rows$cell[i];fold <- rows$fold[i];source_label <- "feature_sensitivity"
    fit <- readRDS(file.path(out,sprintf("cell_%02d.rds",cell)))$fit
  }
  for(k in seq_len(ncol(fit$EF))) {
    variance <- fit$EF2[,k]-fit$EF[,k]^2
    loading_rows[[length(loading_rows)+1L]] <- data.frame(source=source_label,
      cell=cell,fold=fold,p=nrow(fit$EF),factor=k,
      nonzero_posterior_means=sum(fit$EF[,k]!=0),
      mean_abs_loading=mean(abs(fit$EF[,k])),median_abs_loading=stats::median(abs(fit$EF[,k])),
      mean_posterior_variance=mean(variance),min_posterior_variance=min(variance),
      negative_variance_entries=sum(variance<0))
  }
}
write.csv(do.call(rbind,loading_rows),file.path(out,"loading_shrinkage_with_benchmark.csv"),row.names=FALSE)
