# ============================================================
# Script: run_yfb_precision_audit.R
# Purpose: Audit scaled per-gene survival precision on training data only.
# Created: 2026-10-30
# Dependencies: yaml; existing preprocessing and log-space Cox helpers
# ============================================================
source("code/preprocess_desurv.R")
source("code/cox_working_logspace.R")
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
out <- "results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
raw <- lapply(setNames(cfg$pdac$training_cohorts, cfg$pdac$training_cohorts),
              load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
time <- unlist(lapply(raw, `[[`, "time")); event <- unlist(lapply(raw, `[[`, "status"))
pp <- preprocess_merged_cohorts(raw, PLATFORM_LOG_TRANSFORM[names(raw)],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = FALSE, selection_per_cohort = TRUE,
  selection_method = "combined_rank")
# The single-log saved fit uses this same gene set, z-standardized.
x <- readRDS("results/benchmark_sim/outputs/kl_fix_rerun_2026-10-30/fits/full_single.rds")
stopifnot(identical(pp$gene_names, x$genes))
fit <- x$fit
Y <- do.call(rbind,lapply(raw,function(z)z$Y[,match(pp$gene_names,z$gene_names),drop=FALSE]))
for (co in unique(pp$dataset_labels)) {
  ii <- which(pp$dataset_labels == co)
  sdj <- apply(Y[ii,,drop=FALSE], 2, sd)
  if (any(sdj <= 0)) stop("Zero training cohort SD")
  Y[ii,] <- sweep(Y[ii,,drop=FALSE], 2, sdj, "/")
}
Yc <- scale(Y, scale = FALSE); n <- nrow(Y); p <- ncol(Y)
# RMS canonicalization leaves LF′ unchanged; β is in the new normalized projection units.
c <- sqrt(colMeans(fit$EL2)); L <- sweep(fit$EL,2,c,"/")
L2 <- sweep(fit$EL2,2,c^2,"/"); F <- sweep(fit$EF,2,c,"*")
F2 <- sweep(fit$EF2,2,c^2,"*")
s <- 1 / sqrt(colSums(F^2) + p * 1e-10)
# Preserve the saved coefficient's risk contribution on the rescaled expression.
b <- fit$EBeta / (fit$EF_norms * c * s)
b2 <- fit$EBeta2 / (fit$EF_norms * c * s)^2
Z <- sweep(Yc %*% F, 2, s, "*")
tay <- multimodal_yfb_cox_working_logspace(time,event,drop(Z %*% b))
w <- pmax(tay$w,1e-6)
Arec <- outer(fit$Tau,colSums(L2))
Asurv <- outer(colSums(Yc^2*w), b2*s^2)
rows <- list(); z <- 0L
for (alpha in c(.5,.9,.99,.999)) for(k in seq_along(b)) {
 gamma <- alpha/(1-alpha)*p/sum(event)
 # Numerical direct sum independently checks the matrix identity.
 numeric_A <- vapply(seq_len(p),function(j) b2[k]*sum(w*(s[k]*Yc[,j])^2),0.)
 stopifnot(max(abs(numeric_A-Asurv[,k])) < 1e-8)
 z<-z+1L; rows[[z]]<-data.frame(n=n,p=p,events=sum(event),factor=k,alpha=alpha,
 gamma=gamma,median_share=median(gamma*Asurv[,k]/(Arec[,k]+gamma*Asurv[,k])),
 gamma_share_01=.1/.9/median(Asurv[,k]/Arec[,k]),
 gamma_share_05=1/median(Asurv[,k]/Arec[,k]))
}
write.csv(do.call(rbind,rows),file.path(out,"precision_audit.csv"),row.names=FALSE)
write.csv(data.frame(factor=seq_along(b),beta=b,variance=b2-b^2,
 stored_beta=fit$EBeta,stored_variance=fit$EBeta2-fit$EBeta^2),
 file.path(out,"sign_variance_audit.csv"),row.names=FALSE)
saveRDS(list(Y=Y,time=time,event=event,cohort=pp$dataset_labels,genes=pp$gene_names),
 file.path(out,"training_only.rds"))
