# ============================================================
# Script: run_k_selection_real_1se.R
# Purpose: Apply the conditional training-only 1-SE starting-rank rule
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, yaml; run_supervised_train_cv.R loading/fold definitions
# Run: Rscript results/benchmark_sim/run_k_selection_real_1se.R
# ============================================================
# Fixed before execution: K={3,5,7,9,11,13}, original 5 cohort/event-stratified
# folds (seed 2026), fit seed 42; selected shared point-Laplace configuration
# with protected pruning. Highest mean fold C; threshold = that mean minus
# SD(fold C)/sqrt(5); smallest qualifying K. No external cohorts loaded.
# Run ONLY after the simulation adoption table explicitly adopts single modality.

out_dir <- "results/benchmark_sim/outputs/k_selection/round2/attempt2"
adoption <- read.csv(file.path(out_dir,"adoption.csv"))
if (!any(adoption$model == "single" & adoption$adopted))
  stop("Protected pruning not adopted for single modality; the conditional 1-SE run is ineligible.")
# Read exactly the existing training-loading and fold code, ending before fits.
for (e in parse("results/benchmark_sim/run_supervised_train_cv.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("common"))) break
  eval(e)
}
# The prefix sets OUT_DIR to the old CV location; every new output goes here.
OUT_DIR <- out_dir
rows <- list()
for (K in c(3L,5L,7L,9L,11L,13L)) for (f in 1:5) {
  tr <- which(fold!=f); te <- which(fold==f);set.seed(42)
  fit <- fit_cox_on_yf(Y_train[tr,,drop=FALSE],time_train[tr],status_train[tr],K=K,
    max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,prior_LF="point_laplace",
    canonicalize=TRUE,w_S=1,init_supervised=TRUE,beta_prior_mode="shared",
    prior_beta="point_laplace",beta_warmup=10L,prune=TRUE,prune_protect_survival=TRUE,
    verbose=FALSE)
  risk <- predict_cox_on_yf(Y_train[te,,drop=FALSE],fit$EF,fit$EBeta,EF_norms=fit$EF_norms)$risk_scores
  rows[[length(rows)+1L]] <- data.frame(K=K,fold=f,
    C=frozen_reverse_cindex(risk,time_train[te],status_train[te]),
    K_final=ncol(fit$EF),converged=fit$history$converged)
}
r <- do.call(rbind,rows)
write.csv(r,file.path(out_dir,"real_1se_folds.csv"),row.names=FALSE)
s <- do.call(rbind,lapply(split(r,r$K),function(x)
  data.frame(K=x$K[1],mean_C=mean(x$C),SE=sd(x$C)/sqrt(5),converged=sum(x$converged))))
best <- which.max(s$mean_C);cutoff <- s$mean_C[best]-s$SE[best]
chosen <- min(s$K[s$mean_C>=cutoff]);s$selected <- s$K==chosen
write.csv(s,file.path(out_dir,"real_1se_summary.csv"),row.names=FALSE)
set.seed(42)
fit <- fit_cox_on_yf(Y_train,time_train,status_train,K=chosen,
  max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,prior_LF="point_laplace",
  canonicalize=TRUE,w_S=1,init_supervised=TRUE,beta_prior_mode="shared",
  prior_beta="point_laplace",beta_warmup=10L,prune=TRUE,prune_protect_survival=TRUE,verbose=FALSE)
beta_sd <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
active <- abs(fit$EBeta)>=1.96*beta_sd & beta_sd>0
write.csv(data.frame(K_init=chosen,K_final=ncol(fit$EF),program=seq_along(active),
  beta=fit$EBeta,survival_active=active,converged=fit$history$converged),
  file.path(out_dir,"real_1se_full_fit.csv"),row.names=FALSE)
print(s)
