# ============================================================
# Script: k5_joint_vs_twostep_check.R
# Purpose: Compare the rank-5 joint SSBMF fit (w_S = 1) with the w_S = 0 fit behind the
#          two-step EBMF and Cox comparator: loading and score agreement, top-50 overlap,
#          and Efron vs Breslow ties in the second-step Cox model. Source for Chapter 4
#          Appendix C.12.
# Author: Andrew Walther
# Created: 2026-11-15
# Dependencies: survival
# Run from the repository root:
#   Rscript results/benchmark_sim/k5_joint_vs_twostep_check.R \
#     > results/benchmark_sim/outputs/k5_working_2026-11-12/joint_vs_twostep_check.txt
# ============================================================
suppressMessages(library(survival))
z <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
d <- "results/benchmark_sim/outputs/k5_working_2026-11-12/s11/fits/"
f1 <- readRDS(paste0(d,"size_D4_N2218_w1_f0.rds")); f0 <- readRDS(paste0(d,"size_D4_N2218_w0_f0.rds"))
stopifnot(identical(f1$genes,f0$genes))
F1 <- f1$state$EF; F0 <- f0$state$EF
cat("loading cor matrix (rows w1, cols w0):\n"); print(round(cor(F1,F0),7))
Y <- z$Y[, f1$genes]
Z1 <- Y %*% F1; Z0 <- Y %*% F0
cat("score cor diag:", round(diag(cor(Z1,Z0)),7), "\n")
top <- function(F,k,n=50) order(-abs(F[,k]))[1:n]
cat("top50 overlap:", sapply(1:5, function(k) length(intersect(top(F1,k),top(F0,k)))), "\n")
cat("conv iter w1/w0:", length(f1$state$elbo), length(f0$state$elbo), "\n")
S <- Surv(z$time, z$status)
ce <- coxph(S ~ Z0, ties="efron"); cb <- coxph(S ~ Z0, ties="breslow")
print(rbind(efron=coef(ce), breslow=coef(cb)))
cat("risk cor:", cor(Z0%*%coef(ce), Z0%*%coef(cb)), " C efron/breslow:",
  concordance(S ~ I(-Z0%*%coef(ce)))$concordance, concordance(S ~ I(-Z0%*%coef(cb)))$concordance, "\n")
cat("EBeta w1:", round(f1$state$EBeta,4), "\n")
