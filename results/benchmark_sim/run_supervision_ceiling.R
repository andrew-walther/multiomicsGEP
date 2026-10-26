# ============================================================
# Script: run_supervision_ceiling.R
# Purpose: Survival-only reference for the planted-program test bed: how often
#          does a Cox model fit directly on all genes (ridge, lasso; glmnet,
#          lambda.min by 10-fold CV) recover the planted loading direction
#          (|cor(coefficients, F_true[,3])| >= 0.7)? Any supervised
#          factorization uses no more survival information than this.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: glmnet, survival; run_supervision_diagnostics.R's test bed
# Run: Rscript results/benchmark_sim/run_supervision_ceiling.R [n_seeds] [n_train,...]
#      Several training sizes calibrate the test bed (DECISIONS.md 2026-10-15):
#      the chosen n is the smallest at which the lasso recovers the planted
#      program in >= 90% of seeds at both noise levels.
# ============================================================

suppressPackageStartupMessages({ library(survival); library(glmnet) })
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL))
source("code/concordance_ci.R"); source("code/simulate_multimodal_yfb.R")

args    <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 20L
n_trains <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else 120L
# Planted log hazard ratio of program 3 (scenario value 1.10)
beta3s   <- if (length(args) >= 3) as.numeric(strsplit(args[3], ",")[[1]]) else 1.10
out_dir <- "results/benchmark_sim/outputs/supervision_diagnostics/ceiling"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Same data as make_testbed() in run_supervision_diagnostics.R
testbed <- function(seed, noise, n_train = 120L, beta3 = 1.10) {
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor", n_train, 120L, 50L, 50L,
                                    seed = seed, truncate = FALSE, noise_scale = noise,
                                    beta = c(0, 0, beta3))
  list(Y = scale(d$training$Y$expression), time = d$training$time,
       event = d$training$event, f3 = d$truth$F$expression[, 3])
}

rows <- list()
for (beta3 in beta3s) for (n_train in n_trains) for (noise in c(0.25, 1))
    for (seed in 20261015L + seq_len(n_seeds)) {
  tb <- testbed(seed, noise, n_train, beta3)
  y  <- Surv(tb$time, tb$event)
  coef_cor <- function(a) {
    set.seed(1)
    b <- as.vector(coef(cv.glmnet(tb$Y, y, family = "cox", alpha = a), s = "lambda.min"))
    if (sd(b) > 0) abs(cor(b, tb$f3)) else 0
  }
  # Null-model score direction Y'u0: the supervised initialization of C3
  u0 <- calc_cox_taylor_yf(rep(0, nrow(tb$Y)), tb$time, tb$event)$u
  rows[[length(rows) + 1]] <- data.frame(
    beta3 = beta3, n_train = n_train, noise = noise, seed = seed, events = sum(tb$event), cor_ridge = coef_cor(0), cor_lasso = coef_cor(1),
    cor_score_direction = abs(cor(as.vector(crossprod(tb$Y, u0)), tb$f3)))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(split(res, list(res$noise, res$n_train, res$beta3), drop = TRUE),
                               function(d) data.frame(
  beta3 = d$beta3[1], n_train = d$n_train[1], noise = d$noise[1], n = nrow(d), events_median = median(d$events),
  ridge_recovered = sum(d$cor_ridge >= 0.7), lasso_recovered = sum(d$cor_lasso >= 0.7),
  score_dir_recovered = sum(d$cor_score_direction >= 0.7),
  ridge_median = median(d$cor_ridge), lasso_median = median(d$cor_lasso),
  score_dir_median = median(d$cor_score_direction))))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summ, digits = 3)
