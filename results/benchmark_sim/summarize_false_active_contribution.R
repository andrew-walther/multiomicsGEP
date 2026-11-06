# ============================================================
# Script: summarize_false_active_contribution.R
# Purpose: How much the programs other than the planted one contribute to the validation risk
#          score in S-13 (CV-selected fits) and S-14 (all arms): risk-variance share, standardized
#          coefficients and the change in validation C when they are removed.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: survival, ebnm; run_yfb_survival_weight_p2000.R, run_yfb_supervision_strength.R,
#               run_yfb_near_threshold.R, code/concordance_ci.R
# Run: Rscript results/benchmark_sim/summarize_false_active_contribution.R
# ============================================================
# Descriptive, post hoc: reads saved checkpoints, regenerates the synthetic data from their seeds,
# refits nothing and changes no registered verdict. Synthetic data only.
#
# For a fit with loadings F (p x K), coefficients b and validation input Y_v, the risk score is
# eta = Y_v F b = sum_k eta_k with eta_k = (Y_v F_k) b_k. With m the program best correlated with
# the planted direction:
#   other_share = Var(sum_{k != m} eta_k) / Var(eta)
#   std_beta_k  = |b_k| * sd(Y_tr F_k)          (log hazard ratio per training SD of program k)
#   delta_C     = C(eta) - C(eta_m)             (what the other programs add to validation C)

Sys.setenv(OPENBLAS_NUM_THREADS = 1)
source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
source("results/benchmark_sim/run_yfb_supervision_strength.R")
source("results/benchmark_sim/run_yfb_near_threshold.R")
source("code/concordance_ci.R")

#' Contribution of the non-matched programs for one fit
#' @param d Data list with Ytr, Yva, truth, tv, ev. @param EF Loadings (p x K).
#' @param b Coefficients (length K). @param active Logical K-vector of active programs, or NULL.
#' @return One-row data frame.
#' @examples
#' # fa_contribution(d, fit$EF, fit$EBeta, NULL)
fa_contribution <- function(d, EF, b, active = NULL) {
  r <- suppressWarnings(cor(d$truth, EF)); r[!is.finite(r)] <- 0; m <- which.max(abs(r))
  Zv <- d$Yva %*% EF; Zt <- d$Ytr %*% EF
  eta <- drop(Zv %*% b); eta_m <- Zv[, m] * b[m]; other <- eta - eta_m
  sb <- abs(b) * apply(Zt, 2, sd)
  fa <- if (is.null(active)) setdiff(seq_along(b), m) else setdiff(which(active), m)
  cfun <- function(x) if (sd(x) > 1e-8) frozen_reverse_cindex(x, d$tv, d$ev) else NA_real_
  data.frame(other_share = if (var(eta) > 0) var(other) / var(eta) else NA_real_,
             std_beta_matched = sb[m],
             std_beta_other_max = if (length(fa)) max(sb[fa]) else 0,
             n_counted = length(fa), c_full = cfun(eta), c_matched_only = cfun(eta_m),
             delta_c = cfun(eta) - cfun(eta_m))
}

# S-13: CV-selected signal fits (activity from the posterior, |mean| / sd >= 1.96) ----
S13 <- "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05"
pick <- read.csv(file.path(S13, "cv_selected.csv"))
s13 <- do.call(rbind, lapply(seq_len(nrow(pick)), function(i) {
  p <- pick[i, ]
  ck <- readRDS(file.path(S13, "fits", sprintf("sig_%s_s%02d_%s_w%g.rds", p$cell, p$seed, p$arm, p$w_S)))
  d <- s13_data(p$seed, p$cell)
  b <- ck$state$EBeta; act <- abs(b) / sqrt(pmax(ck$state$EBeta2 - b^2, 1e-300)) >= 1.96
  cbind(p[, c("cell", "seed", "arm", "w_S", "false_active")], fa_contribution(d, ck$state$EF, b, act))
}))
write.csv(s13, file.path(S13, "false_active_contribution.csv"), row.names = FALSE)

# S-14: all arms (activity not stored, so every non-matched program is counted) ----
S14 <- "results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05"
m14 <- read.csv(file.path(S14, "metrics.csv"))
s14 <- do.call(rbind, lapply(seq_len(nrow(m14)), function(i) {
  t <- m14[i, ]
  key <- sprintf("%s_sh%g_s%02d_%s", t$strength, t$share, t$seed, t$arm)
  ck <- readRDS(file.path(S14, "fits", paste0(key, ".rds")))
  d <- s14_data(t$seed, t$share, t$strength)
  cbind(t[, c("strength", "share", "seed", "arm", "recovered", "false_active")],
        fa_contribution(d, ck$EF, ck$beta, NULL))
}))
write.csv(s14, file.path(S14, "false_active_contribution.csv"), row.names = FALSE)

# Cell summaries ----
agg <- function(x, by) do.call(rbind, lapply(split(x, x[by], drop = TRUE), function(z)
  cbind(z[1, by, drop = FALSE], n = nrow(z), false_active = mean(z$false_active),
        other_share_median = median(z$other_share, na.rm = TRUE),
        std_beta_matched_median = median(z$std_beta_matched),
        std_beta_other_max_median = median(z$std_beta_other_max),
        delta_c_mean = mean(z$delta_c, na.rm = TRUE))))
cat("\nS-13, CV-selected fits (false actives only)\n")
print(agg(s13, c("arm", "cell")), digits = 3, row.names = FALSE)
cat("\nS-14, all non-matched programs\n")
print(agg(s14, c("strength", "share", "arm")), digits = 3, row.names = FALSE)
