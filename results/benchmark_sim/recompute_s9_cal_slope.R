# ============================================================
# Script: recompute_s9_cal_slope.R
# Purpose: Recompute the S-9 validation calibration slopes from the saved fits with exact event
#          times (coxph timefix = FALSE), plus the true predictor's slope as a reference.
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: survival, parallel; results/benchmark_sim/run_yfb_survival_weight_p2000.R
# Run: Rscript results/benchmark_sim/recompute_s9_cal_slope.R   (then the S-9 summarizer)
# ============================================================
# The grid's saved slopes used coxph's default tie correction, which merges the extreme simulated
# times at p = 2,000 and flattens every slope. Synthetic data only. PSOCK, one BLAS thread each.
Sys.setenv(OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
G <- file.path(S9_OUT, "grid")
design <- expand.grid(seed = 1:10, p = c(50L, 2000L), variance = c(.0025, .01, .05))
one <- function(row) {
  suppressPackageStartupMessages(library(survival)); s9_load()
  d <- yfb_decision_data(row$seed, row$p, row$variance, truth_mode = "shared")
  Yva <- s9_zscore(d$Yv); ctl <- coxph.control(timefix = FALSE)
  sl <- function(x) if (sd(x) > 0) unname(coef(coxph(Surv(d$tv, d$ev) ~ x, control = ctl))) else NA_real_
  pre <- sprintf("p%d_v%s_s%02d_", row$p, format(row$variance), row$seed)
  fs <- list.files(G, paste0("^", gsub("\\.", "\\\\.", pre), ".*\\.rds$"), full.names = TRUE)
  out <- lapply(fs, function(f) { x <- readRDS(f); if (is.null(x$fit)) return(NULL)
    risk <- drop(Yva %*% x$fit$EF %*% x$fit$EBeta)
    data.frame(row, arm = x$metrics$arm, cal_slope_exact = sl(risk), cal_slope_default = x$metrics$cal_slope,
               risk_sd = sd(risk)) })
  rbind(do.call(rbind, out), data.frame(row, arm = "oracle", cal_slope_exact = sl(d$oracle),
        cal_slope_default = NA, risk_sd = sd(d$oracle)))
}
cl <- parallel::makePSOCKcluster(8)
parallel::clusterExport(cl, c("G", "S9_OUT", "s9_load", "s9_zscore", "one"))
res <- do.call(rbind, parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), one))
parallel::stopCluster(cl)
write.csv(res, file.path(G, "cal_slope_exact.csv"), row.names = FALSE)
agg <- aggregate(cbind(cal_slope_exact, cal_slope_default, risk_sd) ~ p + variance + arm, res, median, na.action = na.pass)
print(agg[order(agg$p, agg$variance, agg$arm), ], digits = 3)
