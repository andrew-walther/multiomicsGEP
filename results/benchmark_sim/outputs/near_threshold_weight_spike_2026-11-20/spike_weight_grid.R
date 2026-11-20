# ============================================================
# Script: spike_weight_grid.R
# Purpose: Refit supervised YFB at w_S = 10, 100, 1000 in the near-threshold test bed.
# Author: Andrew Walther
# Created: 2026-11-20
# Dependencies: parallel, yaml; results/benchmark_sim/run_yfb_near_threshold.R
# Run (from the repo root): Rscript results/benchmark_sim/outputs/near_threshold_weight_spike_2026-11-20/spike_weight_grid.R
# Needs the two-step checkpoints of the 2026-11-19 run (git-ignored .rds files).
# Findings: docs/findings/2026-11-20-supervision-influence-brainstorm-findings.md, Section 7.
# ============================================================
# Spike: supervised YFB at w_S in {10, 100, 1000} in the
# real-data-like and moderate near-threshold settings, K = 5. Records the usual metrics
# plus loading norm growth and loading similarity to the saved two-step fit.
# Synthetic data only. Reuses the committed runner's data generator and scoring.
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
OUT <- "results/benchmark_sim/outputs/near_threshold_weight_spike_2026-11-20"  # per-fit .rds checkpoints are git-ignored
REF <- "results/benchmark_sim/outputs/yfb_near_threshold_k5_2026-11-19/fits"
dir.create(OUT, showWarnings = FALSE)
design <- expand.grid(seed = 1:10, share = c(.0005, .001), strength = c("real_like", "realistic"),
                      w = c(10, 100, 1000), stringsAsFactors = FALSE)
task <- function(t) {
  key <- sprintf("%s_sh%g_s%02d_w%g", t$strength, t$share, t$seed, t$w)
  path <- file.path(OUT, paste0(key, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  d <- s14_data(t$seed, t$share, t$strength)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(fit_ssbmf_yfb(d$Ytr, d$time, d$event, K = 5L, max_iter = 300L, tol = cfg$cavi$tol,
                     prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE, w_S = t$w,
                     survival_block_weight = TRUE, beta_posterior = "meanfield",
                     alpha = cfg$benchmark$alpha, verbose = FALSE), error = function(e) e)
  if (inherits(f, "error")) { r <- data.frame(t, error = conditionMessage(f)); saveRDS(r, path); return(r) }
  r <- suppressWarnings(cor(d$truth, f$EF)); r[!is.finite(r)] <- 0; k <- which.max(abs(r))
  b <- f$EBeta; sdb <- sqrt(pmax(f$EBeta2 - b^2, 0)); zs <- ifelse(sdb > 0, b / sdb, 0); act <- abs(zs) >= 1.96
  risk <- drop(d$Yva %*% f$EF %*% b)
  two <- readRDS(file.path(REF, sprintf("%s_sh%g_s%02d_twostep.rds", t$strength, t$share, t$seed)))
  Fs <- f$EF[, colSums(f$EF^2) > 1e-8, drop = FALSE]; Ft <- two$EF[, colSums(two$EF^2) > 1e-8, drop = FALSE]
  sim <- if (ncol(Fs)) mean(apply(abs(suppressWarnings(cor(Fs, Ft))), 1, max, na.rm = TRUE)) else NA
  out <- data.frame(t, error = NA, abs_r = abs(r[k]), recovered = abs(r[k]) >= .7 && act[k] && b[k] * r[k] > 0,
    false_active = sum(act[-k]),
    c_validation = if (sd(risk) > 1e-8) frozen_reverse_cindex(risk, d$tv, d$ev) else NA_real_,
    c_two = two$metrics$c_validation, oracle_c = frozen_reverse_cindex(d$oracle, d$tv, d$ev),
    converged = f$history$converged, sweeps = f$history$n_iter, max_abs_beta = max(abs(b)),
    max_normF = max(sqrt(colSums(f$EF^2))), sim_two = sim, fit_sec = proc.time()[["elapsed"]] - t0)
  saveRDS(out, path); out
}
if (sys.nframe() == 0L) {
  cl <- parallel::makePSOCKcluster(8L)
  parallel::clusterExport(cl, c("OUT", "REF", "task"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
    suppressMessages({ source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
      source("code/concordance_ci.R"); source("results/benchmark_sim/run_yfb_near_threshold.R") })
    cfg <- yaml::read_yaml("config/globals.yml"); NULL })
  res <- parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(t)
    tryCatch(task(t), error = function(e) data.frame(t, error = conditionMessage(e))), chunk.size = 1L)
  parallel::stopCluster(cl)
  all <- do.call(rbind, lapply(res, function(x) { for (nm in setdiff(names(res[[1]]), names(x))) x[[nm]] <- NA; x[names(res[[1]])] }))
  write.csv(all, file.path(OUT, "spike_metrics.csv"), row.names = FALSE)
  cat("done", nrow(all), "rows\n")
}
