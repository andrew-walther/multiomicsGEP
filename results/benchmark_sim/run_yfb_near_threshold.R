# ============================================================
# Script: run_yfb_near_threshold.R
# Purpose: Attempt S-14. Supervised YFB (A1, survival-block weight) against the two-step fit for
#          planted programs near the detection threshold (shares 0.05-0.2% at p = 2,000).
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: survival, ebnm, parallel, yaml; results/benchmark_sim/run_yfb_survival_weight_p2000.R
# Run: Rscript results/benchmark_sim/run_yfb_near_threshold.R [cores]
# Registration: docs/plans/2026-11-05-near-threshold-plan.md
# ============================================================
# Parallelization: PSOCK, one (dataset, arm) fit per item, one BLAS thread, one checkpoint per fit.
# Synthetic data only.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
S14_OUT <- "results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05"

#' One S-14 dataset
#' @param seed Seed. @param share Planted share. @param strength "strong" or "realistic".
#' @return Data list with z-scored Ytr, Yva and outcomes.
#' @examples
#' # s14_data(1, .001, "realistic")
s14_data <- function(seed, share, strength) {
  d <- yfb_decision_data(seed, 2000L, share, truth_mode = "shared")
  if (strength == "realistic") {
    f <- 1.1 / sd(d$training_eta)
    set.seed(seed + 5000L)
    draw <- function(eta) {
      et <- (-log(runif(length(eta))) / (.01 * exp(eta)))^(1 / 1.5); cb <- rexp(length(eta))
      cs <- multimodal_yfb_calibrate_censoring(et, cb, .3)
      list(time = pmin(et, cb * cs), event = as.numeric(et <= cb * cs))
    }
    d$training_eta <- d$training_eta * f; d$oracle <- d$oracle * f
    a <- draw(d$training_eta); b <- draw(d$oracle)
    d$time <- a$time; d$event <- a$event; d$tv <- b$time; d$ev <- b$event
  }
  d$Ytr <- s9_zscore(d$Y); d$Yva <- s9_zscore(d$Yv)
  d
}

#' Fit and score one arm
#' @param task Row: seed, share, strength, arm. @param cfg Config. @return Status string.
#' @examples
#' # s14_task(task, cfg)
s14_task <- function(task, cfg) {
  key <- sprintf("%s_sh%g_s%02d_%s", task$strength, task$share, task$seed, task$arm)
  path <- file.path(S14_OUT, "fits", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  d <- s14_data(task$seed, task$share, task$strength)
  w <- switch(task$arm, twostep = 0, A1_w1 = 1, A1_w5 = 5, A2_w1 = 1)
  t0 <- proc.time()[["elapsed"]]
  f <- fit_ssbmf_yfb(d$Ytr, d$time, d$event, K = cfg$benchmark$k_merged_yfb_desurv,
                     max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol, prior_LF = "point_laplace",
                     prior_beta = "normal", canonicalize = TRUE, w_S = w, survival_block_weight = TRUE,
                     beta_posterior = if (task$arm == "A2_w1") "joint" else "meanfield",
                     alpha = cfg$benchmark$alpha, verbose = FALSE)
  r <- suppressWarnings(cor(d$truth, f$EF)); r[!is.finite(r)] <- 0; k <- which.max(abs(r))
  if (task$arm == "twostep") {
    live <- which(colSums(f$EF^2) > 1e-8); Z <- d$Ytr %*% f$EF[, live, drop = FALSE]
    cf <- survival::coxph(survival::Surv(d$time, d$event) ~ Z, control = survival::coxph.control(timefix = FALSE))
    b <- zs <- rep(0, ncol(f$EF)); b[live] <- ifelse(is.finite(coef(cf)), coef(cf), 0)
    zs[live] <- coef(cf) / sqrt(diag(as.matrix(vcov(cf))))
  } else {
    b <- f$EBeta; sdb <- sqrt(pmax(f$EBeta2 - b^2, 0)); zs <- ifelse(sdb > 0, b / sdb, 0)
  }
  act <- abs(zs) >= 1.96
  risk <- drop(d$Yva %*% f$EF %*% b)
  row <- data.frame(task, abs_r = abs(r[k]), recovered = abs(r[k]) >= .7 && act[k] && b[k] * r[k] > 0,
    false_active = sum(act[-k]),
    c_validation = if (sd(risk) > 1e-8) frozen_reverse_cindex(risk, d$tv, d$ev) else NA_real_,
    cal_slope = if (sd(risk) > 1e-8) unname(coef(survival::coxph(survival::Surv(d$tv, d$ev) ~ risk,
      control = survival::coxph.control(timefix = FALSE)))) else NA_real_,
    oracle_c = frozen_reverse_cindex(d$oracle, d$tv, d$ev),
    converged = f$history$converged, sweeps = f$history$n_iter,
    fit_sec = proc.time()[["elapsed"]] - t0)
  saveRDS(list(metrics = row, EF = f$EF, beta = b), path)
  sprintf("%s ok %.0fs", key, row$fit_sec)
}

if (sys.nframe() == 0L) {
  cores <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 1L
  dir.create(file.path(S14_OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time())),
             file.path(S14_OUT, "provenance.txt"))
  design <- expand.grid(seed = 1:10, share = c(.0005, .001, .002), strength = c("strong", "realistic"),
                        arm = c("twostep", "A1_w1", "A1_w5", "A2_w1"), stringsAsFactors = FALSE)
  design <- design[order(design$arm != "twostep" & design$arm != "A1_w1", design$strength, design$share), ]
  cl <- parallel::makePSOCKcluster(cores)
  parallel::clusterExport(cl, c("S14_OUT", "s14_data", "s14_task"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1); source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
    s9_load(); source("code/concordance_ci.R"); cfg <- yaml::read_yaml("config/globals.yml"); NULL })
  st <- parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(t)
    tryCatch(s14_task(t, cfg), error = function(e) paste("TASK ERROR", conditionMessage(e))), chunk.size = 1L)
  parallel::stopCluster(cl)
  writeLines(unlist(st), file.path(S14_OUT, "task_status.txt"))
  m <- do.call(rbind, lapply(list.files(file.path(S14_OUT, "fits"), full.names = TRUE), function(f) readRDS(f)$metrics))
  write.csv(m, file.path(S14_OUT, "metrics.csv"), row.names = FALSE)
}
