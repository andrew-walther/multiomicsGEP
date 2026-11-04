# ============================================================
# Script: run_yfb_survival_weight_p2000.R
# Purpose: Attempt S-9. The survival-weighted YFB fit, fit_ssbmf_yfb(canonicalize = TRUE, w_S),
#          on the corrected S-7 grid at p = 50 and 2,000, over a registered range of w_S.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: survival, ebnm, parallel, yaml; code/fit_ssbmf_yfb.R,
#               code/predict_ssbmf_yfb.R, code/concordance_ci.R,
#               results/benchmark_sim/run_yfb_supervision_decision.R (data generator,
#               gamma solver; its main block does not run when sourced)
# Run: Rscript results/benchmark_sim/run_yfb_survival_weight_p2000.R [cores] [seeds] [mode]
#      mode = "grid" (default) or "smoke" (one off-grid seed, writes to smoke/)
# Registration: docs/plans/2026-11-03-survival-weight-p2000-plan.md and the S-9
#               entry in docs/plans/2026-10-16-round2-attempts-log.md.
# ============================================================
# Parallelization: PSOCK workers, one dataset per work item, one BLAS thread
# each. A worker holds two 400 x 2,000 matrices and up to 42 fits' scratch;
# about 1 GB per worker is ample. No patient data are loaded.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)

S9_OUT  <- "results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04"
S7_DIR  <- "results/benchmark_sim/outputs/yfb_supervision_2026-10-30/attempt7"
S9_FIXED_WS <- c(1, 10, 100, 1000, 10000)

# Setup ----

#' Load the fitter, prediction, concordance and the S-7 generator
#' @return NULL, invisibly; functions are defined in the calling environment.
#' @examples
#' s9_load()
s9_load <- function() {
  suppressPackageStartupMessages({ library(survival); library(ebnm) })
  source("results/benchmark_sim/run_yfb_supervision_decision.R")  # sources fitters too
  multimodal_yfb_use_cpp()  # compiled sweep, as in the S-7 workers
  invisible(NULL)
}

# Preprocessing ----

#' Signed per-cohort z-score, the input fit_ssbmf_yfb() was written for
#'
#' yfb_decision_data() returns X/sd within the cohort (nonnegative, uncentered).
#' Centering that matrix gives the per-cohort z-score, because its columns
#' already have unit SD. The planted direction (F/sd, unit norm) and the
#' generated risk Yc %*% truth * 3.3 are therefore unchanged.
#' @param Y Nonnegative X/sd matrix for one cohort (patients x genes).
#' @return Column-centered, unit-SD matrix.
#' @examples
#' s9_zscore(matrix(abs(rnorm(20)), 5))
s9_zscore <- function(Y) scale(Y)

#' Planted program's variance share under each preprocessing
#'
#' Centered share: ||C(S3)||^2 / ||C(Y)||^2, the share seen by the z-scored
#' fitter (equal to the S-7 prepared_share). Uncentered share:
#' ||S3||^2 / ||Y||^2, the share of the raw sum of squares that the
#' nonnegative uncentered input carries, where the column means dominate.
#' @param d Output of yfb_decision_data(truth_mode = "shared").
#' @return One-row data frame of the two shares.
#' @examples
#' # s9_shares(yfb_decision_data(1, 50, .01, "shared"))
s9_shares <- function(d) {
  S3 <- d$training_component3
  data.frame(share_centered   = sum(scale(S3, scale = FALSE)^2) / sum(scale(d$Y, scale = FALSE)^2),
             share_uncentered = sum(S3^2) / sum(d$Y^2))
}

# Fitting ----

#' Fit the canonical survival-weighted YFB model at one weight
#'
#' Exactly the S-7 canonical_w1 call (point-Laplace L/F, per-program Normal
#' beta, K, sweeps, tolerance and alpha from globals) except the weight and
#' the input. The F update precision is A_j = tau_j sum_i E[l_ik^2]
#' + w_S E[beta_k^2] sum_i w_i y_ij^2.
#' @param Y Signed z-scored training matrix.
#' @param time,event Training outcomes.
#' @param w_S Nonnegative survival weight.
#' @param cfg Parsed config/globals.yml.
#' @return fit_ssbmf_yfb() result.
#' @examples
#' # s9_fit(Y, time, event, 10, cfg)
s9_fit <- function(Y, time, event, w_S, cfg) {
  fit_ssbmf_yfb(Y, time, event, K = cfg$benchmark$k_merged_yfb_desurv,
                max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
                w_S = w_S, alpha = cfg$benchmark$alpha, verbose = FALSE,
                trace_supervision = TRUE, sign_correction = FALSE)
}

#' Per-gene loading precision terms at a fitted state
#'
#' A_gen[j,k] = tau_j sum_i E[l_ik^2]; A_surv[j,k] = E[beta_k^2] sum_i w_i y_ij^2,
#' with Cox working weights w at eta = Y E[F] E[beta] (canonical: EF_norms = 1).
#' These are the unweighted terms of the canonical F update.
#' @param fit fit_ssbmf_yfb() result with canonicalize = TRUE.
#' @param Y,time,event Training input and outcomes used for the fit.
#' @return list(A_gen, A_surv, live) with p x K matrices and a K logical.
#' @examples
#' # s9_precision_terms(fit, Y, time, event)
s9_precision_terms <- function(fit, Y, time, event) {
  eta <- drop(Y %*% fit$EF %*% fit$EBeta)
  w   <- calc_cox_taylor_yf(eta, time, event)$w
  list(A_gen  = outer(fit$Tau, colSums(fit$EL2)),
       A_surv = outer(colSums(Y^2 * w), fit$EBeta2),
       live   = colSums(fit$EF^2) > 0 & colSums(fit$EL2) > 0)
}

#' Weight giving a median survival share of 0.2 at the paired w_S = 0 fit
#'
#' Solves median_{j,k live}( w A_surv / (A_gen + w A_surv) ) = 0.2 with the
#' S-2 solver (no damping: rho = 1). Computed before any supervised fit and
#' never revised.
#' @param fit0 Canonical w_S = 0 fit on the same training data.
#' @param Y,time,event Its training input and outcomes.
#' @param target Median share target.
#' @return Positive scalar weight.
#' @examples
#' # s9_share_weight(fit0, Y, time, event)
s9_share_weight <- function(fit0, Y, time, event, target = 0.2) {
  pt <- s9_precision_terms(fit0, Y, time, event)
  yfb_calibrate_gamma(pt$A_gen, pt$A_surv, pt$live, target, rho = 1)
}

# Scoring ----

#' Score one canonical fit against the planted program
#'
#' Recovered: best |r(truth, EF_k)| >= .7, survival-active (|E beta| >= 1.96 SD,
#' SD > 0) and correct sign. False active: active programs other than the
#' match. Share: final in-fit w_S A_surv / A_F median over genes, for the match
#' and median over live programs. Calibration: validation Cox slope of the
#' risk score (1 = calibrated, < 1 = overconfident). beta*: the coefficient of
#' the true training predictor projected on the matched Z_k; coverage asks
#' whether E beta_k +- 1.96 SD contains it.
#' @param fit Fit or error condition.
#' @param d Generated data; Ytr, Yva the z-scored inputs.
#' @return One-row data frame; an error gives the registered failure penalty.
#' @examples
#' # s9_score(fit, d, Ytr, Yva)
s9_score <- function(fit, d, Ytr, Yva) {
  if (inherits(fit, "error"))
    return(data.frame(cor_planted = NA_real_, recovered = FALSE, false_active = 7L,
                      converged = FALSE, c_validation = NA_real_, cal_slope = NA_real_,
                      beta_match = NA_real_, beta_sd_match = NA_real_, beta_star = NA_real_,
                      beta_covered = NA, share_match = NA_real_, share_live = NA_real_,
                      sweeps = NA_integer_, decreases = NA_integer_, max_beta = NA_real_,
                      recon_rel = NA_real_, min_tau = NA_real_, lost_variance = NA,
                      error = conditionMessage(fit)))
  r <- suppressWarnings(cor(d$truth, fit$EF)); r[!is.finite(r)] <- 0
  k <- which.max(abs(r)); matched <- abs(r[k]) >= .7
  v <- fit$EBeta2 - fit$EBeta^2
  sd <- sqrt(pmax(v, 0))
  active <- abs(fit$EBeta) >= 1.96 * sd & sd > 0
  live <- colSums(fit$EF^2) > 0
  risk <- predict_ssbmf_yfb(Yva, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
  slope <- if (sd(risk) > 0) unname(coef(survival::coxph(survival::Surv(d$tv, d$ev) ~ risk))) else NA_real_
  Zk <- drop(Ytr %*% fit$EF[, k])
  bstar <- if (var(Zk) > 0) cov(d$training_eta, Zk) / var(Zk) else NA_real_
  tr <- fit$history$trace
  sh <- tail(tr$shareA_median, 1)
  elbo <- fit$history$elbo_full
  data.frame(cor_planted = abs(r[k]),
             recovered = matched && active[k] && fit$EBeta[k] * r[k] > 0,
             false_active = sum(active[setdiff(seq_along(active), if (matched) k else integer())]),
             converged = isTRUE(fit$history$converged),
             c_validation = frozen_reverse_cindex(risk, d$tv, d$ev),
             cal_slope = slope,
             beta_match = fit$EBeta[k], beta_sd_match = sd[k], beta_star = bstar,
             beta_covered = abs(fit$EBeta[k] - bstar) <= 1.96 * sd[k],
             share_match = sh[k], share_live = median(sh[live], na.rm = TRUE),
             sweeps = fit$history$n_iter,
             decreases = sum(diff(elbo) < -1e-8 * abs(head(elbo, -1)), na.rm = TRUE),
             max_beta = max(abs(fit$EBeta)),
             recon_rel = sum((Ytr - fit$EL %*% t(fit$EF))^2) / sum(Ytr^2),
             min_tau = min(fit$Tau),
             lost_variance = any(live & v <= 0),
             error = "")
}

#' Training five-fold cross-validated C for one arm
#'
#' Folds are stratified by event and fixed per dataset. Each fold z-scores its
#' training and held-out parts within themselves (the per-cohort convention).
#' The share0.2 arm is recalibrated inside each fold from that fold's own
#' w_S = 0 fit, so no held-out patient informs its weight.
#' @param Yraw Nonnegative X/sd training matrix (before z-scoring).
#' @param time,event Training outcomes.
#' @param folds Integer fold labels.
#' @param arm_weight Fixed weight, or NA for the share0.2 arm.
#' @param cfg Parsed globals.
#' @return Mean held-out C over folds (NA folds dropped) and the fold count used.
#' @examples
#' # s9_cv_c(Y, time, event, folds, 10, cfg)
s9_cv_c <- function(Yraw, time, event, folds, arm_weight, cfg) {
  cs <- vapply(sort(unique(folds)), function(f) {
    tr <- folds != f
    Yt <- s9_zscore(Yraw[tr, ]); Yh <- s9_zscore(Yraw[!tr, ])
    ws <- arm_weight
    if (is.na(ws)) {
      f0 <- tryCatch(s9_fit(Yt, time[tr], event[tr], 0, cfg), error = function(e) e)
      if (inherits(f0, "error")) return(NA_real_)
      ws <- tryCatch(s9_share_weight(f0, Yt, time[tr], event[tr]), error = function(e) NA_real_)
      if (!is.finite(ws)) return(NA_real_)
    }
    fit <- tryCatch(s9_fit(Yt, time[tr], event[tr], ws, cfg), error = function(e) e)
    if (inherits(fit, "error")) return(NA_real_)
    risk <- predict_ssbmf_yfb(Yh, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
    frozen_reverse_cindex(risk, time[!tr], event[!tr])
  }, numeric(1))
  data.frame(cv_c = mean(cs, na.rm = TRUE), cv_folds_ok = sum(is.finite(cs)))
}

# One dataset ----

#' Fit every registered S-9 arm on one dataset and checkpoint each result
#'
#' Arms: canonical w_S = 0 (same-fitter reference), w_S in S9_FIXED_WS,
#' w_S = p/50 (equal to w_S = 1 at p = 50, so that fit is reused), and share0.2.
#' Each arm also gets its training five-fold CV C. The S-7 unsupervised
#' reference is read from its saved checkpoint, not refitted.
#' @param row Design row (seed, p, variance).
#' @param out Output directory; existing checkpoints are never overwritten.
#' @return Status string.
#' @examples
#' # s9_task(data.frame(seed = 1, p = 50, variance = .01), S9_OUT)
s9_task <- function(row, out) {
  cfg <- yaml::read_yaml("config/globals.yml")
  d <- yfb_decision_data(row$seed, row$p, row$variance, truth_mode = "shared")
  Ytr <- s9_zscore(d$Y); Yva <- s9_zscore(d$Yv)
  prefix <- sprintf("p%d_v%s_s%02d", row$p, row$variance, row$seed)
  set.seed(row$seed + 7000L)
  folds <- integer(length(d$event))
  for (e in 0:1) { ii <- which(d$event == e); folds[ii] <- sample(rep_len(seq_len(cfg$cavi$n_cv_folds), length(ii))) }
  shares <- s9_shares(d)
  run <- function(arm, w) {
    path <- file.path(out, paste0(prefix, "_", arm, ".rds"))
    if (file.exists(path)) return(readRDS(path))
    t0 <- proc.time()[["elapsed"]]
    warn <- character()
    fit <- tryCatch(withCallingHandlers(s9_fit(Ytr, d$time, d$event, w, cfg),
                     warning = function(x) { warn <<- c(warn, conditionMessage(x)); invokeRestart("muffleWarning") }),
                    error = function(e) e)
    m <- s9_score(fit, d, Ytr, Yva)
    cv <- if (arm == "w0") data.frame(cv_c = NA_real_, cv_folds_ok = NA_integer_)
          else s9_cv_c(d$Y, d$time, d$event, folds, if (arm == "share0.2") NA_real_ else w, cfg)
    rec <- list(metrics = cbind(row, arm = arm, w_S = w, shares, m, cv,
                                elapsed = proc.time()[["elapsed"]] - t0,
                                warnings = paste(unique(warn), collapse = " | ")),
                fit = if (inherits(fit, "error")) NULL else fit)
    saveRDS(rec, path)
    rec
  }
  r0 <- run("w0", 0)
  ws02 <- if (is.null(r0$fit)) NA_real_ else
    tryCatch(s9_share_weight(r0$fit, Ytr, d$time, d$event), error = function(e) NA_real_)
  for (w in S9_FIXED_WS) run(paste0("w", w), w)
  if (row$p != 50L) run("p_over_50", row$p / 50)
  if (is.finite(ws02)) run("share0.2", ws02) else
    saveRDS(list(metrics = cbind(row, arm = "share0.2", w_S = NA_real_, shares,
                                 s9_score(simpleError("share0.2 weight could not be calibrated"), d, Ytr, Yva),
                                 cv_c = NA_real_, cv_folds_ok = NA_integer_, elapsed = NA_real_, warnings = ""),
                 fit = NULL), file.path(out, paste0(prefix, "_share0.2.rds")))
  prefix
}

# Main ----
if (sys.nframe() == 0L) {
  args  <- commandArgs(TRUE)
  cores <- if (length(args)) as.integer(args[1]) else 8L
  seeds <- if (length(args) > 1) as.integer(args[2]) else 10L
  mode  <- if (length(args) > 2) args[3] else "grid"
  s9_load()
  out <- file.path(S9_OUT, if (mode == "smoke") "smoke" else "grid")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  design <- if (mode == "smoke") expand.grid(seed = 99L, p = c(50L, 2000L), variance = .01)
            else expand.grid(seed = seq_len(seeds), p = c(50L, 2000L), variance = c(.0025, .01, .05))
  write.csv(design, file.path(out, "design.csv"), row.names = FALSE)
  writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time())), file.path(out, "provenance.txt"))
  cl <- parallel::makePSOCKcluster(cores, outfile = file.path(out, "workers.log"))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, c("S9_OUT", "S7_DIR", "S9_FIXED_WS", "s9_load", "s9_zscore", "s9_shares",
                                "s9_fit", "s9_precision_terms", "s9_share_weight", "s9_score", "s9_cv_c",
                                "s9_task", "out"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
                                          VECLIB_MAXIMUM_THREADS = 1); s9_load(); NULL })
  status <- parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(row)
    tryCatch(s9_task(row, out), error = function(e) paste("TASK ERROR", conditionMessage(e))),
    chunk.size = 1L)
  writeLines(unlist(status), file.path(out, "task_status.txt"))
}
