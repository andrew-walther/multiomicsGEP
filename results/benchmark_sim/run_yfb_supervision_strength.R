# ============================================================
# Script: run_yfb_supervision_strength.R
# Purpose: Attempt S-13. The survival-block-weighted YFB fit (A1: untempered beta) and its
#          joint-q(beta) version (A2) over a survival-weight grid, with the weight chosen by
#          five-fold training CV partial log-likelihood, at p = 2,000 (S-9 cells, a weak-signal
#          cell matched to the real data, and permuted-outcome nulls).
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: survival, ebnm, parallel, yaml; code/fit_ssbmf_yfb.R (via the S-9 runner),
#               results/benchmark_sim/run_yfb_survival_weight_p2000.R (s9_load, s9_zscore,
#               s9_score), results/benchmark_sim/run_yfb_supervision_decision.R (generator)
# Run: Rscript results/benchmark_sim/run_yfb_supervision_strength.R [cores] [mode]
#      mode = "grid" (default), "null", or "smoke" (one task per arm, off-grid seed)
# Registration: docs/plans/2026-11-05-supervision-strength-plan.md and the S-13 entry in
#               docs/plans/2026-10-16-round2-attempts-log.md.
# ============================================================
# Parallelization: PSOCK workers, one (dataset, arm, weight) task per work item (a full fit
# and five fold fits), one BLAS thread each, one checkpoint per task so a rerun skips finished
# work. A worker holds two 400 x 2,000 matrices and one fit at a time (well under 1 GB).
# No patient data are loaded.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)

S13_OUT     <- "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05"
S13_WEIGHTS <- c(1, 2, 5, 10, 20, 50)
S13_ARMS    <- c(A1 = "meanfield", A2 = "joint")
S13_WEAK_SCALE <- 0.045   # oracle validation C ~ 0.68 (calibrated on seeds 101-112, generator only)
S13_NFOLD   <- 5L

# Data ----

#' One S-13 dataset: an S-9 cell, or the weak-signal cell with the true predictor scaled
#'
#' The weak cell keeps the S-9 expression data (variance 1%) and redraws survival from the
#' same Weibull model with eta scaled by S13_WEAK_SCALE (seed + 5000), censoring calibrated
#' to 30% as in the generator. For the null, training outcomes are permuted (seed + 9000).
#' @param seed Dataset seed. @param cell "v0.0025", "v0.01", "v0.05" or "weak".
#' @param null Logical: permute training outcomes.
#' @return yfb_decision_data()-style list with z-scored Ytr, Yva.
#' @examples
#' # s13_data(1, "weak")
s13_data <- function(seed, cell, null = FALSE) {
  v <- if (cell == "weak") .01 else as.numeric(sub("v", "", cell))
  d <- yfb_decision_data(seed, 2000L, v, truth_mode = "shared")
  if (cell == "weak") {
    set.seed(seed + 5000L)
    draw <- function(eta) {
      et <- (-log(runif(length(eta))) / (.01 * exp(eta)))^(1 / 1.5); cb <- rexp(length(eta))
      cs <- multimodal_yfb_calibrate_censoring(et, cb, .3)
      list(time = pmin(et, cb * cs), event = as.numeric(et <= cb * cs))
    }
    d$training_eta <- d$training_eta * S13_WEAK_SCALE; d$oracle <- d$oracle * S13_WEAK_SCALE
    a <- draw(d$training_eta); b <- draw(d$oracle)
    d$time <- a$time; d$event <- a$event; d$tv <- b$time; d$ev <- b$event
  }
  if (null) {
    set.seed(seed + 9000L); perm <- sample(length(d$time))
    d$time <- d$time[perm]; d$event <- d$event[perm]
  }
  d$Ytr <- s9_zscore(d$Y); d$Yva <- s9_zscore(d$Yv)
  d
}

#' Five folds stratified by event, fixed per dataset
#' @param event Event indicator. @param seed Dataset seed. @return Fold vector.
#' @examples
#' s13_folds(rbinom(20, 1, .5), 1)
s13_folds <- function(event, seed) {
  set.seed(seed + 7000L); f <- integer(length(event))
  for (e in unique(event)) { i <- which(event == e); f[i] <- sample(rep(seq_len(S13_NFOLD), length.out = length(i))) }
  f
}

# Fitting ----

#' Fit the S-13 model: canonical YFB with the survival-block weight
#' @param Y z-scored training matrix. @param time,event Outcomes. @param w_S Weight.
#' @param beta_posterior "meanfield" (A1) or "joint" (A2). @param cfg Config.
#' @return fit_ssbmf_yfb() result.
#' @examples
#' # s13_fit(Y, time, event, 5, "joint", cfg)
s13_fit <- function(Y, time, event, w_S, beta_posterior, cfg) {
  fit_ssbmf_yfb(Y, time, event, K = cfg$benchmark$k_merged_yfb_desurv,
                max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
                w_S = w_S, survival_block_weight = TRUE, beta_posterior = beta_posterior,
                alpha = cfg$benchmark$alpha, verbose = FALSE, trace_supervision = TRUE,
                sign_correction = FALSE)
}

#' Per-program roles at a fitted state (Gemini guide, issue 5)
#'
#' PVE_k = ||L_k F_k'||^2 / ||Y||^2; active = |E beta| >= 1.96 SD. Roles: prognostic program
#' (active, PVE >= 1%), survival-only direction (active, PVE < 1%), genomics-only program
#' (inactive, PVE >= 1%), inactive.
#' @param fit Fit. @param Y Training matrix. @return Data frame, one row per program.
#' @examples
#' # s13_roles(fit, Y)
s13_roles <- function(fit, Y) {
  pve <- colSums(fit$EL^2) * colSums(fit$EF^2) / sum(Y^2)
  sd <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  act <- abs(fit$EBeta) >= 1.96 * sd & sd > 0
  role <- ifelse(act, ifelse(pve >= .01, "prognostic_program", "survival_only_direction"),
                 ifelse(pve >= .01, "genomics_only_program", "inactive"))
  data.frame(k = seq_along(pve), pve = pve, beta = fit$EBeta, beta_sd = sd, active = act, role = role)
}

#' One task: full fit (scored against the planted truth) and five-fold CV
#' @param task Row with seed, cell, arm, w_S, null. @param cfg Config.
#' @return Status string; writes a checkpoint.
#' @examples
#' # s13_task(data.frame(seed = 1, cell = "weak", arm = "A1", w_S = 1, null = FALSE), cfg)
s13_task <- function(task, cfg) {
  key <- sprintf("%s_%s_s%02d_%s_w%g", if (task$null) "null" else "sig", task$cell, task$seed,
                 task$arm, task$w_S)
  path <- file.path(S13_OUT, "fits", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  d <- s13_data(task$seed, task$cell, task$null); bp <- S13_ARMS[[task$arm]]
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(s13_fit(d$Ytr, d$time, d$event, task$w_S, bp, cfg),
                    warning = function(w) invokeRestart("muffleWarning")), error = function(e) e)
  sc <- s9_score(fit, d, d$Ytr, d$Yva)
  row <- cbind(task, sc, fit_sec = proc.time()[["elapsed"]] - t0)
  roles <- if (!inherits(fit, "error")) s13_roles(fit, d$Ytr) else NULL
  # matched program's role
  row$match_role <- NA_character_; row$match_pve <- NA_real_
  if (!is.null(roles)) {
    r <- suppressWarnings(abs(cor(d$truth, fit$EF))); r[!is.finite(r)] <- 0; k <- which.max(r)
    row$match_role <- roles$role[k]; row$match_pve <- roles$pve[k]
    row$n_survival_only <- sum(roles$role == "survival_only_direction")
    row$n_prognostic <- sum(roles$role == "prognostic_program")
  }
  # Five-fold CV: Verweij-van Houwelingen partial log-likelihood and held-out C
  cv_pl <- 0; cv_c <- numeric(0); folds_ok <- 0L; folds_conv <- 0L
  if (!task$null) {
    fo <- s13_folds(d$event, task$seed)
    for (f in seq_len(S13_NFOLD)) {
      tr <- fo != f; te <- fo == f
      ff <- tryCatch(withCallingHandlers(s13_fit(d$Ytr[tr, ], d$time[tr], d$event[tr], task$w_S, bp, cfg),
                       warning = function(w) invokeRestart("muffleWarning")), error = function(e) e)
      if (inherits(ff, "error")) { cv_pl <- NA_real_; next }
      eta <- drop(d$Ytr %*% ff$EF %*% ff$EBeta)
      cv_pl <- cv_pl + calc_cox_taylor_yf(eta, d$time, d$event)$logPL -
                       calc_cox_taylor_yf(eta[tr], d$time[tr], d$event[tr])$logPL
      if (sd(eta[te]) > 1e-8) cv_c <- c(cv_c, frozen_reverse_cindex(eta[te], d$time[te], d$event[te]))
      folds_ok <- folds_ok + 1L; folds_conv <- folds_conv + as.integer(ff$history$converged)
    }
  }
  row$cv_pl <- if (task$null) NA_real_ else cv_pl
  row$cv_c <- if (length(cv_c)) mean(cv_c) else NA_real_
  row$cv_folds_ok <- folds_ok; row$cv_folds_converged <- folds_conv
  row$elapsed <- proc.time()[["elapsed"]] - t0
  saveRDS(list(metrics = row, roles = roles,
               state = if (!inherits(fit, "error")) list(EF = fit$EF, EBeta = fit$EBeta,
                 EBeta2 = fit$EBeta2, EBeta_cov = fit$EBeta_cov, L_norm = colSums(fit$EL^2),
                 elbo = fit$history$elbo_full[seq_len(fit$history$n_iter)]) else NULL), path)
  sprintf("%s ok %.0fs", key, row$elapsed)
}

#' Task design: A1 tasks first, then A2; largest weights first within an arm
#' @param mode "grid", "null" or "smoke". @return Data frame of tasks.
#' @examples
#' s13_design("smoke")
s13_design <- function(mode) {
  if (mode == "smoke")
    return(data.frame(seed = 99L, cell = "weak", arm = c("A1", "A2"), w_S = 5, null = FALSE))
  cells <- if (mode == "null") c("v0.0025", "v0.01", "v0.05") else c("weak", "v0.0025", "v0.01", "v0.05")
  d <- expand.grid(seed = 1:10, cell = cells, arm = names(S13_ARMS), w_S = S13_WEIGHTS,
                   stringsAsFactors = FALSE)
  d$null <- mode == "null"
  d[order(d$arm, -d$w_S, d$cell, d$seed), ]
}

# Two-step comparison (registered criterion 5) ----

#' Two-step fit: canonical w_S = 0 loadings, then Cox on the live training program scores
#'
#' Recovered: best |r(truth, F_k)| >= .7, Cox Wald |z| >= 1.96 and correct sign. Validation C and
#' the exact-times calibration slope of Yva E[F] b.
#' @param seed Dataset seed. @param cell Cell. @param cfg Config.
#' @return One-row data frame.
#' @examples
#' # s13_twostep(1, "weak", cfg)
s13_twostep <- function(seed, cell, cfg) {
  d <- s13_data(seed, cell)
  f0 <- fit_ssbmf_yfb(d$Ytr, d$time, d$event, K = cfg$benchmark$k_merged_yfb_desurv,
                      max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol, prior_LF = "point_laplace",
                      prior_beta = "normal", canonicalize = TRUE, w_S = 0,
                      alpha = cfg$benchmark$alpha, verbose = FALSE, sign_correction = FALSE)
  live <- which(colSums(f0$EF^2) > 1e-8)
  Z <- d$Ytr %*% f0$EF[, live, drop = FALSE]
  cf <- survival::coxph(survival::Surv(d$time, d$event) ~ Z, control = survival::coxph.control(timefix = FALSE))
  b <- rep(0, ncol(f0$EF)); b[live] <- ifelse(is.finite(coef(cf)), coef(cf), 0)
  zst <- rep(0, ncol(f0$EF)); zst[live] <- coef(cf) / sqrt(diag(as.matrix(vcov(cf))))
  r <- suppressWarnings(cor(d$truth, f0$EF)); r[!is.finite(r)] <- 0; k <- which.max(abs(r))
  risk <- drop(d$Yva %*% f0$EF %*% b)
  data.frame(seed = seed, cell = cell, cor_planted = abs(r[k]),
             recovered = abs(r[k]) >= .7 && abs(zst[k]) >= 1.96 && b[k] * r[k] > 0,
             false_active = sum(abs(zst[-k]) >= 1.96, na.rm = TRUE),
             c_validation = frozen_reverse_cindex(risk, d$tv, d$ev),
             cal_slope = unname(coef(survival::coxph(survival::Surv(d$tv, d$ev) ~ risk,
                               control = survival::coxph.control(timefix = FALSE)))),
             converged = f0$history$converged)
}

# Main ----
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  cores <- if (length(args) >= 1) as.integer(args[1]) else 9L
  mode  <- if (length(args) >= 2) args[2] else "grid"
  if (mode == "twostep") {
    source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
    source("code/concordance_ci.R"); cfg <- yaml::read_yaml("config/globals.yml")
    des <- expand.grid(seed = 1:10, cell = c("weak", "v0.0025", "v0.01", "v0.05"), stringsAsFactors = FALSE)
    res <- do.call(rbind, lapply(seq_len(nrow(des)), function(i) s13_twostep(des$seed[i], des$cell[i], cfg)))
    write.csv(res, file.path(S13_OUT, "twostep.csv"), row.names = FALSE); print(res)
    quit(save = "no")
  }
  if (mode == "smoke") S13_OUT <- file.path(S13_OUT, "smoke")
  dir.create(file.path(S13_OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()), mode),
             file.path(S13_OUT, sprintf("provenance_%s.txt", mode)))
  design <- s13_design(mode)
  cl <- parallel::makePSOCKcluster(cores, outfile = file.path(S13_OUT, sprintf("workers_%s.log", mode)))
  parallel::clusterExport(cl, c("S13_OUT", "S13_WEIGHTS", "S13_ARMS", "S13_WEAK_SCALE", "S13_NFOLD",
                                "s13_data", "s13_folds", "s13_fit", "s13_roles", "s13_task"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1)
    source("results/benchmark_sim/run_yfb_survival_weight_p2000.R"); s9_load()
    source("code/concordance_ci.R"); cfg <- yaml::read_yaml("config/globals.yml"); NULL })
  st <- parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(task)
    tryCatch(s13_task(task, cfg), error = function(e) paste("TASK ERROR", conditionMessage(e))),
    chunk.size = 1L)
  parallel::stopCluster(cl)
  writeLines(unlist(st), file.path(S13_OUT, sprintf("task_status_%s.txt", mode)))
}
