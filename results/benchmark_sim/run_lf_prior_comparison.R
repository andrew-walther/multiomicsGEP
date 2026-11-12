# ============================================================
# Script: run_lf_prior_comparison.R
# Purpose: Task A. Point-exponential vs point-Laplace priors on L and F for the
#          survival-weighted YFB fit, in simulation and by pooled training CV.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm, yaml, parallel; code/fit_ssbmf_yfb.R,
#               code/predict_ssbmf_yfb.R, code/concordance_ci.R, code/preprocess_desurv.R,
#               results/benchmark_sim/benchmark_helpers.R,
#               results/benchmark_sim/run_yfb_survival_weight_p2000.R (S-9 generator, fit, score),
#               results/benchmark_sim/run_yfb_near_threshold.R (S-14 near-threshold data)
# Run: Rscript results/benchmark_sim/run_lf_prior_comparison.R prepare
#      Rscript results/benchmark_sim/run_lf_prior_comparison.R smoke [cores]
#      Rscript results/benchmark_sim/run_lf_prior_comparison.R grid [cores]
#      Rscript results/benchmark_sim/run_lf_prior_comparison.R summary
# Registration: 2026-11-12 entry "A" in docs/plans/2026-10-16-round2-attempts-log.md;
#               task A of docs/plans/2026-11-09-week-plan.md, pooled per DECISIONS.md 2026-11-12.
# ============================================================
# Arms (all: canonicalize = TRUE, w_S = 1, Normal beta, K = 7, no pruning, 300 sweeps, tol 1e-5):
#   PL   z-score within cohort, point-Laplace L and F (the working model, S-11 call).
#   PE   X/sd within cohort (nonnegative), point-exponential L and F,
#        center_projection = TRUE with one pooled center.
#   RANK real data only (task F fixed arm): within-patient ranks over the D4 genes, each gene
#        divided by its within-cohort SD of the ranks (no centering), as PE otherwise.
# Training is pooled over TCGA and CPTAC with no strata_id. The primary CV metric is
# within-cohort concordance: only TCGA-TCGA and CPTAC-CPTAC pairs are compared, pooled by
# pair count (survival::concordance with strata()). X/sd keeps each cohort's mean level, so a
# score could partly encode cohort membership; within-cohort pairs cannot reward that.
#
# Parallelization: PSOCK workers (forked workers segfaulted on macOS with the compiled sweep),
# one BLAS/OpenMP thread each, one checkpoint per fit, so a rerun skips finished fits.
# Memory: a worker holds two 400 x 2,000 synthetic matrices or the 273 x 2,218 real matrix.
# Patient-level inputs, fits and risk scores go to data/ (git-ignored). Only aggregate CSVs
# are written to OUT.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
OUT       <- "results/benchmark_sim/outputs/lf_prior_2026-11-09"
PRIV      <- "data/lf_prior_2026-11-12"                       # patient-level checkpoints
REAL_IN   <- "data/lf_prior_real_2026-11-12.rds"              # prepared real inputs
S11_IN    <- "data/yfb_supervised_real_2026-11-04.rds"        # S-11 prepared data and folds
S11_FITS  <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/fits"
S9_GRID   <- "results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid"
S14_FITS  <- "results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05/fits"
BETA_ARMS <- c(B_PLAP = "shared point-Laplace", B_PN = "shared point-normal",
               B_NSH = "shared Normal (reference)")   # task B arms; Normal per-program = A's PL
SIM_CELLS <- data.frame(cell = c("s9_0.25pct", "s9_1pct", "s9_5pct", "s14_0.05pct"),
                        share = c(.0025, .01, .05, .0005),
                        source = c("s9", "s9", "s9", "s14"), stringsAsFactors = FALSE)

# Code loading ----

#' Source the fitter, generators, prediction and concordance helpers
#'
#' Also run inside every worker. The S-9 runner sources the fitter and the S-7 generator;
#' the S-14 runner defines s14_data(); neither main block runs when sourced from a function.
#' @return NULL, invisibly.
#' @examples
#' # lfp_load()
lfp_load <- function() {
  suppressPackageStartupMessages({ library(survival); library(ebnm); library(yaml) })
  source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
  s9_load()
  source("results/benchmark_sim/run_yfb_near_threshold.R")
  source("code/predict_ssbmf_yfb.R"); source("code/concordance_ci.R")
  invisible(NULL)
}

# Fitting ----

#' Fit one arm with the working-model settings
#'
#' PL is exactly the S-9/S-11 call. PE and RANK change only prior_LF and turn on
#' center_projection, so the Cox projection uses Y - 1 ybar' while the reconstruction uses
#' the nonnegative Y.
#' @param Y Training matrix in the arm's own input scale (patients x genes).
#' @param time,status Training outcomes.
#' @param arm "PL", "PE" or "RANK".
#' @param cfg Parsed config/globals.yml.
#' @return fit_ssbmf_yfb() result.
#' @examples
#' # lfp_fit(Y, time, status, "PE", cfg)
lfp_fit <- function(Y, time, status, arm, cfg) {
  pe <- arm %in% c("PE", "RANK")
  # Task B arms: the PL call with one beta prior shared across the K programs.
  beta <- switch(arm, B_PLAP = "point_laplace", B_PN = "point_normal", B_NSH = "normal", "normal")
  shared <- arm %in% names(BETA_ARMS)
  fit_ssbmf_yfb(Y, time, status, K = cfg$benchmark$k_merged_yfb_desurv,
                max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = if (pe) "point_exponential" else "point_laplace",
                prior_beta = beta, beta_prior_mode = if (shared) "shared" else "per_program",
                canonicalize = TRUE, w_S = 1,
                alpha = cfg$benchmark$alpha, verbose = FALSE, trace_supervision = TRUE,
                center_projection = pe)
}

#' Posterior inclusion probability of each program's beta under a point-normal prior
#'
#' At the final sweep each program's beta posterior is the EBNM posterior for the
#' pseudo-observation x_k ~ N(beta_k, s_k^2) under g = pi0 delta_0 + (1 - pi0) N(0, sigma^2):
#' PIP_k = (1 - pi0) N(x_k; 0, sigma^2 + s_k^2) /
#'         [pi0 N(x_k; 0, s_k^2) + (1 - pi0) N(x_k; 0, sigma^2 + s_k^2)].
#' Uses the prior in force during that sweep (prior_before). Flat programs (s = Inf) get NA.
#' @param fit fit_ssbmf_yfb() result with beta_prior_mode = "shared", trace_supervision = TRUE.
#' @return K-vector of PIPs, or NULL if the fit has no point-normal prior trace.
#' @examples
#' # lfp_pip(fit)
lfp_pip <- function(fit) {
  tr <- Filter(Negate(is.null), fit$history$prior_trace)
  if (!length(tr)) return(NULL)
  last <- tr[[length(tr)]]; g <- last$prior_before
  if (!inherits(g, "normalmix") || length(g$pi) != 2) return(NULL)
  i0 <- which(g$sd == 0); i1 <- which(g$sd > 0)
  if (length(i0) != 1 || length(i1) != 1) return(NULL)
  x <- last$x; s <- last$s
  out <- rep(NA_real_, length(x)); ok <- is.finite(s)
  l0 <- g$pi[i0] * dnorm(x[ok], 0, s[ok]); l1 <- g$pi[i1] * dnorm(x[ok], 0, sqrt(g$sd[i1]^2 + s[ok]^2))
  out[ok] <- l1 / (l0 + l1)
  out
}

#' Run an expression, capturing warnings and errors instead of stopping
#' @param expr Expression to evaluate.
#' @return list(value, warnings, seconds); value is an error condition on failure.
#' @examples
#' lfp_try(log(-1))
lfp_try <- function(expr) {
  warn <- character(); t0 <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  list(value = value, warnings = paste(unique(warn), collapse = " | "),
       seconds = proc.time()[["elapsed"]] - t0)
}

# Simulation ----

#' One synthetic dataset for a registered cell
#'
#' S-9 cells use yfb_decision_data(seed, 2000, share, "shared"); the S-14 cell uses
#' s14_data(seed, 0.0005, "strong"), which returns the same object for strong survival.
#' Y and Yv are X/sd within each cohort (nonnegative), the PE input; PL z-scores them.
#' @param cell One SIM_CELLS row. @param seed Dataset seed.
#' @return Generator output list.
#' @examples
#' # lfp_sim_data(SIM_CELLS[1, ], 1)
lfp_sim_data <- function(cell, seed) {
  if (cell$source == "s14") s14_data(seed, cell$share, "strong")
  else yfb_decision_data(seed, 2000L, cell$share, truth_mode = "shared")
}

#' Saved S-9 or S-14 PL metrics for the reuse spot check
#' @param cell One SIM_CELLS row. @param seed Dataset seed.
#' @return One-row data frame (recovered, false_active, c_validation) or NULL if absent.
#' @examples
#' # lfp_saved_pl(SIM_CELLS[1, ], 1)
lfp_saved_pl <- function(cell, seed) {
  path <- if (cell$source == "s9") file.path(S9_GRID, sprintf("p2000_v%s_s%02d_w1.rds", cell$share, seed))
          else file.path(S14_FITS, sprintf("strong_sh%s_s%02d_A1_w1.rds", format(cell$share, scientific = FALSE), seed))
  if (!file.exists(path)) return(NULL)
  m <- readRDS(path)$metrics
  data.frame(saved_recovered = m$recovered, saved_false_active = m$false_active,
             saved_c_validation = m$c_validation)
}

#' Fit and score one (cell, seed, arm) simulation task
#'
#' Scoring is s9_score() for both arms: recovery is best |r(truth, E[F_k])| >= 0.7 with the
#' matched program survival-active (|E beta| >= 1.96 SD) and the correct sign; false-active
#' counts other active programs; calibration slope uses exact event times (timefix = FALSE).
#' The validation matrix is scaled within itself, in the arm's input scale.
#' @param task Row with cell, share, source, seed, arm.
#' @param cfg Parsed config.
#' @return Status string; the checkpoint holds the metric row and a compact fit state.
#' @examples
#' # lfp_sim_task(task, cfg)
lfp_sim_task <- function(task, cfg) {
  key <- sprintf("sim_%s_s%02d_%s", task$cell, task$seed, task$arm)
  path <- file.path(PRIV, "sim", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  d <- lfp_sim_data(task, task$seed)
  pl <- task$arm != "PE"   # PL and every task B arm use the z-scored input
  Ytr <- if (pl) s9_zscore(d$Y) else d$Y
  Yva <- if (pl) s9_zscore(d$Yv) else d$Yv
  set.seed(task$seed + 42L)
  run <- lfp_try(lfp_fit(Ytr, d$time, d$event, task$arm, cfg))
  m <- s9_score(run$value, d, Ytr, Yva)
  live <- if (inherits(run$value, "error")) NA_integer_ else sum(colSums(run$value$EF^2) > 1e-8)
  row <- cbind(task[c("cell", "share", "seed", "arm")], m, n_alive = live,
               fit_sec = run$seconds, warnings = run$warnings)
  if (task$arm == "PL") { sv <- lfp_saved_pl(task, task$seed); if (!is.null(sv)) row <- cbind(row, sv) }
  st <- if (inherits(run$value, "error")) NULL else
    run$value[c("EF", "EBeta", "EBeta2", "EF_norms", "Tau")]
  pip <- if (inherits(run$value, "error")) NULL else lfp_pip(run$value)
  if (!is.null(pip)) row$pip <- paste(signif(pip, 4), collapse = ";")
  saveRDS(list(metrics = row, state = st, pip = pip), path)
  sprintf("%s ok %.0fs", key, run$seconds)
}

# Real data ----

#' Prepare the three real-data inputs on the D4 genes, in the S-11 patient order
#'
#' Rebuilds the D4 log-scale matrix without standardization through the same
#' preprocess_merged_cohorts() call S-11 used (single log, combined mean/variance rank,
#' top 3,000 per cohort, intersection), checks that its genes equal S-11's D4 and that its
#' per-cohort z-score reproduces S-11's matrix, then forms the X/sd and rank inputs.
#' Uses the S-11 folds, so every arm sees identical patients and folds.
#' @return NULL; writes REAL_IN (patient-level, git-ignored) and an aggregate audit CSV.
#' @examples
#' # lfp_prepare()
lfp_prepare <- function() {
  if (file.exists(REAL_IN)) stop("Prepared input exists; do not overwrite a registered design.")
  cfg <<- yaml::read_yaml("config/globals.yml")   # benchmark_helpers.R reads cfg globally
  source("results/benchmark_sim/benchmark_helpers.R"); source("code/preprocess_desurv.R")
  z <- readRDS(S11_IN)
  train <- cfg$pdac$training_cohorts
  raw <- lapply(setNames(train, train), load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
  x <- preprocess_merged_cohorts(raw, PLATFORM_LOG_TRANSFORM[train],
    top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
    per_platform_standardize = FALSE, normalize_method = "none",
    selection_per_cohort = TRUE, selection_method = "combined_rank")
  stopifnot(setequal(x$gene_names, z$d4), nrow(x$Y) == length(z$time))
  X <- x$Y[, match(z$d4, x$gene_names), drop = FALSE]; colnames(X) <- z$d4
  cohort <- z$cohort
  by_cohort <- function(M, f) {           # apply f to each cohort's rows
    out <- M
    for (co in unique(cohort)) out[cohort == co, ] <- f(M[cohort == co, , drop = FALSE])
    out
  }
  sd_scale <- function(M) sweep(M, 2, pmax(apply(M, 2, sd), 1e-8), "/")
  Yz <- by_cohort(X, function(M) scale(M))
  ref <- z$Y[, match(z$d4, z$genes), drop = FALSE]
  zdiff <- max(abs(unname(Yz) - unname(ref)))
  if (zdiff > 1e-8) stop(sprintf("Per-cohort z-score does not reproduce S-11 (max diff %.3g)", zdiff))
  Ysd <- by_cohort(X, sd_scale)
  R <- rank_transform_subjects(X)
  Yrank <- by_cohort(R, sd_scale)
  dir.create(dirname(REAL_IN), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(Y = list(PL = ref, PE = Ysd, RANK = Yrank), X = X, genes = z$d4,
               time = z$time, status = z$status, cohort = cohort, fold = z$fold), REAL_IN)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  audit <- do.call(rbind, lapply(c("PL", "PE", "RANK"), function(a) {
    M <- list(PL = ref, PE = Ysd, RANK = Yrank)[[a]]
    do.call(rbind, lapply(unique(cohort), function(co) data.frame(arm = a, cohort = co,
      n = sum(cohort == co), min = min(M[cohort == co, ]), max = max(M[cohort == co, ]),
      median_gene_mean = median(colMeans(M[cohort == co, ])),
      median_gene_sd = median(apply(M[cohort == co, ], 2, sd)))))
  }))
  audit$raw_log_min <- min(X); audit$z_reproduction_max_diff <- zdiff
  write.csv(audit, file.path(OUT, "real_input_audit.csv"), row.names = FALSE)
  writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()),
               sprintf("n = %d, events = %d, genes = %d", length(z$time), sum(z$status), ncol(X))),
             file.path(OUT, "real_provenance.txt"))
  print(audit)
  invisible(NULL)
}

#' Pooled, within-cohort and per-cohort concordance of a fixed risk score
#'
#' Orientation is fixed by the model (higher risk = higher hazard, reverse = TRUE), never
#' re-chosen on the scored patients, as in frozen_reverse_cindex().
#' @param risk Held-out risk scores. @param time,status,cohort Held-out outcomes and labels.
#' @return One-row data frame of concordances (NA when the score is constant).
#' @examples
#' # lfp_concordance(rnorm(20), rexp(20), rbinom(20, 1, .7), rep(c("a", "b"), 10))
lfp_concordance <- function(risk, time, status, cohort) {
  if (sd(risk) < 1e-12) return(data.frame(C_within = NA_real_, C_pooled = NA_real_,
                                          C_TCGA = NA_real_, C_CPTAC = NA_real_))
  w <- survival::concordance(survival::Surv(time, status) ~ risk + strata(cohort), reverse = TRUE)
  one <- function(co) { i <- cohort == co
    if (sum(i) < 3 || sd(risk[i]) < 1e-12) NA_real_ else frozen_reverse_cindex(risk[i], time[i], status[i]) }
  cos <- unique(cohort)
  out <- data.frame(C_within = as.numeric(w$concordance),
                    C_pooled = frozen_reverse_cindex(risk, time, status))
  out$C_TCGA  <- one(cos[grepl("TCGA", cos)][1])
  out$C_CPTAC <- one(cos[grepl("CPTAC", cos)][1])
  out
}

#' Fit and score one (arm, fold) real-data task; fold 0 is the full-data fit
#'
#' Records convergence, programs with nonzero loadings, survival-active programs, held-out
#' concordances and partial log-likelihood (pooled and within cohort), and per program the
#' correlation of its training scores YF_k with each patient's mean log expression and with
#' the TCGA indicator.
#' @param task Row with arm and fold. @param z Prepared real inputs. @param cfg Config.
#' @return Status string; checkpoint (data/) holds metrics, program table, state and risk.
#' @examples
#' # lfp_real_task(task, z, cfg)
lfp_real_task <- function(task, z, cfg) {
  key <- sprintf("real_%s_f%d", task$arm, task$fold)
  path <- file.path(PRIV, "real", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  Y <- z$Y[[if (task$arm %in% names(BETA_ARMS)) "PL" else task$arm]]
  te <- which(z$fold == task$fold)
  tr <- if (task$fold == 0L) seq_along(z$time) else which(z$fold != task$fold)
  set.seed(42L)
  run <- lfp_try(lfp_fit(Y[tr, , drop = FALSE], z$time[tr], z$status[tr], task$arm, cfg))
  row <- data.frame(arm = task$arm, fold = task$fold, n_train = length(tr), n_test = length(te),
                    fit_sec = run$seconds, warnings = run$warnings, error = NA_character_)
  if (inherits(run$value, "error")) {
    row$error <- conditionMessage(run$value); row$converged <- FALSE
    saveRDS(list(metrics = row), path); return(paste(key, "ERROR", row$error))
  }
  fit <- run$value
  sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  alive <- colSums(fit$EF^2) > 1e-8
  active <- abs(fit$EBeta) >= 1.96 * sdb & sdb > 0 & alive
  Ztr <- Y[tr, , drop = FALSE] %*% fit$EF
  meanexpr <- rowMeans(z$X[tr, , drop = FALSE])
  tcga <- as.numeric(grepl("TCGA", z$cohort[tr]))
  safe_cor <- function(a, b) if (sd(a) < 1e-12) NA_real_ else cor(a, b)
  prog <- data.frame(arm = task$arm, fold = task$fold, program = seq_len(ncol(fit$EF)),
    alive = alive, active = active, beta = fit$EBeta, beta_sd = sdb,
    n_nonzero_genes = colSums(abs(fit$EF) > 1e-8),
    cor_mean_expr = apply(Ztr, 2, safe_cor, b = meanexpr),
    cor_cohort = apply(Ztr, 2, safe_cor, b = tcga))
  pip <- lfp_pip(fit); prog$pip <- if (is.null(pip)) NA_real_ else pip
  row <- cbind(row, data.frame(converged = isTRUE(fit$history$converged), n_iter = fit$history$n_iter,
    n_alive = sum(alive), n_active = sum(active), max_abs_beta = max(abs(fit$EBeta)),
    max_abs_cor_mean_expr = suppressWarnings(max(abs(prog$cor_mean_expr[alive]), na.rm = TRUE)),
    max_abs_cor_cohort = suppressWarnings(max(abs(prog$cor_cohort[alive]), na.rm = TRUE)),
    recon_train = sum((Y[tr, ] - fit$EL %*% t(fit$EF))^2) / sum(Y[tr, ]^2)))
  risk <- NULL
  if (task$fold > 0L) {
    risk <- drop(Y[te, , drop = FALSE] %*% fit$EF %*% fit$EBeta)
    row <- cbind(row, lfp_concordance(risk, z$time[te], z$status[te], z$cohort[te]))
    eta <- drop(Y %*% fit$EF %*% fit$EBeta)
    row$cvpl_pooled <- calc_cox_taylor_yf(eta, z$time, z$status)$logPL -
                       calc_cox_taylor_yf(eta[tr], z$time[tr], z$status[tr])$logPL
    row$cvpl_within <- calc_cox_taylor_yf(eta, z$time, z$status, strata = z$cohort)$logPL -
                       calc_cox_taylor_yf(eta[tr], z$time[tr], z$status[tr], strata = z$cohort[tr])$logPL
    row$risk_sd <- sd(risk)
  }
  saveRDS(list(metrics = row, programs = prog, te = te, risk = risk,
               state = list(EF = fit$EF, EBeta = fit$EBeta, EBeta2 = fit$EBeta2,
                            EL = if (task$fold == 0L) fit$EL else NULL, Tau = fit$Tau,
                            elbo = fit$history$elbo_full)), path)
  sprintf("%s ok %.0fs", key, run$seconds)
}

# Design and execution ----

#' Task list: real-data fits first (longest), then simulation fits
#' @param smoke If TRUE, one dataset per arm (seed 1, 1% share) and real fold 1 per arm.
#' @return list(real, sim) of design data frames.
#' @examples
#' lfp_design(TRUE)
lfp_design <- function(smoke = FALSE, part = "lf") {
  if (part == "beta") {
    b <- names(BETA_ARMS)
    if (smoke) return(list(real = data.frame(arm = b, fold = 1L, stringsAsFactors = FALSE),
      sim = merge(cbind(SIM_CELLS[c(2, 4), ], seed = 1L), data.frame(arm = b))))
    return(list(real = expand.grid(arm = b, fold = 0:5, stringsAsFactors = FALSE),
                sim = merge(merge(SIM_CELLS, data.frame(seed = 1:10)), data.frame(arm = b))))
  }
  arms_real <- c("PL", "PE", "RANK")
  if (smoke) return(list(real = data.frame(arm = arms_real, fold = 1L, stringsAsFactors = FALSE),
    sim = merge(cbind(SIM_CELLS[2, ], seed = 1L), data.frame(arm = c("PL", "PE")))))
  list(real = expand.grid(arm = arms_real, fold = 0:5, stringsAsFactors = FALSE),
       sim = merge(merge(SIM_CELLS, data.frame(seed = 1:10)), data.frame(arm = c("PL", "PE"))))
}

#' Run every task on a PSOCK cluster, keeping failures as status lines
#' @param design lfp_design() output. @param cores Worker count.
#' @return Character status vector.
#' @examples
#' # lfp_run(lfp_design(TRUE), 5)
lfp_run <- function(design, cores) {
  dir.create(file.path(PRIV, "sim"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PRIV, "real"), recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  tasks <- c(lapply(split(design$real, seq_len(nrow(design$real))), function(t) c(type = "real", t)),
             lapply(split(design$sim, seq_len(nrow(design$sim))), function(t) c(type = "sim", t)))
  cl <- parallel::makePSOCKcluster(cores, outfile = file.path(PRIV, "workers.log"))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, c("OUT", "PRIV", "REAL_IN", "S9_GRID", "S14_FITS", "SIM_CELLS",
    "BETA_ARMS", "lfp_load", "lfp_fit", "lfp_pip", "lfp_try", "lfp_sim_data", "lfp_saved_pl", "lfp_sim_task",
    "lfp_concordance", "lfp_real_task"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1); lfp_load(); cfg <- yaml::read_yaml("config/globals.yml")
    z <- readRDS(REAL_IN); NULL })
  st <- parallel::parLapplyLB(cl, tasks, function(t) tryCatch({
    t <- as.data.frame(t, stringsAsFactors = FALSE)
    if (t$type == "real") { t$fold <- as.integer(t$fold); lfp_real_task(t, z, cfg) }
    else { t$seed <- as.integer(t$seed); t$share <- as.numeric(t$share); lfp_sim_task(t, cfg) }
  }, error = function(e) paste("TASK ERROR", conditionMessage(e))), chunk.size = 1L)
  unlist(st)
}

# Summary ----

#' Aggregate checkpoints into CSVs and apply the four registered criteria
#'
#' Criteria (PE adopted if all hold): (1) recovery in every cell >= PL's minus 1 of 10;
#' (2) mean false-active programs per fit <= 1 in every cell; (3) >= 9 of 10 converged in
#' every cell; (4) mean within-cohort CV concordance not lower than PL's by more than one
#' SE of the paired fold difference, and every PE fold converges. A missing fit counts as
#' not recovered, not converged and 7 false-active programs (s9_score's failure row).
#' @return NULL; writes sim_fits.csv, sim_cells.csv, real_fits.csv, real_programs.csv,
#'   real_paired.csv, real_gene_overlap.csv and criteria.csv to OUT.
#' @examples
#' # lfp_summary()
lfp_summary <- function() {
  rd <- function(dir) lapply(list.files(file.path(PRIV, dir), full.names = TRUE, pattern = "rds$"), readRDS)
  # lfp_bind: task B checkpoints in the same folder carry an extra pip column; the merge
  # with the A design below keeps only A's arms.
  sim <- lfp_bind(lapply(rd("sim"), function(x) {
    m <- x$metrics; for (nm in c("saved_recovered", "saved_false_active", "saved_c_validation"))
      if (is.null(m[[nm]])) m[[nm]] <- NA; m }))
  des <- lfp_design()$sim
  sim <- merge(des[c("cell", "seed", "arm")], sim, all.x = TRUE)
  sim$missing <- is.na(sim$converged)
  sim$recovered[sim$missing] <- FALSE; sim$converged[sim$missing] <- FALSE
  sim$false_active[sim$missing] <- 7L
  write.csv(sim, file.path(OUT, "sim_fits.csv"), row.names = FALSE)
  cells <- do.call(rbind, lapply(split(sim, list(sim$cell, sim$arm), drop = TRUE), function(s)
    data.frame(cell = s$cell[1], arm = s$arm[1], n = nrow(s), recovered = sum(s$recovered),
      mean_false_active = mean(s$false_active), converged = sum(s$converged),
      median_cor = median(s$cor_planted, na.rm = TRUE),
      median_c_validation = median(s$c_validation, na.rm = TRUE),
      median_cal_slope = median(s$cal_slope, na.rm = TRUE),
      median_alive = median(s$n_alive, na.rm = TRUE), missing = sum(s$missing),
      median_fit_sec = median(s$fit_sec, na.rm = TRUE))))
  write.csv(cells, file.path(OUT, "sim_cells.csv"), row.names = FALSE)
  pl <- sim[sim$arm == "PL" & !is.na(sim$saved_recovered), ]
  if (nrow(pl)) write.csv(data.frame(cell = pl$cell, seed = pl$seed,
    same_recovered = pl$recovered == pl$saved_recovered,
    same_false_active = pl$false_active == pl$saved_false_active,
    same_c_3dp = round(pl$c_validation, 3) == round(pl$saved_c_validation, 3),
    c_new = pl$c_validation, c_saved = pl$saved_c_validation),
    file.path(OUT, "pl_reuse_check.csv"), row.names = FALSE)

  real <- Filter(function(x) x$metrics$arm %in% c("PL", "PE", "RANK"), rd("real"))
  rf <- do.call(rbind, lapply(real, function(x) { m <- x$metrics
    for (nm in c("C_within", "C_pooled", "C_TCGA", "C_CPTAC", "cvpl_pooled", "cvpl_within", "risk_sd",
                 "n_iter", "n_alive", "n_active", "max_abs_beta", "max_abs_cor_mean_expr",
                 "max_abs_cor_cohort", "recon_train")) if (is.null(m[[nm]])) m[[nm]] <- NA
    m[c("arm", "fold", "n_train", "n_test", "converged", "n_iter", "n_alive", "n_active",
        "max_abs_beta", "max_abs_cor_mean_expr", "max_abs_cor_cohort", "recon_train",
        "C_within", "C_pooled", "C_TCGA", "C_CPTAC", "cvpl_pooled", "cvpl_within", "risk_sd",
        "fit_sec", "warnings", "error")] }))
  rf <- merge(lfp_design()$real, rf, all.x = TRUE)
  rf$converged[is.na(rf$converged)] <- FALSE
  write.csv(rf, file.path(OUT, "real_fits.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(real, `[[`, "programs")), file.path(OUT, "real_programs.csv"),
            row.names = FALSE)
  # S-11 saved PL fold concordance (pooled) for the real-data reuse comparison
  s11 <- do.call(rbind, lapply(0:5, function(f) {
    p <- file.path(S11_FITS, sprintf("size_D4_N2218_w1_f%d.rds", f))
    if (!file.exists(p)) return(NULL)
    m <- readRDS(p)$metrics; data.frame(fold = f, s11_C_pooled = m$C, s11_converged = m$converged,
                                        s11_n_active = m$n_active, s11_n_alive = m$n_alive) }))
  cv <- rf[rf$fold > 0, ]
  paired <- do.call(rbind, lapply(c("PE", "RANK"), function(a) {
    x <- merge(cv[cv$arm == a, c("fold", "C_within", "C_pooled", "converged")],
               cv[cv$arm == "PL", c("fold", "C_within", "C_pooled", "converged")],
               by = "fold", suffixes = c("", "_PL"))
    dd <- x$C_within - x$C_within_PL; dp <- x$C_pooled - x$C_pooled_PL
    data.frame(arm = a, n_folds = sum(is.finite(dd)), mean_C_within = mean(x$C_within),
      mean_C_within_PL = mean(x$C_within_PL), diff_within = mean(dd),
      se_within = sd(dd) / sqrt(sum(is.finite(dd))), mean_C_pooled = mean(x$C_pooled),
      mean_C_pooled_PL = mean(x$C_pooled_PL), diff_pooled = mean(dp),
      se_pooled = sd(dp) / sqrt(sum(is.finite(dp))), folds_converged = sum(x$converged),
      folds_converged_PL = sum(x$converged_PL))
  }))
  write.csv(paired, file.path(OUT, "real_paired.csv"), row.names = FALSE)
  if (!is.null(s11)) write.csv(merge(rf[rf$arm == "PL", c("fold", "C_pooled", "converged", "n_active", "n_alive")],
    s11, all = TRUE), file.path(OUT, "pl_s11_reuse_check.csv"), row.names = FALSE)
  # Top-100 gene overlap of each full-data program with every program of the current S-11 fit
  cur <- readRDS(file.path(S11_FITS, "size_D4_N2218_w1_f0.rds"))
  zin <- readRDS(REAL_IN)
  top <- function(EF, k) order(-abs(EF[, k]))[1:100]
  cur_EF <- cur$state$EF[match(zin$genes, cur$genes), , drop = FALSE]
  cur_alive <- which(colSums(cur_EF^2) > 1e-8)
  ov <- do.call(rbind, lapply(real, function(x) {
    if (x$metrics$fold != 0 || is.null(x$state)) return(NULL)
    EF <- x$state$EF; al <- which(colSums(EF^2) > 1e-8)
    if (!length(al)) return(NULL)
    do.call(rbind, lapply(al, function(k) data.frame(arm = x$metrics$arm, program = k,
      current_program = cur_alive,
      overlap_top100 = vapply(cur_alive, function(j) length(intersect(top(EF, k), top(cur_EF, j))), 1L))))
  }))
  if (!is.null(ov)) write.csv(ov, file.path(OUT, "real_gene_overlap.csv"), row.names = FALSE)

  ce <- merge(cells[cells$arm == "PE", ], cells[cells$arm == "PL", ], by = "cell", suffixes = c("_PE", "_PL"))
  pe <- paired[paired$arm == "PE", ]
  crit <- data.frame(
    criterion = c("1 recovery >= PL - 1 in every cell", "2 mean false-active <= 1 in every cell",
                  "3 >= 9/10 converged in every cell",
                  "4 CV within-cohort C not lower than PL by > 1 SE; all folds converge"),
    detail = c(paste(sprintf("%s %d vs %d", ce$cell, ce$recovered_PE, ce$recovered_PL), collapse = "; "),
               paste(sprintf("%s %.2f", ce$cell, ce$mean_false_active_PE), collapse = "; "),
               paste(sprintf("%s %d/10", ce$cell, ce$converged_PE), collapse = "; "),
               sprintf("diff %.4f, SE %.4f, PE folds converged %d/5", pe$diff_within, pe$se_within,
                       pe$folds_converged)),
    pass = c(all(ce$recovered_PE >= ce$recovered_PL - 1), all(ce$mean_false_active_PE <= 1),
             all(ce$converged_PE >= 9),
             isTRUE(pe$diff_within >= -pe$se_within && pe$folds_converged == 5)))
  write.csv(crit, file.path(OUT, "criteria.csv"), row.names = FALSE)
  print(cells); print(paired); print(crit)
  invisible(NULL)
}

#' Bind data frames whose columns differ, filling missing columns with NA
#' @param xs List of data frames. @return One data frame with the union of columns.
#' @examples
#' lfp_bind(list(data.frame(a = 1), data.frame(a = 2, b = 3)))
lfp_bind <- function(xs) {
  xs <- Filter(Negate(is.null), xs); cols <- unique(unlist(lapply(xs, names)))
  do.call(rbind, lapply(xs, function(x) { for (nm in setdiff(cols, names(x))) x[[nm]] <- NA; x[cols] }))
}

#' Task B summary: beta-prior arms against the per-program Normal working model
#'
#' Normal is A's PL arm (same checkpoints). Rule (registered): a sparse prior replaces Normal
#' if (a) recovery >= Normal's minus 1, mean false-active <= 1 and >= 9/10 converged in every
#' cell; (b) median calibration slope closer to 1 than Normal's in >= 2 of 4 cells; (c) mean
#' within-cohort CV C not lower than Normal's by more than one SE of the paired difference.
#' Missing fits count as failures. Also reports how the S-14 cell alone affects (a), and
#' whether Normal itself would meet A's false-activity limit there.
#' @return NULL; writes beta_*.csv to OUT.
#' @examples
#' # lfp_summary_beta()
lfp_summary_beta <- function() {
  arms <- c("PL", names(BETA_ARMS))
  label <- c(PL = "Normal (per program)", BETA_ARMS)
  rd <- function(dir) lapply(list.files(file.path(PRIV, dir), full.names = TRUE, pattern = "rds$"), readRDS)
  sim <- lfp_bind(lapply(rd("sim"), `[[`, "metrics"))
  des <- merge(merge(SIM_CELLS[c("cell")], data.frame(seed = 1:10)), data.frame(arm = arms))
  sim <- merge(des, sim[sim$arm %in% arms, ], all.x = TRUE)
  sim$missing <- is.na(sim$converged)
  sim$recovered[sim$missing] <- FALSE; sim$converged[sim$missing] <- FALSE
  sim$false_active[sim$missing] <- 7L
  sim$prior <- label[sim$arm]
  write.csv(sim, file.path(OUT, "beta_sim_fits.csv"), row.names = FALSE)
  cells <- do.call(rbind, lapply(split(sim, list(sim$cell, sim$arm), drop = TRUE), function(s)
    data.frame(cell = s$cell[1], arm = s$arm[1], prior = s$prior[1], n = nrow(s),
      recovered = sum(s$recovered), mean_false_active = mean(s$false_active),
      converged = sum(s$converged), median_cor = median(s$cor_planted, na.rm = TRUE),
      median_c_validation = median(s$c_validation, na.rm = TRUE),
      median_cal_slope = median(s$cal_slope, na.rm = TRUE),
      all_beta_zero = sum(s$max_beta < 1e-8, na.rm = TRUE), missing = sum(s$missing))))
  write.csv(cells, file.path(OUT, "beta_sim_cells.csv"), row.names = FALSE)

  real <- Filter(function(x) x$metrics$arm %in% arms, rd("real"))
  rf <- lfp_bind(lapply(real, `[[`, "metrics"))
  rf <- merge(expand.grid(arm = arms, fold = 0:5, stringsAsFactors = FALSE), rf, all.x = TRUE)
  rf$converged[is.na(rf$converged)] <- FALSE
  write.csv(rf[setdiff(names(rf), "warnings")], file.path(OUT, "beta_real_fits.csv"), row.names = FALSE)
  progs <- lfp_bind(lapply(real, `[[`, "programs"))
  write.csv(progs, file.path(OUT, "beta_real_programs.csv"), row.names = FALSE)
  cv <- rf[rf$fold > 0, ]
  paired <- do.call(rbind, lapply(names(BETA_ARMS), function(a) {
    x <- merge(cv[cv$arm == a, c("fold", "C_within", "C_pooled", "converged", "n_active", "max_abs_beta")],
               cv[cv$arm == "PL", c("fold", "C_within", "C_pooled", "converged", "n_active")],
               by = "fold", suffixes = c("", "_N"))
    dd <- x$C_within - x$C_within_N; dp <- x$C_pooled - x$C_pooled_N
    data.frame(arm = a, prior = label[a], n_folds = sum(is.finite(dd)),
      mean_C_within = mean(x$C_within), mean_C_within_Normal = mean(x$C_within_N),
      diff_within = mean(dd), se_within = sd(dd) / sqrt(sum(is.finite(dd))),
      mean_C_pooled = mean(x$C_pooled), diff_pooled = mean(dp),
      se_pooled = sd(dp) / sqrt(sum(is.finite(dp))), folds_converged = sum(x$converged),
      mean_active = mean(x$n_active), mean_active_Normal = mean(x$n_active_N),
      folds_no_active = sum(x$n_active == 0),
      folds_all_beta_zero = sum(x$max_abs_beta < 1e-8))   # every E[beta_k] numerically zero
  }))
  write.csv(paired, file.path(OUT, "beta_real_paired.csv"), row.names = FALSE)

  nor <- cells[cells$arm == "PL", ]
  rule <- do.call(rbind, lapply(names(BETA_ARMS)[1:2], function(a) {
    ce <- merge(cells[cells$arm == a, ], nor, by = "cell", suffixes = c("", "_N"))
    p <- paired[paired$arm == a, ]
    a_cell <- ce$recovered >= ce$recovered_N - 1 & ce$mean_false_active <= 1 & ce$converged >= 9
    slope_closer <- abs(ce$median_cal_slope - 1) < abs(ce$median_cal_slope_N - 1)
    non_s14 <- ce$cell != "s14_0.05pct"
    data.frame(arm = a, prior = label[a],
      a_all_cells = all(a_cell), a_cells_failed = paste(ce$cell[!a_cell], collapse = ";"),
      a_without_s14 = all(a_cell[non_s14]),
      b_slope_closer_cells = sum(slope_closer), b_pass = sum(slope_closer) >= 2,
      c_diff = p$diff_within, c_se = p$se_within, c_pass = isTRUE(p$diff_within >= -p$se_within),
      all_beta_zero_real = p$folds_all_beta_zero == p$n_folds,
      mean_false_active_all_cells = mean(ce$mean_false_active),
      qualifies = all(a_cell) && sum(slope_closer) >= 2 && isTRUE(p$diff_within >= -p$se_within) &&
                  p$folds_all_beta_zero < p$n_folds,
      b_slope_closer_cells_without_s14 = sum(slope_closer[non_s14]),
      # Sensitivity 1: remove the S-14 cell from parts (a) and (b), keeping the two-cell threshold.
      qualifies_without_s14 = all(a_cell[non_s14]) && sum(slope_closer[non_s14]) >= 2 &&
                              isTRUE(p$diff_within >= -p$se_within) && p$folds_all_beta_zero < p$n_folds,
      # Sensitivity 2: exempt S-14 from part (a) only, keeping its calibration result in (b).
      qualifies_s14_exempt_a_only = all(a_cell[non_s14]) && sum(slope_closer) >= 2 &&
                              isTRUE(p$diff_within >= -p$se_within) && p$folds_all_beta_zero < p$n_folds)
  }))
  rule$normal_meets_false_active_all_cells <- all(nor$mean_false_active <= 1)
  rule$normal_s14_false_active <- nor$mean_false_active[nor$cell == "s14_0.05pct"]
  write.csv(rule, file.path(OUT, "beta_rule.csv"), row.names = FALSE)
  print(cells[order(cells$cell, cells$arm), ]); print(paired); print(rule)
  invisible(NULL)
}

# Main ----
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE); mode <- args[1]
  cores <- if (length(args) > 1) as.integer(args[2]) else 5L
  if (mode == "prepare") { lfp_prepare(); quit(save = "no") }
  if (mode %in% c("smoke", "grid")) {
    dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
    writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()), mode),
               file.path(OUT, sprintf("provenance_%s.txt", mode)))
    st <- lfp_run(lfp_design(mode == "smoke"), cores)
    writeLines(st, file.path(OUT, sprintf("task_status_%s.txt", mode))); print(st)
  }
  if (mode %in% c("beta_smoke", "beta_grid")) {
    dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
    writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()), mode),
               file.path(OUT, sprintf("provenance_%s.txt", mode)))
    st <- lfp_run(lfp_design(mode == "beta_smoke", "beta"), cores)
    writeLines(st, file.path(OUT, sprintf("task_status_%s.txt", mode))); print(st)
  }
  if (mode == "summary") { lfp_load(); lfp_summary() }
  if (mode == "beta_summary") { lfp_load(); lfp_summary_beta() }
}
