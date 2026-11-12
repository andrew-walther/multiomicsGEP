# ============================================================
# Script: run_lf_prior_followup.R
# Purpose: Follow-up A2. Point-exponential L and F with the TCGA vs CPTAC mean difference
#          removed, by a shifted per-cohort z-score (pooled) or by cohort strata on X/sd.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm, yaml, parallel, tools;
#               results/benchmark_sim/run_lf_prior_comparison.R (task A: lfp_load, lfp_try,
#               lfp_concordance, lfp_sim_data, SIM_CELLS, REAL_IN, S11_FITS) and, through
#               lfp_load(), the S-9 scorer s9_score()/s9_zscore() and calc_cox_taylor_yf().
# Run: Rscript results/benchmark_sim/run_lf_prior_followup.R smoke [cores]
#      Rscript results/benchmark_sim/run_lf_prior_followup.R grid [cores]
#      Rscript results/benchmark_sim/run_lf_prior_followup.R summary
# Registration: 2026-11-12 entry "A2" in docs/plans/2026-10-16-round2-attempts-log.md.
# ============================================================
# Setups (both: point-exponential L and F, per-program Normal beta, K = 7, canonicalize = TRUE,
# w_S = 1, alpha, max_iter and tol from config, trace_supervision = TRUE, no pruning):
#   SZ   Shifted z-score. The per-cohort z-score plus one constant per gene,
#        c_j = -min_i z_ij over all 273 training patients (outcome-free, computed once), so the
#        input is nonnegative and both cohorts keep mean 0 + c_j and SD 1. Pooled,
#        center_projection = TRUE, no strata. Real data (five S-11 folds + full fit) and A's four
#        simulation cells (10 seeds); there the synthetic training z-score is shifted by its own
#        per-gene minimum and the validation z-score by the training constants.
#   PES  X/sd within cohort (A's PE input) with strata_id = cohort and center_projection = TRUE,
#        which centers the projection within cohort and gives each cohort its own baseline
#        hazard. Real data only. In A's one-cohort simulation a single stratum is identical to
#        A's PE arm, so A's PE simulation results stand for PES there.
# Held-out risk is Y[te, ] %*% E[F] %*% E[beta] in the setup's own input. A per-gene constant
# shifts every held-out score by the same amount, and a per-cohort shift does not change
# within-cohort pairs, so within-cohort concordance (primary) is unaffected by either. Pooled
# concordance is not interpretable for PES (its scores keep the cohort offset).
#
# Parallelization: PSOCK workers (forked workers segfaulted on macOS with the compiled sweep),
# one BLAS/OpenMP thread each, default 2 workers (task B and other sessions share the machine).
# Task A's runner is sourced once in the main process and its functions exported, so edits made
# to it by another session during the run cannot reach the workers. One checkpoint per fit.
# Memory: a worker holds two 400 x 2,000 synthetic matrices or the 273 x 2,218 real matrix.
# Patient-level fits and risk scores go to data/ (git-ignored); only aggregate CSVs to OUT2.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_lf_prior_comparison.R")   # task A; its main block does not run
OUT2  <- "results/benchmark_sim/outputs/lf_prior_2026-11-09/a2"
PRIV2 <- "data/lf_prior_2026-11-12/a2"
A_OUT <- "results/benchmark_sim/outputs/lf_prior_2026-11-09"   # task A aggregate outputs (read only)

# Inputs ----

#' Shift each column so its minimum is 0
#'
#' Returns the shifted matrix and the constants c_j = -min_i Y_ij, so the same constants can be
#' applied to another matrix (the synthetic validation set).
#' @param Y Numeric matrix (patients x genes).
#' @return list(Y = Y + 1 c', c = c).
#' @examples
#' a2_shift(matrix(rnorm(6), 3))
a2_shift <- function(Y) {
  cc <- -apply(Y, 2, min)
  list(Y = sweep(Y, 2, cc, "+"), c = cc)
}

# Fitting ----

#' Fit one A2 setup with the working-model settings and point-exponential L and F
#'
#' Differs from task A's PE call only by strata_id for PES. center_projection = TRUE makes the
#' Cox projection use Y - 1 ybar' (SZ) or Y minus each cohort's training means (PES), while the
#' reconstruction Y ~ L F' uses the nonnegative Y.
#' @param Y Nonnegative training matrix in the setup's input.
#' @param time,status Training outcomes.
#' @param setup "SZ" or "PES".
#' @param cfg Parsed config/globals.yml.
#' @param strata Training cohort labels (used only by PES).
#' @return fit_ssbmf_yfb() result.
#' @examples
#' # a2_fit(Y, time, status, "PES", cfg, cohort)
a2_fit <- function(Y, time, status, setup, cfg, strata = NULL) {
  fit_ssbmf_yfb(Y, time, status, K = cfg$benchmark$k_merged_yfb_desurv,
                max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = "point_exponential", prior_beta = "normal",
                beta_prior_mode = "per_program", canonicalize = TRUE, w_S = 1,
                alpha = cfg$benchmark$alpha, verbose = FALSE, trace_supervision = TRUE,
                center_projection = TRUE,
                strata_id = if (setup == "PES") strata else NULL)
}

# Simulation ----

#' Fit and score one (cell, seed) simulation task for SZ
#'
#' Scoring is s9_score() as in task A. Also records, for every false-active program j, the
#' ratio |E beta_j| / |E beta_k| to the matched program k (NA when no program matches).
#' @param task Row with cell, share, source, seed.
#' @param cfg Parsed config.
#' @return Status string; the checkpoint holds the metric row and a compact fit state.
#' @examples
#' # a2_sim_task(task, cfg)
a2_sim_task <- function(task, cfg) {
  key <- sprintf("sim_%s_s%02d_SZ", task$cell, task$seed)
  path <- file.path(PRIV2, "sim", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  d <- lfp_sim_data(task, task$seed)
  sh <- a2_shift(s9_zscore(d$Y))
  Ytr <- sh$Y
  Yva <- sweep(s9_zscore(d$Yv), 2, sh$c, "+")    # validation shifted by the training constants
  set.seed(task$seed + 42L)
  run <- lfp_try(a2_fit(Ytr, d$time, d$event, "SZ", cfg))
  m <- s9_score(run$value, d, Ytr, Yva)
  fa <- data.frame(n_false_active_beta = NA_integer_, max_rel_beta_false = NA_real_,
                   mean_rel_beta_false = NA_real_, abs_beta_match = NA_real_)
  live <- NA_integer_
  if (!inherits(run$value, "error")) {
    fit <- run$value
    live <- sum(colSums(fit$EF^2) > 1e-8)
    r <- suppressWarnings(cor(d$truth, fit$EF)); r[!is.finite(r)] <- 0
    k <- which.max(abs(r)); matched <- abs(r[k]) >= .7
    sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
    active <- abs(fit$EBeta) >= 1.96 * sdb & sdb > 0
    fal <- setdiff(which(active), if (matched) k else integer())
    rel <- if (matched && length(fal)) abs(fit$EBeta[fal]) / abs(fit$EBeta[k]) else numeric()
    fa <- data.frame(n_false_active_beta = length(fal),
                     max_rel_beta_false = if (length(rel)) max(rel) else NA_real_,
                     mean_rel_beta_false = if (length(rel)) mean(rel) else NA_real_,
                     abs_beta_match = if (matched) abs(fit$EBeta[k]) else NA_real_)
  }
  row <- cbind(task[c("cell", "share", "seed")], setup = "SZ", m, fa, n_alive = live,
               fit_sec = run$seconds, warnings = run$warnings)
  st <- if (inherits(run$value, "error")) NULL else run$value[c("EF", "EBeta", "EBeta2", "EF_norms", "Tau")]
  saveRDS(list(metrics = row, state = st), path)
  sprintf("%s ok %.0fs", key, run$seconds)
}

# Real data ----

#' Fit and score one (setup, fold) real-data task; fold 0 is the full-data fit
#'
#' Records what task A recorded: convergence, programs with nonzero loadings, survival-active
#' programs, held-out concordances (within-cohort, pooled, per cohort), cross-validated partial
#' log-likelihood, and per program the correlation of its training scores Y F_k with each
#' patient's mean log expression and with the TCGA indicator.
#' @param task Row with setup and fold. @param z Prepared real inputs (REAL_IN).
#' @param Ysz Shifted z-score matrix (273 x genes). @param cfg Config.
#' @return Status string; checkpoint (data/) holds metrics, program table, state and risk.
#' @examples
#' # a2_real_task(task, z, Ysz, cfg)
a2_real_task <- function(task, z, Ysz, cfg) {
  key <- sprintf("real_%s_f%d", task$setup, task$fold)
  path <- file.path(PRIV2, "real", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  Y <- if (task$setup == "SZ") Ysz else z$Y$PE
  te <- which(z$fold == task$fold)
  tr <- if (task$fold == 0L) seq_along(z$time) else which(z$fold != task$fold)
  set.seed(42L)
  run <- lfp_try(a2_fit(Y[tr, , drop = FALSE], z$time[tr], z$status[tr], task$setup, cfg,
                        strata = z$cohort[tr]))
  row <- data.frame(setup = task$setup, fold = task$fold, n_train = length(tr), n_test = length(te),
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
  prog <- data.frame(setup = task$setup, fold = task$fold, program = seq_len(ncol(fit$EF)),
    alive = alive, active = active, beta = fit$EBeta, beta_sd = sdb,
    n_nonzero_genes = colSums(abs(fit$EF) > 1e-8),
    cor_mean_expr = apply(Ztr, 2, safe_cor, b = meanexpr),
    cor_cohort = apply(Ztr, 2, safe_cor, b = tcga))
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
#' @param smoke If TRUE, real fold 1 for each setup and one synthetic dataset (1% cell, seed 1).
#' @return list(real, sim) of design data frames.
#' @examples
#' a2_design(TRUE)
a2_design <- function(smoke = FALSE) {
  if (smoke) return(list(real = data.frame(setup = c("SZ", "PES"), fold = 1L, stringsAsFactors = FALSE),
                         sim = cbind(SIM_CELLS[2, ], seed = 1L)))
  list(real = expand.grid(setup = c("SZ", "PES"), fold = 0:5, stringsAsFactors = FALSE),
       sim = merge(SIM_CELLS, data.frame(seed = 1:10)))
}

#' Run every task on a PSOCK cluster, keeping failures as status lines
#' @param design a2_design() output. @param cores Worker count (default 2).
#' @return Character status vector.
#' @examples
#' # a2_run(a2_design(TRUE), 2)
a2_run <- function(design, cores) {
  dir.create(file.path(PRIV2, "sim"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PRIV2, "real"), recursive = TRUE, showWarnings = FALSE)
  tasks <- c(lapply(split(design$real, seq_len(nrow(design$real))), function(t) c(type = "real", t)),
             lapply(split(design$sim, seq_len(nrow(design$sim))), function(t) c(type = "sim", t)))
  cl <- parallel::makePSOCKcluster(cores, outfile = file.path(PRIV2, "workers.log"))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, c("PRIV2", "REAL_IN", "SIM_CELLS", "lfp_load", "lfp_try",
    "lfp_sim_data", "lfp_concordance", "a2_shift", "a2_fit", "a2_sim_task", "a2_real_task"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1); lfp_load(); cfg <- yaml::read_yaml("config/globals.yml")
    z <- readRDS(REAL_IN); Ysz <- a2_shift(z$Y$PL)$Y; NULL })
  st <- parallel::parLapplyLB(cl, tasks, function(t) tryCatch({
    t <- as.data.frame(t, stringsAsFactors = FALSE)
    if (t$type == "real") { t$fold <- as.integer(t$fold); a2_real_task(t, z, Ysz, cfg) }
    else { t$seed <- as.integer(t$seed); t$share <- as.numeric(t$share); a2_sim_task(t, cfg) }
  }, error = function(e) paste("TASK ERROR", conditionMessage(e))), chunk.size = 1L)
  unlist(st)
}

# Summary ----

#' Paired fold differences of a setup's held-out concordance against a task A arm
#' @param cv A2 fold rows (fold > 0) for one setup. @param ref Task A fold rows for one arm.
#' @param ref_arm Label of the reference arm.
#' @return One-row data frame: means, mean paired difference and its SE (SD / sqrt(folds)).
#' @examples
#' # a2_paired(cv, ref, "PL")
a2_paired <- function(cv, ref, ref_arm) {
  x <- merge(cv[c("fold", "C_within", "C_pooled", "converged")],
             ref[c("fold", "C_within", "C_pooled", "converged")], by = "fold", suffixes = c("", "_ref"))
  dd <- x$C_within - x$C_within_ref; dp <- x$C_pooled - x$C_pooled_ref
  data.frame(setup = cv$setup[1], reference = ref_arm, n_folds = sum(is.finite(dd)),
    mean_C_within = mean(x$C_within), mean_C_within_ref = mean(x$C_within_ref),
    diff_within = mean(dd), se_within = sd(dd) / sqrt(sum(is.finite(dd))),
    mean_C_pooled = mean(x$C_pooled), mean_C_pooled_ref = mean(x$C_pooled_ref),
    diff_pooled = mean(dp), se_pooled = sd(dp) / sqrt(sum(is.finite(dp))),
    folds_converged = sum(x$converged), folds_converged_ref = sum(x$converged_ref))
}

#' Within-cohort correlation of each program's training scores with mean expression
#'
#' Added after the smoke test: under PES the uncentered scores Y F_k keep each cohort's offset,
#' which the strata absorb, so A's correlation with the TCGA indicator cannot show whether a
#' program carries anything beyond that offset. Here both the score Y F_k and the patient's mean
#' log expression are centered within cohort before correlating. Computed from the saved E[F]
#' for SZ and PES and for task A's PL and PE checkpoints (read only).
#' @return Data frame: setup, fold, program, cor_mean_expr_within.
#' @examples
#' # a2_within_cor()
a2_within_cor <- function() {
  z <- readRDS(REAL_IN)
  Ys <- list(SZ = a2_shift(z$Y$PL)$Y, PES = z$Y$PE, PL = z$Y$PL, PE = z$Y$PE)
  files <- c(list.files(file.path(PRIV2, "real"), full.names = TRUE, pattern = "rds$"),
             list.files(file.path(PRIV, "real"), full.names = TRUE, pattern = "^real_(PL|PE)_f[0-9]\\.rds$"))
  cw <- function(v, g) v - ave(v, g)
  do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f); if (is.null(x$state)) return(NULL)
    s <- if (!is.null(x$metrics$setup)) x$metrics$setup else x$metrics$arm
    fold <- x$metrics$fold
    tr <- if (fold == 0L) seq_along(z$time) else which(z$fold != fold)
    g <- z$cohort[tr]
    Z <- Ys[[s]][tr, , drop = FALSE] %*% x$state$EF
    me <- cw(rowMeans(z$X[tr, , drop = FALSE]), g)
    data.frame(setup = s, fold = fold, program = seq_len(ncol(Z)),
      cor_mean_expr_within = apply(Z, 2, function(v) { v <- cw(v, g)
        if (sd(v) < 1e-12) NA_real_ else cor(v, me) }))
  }))
}

#' Aggregate checkpoints into CSVs and apply task A's four criteria to each setup
#'
#' Criteria (as registered for A, applied literally): (1) recovery in every cell >= PL's minus
#' 1 of 10; (2) mean false-active programs per fit <= 1 in every cell; (3) >= 9 of 10 converged
#' in every cell; (4) mean within-cohort CV concordance not lower than PL's by more than one SE
#' of the paired fold difference, and every fold converges. A missing fit counts as not
#' recovered, not converged and 7 false-active programs. PES has no simulation fits of its own:
#' with one cohort its fit is task A's PE fit, so criteria 1-3 use A's PE cells.
#' @return NULL; writes a2_*.csv to OUT2.
#' @examples
#' # a2_summary()
a2_summary <- function() {
  dir.create(OUT2, recursive = TRUE, showWarnings = FALSE)
  rd <- function(dir) lapply(list.files(file.path(PRIV2, dir), full.names = TRUE, pattern = "rds$"), readRDS)
  # Simulation
  sim <- lfp_bind(lapply(rd("sim"), `[[`, "metrics"))
  sim <- merge(a2_design()$sim[c("cell", "seed")], sim, all.x = TRUE)
  sim$missing <- is.na(sim$converged)
  sim$recovered[sim$missing] <- FALSE; sim$converged[sim$missing] <- FALSE
  sim$false_active[sim$missing] <- 7L
  write.csv(sim[setdiff(names(sim), "warnings")], file.path(OUT2, "a2_sim_fits.csv"), row.names = FALSE)
  cells <- do.call(rbind, lapply(split(sim, sim$cell), function(s)
    data.frame(cell = s$cell[1], arm = "SZ", n = nrow(s), recovered = sum(s$recovered),
      mean_false_active = mean(s$false_active), converged = sum(s$converged),
      median_cor = median(s$cor_planted, na.rm = TRUE),
      median_c_validation = median(s$c_validation, na.rm = TRUE),
      median_cal_slope = median(s$cal_slope, na.rm = TRUE),
      median_alive = median(s$n_alive, na.rm = TRUE), missing = sum(s$missing),
      median_fit_sec = median(s$fit_sec, na.rm = TRUE),
      fits_with_false_active = sum(s$false_active > 0),
      median_max_rel_beta_false = median(s$max_rel_beta_false, na.rm = TRUE),
      max_rel_beta_false = suppressWarnings(max(s$max_rel_beta_false, na.rm = TRUE)),
      mean_rel_beta_false = mean(s$mean_rel_beta_false, na.rm = TRUE))))
  a_cells <- read.csv(file.path(A_OUT, "sim_cells.csv"))
  all_cells <- lfp_bind(list(cells, a_cells[a_cells$arm %in% c("PL", "PE"), ]))
  write.csv(all_cells[order(all_cells$cell, all_cells$arm), ], file.path(OUT2, "a2_sim_cells.csv"),
            row.names = FALSE)
  # Real data
  real <- rd("real")
  rf <- lfp_bind(lapply(real, `[[`, "metrics"))
  rf <- merge(a2_design()$real, rf, all.x = TRUE)
  rf$converged[is.na(rf$converged)] <- FALSE
  write.csv(rf[setdiff(names(rf), "warnings")], file.path(OUT2, "a2_real_fits.csv"), row.names = FALSE)
  progs <- lfp_bind(lapply(real, `[[`, "programs"))
  wc <- a2_within_cor()
  progs <- merge(progs, wc[wc$setup %in% c("SZ", "PES"), ], all.x = TRUE)
  write.csv(progs, file.path(OUT2, "a2_real_programs.csv"), row.names = FALSE)
  write.csv(wc, file.path(OUT2, "a2_real_within_cohort_cor.csv"), row.names = FALSE)
  a_real <- read.csv(file.path(A_OUT, "real_fits.csv"))
  a_cv <- a_real[a_real$fold > 0, ]
  cv <- rf[rf$fold > 0, ]
  paired <- do.call(rbind, lapply(c("SZ", "PES"), function(s) do.call(rbind, lapply(c("PL", "PE"),
    function(a) a2_paired(cv[cv$setup == s, ], a_cv[a_cv$arm == a, ], a)))))
  write.csv(paired, file.path(OUT2, "a2_real_paired.csv"), row.names = FALSE)
  # Top-100 gene overlap of each full-data program with the current S-11 fit's programs
  cur <- readRDS(file.path(S11_FITS, "size_D4_N2218_w1_f0.rds"))
  zin <- readRDS(REAL_IN)
  top <- function(EF, k) order(-abs(EF[, k]))[1:100]
  cur_EF <- cur$state$EF[match(zin$genes, cur$genes), , drop = FALSE]
  cur_alive <- which(colSums(cur_EF^2) > 1e-8)
  ov <- do.call(rbind, lapply(real, function(x) {
    if (x$metrics$fold != 0 || is.null(x$state)) return(NULL)
    EF <- x$state$EF; al <- which(colSums(EF^2) > 1e-8)
    if (!length(al)) return(NULL)
    do.call(rbind, lapply(al, function(k) data.frame(setup = x$metrics$setup, program = k,
      current_program = cur_alive,
      overlap_top100 = vapply(cur_alive, function(j) length(intersect(top(EF, k), top(cur_EF, j))), 1L))))
  }))
  if (!is.null(ov)) write.csv(ov, file.path(OUT2, "a2_real_gene_overlap.csv"), row.names = FALSE)
  # Criteria
  pl <- a_cells[a_cells$arm == "PL", ]
  crit <- do.call(rbind, lapply(c("SZ", "PES"), function(s) {
    ce <- if (s == "SZ") cells else transform(a_cells[a_cells$arm == "PE", ], arm = "PES")
    ce <- merge(ce, pl, by = "cell", suffixes = c("", "_PL"))
    p <- paired[paired$setup == s & paired$reference == "PL", ]
    src <- if (s == "SZ") "" else " (A's PE fits: one stratum)"
    data.frame(setup = s,
      criterion = c("1 recovery >= PL - 1 in every cell", "2 mean false-active <= 1 in every cell",
                    "3 >= 9/10 converged in every cell",
                    "4 CV within-cohort C not lower than PL by > 1 SE; all folds converge"),
      detail = c(paste0(paste(sprintf("%s %d vs %d", ce$cell, ce$recovered, ce$recovered_PL), collapse = "; "), src),
                 paste0(paste(sprintf("%s %.2f", ce$cell, ce$mean_false_active), collapse = "; "), src),
                 paste0(paste(sprintf("%s %d/10", ce$cell, ce$converged), collapse = "; "), src),
                 sprintf("diff %.4f, SE %.4f, folds converged %d/5", p$diff_within, p$se_within,
                         p$folds_converged)),
      pass = c(all(ce$recovered >= ce$recovered_PL - 1), all(ce$mean_false_active <= 1),
               all(ce$converged >= 9),
               isTRUE(p$diff_within >= -p$se_within && p$folds_converged == 5)))
  }))
  write.csv(crit, file.path(OUT2, "a2_criteria.csv"), row.names = FALSE)
  print(all_cells[order(all_cells$cell, all_cells$arm), ]); print(rf[c("setup", "fold", "converged",
    "n_iter", "n_alive", "n_active", "max_abs_cor_cohort", "max_abs_cor_mean_expr", "C_within", "C_pooled")])
  print(paired); print(crit)
  invisible(NULL)
}

# Main ----
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE); mode <- args[1]
  cores <- if (length(args) > 1) as.integer(args[2]) else 2L
  if (mode %in% c("smoke", "grid")) {
    dir.create(OUT2, recursive = TRUE, showWarnings = FALSE)
    writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()), mode,
                 paste("run_lf_prior_comparison.R md5 (working copy):",
                       tools::md5sum("results/benchmark_sim/run_lf_prior_comparison.R"))),
               file.path(OUT2, sprintf("provenance_%s.txt", mode)))
    st <- a2_run(a2_design(mode == "smoke"), cores)
    writeLines(st, file.path(OUT2, sprintf("task_status_%s.txt", mode))); print(st)
  }
  if (mode == "summary") { lfp_load(); a2_summary() }
}
