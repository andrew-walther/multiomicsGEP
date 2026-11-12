# ============================================================
# Script: run_training_cv_redo.R
# Purpose: Task F. Redo by pooled five-fold training CV the three choices that external
#          cohorts once informed: gene screen (D2 vs D4), input scaling, and rank K.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm, yaml, parallel; results/benchmark_sim/run_lf_prior_comparison.R
#               (helpers and task A's point-Laplace checkpoints), code/preprocess_desurv.R,
#               results/benchmark_sim/benchmark_helpers.R
# Run: Rscript results/benchmark_sim/run_training_cv_redo.R prepare
#      Rscript results/benchmark_sim/run_training_cv_redo.R smoke|grid [cores]
#      Rscript results/benchmark_sim/run_training_cv_redo.R brerun 3 5
#      Rscript results/benchmark_sim/run_training_cv_redo.R summary
# Registration: 2026-11-12 entry "F" in docs/plans/2026-10-16-round2-attempts-log.md.
# ============================================================
# Settled model: fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1), point-Laplace L and F,
# per-program Normal beta, no pruning, pooled TCGA + CPTAC, no strata. Every arm uses the
# S-11 folds (cohort x event, seed 2026), so it is paired with task A's point-Laplace fits,
# which are the reference (D4 screen, per-cohort z-score, K = 7) and are not refitted.
# Primary metric: within-cohort concordance (pairs only within TCGA or within CPTAC).
#
# Parallelization: PSOCK workers, one BLAS/OpenMP thread each, one checkpoint per fit.
# Patient-level inputs and fits go to data/ (git-ignored); OUT holds aggregate CSVs only.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_lf_prior_comparison.R")   # main block does not run
F_OUT   <- "results/benchmark_sim/outputs/training_cv_redo_2026-11-12"
F_PRIV  <- "data/training_cv_redo_2026-11-12"
F_IN    <- "data/training_cv_redo_inputs_2026-11-12.rds"
K_GRID  <- 2:12

# Inputs ----

#' Prepare the screen and scaling inputs in the S-11 patient order
#'
#' D4_Z is task A's PL input (per-cohort z-score of the 2,218 D4 genes). D2_Z is the
#' 2026-05-27 D2 preprocessing: per-cohort z-score, quantile normalization of the merged
#' matrix, top 2,000 genes by variance after normalization. D4_XSD is task A's X/sd input.
#' D4_RANKZ ranks the D4 genes within each patient and then z-scores each gene within
#' cohort. All are outcome-free and computed once on the 273 training patients.
#' @return NULL; writes F_IN (patient-level, git-ignored) and an aggregate audit CSV.
#' @examples
#' # f_prepare()
f_prepare <- function() {
  if (file.exists(F_IN)) stop("Prepared input exists; do not overwrite a registered design.")
  z <- readRDS(REAL_IN)
  cfg <<- yaml::read_yaml("config/globals.yml")
  source("results/benchmark_sim/benchmark_helpers.R"); source("code/preprocess_desurv.R")
  train <- cfg$pdac$training_cohorts
  raw <- lapply(setNames(train, train), load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
  stopifnot(identical(as.numeric(unlist(lapply(train, function(ds) raw[[ds]]$time))), as.numeric(z$time)))
  d2 <- preprocess_merged_cohorts(raw, PLATFORM_LOG_TRANSFORM[train], top_n = cfg$preprocessing$top_n_genes,
    rank_transform = FALSE, per_platform_standardize = TRUE, normalize_method = "quantile",
    selection_per_cohort = FALSE, selection_method = "variance")
  by_cohort_z <- function(M) { for (co in unique(z$cohort)) { i <- z$cohort == co
    M[i, ] <- scale(M[i, , drop = FALSE]) }; M[!is.finite(M)] <- 0; M }
  Y <- list(D4_Z = z$Y$PL, D2_Z = unname(d2$Y), D4_XSD = z$Y$PE,
            D4_RANKZ = by_cohort_z(rank_transform_subjects(z$X)))
  colnames(Y$D2_Z) <- d2$gene_names
  saveRDS(list(Y = Y, genes = list(D4 = z$genes, D2 = d2$gene_names), X = z$X, time = z$time,
               status = z$status, cohort = z$cohort, fold = z$fold), F_IN)
  dir.create(F_OUT, recursive = TRUE, showWarnings = FALSE)
  audit <- data.frame(input = names(Y), n_genes = vapply(Y, ncol, 1L),
    min = vapply(Y, min, 1), max = vapply(Y, max, 1),
    overlap_with_D4 = vapply(names(Y), function(n) length(intersect(colnames(Y[[n]]), z$genes)), 1L))
  write.csv(audit, file.path(F_OUT, "input_audit.csv"), row.names = FALSE)
  writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time())),
             file.path(F_OUT, "provenance_prepare.txt"))
  print(audit)
  invisible(NULL)
}

# Fitting ----

#' Fit the settled model at rank K on one input
#' @param Y Training matrix. @param time,status Outcomes. @param K Rank.
#' @param center TRUE for uncentered input (X/sd), which needs center_projection.
#' @param cfg Parsed config.
#' @param prior_LF L/F prior; "point_exponential" only for the task A rerun at a new K.
#' @param prior_beta Coefficient prior, default Normal.
#' @param beta_prior_mode Per-program by default; shared for the task B rerun.
#' @return fit_ssbmf_yfb() result.
#' @examples
#' # f_fit(Y, time, status, 7, FALSE, cfg)
f_fit <- function(Y, time, status, K, center, cfg, prior_LF = "point_laplace",
                  prior_beta = "normal", beta_prior_mode = "per_program") {
  fit_ssbmf_yfb(Y, time, status, K = K, max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = prior_LF, prior_beta = prior_beta, beta_prior_mode = beta_prior_mode,
                canonicalize = TRUE, w_S = 1,
                alpha = cfg$benchmark$alpha, verbose = FALSE, trace_supervision = TRUE,
                center_projection = center)
}

#' Fit and score one (part, input, K, fold) task; fold 0 is the full-data fit
#'
#' Same metrics as task A: convergence, programs with nonzero loadings, survival-active
#' programs (|E beta| >= 1.96 SD), held-out within-cohort, pooled and per-cohort
#' concordance, held-out partial log-likelihood, and program correlations with the TCGA
#' indicator and with each patient's mean log expression.
#' Task B parts b_rerun_PLAP, b_rerun_PN and b_rerun_NSH use shared beta priors.
#' Point-normal inclusion probabilities are saved with each checkpoint.
#' @param task Row with part, input, K, fold. @param z Prepared inputs. @param cfg Config.
#' @return Status string.
#' @examples
#' # f_task(task, z, cfg)
f_task <- function(task, z, cfg) {
  key <- sprintf("%s_%s_K%02d_f%d", task$part, task$input, task$K, task$fold)
  path <- file.path(F_PRIV, paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  Y <- z$Y[[task$input]]
  te <- which(z$fold == task$fold)
  tr <- if (task$fold == 0L) seq_along(z$time) else which(z$fold != task$fold)
  set.seed(42L)
  beta <- switch(task$part, b_rerun_PLAP = "point_laplace",
                 b_rerun_PN = "point_normal", b_rerun_NSH = "normal", "normal")
  shared <- startsWith(task$part, "b_rerun_")
  run <- lfp_try(f_fit(Y[tr, , drop = FALSE], z$time[tr], z$status[tr], task$K,
                       task$input == "D4_XSD", cfg,
                       if (task$part == "a_rerun") "point_exponential" else "point_laplace",
                       prior_beta = beta, beta_prior_mode = if (shared) "shared" else "per_program"))
  row <- data.frame(part = task$part, input = task$input, K = task$K, fold = task$fold,
                    n_genes = ncol(Y), fit_sec = run$seconds, warnings = run$warnings,
                    error = NA_character_)
  if (inherits(run$value, "error")) {
    row$error <- conditionMessage(run$value); row$converged <- FALSE
    saveRDS(list(metrics = row), path); return(paste(key, "ERROR", row$error))
  }
  fit <- run$value
  pip <- if (task$part == "b_rerun_PN") lfp_pip(fit) else NULL
  sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  alive <- colSums(fit$EF^2) > 1e-8
  active <- abs(fit$EBeta) >= 1.96 * sdb & sdb > 0 & alive
  Ztr <- Y[tr, , drop = FALSE] %*% fit$EF
  tcga <- as.numeric(grepl("TCGA", z$cohort[tr])); meanexpr <- rowMeans(z$X[tr, , drop = FALSE])
  sc <- function(a, b) if (sd(a) < 1e-12) NA_real_ else cor(a, b)
  cc <- apply(Ztr, 2, sc, b = tcga); cm <- apply(Ztr, 2, sc, b = meanexpr)
  row <- cbind(row, data.frame(converged = isTRUE(fit$history$converged), n_iter = fit$history$n_iter,
    n_alive = sum(alive), n_active = sum(active), max_abs_beta = max(abs(fit$EBeta)),
    max_abs_cor_cohort = suppressWarnings(max(abs(cc[alive]), na.rm = TRUE)),
    max_abs_cor_mean_expr = suppressWarnings(max(abs(cm[alive]), na.rm = TRUE))))
  risk <- NULL
  if (task$fold > 0L) {
    risk <- drop(Y[te, , drop = FALSE] %*% fit$EF %*% fit$EBeta)
    row <- cbind(row, lfp_concordance(risk, z$time[te], z$status[te], z$cohort[te]))
    eta <- drop(Y %*% fit$EF %*% fit$EBeta)
    row$cvpl_pooled <- calc_cox_taylor_yf(eta, z$time, z$status)$logPL -
                       calc_cox_taylor_yf(eta[tr], z$time[tr], z$status[tr])$logPL
  }
  saveRDS(list(metrics = row, te = te, risk = risk,
               state = list(EF = fit$EF, EBeta = fit$EBeta, EBeta2 = fit$EBeta2), pip = pip), path)
  sprintf("%s ok %.0fs", key, run$seconds)
}

#' Task list for the three parts (task A's PL fits are the D4_Z, K = 7 reference)
#' @param smoke If TRUE, fold 1 only of D2_Z, D4_RANKZ and K = 3.
#' @return Data frame of tasks, larger fits first.
#' @examples
#' f_design(TRUE)
f_design <- function(smoke = FALSE) {
  if (smoke) return(data.frame(part = c("screen", "scaling", "rank"),
    input = c("D2_Z", "D4_RANKZ", "D4_Z"), K = c(7L, 7L, 3L), fold = 1L, stringsAsFactors = FALSE))
  d <- rbind(
    expand.grid(part = "screen", input = "D2_Z", K = 7L, fold = 0:5, stringsAsFactors = FALSE),
    expand.grid(part = "scaling", input = c("D4_XSD", "D4_RANKZ"), K = 7L, fold = 0:5, stringsAsFactors = FALSE),
    expand.grid(part = "rank", input = "D4_Z", K = setdiff(K_GRID, 7L), fold = 0:5, stringsAsFactors = FALSE))
  d[order(-d$K), ]
}

#' Run tasks on a PSOCK cluster, keeping failures as status lines
#' @param design f_design() output. @param cores Worker count.
#' @return Character status vector.
#' @examples
#' # f_run(f_design(TRUE), 5)
f_run <- function(design, cores) {
  dir.create(F_PRIV, recursive = TRUE, showWarnings = FALSE)
  dir.create(F_OUT, recursive = TRUE, showWarnings = FALSE)
  cl <- parallel::makePSOCKcluster(cores, outfile = file.path(F_PRIV, "workers.log"))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, c("F_PRIV", "F_IN", "f_fit", "f_task"))
  parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1); source("results/benchmark_sim/run_lf_prior_comparison.R")
    lfp_load(); cfg <- yaml::read_yaml("config/globals.yml"); z <- readRDS(F_IN); NULL })
  unlist(parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(t)
    tryCatch(f_task(t, z, cfg), error = function(e) paste("TASK ERROR", conditionMessage(e))),
    chunk.size = 1L))
}

# Summary ----

#' Aggregate F fits with task A's PL reference and apply the registered rules
#'
#' Screen and scaling: the current choice is kept unless an alternative's mean
#' within-cohort C exceeds it by more than one SE of the paired fold difference.
#' K: one-SE rule, smallest K with mean within-cohort C >= best mean - SE(fold C at best K).
#' @return NULL; writes fits.csv, curves.csv, paired.csv, k_curve.csv and choices.csv.
#' @examples
#' # f_summary()
f_summary <- function() {
  rd <- function(files) lapply(files, function(f) readRDS(f)$metrics)
  fr <- lfp_bind(rd(list.files(F_PRIV, full.names = TRUE, pattern = "rds$")))
  fr <- fr[!startsWith(fr$part, "b_rerun_"), , drop = FALSE]
  ref <- lfp_bind(rd(file.path(PRIV, "real", sprintf("real_PL_f%d.rds", 0:5))))
  ref <- data.frame(part = "reference", input = "D4_Z", K = 7L, fold = ref$fold, n_genes = 2218L,
                    ref[setdiff(names(ref), c("arm", "fold", "n_train", "n_test"))])
  all <- lfp_bind(list(fr, ref))
  all$converged[is.na(all$converged)] <- FALSE
  keep <- c("part", "input", "K", "fold", "n_genes", "converged", "n_iter", "n_alive", "n_active",
            "max_abs_beta", "max_abs_cor_cohort", "max_abs_cor_mean_expr", "C_within", "C_pooled",
            "C_TCGA", "C_CPTAC", "cvpl_pooled", "fit_sec", "error")
  write.csv(all[intersect(keep, names(all))], file.path(F_OUT, "fits.csv"), row.names = FALSE)
  cv <- all[all$fold > 0, ]; full <- all[all$fold == 0, ]
  curve <- do.call(rbind, lapply(split(cv, list(cv$input, cv$K), drop = TRUE), function(s)
    data.frame(input = s$input[1], K = s$K[1], n_folds = sum(is.finite(s$C_within)),
      mean_C_within = mean(s$C_within), se_C_within = sd(s$C_within) / sqrt(nrow(s)),
      mean_C_pooled = mean(s$C_pooled), mean_cvpl = mean(s$cvpl_pooled),
      folds_converged = sum(s$converged), mean_alive_cv = mean(s$n_alive),
      mean_active_cv = mean(s$n_active))))
  curve <- merge(curve, data.frame(input = full$input, K = full$K, full_alive = full$n_alive,
                                   full_active = full$n_active, full_converged = full$converged), all.x = TRUE)
  write.csv(curve[order(curve$input, curve$K), ], file.path(F_OUT, "curves.csv"), row.names = FALSE)
  base <- cv[cv$input == "D4_Z" & cv$K == 7L, c("fold", "C_within", "C_pooled")]
  paired <- do.call(rbind, lapply(c("D2_Z", "D4_XSD", "D4_RANKZ"), function(a) {
    x <- merge(cv[cv$input == a & cv$K == 7L, c("fold", "C_within", "C_pooled")], base,
               by = "fold", suffixes = c("", "_ref"))
    d <- x$C_within - x$C_within_ref; dp <- x$C_pooled - x$C_pooled_ref
    data.frame(input = a, part = if (a == "D2_Z") "screen" else "scaling", n_folds = nrow(x),
      diff_within = mean(d), se_within = sd(d) / sqrt(nrow(x)), diff_pooled = mean(dp),
      se_pooled = sd(dp) / sqrt(nrow(x)), replaces_current = mean(d) > sd(d) / sqrt(nrow(x)))
  }))
  write.csv(paired, file.path(F_OUT, "paired.csv"), row.names = FALSE)
  kc <- curve[curve$input == "D4_Z", ]; kc <- kc[order(kc$K), ]
  best <- which.max(kc$mean_C_within)
  thr <- kc$mean_C_within[best] - kc$se_C_within[best]
  k1se <- min(kc$K[kc$mean_C_within >= thr])
  # Paired differences of each K against K = 7 (same folds)
  kp <- do.call(rbind, lapply(kc$K, function(k) {
    x <- merge(cv[cv$input == "D4_Z" & cv$K == k, c("fold", "C_within")], base[c("fold", "C_within")],
               by = "fold", suffixes = c("", "_7"))
    data.frame(K = k, diff_vs_K7 = mean(x$C_within - x$C_within_7),
               se_diff_vs_K7 = sd(x$C_within - x$C_within_7) / sqrt(nrow(x)))
  }))
  kc <- merge(kc, kp); kc$threshold_1se <- thr; kc$selected <- kc$K == k1se; kc$best <- kc$K == kc$K[best]
  write.csv(kc, file.path(F_OUT, "k_curve.csv"), row.names = FALSE)
  choices <- data.frame(choice = c("gene screen", "input scaling", "rank K"),
    current = c("D4", "per-cohort z-score", "7"),
    selected = c(if (paired$replaces_current[paired$input == "D2_Z"]) "D2" else "D4",
                 { s <- paired[paired$part == "scaling" & paired$replaces_current, ]
                   if (nrow(s)) s$input[which.max(s$diff_within)] else "per-cohort z-score" },
                 as.character(k1se)),
    rule = c("alternative replaces current only if better by > 1 SE (paired)",
             "alternative replaces current only if better by > 1 SE (paired)",
             sprintf("one-SE: best K = %d (mean %.4f, SE %.4f), threshold %.4f",
                     kc$K[kc$best], max(kc$mean_C_within), kc$se_C_within[kc$best], thr)))
  write.csv(choices, file.path(F_OUT, "choices.csv"), row.names = FALSE)
  # Task A rerun at the selected K: point-exponential (X/sd) vs point-Laplace (z-score)
  ar <- cv[cv$part == "a_rerun", ]
  if (nrow(ar)) {
    k <- ar$K[1]
    x <- merge(ar[c("fold", "C_within", "converged", "n_active", "max_abs_cor_cohort")],
               cv[cv$input == "D4_Z" & cv$K == k, c("fold", "C_within", "converged", "n_active")],
               by = "fold", suffixes = c("_PE", "_PL"))
    d <- x$C_within_PE - x$C_within_PL
    write.csv(data.frame(K = k, mean_C_within_PE = mean(x$C_within_PE),
      mean_C_within_PL = mean(x$C_within_PL), diff = mean(d), se = sd(d) / sqrt(nrow(x)),
      folds_converged_PE = sum(x$converged_PE), folds_converged_PL = sum(x$converged_PL),
      mean_active_PE = mean(x$n_active_PE), mean_active_PL = mean(x$n_active_PL),
      max_abs_cor_cohort_PE = max(x$max_abs_cor_cohort), criterion4_pass = mean(d) >= -sd(d) / sqrt(nrow(x)) &&
      all(x$converged_PE)), file.path(F_OUT, "a_rerun_paired.csv"), row.names = FALSE)
  }
  print(curve[order(curve$input, curve$K), ]); print(paired); print(kc); print(choices)
  invisible(NULL)
}

#' Summarize the real-data beta-prior rerun against saved rank-sweep Normal fits
#'
#' Pairs folds explicitly. Rule (c) passes if mean(C_arm - C_Normal) >= -SE of
#' that paired difference. Simulation rules (a) and (b) are not evaluated here.
#' Missing, failed or nonfinite fold results stop the summary rather than dropping folds.
#' @param K Rerun rank.
#' @return Data frame, also written to b_rerun_paired.csv. PIPs are full-fit values.
#' @examples
#' # f_b_summary(5)
f_b_summary <- function(K) {
  parts <- c("rank", "b_rerun_PLAP", "b_rerun_PN", "b_rerun_NSH")
  runs <- lapply(parts, function(part) lapply(0:5, function(fold) {
    path <- file.path(F_PRIV, sprintf("%s_D4_Z_K%02d_f%d.rds", part, K, fold))
    x <- readRDS(path)
    if (!is.na(x$metrics$error) || is.null(x$state)) stop("Failed fit: ", path)
    stopifnot(x$metrics$fold == fold, x$metrics$K == K)
    if (fold > 0L && !all(is.finite(unlist(x$metrics[c("C_within", "C_pooled")]))))
      stop("Nonfinite concordance: ", path)
    x
  }))
  ref <- do.call(rbind, lapply(runs[[1]][-1], function(x) x$metrics))
  out <- do.call(rbind, lapply(seq_along(parts), function(i) {
    cv <- do.call(rbind, lapply(runs[[i]][-1], function(x) x$metrics))
    stopifnot(identical(cv$fold, ref$fold))
    d <- cv$C_within - ref$C_within; dp <- cv$C_pooled - ref$C_pooled
    se <- sd(d) / sqrt(5)
    pip <- runs[[i]][[1]]$pip
    if (parts[i] == "b_rerun_PN" && (length(pip) != K || any(!is.finite(pip))))
      stop("Missing full-fit point-normal inclusion probabilities")
    row <- data.frame(arm = parts[i], K = K, n_folds = nrow(cv),
      mean_C_within = mean(cv$C_within), mean_C_pooled = mean(cv$C_pooled),
      diff_within = mean(d), se_within = se, diff_pooled = mean(dp),
      se_pooled = sd(dp) / sqrt(5), folds_converged = sum(cv$converged),
      mean_active_cv = mean(cv$n_active), folds_all_beta_zero = sum(vapply(
        runs[[i]][-1], function(x) all(x$state$EBeta == 0), TRUE)),
      full_converged = runs[[i]][[1]]$metrics$converged,
      full_active = runs[[i]][[1]]$metrics$n_active,
      fit_warnings = sum(vapply(runs[[i]], function(x)
        !is.na(x$metrics$warnings) && nzchar(x$metrics$warnings), TRUE)),
      criterion_c_pass = mean(d) >= -se)
    for (k in seq_len(K)) row[[paste0("full_pip_", k)]] <- if (is.null(pip)) NA_real_ else pip[k]
    row
  }))
  write.csv(out, file.path(F_OUT, "b_rerun_paired.csv"), row.names = FALSE)
  print(out)
  invisible(out)
}

# Main ----
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE); mode <- args[1]
  cores <- if (length(args) > 1) as.integer(args[2]) else 5L
  if (mode == "prepare") { f_prepare(); quit(save = "no") }
  if (mode %in% c("smoke", "grid")) {
    dir.create(F_OUT, recursive = TRUE, showWarnings = FALSE)
    writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()), mode),
               file.path(F_OUT, sprintf("provenance_%s.txt", mode)))
    st <- f_run(f_design(mode == "smoke"), cores)
    writeLines(st, file.path(F_OUT, sprintf("task_status_%s.txt", mode))); print(st)
  }
  if (mode == "arerun") {
    # Task A's real-data CV at the newly selected K (registered: rerun once if K changes).
    # PE (X/sd, point-exponential, centered projection) at K; PL at K is the rank-sweep fit.
    k <- as.integer(args[3])
    st <- f_run(data.frame(part = "a_rerun", input = "D4_XSD", K = k, fold = 0:5,
                           stringsAsFactors = FALSE), cores)
    writeLines(st, file.path(F_OUT, sprintf("task_status_arerun_K%d.txt", k))); print(st)
  }
  if (mode == "brerun") {
    k <- if (length(args) > 2) as.integer(args[3]) else
      yaml::read_yaml("config/globals.yml")$benchmark$k_working_yfb
    stopifnot(length(k) == 1L, is.finite(k), k > 0L, cores >= 1L, cores <= 3L)
    st <- f_run(expand.grid(part = c("b_rerun_PLAP", "b_rerun_PN", "b_rerun_NSH"),
      input = "D4_Z", K = k, fold = 0:5, stringsAsFactors = FALSE), cores)
    writeLines(st, file.path(F_OUT, sprintf("task_status_brerun_K%d.txt", k))); print(st)
    if (any(grepl("ERROR", st, fixed = TRUE))) stop("Task B rerun contains failed fits")
    f_b_summary(k)
  }
  if (mode == "summary") { lfp_load(); f_summary() }
}
