# ============================================================
# Script: run_yfb_supervised_real.R
# Purpose: Attempt S-11: the survival-weighted YFB fit (w_S = 1) on the real training data:
#          step S3, the feature-set-size curve, survival screening, and a null check.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: survival, ebnm, yaml, parallel; code/fit_ssbmf_yfb.R, code/preprocess_desurv.R,
#               code/feature_sensitivity.R, code/concordance_ci.R, code/predict_ssbmf_yfb.R,
#               results/benchmark_sim/benchmark_helpers.R,
#               results/feature_selection/feature_selection_helpers.R,
#               results/benchmark_sim/run_yfb_survival_weight_p2000.R (null mode)
# Run: Rscript results/benchmark_sim/run_yfb_supervised_real.R prepare
#      Rscript results/benchmark_sim/run_yfb_supervised_real.R size|screen|null [cores]
#      Rscript results/benchmark_sim/run_yfb_supervised_real.R summary
# ============================================================
# Registration: docs/plans/2026-11-04-supervised-real-data-plan.md. Training patients only
# (TCGA_PAAD, CPTAC); no external cohort is loaded. The prepared matrix is patient-level
# and is written to data/ (git-ignored), never to the output folder.
#
# Parallelization: PSOCK workers (forked workers segfaulted on macOS after the compiled
# sweep and BLAS were loaded), one BLAS thread each, one checkpoint per fit, so a rerun
# skips finished fits. Memory: the all-genes matrix is 273 x ~17,900 doubles (~40 MB)
# per worker.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
# EXHIBIT_K (results/benchmark_sim/exhibit_rank.R) reruns at another rank into a separate
# folder; unset, the rank and paths are the historical K = 7 ones. The prepared input LOCAL
# does not depend on K and is shared.
source("results/benchmark_sim/exhibit_rank.R")
OUT   <- exhibit_out("results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04", "s11")
LOCAL <- "data/yfb_supervised_real_2026-11-04.rds"
SIZE_TARGETS <- c(3000L, 5000L, 8000L, 12000L)   # between the D4 size and all shared genes
SCREEN_N     <- c(250L, 500L, 1000L)
ARMS_WS      <- c(w1 = 1, w0 = 0)
NFOLD        <- 5L

# Code loading ----

#' Source the fitter and helpers (also run inside each worker)
#' @return NULL, invisibly.
#' @examples
#' # s11_load()
s11_load <- function() {
  suppressPackageStartupMessages({ library(survival); library(ebnm); library(yaml) })
  suppressMessages(tryCatch(source("code/fit_ssbmf_yfb.R"), error = function(e) {
    if (conditionMessage(e) != "!is.null(real_Y) is not TRUE") stop(e)
  }))
  source("code/predict_ssbmf_yfb.R"); source("code/concordance_ci.R")
  source("code/feature_sensitivity.R")
  source("results/feature_selection/feature_selection_helpers.R")
  multimodal_yfb_use_cpp()
  invisible(NULL)
}

#' The S-9 fixed call: canonical survival-weighted YFB fit
#' @param Y Centered patients x genes matrix. @param time,status Survival outcome.
#' @param w_S Survival weight in the F update (1 = the S-9 passing weight; 0 = unsupervised loadings).
#' @param cfg Parsed config/globals.yml.
#' @return fit_ssbmf_yfb() result.
#' @examples
#' # s11_fit(Y, time, status, 1, cfg)
s11_fit <- function(Y, time, status, w_S, cfg) {
  fit_ssbmf_yfb(Y, time, status, K = exhibit_fit_k(cfg),
                max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
                w_S = w_S, alpha = cfg$benchmark$alpha, verbose = FALSE,
                trace_supervision = TRUE)
}

# Prepare ----

#' Build and save the single-log, per-cohort z-scored training matrix and the design
#'
#' Gene ranks are the outcome-free combined mean/variance ranks per cohort; folds are
#' stratified by cohort x event (set.seed(2026)), as in run_a1_real_ws_cv.R.
#' @return NULL; writes LOCAL (patient-level, git-ignored) and design CSVs to OUT.
#' @examples
#' # s11_prepare()
s11_prepare <- function() {
  if (file.exists(LOCAL)) stop("Prepared input exists; do not overwrite a registered design.")
  cfg <<- yaml::read_yaml("config/globals.yml")   # benchmark_helpers.R reads cfg globally
  source("results/benchmark_sim/benchmark_helpers.R"); source("code/preprocess_desurv.R")
  source("code/feature_sensitivity.R")
  train <- cfg$pdac$training_cohorts
  stopifnot(!any(unlist(PLATFORM_LOG_TRANSFORM[train])))   # single log: loaders return log scale
  raw <- lapply(setNames(train, train), load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
  common <- Reduce(intersect, lapply(raw, function(z) z$gene_names))
  ranked <- lapply(train, function(ds) {
    x <- raw[[ds]]$Y[, match(common, raw[[ds]]$gene_names), drop = FALSE]
    colnames(x) <- common; x })
  names(ranked) <- train
  d4 <- preprocess_merged_cohorts(raw, PLATFORM_LOG_TRANSFORM[train],
    top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
    per_platform_standardize = TRUE, normalize_method = "none",
    selection_per_cohort = TRUE, selection_method = "combined_rank")
  all <- preprocess_merged_cohorts(raw, PLATFORM_LOG_TRANSFORM[train], top_n = length(common),
    rank_transform = FALSE, per_platform_standardize = TRUE, normalize_method = "none",
    selection_per_cohort = TRUE, selection_method = "combined_rank")
  targets <- c(length(d4$gene_names), SIZE_TARGETS, length(common))
  man <- feature_count_manifest(ranked, common, targets, cfg$preprocessing$top_n_genes_desurv)
  stopifnot(setequal(man$genes[[1]], d4$gene_names), all(is.finite(all$Y)),
            identical(unname(d4$Y), unname(all$Y[, match(d4$gene_names, all$gene_names)])))
  time   <- unlist(lapply(train, function(ds) raw[[ds]]$time))
  status <- unlist(lapply(train, function(ds) raw[[ds]]$status))
  cohort <- rep(train, vapply(train, function(ds) raw[[ds]]$n, 1L))
  set.seed(2026)
  fold <- integer(length(time))
  for (co in unique(cohort)) for (st in unique(status)) {
    idx <- which(cohort == co & status == st)
    if (length(idx)) fold[idx] <- sample(rep(seq_len(NFOLD), length.out = length(idx)))
  }
  sets <- setNames(lapply(man$genes, function(g) g[order(match(g, all$gene_names))]),
                   sprintf("size%d", man$design$actual_shared))
  names(sets)[1] <- "D4"
  saveRDS(list(Y = all$Y, genes = all$gene_names, time = time, status = status, cohort = cohort,
               fold = fold, sets = sets, entry_rank = man$entry_rank, d4 = d4$gene_names), LOCAL)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(set = names(sets), n_genes = lengths(sets)), file.path(OUT, "sizes.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(cohort = cohort, event = status, fold = fold)),
            file.path(OUT, "folds.csv"), row.names = FALSE)
  writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()),
               sprintf("n = %d, events = %d, shared genes = %d", length(time), sum(status),
                       length(common))), file.path(OUT, "provenance.txt"))
  print(data.frame(set = names(sets), n_genes = lengths(sets)))
  invisible(NULL)
}

# One fit ----

#' Genes for one task; S(N) is screened on the fold-training patients only
#' @param z Prepared data. @param task One design row (set, N, fold).
#' @return Character gene vector.
#' @examples
#' # s11_genes(z, task)
s11_genes <- function(z, task) {
  if (task$part == "size") return(z$sets[[task$set]])
  if (task$set == "U") return(unsupervised_gene_set(z$entry_rank[z$d4], task$N)$genes)
  # S(N): D4 genes ranked by the cohort-stratified Cox Wald p, top N matched to U(N)'s count
  nU <- unsupervised_gene_set(z$entry_rank[z$d4], task$N)$n_shared_actual
  tr <- if (task$fold == 0L) seq_along(z$time) else which(z$fold != task$fold)
  idx <- match(z$d4, z$genes)
  s <- survival_screen_genes(z$Y[tr, idx, drop = FALSE], z$time[tr], z$status[tr], z$cohort[tr], nU)
  z$d4[s]
}

#' Fit one (part, set, N, arm, fold) task and save its checkpoint
#' @param task One design row. @param z Prepared data. @param cfg Config.
#' @return Status string.
#' @examples
#' # s11_task(task, z, cfg)
s11_task <- function(task, z, cfg) {
  key <- sprintf("%s_%s_N%d_%s_f%d", task$part, task$set, task$N, task$arm, task$fold)
  path <- file.path(OUT, "fits", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  genes <- s11_genes(z, task)
  Y <- z$Y[, match(genes, z$genes), drop = FALSE]
  te <- which(z$fold == task$fold); tr <- if (task$fold == 0L) seq_along(z$time) else which(z$fold != task$fold)
  set.seed(42L); t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(s11_fit(Y[tr, , drop = FALSE], z$time[tr], z$status[tr], ARMS_WS[[task$arm]], cfg),
                  error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  row <- data.frame(part = task$part, set = task$set, N = task$N, arm = task$arm, fold = task$fold,
                    n_genes = length(genes), error = NA_character_)
  if (inherits(fit, "error")) {
    row$error <- conditionMessage(fit)
    saveRDS(list(metrics = row), path); return(paste(key, "ERROR", row$error))
  }
  sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  alive <- colSums(fit$EF^2) > 1e-8
  active <- abs(fit$EBeta) >= 1.96 * sdb & sdb > 0 & alive
  shareA <- tail(fit$history$trace$shareA_median, 1)[1, ]
  Ytr <- Y[tr, , drop = FALSE]
  row <- cbind(row, data.frame(converged = fit$history$converged, n_iter = fit$history$n_iter,
    n_alive = sum(alive), n_active = sum(active), max_abs_beta = max(abs(fit$EBeta)),
    share_active_median = if (any(active)) median(shareA[active]) else NA_real_,
    share_alive_median = median(shareA[alive]),
    recon_train = sum((Ytr - fit$EL %*% t(fit$EF))^2) / sum(Ytr^2), fit_sec = secs,
    C = NA_real_, cvpl = NA_real_, risk_sd = NA_real_))
  risk <- NULL
  if (task$fold > 0L) {
    risk <- predict_ssbmf_yfb(Y[te, , drop = FALSE], fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
    eta <- as.numeric(predict_ssbmf_yfb(Y, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores)
    row$risk_sd <- sd(risk)
    row$C <- if (sd(risk) > 1e-8) frozen_reverse_cindex(risk, z$time[te], z$status[te]) else NA_real_
    row$cvpl <- calc_cox_taylor_yf(eta, z$time, z$status)$logPL -
                calc_cox_taylor_yf(eta[tr], z$time[tr], z$status[tr])$logPL
  }
  saveRDS(list(metrics = row, genes = genes, te = te, risk = risk,
               state = list(EF = fit$EF, EBeta = fit$EBeta, EBeta2 = fit$EBeta2, EF_norms = fit$EF_norms,
                            EL = if (task$fold == 0L) fit$EL else NULL, Tau = fit$Tau,
                            elbo = fit$history$elbo_full)), path)
  sprintf("%s ok %.0fs", key, secs)
}

#' Design rows for the size or screen part
#' @param part "size", "screen" or "d4" (the size part restricted to the D4 set, used to
#'   fit the reported 2,218-gene model first when rerunning at another rank).
#' @param z Prepared data.
#' @return Data frame of tasks, largest gene sets first.
#' @examples
#' # s11_design("size", z)
s11_design <- function(part, z) {
  if (part == "d4") return(subset(s11_design("size", z), set == "D4"))
  sets <- if (part == "size") data.frame(set = names(z$sets), N = lengths(z$sets))
          else expand.grid(set = c("U", "S"), N = SCREEN_N, stringsAsFactors = FALSE)
  d <- merge(merge(sets, data.frame(arm = names(ARMS_WS))), data.frame(fold = 0:NFOLD))
  d$part <- part
  d[order(-d$N, d$arm, d$fold), ]
}

# Null check (simulation) ----

#' Permuted-outcome null fit on one S-7 dataset at w_S = 1
#' @param row seed, p, variance. @return One-row data frame.
#' @examples
#' # s11_null_task(data.frame(seed = 1, p = 50, variance = .01))
s11_null_task <- function(row) {
  d <- yfb_decision_data(row$seed, row$p, row$variance, truth_mode = "shared")
  set.seed(row$seed + 9000L); perm <- sample(length(d$time))
  cfg <- yaml::read_yaml("config/globals.yml")
  fit <- tryCatch(s9_fit(s9_zscore(d$Y), d$time[perm], d$event[perm], 1, cfg), error = function(e) e)
  if (inherits(fit, "error"))
    return(data.frame(row, n_active = 7L, converged = FALSE, n_iter = NA, max_abs_beta = NA,
                      error = conditionMessage(fit)))
  sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  alive <- colSums(fit$EF^2) > 1e-8
  data.frame(row, n_active = sum(abs(fit$EBeta) >= 1.96 * sdb & sdb > 0 & alive),
             converged = fit$history$converged, n_iter = fit$history$n_iter,
             max_abs_beta = max(abs(fit$EBeta)), error = NA_character_)
}

# Main ----
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  mode <- args[1]; cores <- if (length(args) > 1) as.integer(args[2]) else 8L
  dir.create(file.path(OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
  if (mode == "prepare") { s11_prepare(); quit(save = "no") }
  if (mode %in% c("size", "screen", "d4")) {
    z <- readRDS(LOCAL); cfg <- yaml::read_yaml("config/globals.yml")
    writeLines(c(sprintf("K = %d (EXHIBIT_K = '%s')", exhibit_fit_k(cfg), Sys.getenv("EXHIBIT_K")),
                 system("git rev-parse HEAD", intern = TRUE), format(Sys.time())),
               file.path(OUT, sprintf("settings_%s.txt", mode)))
    design <- s11_design(mode, z)
    cl <- parallel::makePSOCKcluster(cores, outfile = file.path(OUT, sprintf("workers_%s.log", mode)))
    parallel::clusterExport(cl, c("OUT", "LOCAL", "ARMS_WS", "NFOLD", "s11_load", "s11_fit",
                                  "s11_genes", "s11_task", "exhibit_k", "exhibit_fit_k"))
    parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1); s11_load(); z <- readRDS(LOCAL)
      cfg <- yaml::read_yaml("config/globals.yml"); NULL })
    st <- parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(task)
      tryCatch(s11_task(task, z, cfg), error = function(e) paste("TASK ERROR", conditionMessage(e))),
      chunk.size = 1L)
    parallel::stopCluster(cl)
    writeLines(unlist(st), file.path(OUT, sprintf("task_status_%s.txt", mode)))
  }
  if (mode == "null") {
    source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
    design <- expand.grid(seed = 1:10, p = c(2000L, 50L), variance = c(.0025, .01, .05))
    cl <- parallel::makePSOCKcluster(cores)
    parallel::clusterExport(cl, c("s9_load", "s9_zscore", "s9_fit", "s11_null_task"))
    parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1); s9_load(); NULL })
    res <- do.call(rbind, parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))),
                                                s11_null_task, chunk.size = 1L))
    parallel::stopCluster(cl)
    write.csv(res, file.path(OUT, "null_fits.csv"), row.names = FALSE)
    cells <- aggregate(cbind(n_active, converged) ~ p + variance, res, mean)
    cells$passes <- cells$n_active <= 1
    write.csv(cells, file.path(OUT, "null_cells.csv"), row.names = FALSE); print(cells)
  }
}
