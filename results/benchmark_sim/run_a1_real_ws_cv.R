# ============================================================
# Script: run_a1_real_ws_cv.R
# Purpose: Survival weight w_S in the F update of fit_cox_on_yf() on the real
#          training data (TCGA_PAAD + CPTAC), 5-fold CV only: survival share of
#          the F precision, convergence, held-out C, reconstruction loss.
# Author: Andrew Walther
# Created: 2026-10-30
# Dependencies: survival, yaml, parallel; code/fit_cox_on_yf.R,
#               code/preprocess_desurv.R, code/predict_cox_on_yf.R,
#               code/concordance_ci.R, results/benchmark_sim/benchmark_helpers.R
# Run: Rscript results/benchmark_sim/run_a1_real_ws_cv.R <arm> [<arm> ...]
#      each <arm> is w_S:init_supervised, e.g. 1:0 100:1 10000:0
# ============================================================
#
# Phase A1 of docs/plans/2026-10-26-parallel-work-sessions-plan.md. External
# cohorts are NOT loaded. Data, preprocessing, folds and the CV partial
# log-likelihood are those of run_supervised_train_cv.R (D4 gene set, D4
# double-log preprocessing done once on all training patients, so the CV is
# conditional on screening; folds stratified by cohort x event, set.seed(2026)).
#
# Configuration per arm: canonicalize = TRUE (exact sequential F update with
# the within-program coupling term), K = 7, point-Laplace L/F, per-program
# Normal beta prior, max_iter = 300, the given w_S and supervised start. This
# is the configuration of the simulation arms (run_supervision_diagnostics.R
# arm set a1_ws), not the shared point-Laplace beta prior of the 10/15 check.
#
# Parallelization: socket cluster, one task per (fold, arm), single-threaded
# fits. Forked workers (mclapply) segfaulted on macOS after the parent had
# loaded the compiled sweep and BLAS, so each worker sources the code itself.

suppressPackageStartupMessages({ library(survival); library(yaml); library(parallel) })
source("code/predict_cox_on_yf.R"); source("code/concordance_ci.R")
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL))
source("code/preprocess_desurv.R")
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
OUT_DIR <- "results/benchmark_sim/outputs/a1_real_ws_cv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
K <- 7L; NFOLD <- 5L

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Give at least one arm as w_S:init_supervised, e.g. 1:0")
ARMS <- lapply(strsplit(args, ":"), function(a) list(w_S = as.numeric(a[1]),
                                                     init_supervised = a[2] == "1"))
names(ARMS) <- sprintf("wS_%g_init%s", sapply(ARMS, `[[`, "w_S"),
                       sapply(strsplit(args, ":"), `[`, 2))

# 1. Training data (D4 preprocessing) ----
TRAIN <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN, TRAIN), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
time_train   <- unlist(lapply(TRAIN, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN, function(ds) train_raw[[ds]]$status))
pp <- preprocess_merged_cohorts(
  cohort_raw_list = train_raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank")
# The D4 gene list is git-ignored; a worktree reads it from the main checkout
d4_path <- Sys.getenv("D4_GENES_RDS",
                      "results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")
stopifnot(identical(pp$gene_names, readRDS(d4_path)))   # same genes as the reported fit
Y_train <- pp$Y
n <- nrow(Y_train)
cat(sprintf("Training: n = %d, p = %d, events = %d\n", n, ncol(Y_train), sum(status_train)))

# 2. Folds (identical to run_supervised_train_cv.R) ----
n_tcga <- train_raw$TCGA_PAAD$n
cohort <- ifelse(seq_len(n) <= n_tcga, "TCGA", "CPTAC")
set.seed(2026)
fold <- integer(n)
for (co in unique(cohort)) for (st in unique(status_train)) {
  idx <- which(cohort == co & status_train == st)
  if (length(idx)) fold[idx] <- sample(rep(seq_len(NFOLD), length.out = length(idx)))
}

common <- list(K = K, max_iter = 300L, prior_LF = "point_laplace", prior_beta = "normal",
               canonicalize = TRUE, verbose = FALSE, trace_supervision = TRUE)

#' Relative reconstruction loss of Y by a least-squares projection on EF
#'
#' ||Y - Y EF (EF'EF)^+ EF'||^2 / ||Y||^2, using only the alive programs, so
#' the same formula scores training and held-out patients.
#' @param Y n x p matrix. @param EF p x K posterior loading means.
#' @return Scalar in [0, 1].
recon_proj <- function(Y, EF) {
  EF <- EF[, colSums(EF^2) > 1e-8, drop = FALSE]
  if (!ncol(EF)) return(1)
  Lhat <- Y %*% EF %*% MASS::ginv(crossprod(EF))
  sum((Y - Lhat %*% t(EF))^2) / sum(Y^2)
}

#' Fit one (fold, arm) cell and score it on the held-out fold
run_cell <- function(f, an) {
  te <- which(fold == f); tr <- which(fold != f)
  set.seed(42L); t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(do.call(fit_cox_on_yf, c(list(Y_train[tr, , drop = FALSE],
                    time_train[tr], status_train[tr]), common, ARMS[[an]])),
                  error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error"))
    return(list(row = data.frame(fold = f, arm = an, error = conditionMessage(fit)),
                risk = NULL, te = te))
  risk <- predict_cox_on_yf(Y_train[te, , drop = FALSE], fit$EF, fit$EBeta,
                            EF_norms = fit$EF_norms)$risk_scores
  eta <- as.numeric(Y_train %*% fit$EF %*% fit$EBeta)   # EF_norms = 1 when canonical
  pl_all <- calc_cox_taylor_yf(eta, time_train, status_train)$logPL
  pl_tr  <- calc_cox_taylor_yf(eta[tr], time_train[tr], status_train[tr])$logPL
  sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  active <- abs(fit$EBeta) >= 1.96 * sdb & sdb > 0
  alive  <- colSums(fit$EF^2) > 1e-8
  shareA <- tail(fit$history$trace$shareA_median, 1)[1, ]   # w_S A_surv / A_F, final sweep
  elbo <- fit$history$elbo_full
  Ytr <- Y_train[tr, , drop = FALSE]
  row <- data.frame(fold = f, arm = an, error = NA_character_, n_test = length(te),
    C = frozen_reverse_cindex(risk, time_train[te], status_train[te]),
    cvpl = pl_all - pl_tr, converged = fit$history$converged, n_iter = fit$history$n_iter,
    elbo_decreases = sum(diff(elbo) < -1e-8 * abs(head(elbo, -1))),
    n_alive = sum(alive), n_active = sum(active),
    shareA_active_median = if (any(active & alive)) median(shareA[active & alive]) else NA_real_,
    shareA_alive_max = max(c(shareA[alive], 0), na.rm = TRUE),
    shareA_alive_median = median(shareA[alive], na.rm = TRUE),
    recon_train_fit = sum((Ytr - fit$EL %*% t(fit$EF))^2) / sum(Ytr^2),
    recon_train_proj = recon_proj(Ytr, fit$EF),
    recon_test_proj  = recon_proj(Y_train[te, , drop = FALSE], fit$EF),
    max_abs_beta = max(abs(fit$EBeta)), fit_sec = secs)
  list(row = row, risk = risk, te = te,
       state = list(EF = fit$EF, EBeta = fit$EBeta, EBeta2 = fit$EBeta2,
                    trace = fit$history$trace, elbo = elbo))
}

# 3. Cross-validation ----
grid <- expand.grid(fold = seq_len(NFOLD), arm = names(ARMS), stringsAsFactors = FALSE)
load_code <- function() {
  suppressPackageStartupMessages({ library(survival); library(ebnm) })
  source("code/predict_cox_on_yf.R"); source("code/concordance_ci.R")
  suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL))
}
cl <- makeCluster(min(nrow(grid), 8L))
clusterExport(cl, c("load_code", "run_cell", "recon_proj", "Y_train", "time_train",
                    "status_train", "fold", "common", "ARMS", "grid"))
invisible(clusterEvalQ(cl, load_code()))
res <- parLapplyLB(cl, seq_len(nrow(grid)), function(i) run_cell(grid$fold[i], grid$arm[i]))
stopCluster(cl)
for (i in seq_along(res))
  if (!is.null(res[[i]]$state))
    saveRDS(res[[i]]$state, file.path(OUT_DIR, sprintf("state_fold%d_%s.rds",
                                                       grid$fold[i], grid$arm[i])))
all_cols <- unique(unlist(lapply(res, function(r) names(r$row))))
folds_df <- do.call(rbind, lapply(res, function(r) {
  x <- r$row; x[setdiff(all_cols, names(x))] <- NA; x[all_cols] }))
if (any(!is.na(folds_df$error))) warning(sum(!is.na(folds_df$error)), " fits failed")
fname <- function(stem) file.path(OUT_DIR, sprintf("%s_%s.csv", stem, paste(names(ARMS), collapse = "-")))
write.csv(folds_df, fname("folds"), row.names = FALSE)

# 4. Summary: mean fold C, pooled out-of-fold C, CV partial log-likelihood ----
summ <- do.call(rbind, lapply(names(ARMS), function(an) {
  d <- folds_df[folds_df$arm == an, ]
  oof <- rep(NA_real_, n)
  for (i in which(grid$arm == an)) if (!is.null(res[[i]]$risk)) oof[res[[i]]$te] <- res[[i]]$risk
  data.frame(arm = an, errors = sum(!is.na(d$error)),
             mean_fold_C = mean(d$C, na.rm = TRUE), se_fold_C = sd(d$C, na.rm = TRUE) / sqrt(sum(!is.na(d$C))),
             pooled_oof_C = if (anyNA(oof)) NA_real_ else frozen_reverse_cindex(oof, time_train, status_train),
             cv_partial_loglik = sum(d$cvpl), converged = sum(d$converged, na.rm = TRUE),
             elbo_decreases = sum(d$elbo_decreases), n_active_mean = mean(d$n_active),
             shareA_active_median = median(d$shareA_active_median, na.rm = TRUE),
             shareA_alive_max_median = median(d$shareA_alive_max),
             recon_train_fit = mean(d$recon_train_fit), recon_test_proj = mean(d$recon_test_proj),
             fit_sec_median = median(d$fit_sec))
}))
write.csv(summ, fname("summary"), row.names = FALSE)
print(folds_df, digits = 4, row.names = FALSE); print(summ, digits = 4, row.names = FALSE)
