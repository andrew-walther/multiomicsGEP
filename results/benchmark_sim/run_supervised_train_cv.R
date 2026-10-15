# ============================================================
# Script: run_supervised_train_cv.R
# Purpose: 5-fold cross-validated C-index and partial log-likelihood on the
#          real training data (TCGA_PAAD + CPTAC) for two survival-supervised
#          fit_cox_on_yf() configurations: point-Laplace vs normal beta prior.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, yaml; code/fit_cox_on_yf.R, code/preprocess_desurv.R,
#               code/predict_cox_on_yf.R, code/concordance_ci.R
# Run: Rscript results/benchmark_sim/run_supervised_train_cv.R
# ============================================================
#
# Known simplification: gene selection and per-platform z-standardization were
# done once on ALL training data (rows of Y_train are not re-preprocessed per
# fold), so the CV estimates are mildly optimistic in absolute terms; both
# configurations share the same folds, so the comparison is unaffected in kind.
# External validation cohorts are NOT loaded.
#
# Folds: 5, stratified by cohort x event status, set.seed(2026), shared by both
# configurations. CV partial log-likelihood (Verweij-van Houwelingen):
#   cvpl_f = PL_all(beta_{-f}) - PL_{-f}(beta_{-f}),  summed over folds,
# with eta = Y %*% EF %*% EBeta from the fit on the training folds.

suppressPackageStartupMessages({ library(survival); library(yaml) })
source("code/predict_cox_on_yf.R"); source("code/concordance_ci.R")
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL))
source("code/preprocess_desurv.R")
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
OUT_DIR <- "results/benchmark_sim/outputs/supervised_train_cv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
K <- 7L; NFOLD <- 5L

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
d4_genes <- readRDS("results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")
stopifnot(identical(pp$gene_names, d4_genes))   # same gene set as the reported fit
Y_train <- pp$Y
n <- nrow(Y_train)
cat(sprintf("Training: n = %d, p = %d, events = %d\n", n, ncol(Y_train), sum(status_train)))

# 2. Folds ----
n_tcga <- train_raw$TCGA_PAAD$n
cohort <- ifelse(seq_len(n) <= n_tcga, "TCGA", "CPTAC")
set.seed(2026)
fold <- integer(n)
for (co in unique(cohort)) for (st in unique(status_train)) {
  idx <- which(cohort == co & status_train == st)
  if (length(idx)) fold[idx] <- sample(rep(seq_len(NFOLD), length.out = length(idx)))
}
print(table(fold, cohort)); print(table(fold, status_train))

# 3. Configurations ----
common <- list(K = K, max_iter = 300L, prior_LF = "point_laplace", canonicalize = TRUE,
               w_S = 1, init_supervised = TRUE, beta_prior_mode = "shared",
               beta_warmup = 10L, verbose = FALSE)
configs <- list(pl = c(common, list(prior_beta = "point_laplace")),
                normal = c(common, list(prior_beta = "normal")))

#' Linear predictor eta = Y (EF / EF_norms) EBeta
#' @param Y n x p matrix. @param fit fit_cox_on_yf() result.
#' @return numeric n-vector eta (EF_norms = 1 under canonicalize = TRUE).
eta_of <- function(Y, fit) {
  as.numeric(Y %*% sweep(fit$EF, 2, fit$EF_norms, "/") %*% fit$EBeta)
}

# 4. Cross-validation ----
rows <- list(); oof <- list(pl = numeric(n), normal = numeric(n))
for (f in seq_len(NFOLD)) {
  te <- which(fold == f); tr <- which(fold != f)
  for (cn in names(configs)) {
    set.seed(42L); t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(do.call(fit_cox_on_yf, c(list(Y_train[tr, , drop = FALSE],
                      time_train[tr], status_train[tr]), configs[[cn]])),
                    error = function(e) { message("FIT ERROR fold ", f, " ", cn, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) {
      rows[[length(rows) + 1]] <- data.frame(fold = f, config = cn, n_test = length(te),
        C = NA, cvpl = NA, converged = NA, n_iter = NA, n_active = NA, fit_sec = NA)
      next
    }
    secs <- proc.time()[["elapsed"]] - t0
    risk <- predict_cox_on_yf(Y_train[te, , drop = FALSE], fit$EF, fit$EBeta,
                              EF_norms = fit$EF_norms)$risk_scores
    oof[[cn]][te] <- risk
    C <- frozen_reverse_cindex(risk, time_train[te], status_train[te])
    eta <- eta_of(Y_train, fit)
    pl_all <- calc_cox_taylor_yf(eta, time_train, status_train)$logPL
    pl_tr  <- calc_cox_taylor_yf(eta[tr], time_train[tr], status_train[tr])$logPL
    sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
    n_act <- sum(abs(fit$EBeta) >= 1.96 * sdb)
    cat(sprintf("fold %d %-6s C = %.4f cvpl = %.2f conv = %s iter = %d active = %d (%.0fs)\n",
                f, cn, C, pl_all - pl_tr, fit$history$converged, fit$history$n_iter, n_act, secs))
    rows[[length(rows) + 1]] <- data.frame(fold = f, config = cn, n_test = length(te),
      C = C, cvpl = pl_all - pl_tr, converged = fit$history$converged,
      n_iter = fit$history$n_iter, n_active = n_act, fit_sec = secs)
  }
}
folds_df <- do.call(rbind, rows)
write.csv(folds_df, file.path(OUT_DIR, "folds.csv"), row.names = FALSE)

# 5. Summary ----
summ <- do.call(rbind, lapply(names(configs), function(cn) {
  d <- folds_df[folds_df$config == cn, ]
  data.frame(config = cn, mean_fold_C = mean(d$C, na.rm = TRUE), sd_fold_C = sd(d$C, na.rm = TRUE),
             pooled_oof_C = frozen_reverse_cindex(oof[[cn]], time_train, status_train),
             cv_partial_loglik = sum(d$cvpl), n_nonconverged = sum(!d$converged | is.na(d$converged)),
             paired_diff_pl_minus_normal = NA_real_, paired_diff_se = NA_real_)
}))
dd <- folds_df$C[folds_df$config == "pl"] - folds_df$C[folds_df$config == "normal"]
summ$paired_diff_pl_minus_normal[summ$config == "pl"] <- mean(dd)
summ$paired_diff_se[summ$config == "pl"] <- sd(dd) / sqrt(length(dd))
write.csv(summ, file.path(OUT_DIR, "summary.csv"), row.names = FALSE)
print(folds_df, digits = 4, row.names = FALSE); print(summ, digits = 4, row.names = FALSE)
