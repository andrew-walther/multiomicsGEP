# ============================================================
# Script: run_single_log_check.R
# Purpose: Measure the effect of the historical double log2(x+1) on the D4 benchmark.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: yaml, survival, existing YFB fitter and preprocessing
# Run from repository root: Rscript results/benchmark_sim/run_single_log_check.R
# ============================================================
# load_data_internal.R (in PDAC_DATA_ROOT) already applies log2(x+1) to TCGA_PAAD,
# Dijk and PACA_AU_seq; cfg$pdac$platform_log_transform then applies it a second
# time. The "single" arm sets every platform flag to FALSE, so each cohort is
# logged exactly once. The "double" arm is the historical pipeline, refitted with
# identical code so the comparison isolates the transform.
#
# Both arms: D4 preprocessing (combined-rank top 3,000 per cohort, per-platform
# z-standardization), K = 7, alpha_F = 0, seed 42, the same 5 stratified folds.
# Risk direction for every fit comes from its own training patients only.
# External cohorts are scored once for each full-training fit; nothing is chosen
# from them (approved with the 2026-10-26 parallel-work plan).
# Local and serial: 12 fits of about a minute each; no parallel workers.

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")
source("code/train_test_split.R")
# The legacy source defines functions before its demonstration runner stops.
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) {
  if (conditionMessage(e) != "!is.null(real_Y) is not TRUE") stop(e)
})
source("code/predict_cox_on_yf.R")
source("code/select_K.R")
source("code/feature_sensitivity.R")
library(survival)

out <- "results/benchmark_sim/outputs/single_log_check"
dir.create(file.path(out, "fits"), recursive = TRUE, showWarnings = FALSE)
snap <- "results/benchmark_snapshot_2026-10-14/single_modality"
train_cohorts <- cfg$pdac$training_cohorts
ext_cohorts <- cfg$pdac$external_cohorts

#' Harrell's C with higher risk meaning shorter survival; NA for constant risk.
#' @param risk numeric risk score.
#' @param time,status survival outcome.
#' @return Scalar C, never direction-flipped.
harrell_c <- function(risk, time, status) {
  if (stats::sd(risk) == 0) return(NA_real_)
  as.numeric(concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance)
}

# Flags: historical (double) versus each cohort logged once (single) ----
flags <- list(double = unlist(cfg$pdac$platform_log_transform),
              single = setNames(rep(FALSE, length(cfg$pdac$platform_log_transform)),
                                names(cfg$pdac$platform_log_transform)))

raw <- lapply(setNames(train_cohorts, train_cohorts), load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
time <- unlist(lapply(train_cohorts, function(ds) raw[[ds]]$time))
status <- unlist(lapply(train_cohorts, function(ds) raw[[ds]]$status))
folds <- create_stratified_folds(status, cfg$cavi$n_cv_folds, seed = 42L)
frozen_genes <- readRDS(file.path(snap, "d4_gene_names.rds"))

controls <- list(K = cfg$benchmark$k_merged_yfb_desurv, max_iter = cfg$cavi$max_iter,
  tol = cfg$cavi$tol, alpha = cfg$benchmark$alpha, alpha_F = 0,
  prior_LF = "point_exponential", prior_beta = "normal", prune = FALSE,
  canonicalize = FALSE, init_supervised = FALSE, verbose = FALSE)

prep <- lapply(flags, function(fl) preprocess_merged_cohorts(raw, fl[train_cohorts],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank"))
stopifnot(identical(prep$double$gene_names, frozen_genes))

gene_rows <- data.frame(arm = names(prep),
  n_genes = sapply(prep, function(p) length(p$gene_names)),
  overlap_with_frozen = sapply(prep, function(p) length(intersect(p$gene_names, frozen_genes))),
  share_negative_entries = sapply(prep, function(p) mean(p$Y < 0)))
write.csv(gene_rows, file.path(out, "gene_sets.csv"), row.names = FALSE)
print(gene_rows)

# Fit: 5 CV folds + full training per arm ----
cv_rows <- list(); full_fits <- list()
for (arm in names(prep)) for (fold in 0:cfg$cavi$n_cv_folds) {
  Y <- prep[[arm]]$Y
  test <- if (fold == 0L) seq_len(nrow(Y)) else folds$folds[[fold]]
  tr <- if (fold == 0L) test else setdiff(seq_len(nrow(Y)), test)
  started <- proc.time()[3]
  set.seed(42L)
  ctl <- controls; ctl$sign_correction <- fold == 0L
  fit <- do.call(fit_cox_on_yf, c(list(Y = Y[tr, , drop = FALSE], time = time[tr],
                                      status = status[tr]), ctl))
  eta_tr <- predict_cox_on_yf(Y[tr, , drop = FALSE], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
  c_tr <- harrell_c(eta_tr, time[tr], status[tr])
  d <- if (is.na(c_tr)) NA_integer_ else if (c_tr >= 0.5) 1L else -1L
  eta_te <- predict_cox_on_yf(Y[test, , drop = FALSE], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
  cv_rows[[length(cv_rows) + 1]] <- data.frame(arm = arm, fold = fold,
    evaluation = if (fold == 0L) "training" else "CV",
    C = if (is.na(d)) NA_real_ else harrell_c(d * eta_te, time[test], status[test]),
    direction = d, risk_sd = sd(eta_te), converged = fit$history$converged,
    iterations = fit$history$n_iter,
    k_eff = sum(abs(fit$EBeta) > cfg$k_selection$beta_threshold),
    runtime_seconds = proc.time()[3] - started)
  print(cv_rows[[length(cv_rows)]])
  if (fold == 0L) {
    full_fits[[arm]] <- list(fit = fit, genes = prep[[arm]]$gene_names, direction = d)
    saveRDS(full_fits[[arm]], file.path(out, "fits", paste0("full_", arm, ".rds")))
  }
}
cv <- do.call(rbind, cv_rows)
write.csv(cv, file.path(out, "cv_results.csv"), row.names = FALSE)

cv_summary <- do.call(rbind, lapply(split(cv[cv$evaluation == "CV", ], ~arm), function(x)
  data.frame(arm = x$arm[1], mean_C = mean(x$C, na.rm = TRUE),
             se_C = sd(x$C, na.rm = TRUE) / sqrt(sum(!is.na(x$C))),
             n_na = sum(is.na(x$C)), n_converged = sum(x$converged))))
paired <- merge(cv[cv$arm == "single" & cv$evaluation == "CV", c("fold", "C")],
                cv[cv$arm == "double" & cv$evaluation == "CV", c("fold", "C")],
                by = "fold", suffixes = c("_single", "_double"))
cv_summary$paired_diff_single_minus_double <- mean(paired$C_single - paired$C_double, na.rm = TRUE)
write.csv(cv_summary, file.path(out, "cv_summary.csv"), row.names = FALSE)
print(cv_summary)

# External: score once; paired donor-level bootstrap against frozen risk ----
frozen_risk <- readRDS(file.path(snap, "desurv_comparison_riskscores.rds"))$D4
ext_rows <- list(); set.seed(20261026L); B <- 2000L
for (ds in ext_cohorts) {
  r <- load_pdac_raw(ds, PDAC_DATA_ROOT)
  scores <- list()
  for (arm in names(full_fits)) {
    pre <- preprocess_desurv_cohort(r$Y, r$gene_names, top_n = NULL,
      log_transform = flags[[arm]][[ds]], cohort_name = ds,
      rank_transform = FALSE, per_platform_standardize = TRUE)
    ff <- full_fits[[arm]]
    common <- intersect(ff$genes, pre$gene_names)
    if (length(common) < 100) stop(sprintf("%s: only %d shared genes", ds, length(common)))
    idx <- match(common, ff$genes)
    eta <- predict_cox_on_yf(pre$Y[, match(common, pre$gene_names), drop = FALSE],
      ff$fit$EF[idx, , drop = FALSE], ff$fit$EBeta, EF_norms = ff$fit$EF_norms)$risk_scores
    scores[[arm]] <- list(risk = ff$direction * eta, n_genes = length(common))
  }
  # Historical D4 predates the sign fix: frozen risk = -eta (DECISIONS 2026-10-02).
  fr <- frozen_risk[[ds]]
  stopifnot(identical(fr$time, r$time), identical(as.integer(fr$status), as.integer(r$status)))
  scores$frozen <- list(risk = -fr$risk, n_genes = NA_integer_)
  cs <- sapply(scores, function(s) harrell_c(s$risk, r$time, r$status))
  boot <- replicate(B, {
    i <- sample.int(length(r$time), replace = TRUE)
    harrell_c(scores$single$risk[i], r$time[i], r$status[i]) -
      harrell_c(scores$frozen$risk[i], r$time[i], r$status[i])
  })
  ext_rows[[ds]] <- data.frame(cohort = ds, n = length(r$time),
    C_frozen = cs[["frozen"]], C_double_refit = cs[["double"]], C_single = cs[["single"]],
    genes_single = scores$single$n_genes, genes_double = scores$double$n_genes,
    diff_single_minus_frozen = cs[["single"]] - cs[["frozen"]],
    diff_ci_lo = unname(quantile(boot, 0.025, na.rm = TRUE)),
    diff_ci_hi = unname(quantile(boot, 0.975, na.rm = TRUE)),
    risk_cor_single_frozen = cor(scores$single$risk, scores$frozen$risk))
  print(ext_rows[[ds]])
}
ext <- do.call(rbind, ext_rows)
ext <- rbind(ext, data.frame(cohort = "mean", n = sum(ext$n),
  C_frozen = mean(ext$C_frozen), C_double_refit = mean(ext$C_double_refit),
  C_single = mean(ext$C_single), genes_single = NA, genes_double = NA,
  diff_single_minus_frozen = mean(ext$diff_single_minus_frozen),
  diff_ci_lo = NA, diff_ci_hi = NA, risk_cor_single_frozen = NA))
write.csv(ext, file.path(out, "external_results.csv"), row.names = FALSE)
writeLines(c("Single-log check of the D4 benchmark (2026-10-26).",
  "CV is conditional on screening and per-cohort scaling done once on all 273 training patients.",
  "External scoring: once per full-training fit, direction from training patients; no selection.",
  "PACA_AU_array and PACA_AU_seq share 47 donors; the bootstrap is per cohort.",
  capture.output(str(controls)), capture.output(sessionInfo())),
  file.path(out, "design.txt"))
