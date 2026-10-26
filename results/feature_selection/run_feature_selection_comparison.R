# ============================================================
# Script: run_feature_selection_comparison.R
# Purpose: Compare unsupervised versus survival-screened gene sets at matched sizes (training CV).
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: yaml, survival, existing YFB fitter, feature-sensitivity helpers
# Run from repository root:
#   Rscript results/feature_selection/run_feature_selection_comparison.R manifest|cv|summary|agreement
# ============================================================
# Design: docs/plans/2026-10-16-feature-selection-experiment-plan.md, Tasks E2-E4, with the
# Gate E2 defaults approved on 2026-10-26: N in {250, 500, 1000} matched to U(N)'s actual
# shared counts; screen = Wald p of coxph(Surv ~ y + strata(cohort)) on fold-training
# patients; no survival-only arm.
#
# Arms: U2064 (D4 baseline genes), U(N) (smaller nested unsupervised mean/variance screen),
# S(N) (the 2,064 baseline genes ranked by the stratified Cox screen; top N).
# Input: data/round2_feature_cold_primary.rds, the frozen 273-patient training matrix of the
# feature-count study (historical double-log preprocessing, per-cohort z-scores of all
# 17,900 shared genes; z-scores are per gene, so any gene subset equals the subset's own
# standardization). Scaling precedes the folds, so CV C is conditional on it.
# External cohorts are not loaded (Task E5 waits for Andrew's review).
# Local and serial: 42 fits of a few seconds each.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args[1] %in% c("manifest", "cv", "summary", "agreement"))
mode <- args[1]
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) {
  if (conditionMessage(e) != "!is.null(real_Y) is not TRUE") stop(e)
})
source("code/predict_cox_on_yf.R")
source("code/select_K.R")
source("code/feature_sensitivity.R")
source("results/feature_selection/feature_selection_helpers.R")
library(survival)

out <- "results/feature_selection/outputs"
dir.create(file.path(out, "fits"), recursive = TRUE, showWarnings = FALSE)
z <- readRDS("data/round2_feature_cold_primary.rds")
frozen_genes <- readRDS("results/benchmark_snapshot_2026-10-14/single_modality/d4_gene_names.rds")
stopifnot(identical(z$manifest$genes[[1]], frozen_genes), nrow(z$Y) == 273L)

# Cohort labels: rows are TCGA_PAAD then CPTAC, as built by run_feature_expansion.R prepare.
n_tcga <- load_pdac_raw("TCGA_PAAD", PDAC_DATA_ROOT)$n
cohort <- rep(c("TCGA_PAAD", "CPTAC"), c(n_tcga, nrow(z$Y) - n_tcga))
stopifnot(sum(cohort == "CPTAC") > 0)

N_GRID <- c(250L, 500L, 1000L)
n_folds <- length(z$folds$folds)

# Gene sets per (arm, N, fold); fold 0 = all 273 training patients ----
gene_sets <- function() {
  sets <- list(); rows <- list()
  add <- function(arm, N, fold, genes, n_failed = NA_integer_, extra = list()) {
    key <- sprintf("%s_%s_f%d", arm, N, fold)
    sets[[key]] <<- genes
    rows[[key]] <<- data.frame(arm = arm, N_requested = N, fold = fold,
      n_genes = length(genes), overlap_U2064 = length(intersect(genes, frozen_genes)),
      n_failed = n_failed, top_n_per_cohort = if (is.null(extra$q)) NA else extra$q)
  }
  for (fold in 0:n_folds) add("U2064", 2064L, fold, frozen_genes, extra = list(q = 3000L))
  base_idx <- match(frozen_genes, z$genes)
  for (N in N_GRID) {
    u <- unsupervised_gene_set(z$manifest$entry_rank, N)
    stopifnot(all(u$genes %in% frozen_genes))
    for (fold in 0:n_folds) {
      add("U", N, fold, u$genes, extra = list(q = u$top_n_per_cohort))
      tr <- if (fold == 0L) seq_len(nrow(z$Y)) else setdiff(seq_len(nrow(z$Y)), z$folds$folds[[fold]])
      s <- survival_screen_genes(z$Y[tr, base_idx, drop = FALSE], z$time[tr], z$status[tr],
                                 cohort[tr], u$n_shared_actual)
      add("S", N, fold, frozen_genes[s], n_failed = attr(s, "n_failed"))
    }
  }
  man <- do.call(rbind, rows)
  man$overlap_S_U <- NA_integer_
  for (i in which(man$arm == "S")) {
    ukey <- sprintf("U_%s_f%d", man$N_requested[i], man$fold[i])
    man$overlap_S_U[i] <- length(intersect(sets[[rownames(man)[i]]], sets[[ukey]]))
  }
  list(sets = sets, manifest = man)
}

# CV summaries ----

#' Summarize CV C by arm and size, with paired S(N) - U(N) differences
#'
#' A fold whose held-out risk is numerically constant (SD < RISK_SD_MIN; this happens
#' when no program is survival-active and η is rounding noise) has C undefined (NA)
#' in the primary summary. The as-computed summary keeps those C values for comparison.
#' @param cv Data frame from cv mode (one row per arm, N, fold).
#' @return Invisibly, the primary summary; writes cv_summary.csv and
#'   cv_summary_as_computed.csv.
summarize_cv <- function(cv) {
  cvf <- cv[cv$evaluation == "CV", ]
  cvf$C_defined <- ifelse(cvf$risk_sd < RISK_SD_MIN, NA_real_, cvf$C)
  one <- function(col) {
    summ <- do.call(rbind, lapply(split(cvf, list(cvf$arm, cvf$N), drop = TRUE), function(x)
      data.frame(arm = x$arm[1], N = x$N[1], n_genes = x$n_genes[1],
        mean_C = mean(x[[col]], na.rm = TRUE),
        se_C = sd(x[[col]], na.rm = TRUE) / sqrt(sum(!is.na(x[[col]]))),
        n_C_na = sum(is.na(x[[col]])), n_converged = sum(x$converged, na.rm = TRUE),
        n_errors = sum(!is.na(x$error) & nzchar(x$error)), mean_k_eff = mean(x$k_eff, na.rm = TRUE))))
    # Paired S(N) - U(N) on folds where both are defined.
    summ$paired_diff_S_minus_U <- NA_real_; summ$paired_diff_se <- NA_real_
    summ$paired_n_folds <- NA_integer_
    for (N in N_GRID) {
      dd <- merge(cvf[cvf$arm == "S" & cvf$N == N, c("fold", col)],
                  cvf[cvf$arm == "U" & cvf$N == N, c("fold", col)], by = "fold")
      dif <- dd[[2]] - dd[[3]]
      i <- summ$arm == "S" & summ$N == N
      summ$paired_diff_S_minus_U[i] <- mean(dif, na.rm = TRUE)
      summ$paired_diff_se[i] <- sd(dif, na.rm = TRUE) / sqrt(sum(!is.na(dif)))
      summ$paired_n_folds[i] <- sum(!is.na(dif))
    }
    summ[order(summ$N, summ$arm), ]
  }
  primary <- one("C_defined"); raw <- one("C")
  write.csv(primary, file.path(out, "cv_summary.csv"), row.names = FALSE)
  write.csv(raw, file.path(out, "cv_summary_as_computed.csv"), row.names = FALSE)
  cat("Primary (near-constant risk -> C = NA):\n"); print(primary)
  cat("As computed:\n"); print(raw)
  invisible(primary)
}
RISK_SD_MIN <- 1e-6

gs_path <- file.path(out, "fits", "gene_sets.rds")
if (mode == "manifest" || !file.exists(gs_path)) {
  gs <- gene_sets()
  saveRDS(gs, gs_path)
  write.csv(gs$manifest, file.path(out, "manifest.csv"), row.names = FALSE)
  print(gs$manifest[gs$manifest$fold %in% c(0, 1), ])
  if (mode == "manifest") quit(status = 0)
}
gs <- readRDS(gs_path)

# Cross-validated fits (Task E3) ----
if (mode == "cv") {
  rows <- list()
  for (key in rownames(gs$manifest)) {
    m <- gs$manifest[key, ]
    genes <- gs$sets[[key]]
    Y <- z$Y[, match(genes, z$genes), drop = FALSE]
    test <- if (m$fold == 0L) seq_len(nrow(Y)) else z$folds$folds[[m$fold]]
    tr <- if (m$fold == 0L) test else setdiff(seq_len(nrow(Y)), test)
    started <- proc.time()[3]; warn <- character(); err <- ""
    row <- data.frame(arm = m$arm, N = m$N_requested, n_genes = length(genes), fold = m$fold,
      evaluation = if (m$fold == 0L) "training" else "CV", direction = NA_integer_,
      C = NA_real_, risk_sd = NA_real_, converged = NA, n_iter = NA_integer_,
      k_eff = NA_integer_, surv_active = NA_integer_, retained = NA_integer_,
      runtime_seconds = NA_real_, warnings = "", error = "")
    tryCatch(withCallingHandlers({
      set.seed(42L)
      ctl <- z$controls; ctl$sign_correction <- m$fold == 0L
      fit <- do.call(fit_cox_on_yf, c(list(Y = Y[tr, , drop = FALSE], time = z$time[tr],
                                          status = z$status[tr]), ctl))
      eta_tr <- predict_cox_on_yf(Y[tr, , drop = FALSE], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
      d <- training_direction(eta_tr, z$time[tr], z$status[tr])
      eta_te <- predict_cox_on_yf(Y[test, , drop = FALSE], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
      classes <- classify_factors(fit, Y[tr, , drop = FALSE],
        beta_thresh = cfg$k_selection$beta_threshold, pve_thresh = cfg$k_selection$pve_threshold)
      row$direction <- d
      row$C <- if (is.na(d)) NA_real_ else harrell_c(d * eta_te, z$time[test], z$status[test])
      row$risk_sd <- sd(eta_te)
      row$converged <- fit$history$converged; row$n_iter <- fit$history$n_iter
      row$k_eff <- sum(abs(fit$EBeta) > cfg$k_selection$beta_threshold)
      row$surv_active <- sum(classes$surv_active); row$retained <- sum(classes$category != "dead")
      if (m$fold == 0L) saveRDS(list(fit = fit, genes = genes, direction = d),
                               file.path(out, "fits", paste0(key, ".rds")))
    }, warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) err <<- conditionMessage(e))
    row$warnings <- paste(unique(warn), collapse = " | "); row$error <- err
    row$runtime_seconds <- proc.time()[3] - started
    rows[[key]] <- row
    print(row[, 1:13])
  }
  cv <- do.call(rbind, rows)
  write.csv(cv, file.path(out, "cv_results.csv"), row.names = FALSE)

  summarize_cv(cv)
  if (any(nzchar(cv$error))) stop("At least one fit failed; all rows were written.")
}

if (mode == "summary") summarize_cv(read.csv(file.path(out, "cv_results.csv"),
  stringsAsFactors = FALSE, na.strings = c("NA")))

# Program agreement with the U2064 full-training fit (Task E4) ----
if (mode == "agreement") {
  ref <- readRDS(file.path(out, "fits", "U2064_2064_f0.rds"))
  ref_Y <- z$Y[, match(ref$genes, z$genes), drop = FALSE]
  rows <- list()
  for (key in grep("_f0$", rownames(gs$manifest), value = TRUE)) {
    if (key == "U2064_2064_f0") next
    f <- readRDS(file.path(out, "fits", paste0(key, ".rds")))
    new_Y <- z$Y[, match(f$genes, z$genes), drop = FALSE]
    ag <- feature_program_agreement(ref$fit, ref$genes, f$fit, f$genes, ref_Y, new_Y)
    risk_ref <- ref$direction * predict_cox_on_yf(ref_Y, ref$fit$EF, ref$fit$EBeta, ref$fit$EF_norms)$risk_scores
    risk_new <- f$direction * predict_cox_on_yf(new_Y, f$fit$EF, f$fit$EBeta, f$fit$EF_norms)$risk_scores
    m <- gs$manifest[key, ]
    tab <- ag$pairs
    tab$arm <- m$arm; tab$N <- m$N_requested; tab$n_genes <- m$n_genes
    tab$risk_cor_with_U2064 <- if (sd(risk_new) == 0) NA_real_ else cor(risk_ref, risk_new)
    rows[[key]] <- tab
  }
  cols <- unique(unlist(lapply(rows, names)))
  for (k in names(rows)) for (nm in setdiff(cols, names(rows[[k]]))) rows[[k]][[nm]] <- NA
  agree <- do.call(rbind, lapply(rows, function(x) x[, cols, drop = FALSE]))
  write.csv(agree, file.path(out, "agreement.csv"), row.names = FALSE)
  print(agree)
}
