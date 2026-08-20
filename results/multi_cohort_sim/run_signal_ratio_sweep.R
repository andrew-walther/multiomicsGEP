# ============================================================
# Script:  results/multi_cohort_sim/run_signal_ratio_sweep.R
# Purpose: Sweep the study-specific signal amplitude (a_specific) in the
#          hybrid scenario to test how the beta false-positive rate changes
#          as study-specific variance dominates shared variance.
#
#          This directly addresses Open Question 5 (Q5) from the proposal:
#          should the cohort indicator be part of the default recommended
#          configuration, and is YFB robust to specific-dominant regimes?
#
#          Only the hybrid scenario and the four supervised arms are run
#          (EBMF has no beta, so it is irrelevant for the FP question).
#          a_shared is held fixed; a_specific varies over a grid so the
#          variance ratio (a_specific/a_shared)^2 spans 1x to 16x.
#
#          Output: results/multi_cohort_sim/outputs/
#            signal_ratio_sweep_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-14
# Usage:   Rscript results/multi_cohort_sim/run_signal_ratio_sweep.R [--quick]
# ============================================================

# --------------------------------------------------------------------------
# 0. Setup (mirrors run_multicohort_sim.R)
# --------------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("code/update_beta.R")
source("code/update_L.R")
source("code/update_F.R")
source("code/update_tau.R")
source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R")
source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))

source("code/preprocess_desurv.R")
source("results/benchmark_sim/benchmark_helpers.R")

source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/build_ebmf_templates.R")
source("results/multi_cohort_sim/sim_scoring.R")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --------------------------------------------------------------------------
# 1. Configuration
# --------------------------------------------------------------------------
mc  <- cfg$synthetic_multicohort
srs <- mc$signal_ratio_sweep

C           <- mc$C
N_PER       <- unlist(mc$n_per)
P           <- mc$p
K_FIT       <- mc$k_fit   # true K for the hybrid scenario (K_shared=2, K_specific=[2,2]); kept
                          # as the reference/"K_true" value, no longer passed to model fits below
BETA_THRESH <- cfg$k_selection$beta_threshold

# K_INIT: Analysis C (docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md) --
# each arm's model fit uses an over-specified K_INIT >> K_FIT=K_true=6, letting ARD prune,
# instead of fitting at K_FIT directly. Confirms the original signal-ratio finding (YFB beats
# EBMF->Cox when survival signal is strong) holds under ARD pruning from over-specified K.
# --k-init N overrides the default of 20.
K_INIT <- 20L
if ("--k-init" %in% args) {
  k_init_val <- suppressWarnings(as.integer(args[which(args == "--k-init") + 1]))
  if (!is.na(k_init_val) && k_init_val >= 1) K_INIT <- k_init_val
}
ALPHA       <- cfg$benchmark$alpha
PRIOR_BETA  <- "normal"
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
SEEDS       <- if (QUICK_MODE) unlist(mc$seeds)[1:2] else unlist(mc$seeds)

# Hybrid scenario only (the one where both factor types coexist)
K_SHARED   <- mc$scenarios$hybrid$K_shared
K_SPECIFIC <- unlist(mc$scenarios$hybrid$K_specific)

A_SHARED        <- srs$a_shared
A_SPECIFIC_GRID <- unlist(srs$a_specific_values)

# Only supervised arms — EBMF has no beta so FP is undefined
ARMS <- c("YFB_base", "YFB_cohort", "LB_base", "LB_cohort")

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" Signal-Ratio Sweep — Hybrid Scenario\n")
cat(sprintf(" K_true=%d | K_INIT (over-specified, ARD-pruned)=%d\n", K_FIT, K_INIT))
cat(sprintf(" a_shared=%.0f | a_specific grid: %s\n",
            A_SHARED, paste(sprintf("%.0f", A_SPECIFIC_GRID), collapse=", ")))
cat(sprintf(" ratios (spec/shared): %s\n",
            paste(sprintf("%.1fx", A_SPECIFIC_GRID / A_SHARED), collapse=", ")))
cat(sprintf(" mode=%s | seeds=%s\n",
            if (QUICK_MODE) "QUICK" else "FULL", paste(SEEDS, collapse=",")))
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# 2. EBMF templates (real-data gene programs; NULL -> synthetic)
# --------------------------------------------------------------------------
templates  <- tryCatch(build_ebmf_templates(Kmax = mc$ebmf_kmax),
                       error = function(e) { message("EBMF templates: ", e$message); NULL })
F_TEMPLATES <- if (!is.null(templates)) {
  cat(sprintf("Using EBMF templates: %d programs.\n\n", templates$K_ebmf))
  Ft <- templates$F
  Ft <- if (nrow(Ft) >= P) Ft[seq_len(P), , drop = FALSE]
        else               Ft[((seq_len(P) - 1L) %% nrow(Ft)) + 1L, , drop = FALSE]
  dimnames(Ft) <- NULL
  abs(Ft)
} else {
  cat("Using synthetic sparse-F fallback.\n\n")
  NULL
}

# --------------------------------------------------------------------------
# 3. Fitting helper (same logic as run_multicohort_sim.R)
# --------------------------------------------------------------------------
fit_arm <- function(arm, d, spl) {
  Ytr    <- d$Y[spl$train_idx, , drop = FALSE]
  ttr    <- d$time[spl$train_idx]
  str_   <- d$status[spl$train_idx]
  cid_tr <- droplevels(d$cohort_id[spl$train_idx])

  if (arm %in% c("YFB_base", "YFB_cohort")) {
    cid <- if (arm == "YFB_cohort") cid_tr else NULL
    f <- suppressMessages(fit_cox_on_yf(
      Ytr, ttr, str_, K = K_INIT, max_iter = MAX_ITER, alpha = ALPHA,
      prior_beta = PRIOR_BETA, verbose = FALSE, cohort_id = cid))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = f$EF_norms,
         EF_cohort = f$EF_cohort, n_iter = f$history$n_iter)

  } else {
    cid <- if (arm == "LB_cohort") cid_tr else NULL
    f <- suppressMessages(fit_supervised_mf_modular(
      Ytr, ttr, str_, K = K_INIT, max_iter = MAX_ITER, alpha = ALPHA,
      prior_beta = PRIOR_BETA, verbose = FALSE, cohort_id = cid))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = NULL,
         EF_cohort = f$EF_cohort, n_iter = f$history$n_iter)
  }
}

predict_test <- function(arm, fit, d, spl) {
  Yte <- d$Y[spl$test_idx, , drop = FALSE]
  tte <- d$time[spl$test_idx]
  ste <- d$status[spl$test_idx]

  if (!is.null(fit$EF_cohort) && ncol(fit$EF_cohort) > 0) {
    lev    <- levels(droplevels(d$cohort_id[spl$train_idx]))
    cid_te <- factor(d$cohort_id[spl$test_idx], levels = lev)
    Lct    <- model.matrix(~ cid_te)[, -1, drop = FALSE]
    Yte    <- Yte - Lct %*% t(fit$EF_cohort)
  }

  pr <- if (grepl("^YFB", arm))
    predict_cox_on_yf(Yte, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)
  else
    predict_supervised_mf(Yte, fit$EF, fit$EBeta)

  list(risk = pr$risk_scores, time = tte, status = ste)
}

# --------------------------------------------------------------------------
# 4. Main loop: a_specific × seed × arm
# --------------------------------------------------------------------------
rows <- list()

for (a_sp in A_SPECIFIC_GRID) {
  ratio_label <- sprintf("%.1fx", a_sp / A_SHARED)
  cat(sprintf("--- a_specific=%.0f (ratio %.1fx specific variance) ---\n",
              a_sp, (a_sp / A_SHARED)^2))

  for (s in SEEDS) {
    d <- generate_multicohort_data(
      C = C, n_per = N_PER, p = P,
      K_shared = K_SHARED, K_specific = K_SPECIFIC,
      F_templates = F_TEMPLATES,
      a_shared = A_SHARED, a_specific = a_sp,
      target_censoring = mc$target_censoring, seed = s)

    spl    <- stratified_split(d$status, test_frac = 0.25, seed = s)
    cid_tr <- droplevels(d$cohort_id[spl$train_idx])

    for (arm in ARMS) {
      fit <- tryCatch(fit_arm(arm, d, spl),
                      error = function(e) { message(arm, " failed: ", e$message); NULL })
      if (is.null(fit)) next

      mf   <- match_factors(fit$EF, d$F_true)
      est  <- classify_specificity(fit$EL, cid_tr)
      sa   <- specificity_accuracy(est, mf$match, d$factor_labels)
      br   <- beta_recovery(fit$EBeta, mf$match, d$factor_labels, BETA_THRESH)
      pred <- predict_test(arm, fit, d, spl)
      ci   <- oriented_cindex(pred$risk, pred$time, pred$status)

      is_sh <- d$factor_labels == "shared"
      rows[[length(rows) + 1]] <- data.frame(
        a_specific       = a_sp,
        var_ratio        = (a_sp / A_SHARED)^2,   # variance ratio: specific / shared
        ratio_label      = ratio_label,
        arm              = arm,
        seed             = s,
        rec_shared       = mean(mf$best_cor[is_sh]),
        rec_specific     = mean(mf$best_cor[!is_sh]),
        spec_acc         = sa$accuracy,
        beta_tp_rate     = br$tp_rate,
        beta_fp_rate     = br$fp_rate,
        mean_abs_shared  = br$mean_abs_shared,
        mean_abs_specific= br$mean_abs_specific,
        c_index          = ci,
        K_init           = K_INIT,
        n_iter           = fit$n_iter %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    }
    cat(sprintf("  seed %d done.\n", s))
  }
  cat("\n")
}

results <- do.call(rbind, rows)

# --------------------------------------------------------------------------
# 5. Save — new filename (k_init-tagged) so the original K_FIT=K_true run
#    (signal_ratio_sweep_results.csv) is preserved for comparison.
# --------------------------------------------------------------------------
out_file <- sprintf("signal_ratio_sweep_results_kinit%d.csv", K_INIT)
write.csv(results, file.path(OUT_DIR, out_file), row.names = FALSE)
cat(sprintf("Results: %s\n", file.path(OUT_DIR, out_file)))

# --------------------------------------------------------------------------
# 6. Console summary
# --------------------------------------------------------------------------
cat("============================================================\n")
cat(" Mean over seeds: FP / TP / spec_acc / C (YFB arms only)\n")
cat("============================================================\n")
yfb <- results[grepl("^YFB", results$arm), ]
agg <- aggregate(cbind(beta_fp_rate, beta_tp_rate, spec_acc, c_index) ~ var_ratio + arm,
                 data = yfb, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
agg <- agg[order(agg$var_ratio, agg$arm), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  var_ratio=%.0fx  %-11s  FP=%.2f  TP=%.2f  specAcc=%.2f  C=%.3f\n",
              agg$var_ratio[i], agg$arm[i],
              agg$beta_fp_rate[i], agg$beta_tp_rate[i],
              agg$spec_acc[i], agg$c_index[i]))
}
cat("============================================================\n")
