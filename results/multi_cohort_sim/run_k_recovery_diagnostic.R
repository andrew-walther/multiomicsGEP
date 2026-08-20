# ============================================================
# Script:  results/multi_cohort_sim/run_k_recovery_diagnostic.R
# Purpose: Follow-up diagnostic to Analysis B (run_k_recovery_sim.R), which
#          found ARD substantially over-counts survival-active factors
#          (K_eff_survival) relative to the KNOWN true count, and Analysis C,
#          which found cohort_id does NOT reliably fix the false-positive
#          rate on non-prognostic factors.
#
#          This asks: are the "extra" survival-active factors (a) fragments
#          of the one true shared signal, split across several estimated
#          columns by CAVI's non-identifiability, (b) driven by cohort
#          membership (a confound cohort_id should absorb but apparently
#          doesn't fully), or (c) neither -- i.e. genuinely spurious noise?
#
#          Uses Condition A (K_shared=1, K_specific=[2,2], K_true_total=5,
#          the simplest case -- exactly one true prognostic factor) at
#          K_init=10, the same 5 seeds as Analysis B, for direct comparability.
#          For every factor the fitted model classifies "survival_active",
#          computes: correlation of its estimated L column with (i) the one
#          true shared L column, (ii) each true specific L column, (iii) a
#          cohort-membership dummy (0/1 for cohort 2).
#
#   Output: results/multi_cohort_sim/outputs/k_recovery_diagnostic_factors.csv (per-factor)
#           results/multi_cohort_sim/outputs/k_recovery_diagnostic_summary.csv (per-seed base vs. cohort_id counts)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-19
# Usage:   Rscript results/multi_cohort_sim/run_k_recovery_diagnostic.R
# ============================================================

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({ library(yaml); library(survival) })
cfg <- yaml::read_yaml("config/globals.yml")

source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/select_K.R")   # classify_factors()

source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/build_ebmf_templates.R")
source("results/multi_cohort_sim/sim_scoring.R")

mc <- cfg$synthetic_multicohort
C           <- mc$C
N_PER       <- unlist(mc$n_per)
P           <- mc$p
A_SHARED    <- mc$a_shared
A_SPECIFIC  <- mc$a_specific
TARGET_CENS <- mc$target_censoring
BETA_THRESH <- cfg$k_selection$beta_threshold
PVE_THRESH  <- cfg$k_selection$pve_threshold
ALPHA       <- cfg$benchmark$alpha
PRIOR_BETA  <- "normal"
MAX_ITER    <- cfg$cavi$max_iter
SEEDS       <- unlist(mc$seeds)
K_INIT      <- 10L   # matches the K_init=10 row of Analysis B's Condition A

K_shared   <- 1L
K_specific <- c(2L, 2L)

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

templates  <- tryCatch(build_ebmf_templates(Kmax = mc$ebmf_kmax),
                       error = function(e) NULL)
F_TEMPLATES <- if (!is.null(templates)) {
  Ft <- templates$F
  Ft <- if (nrow(Ft) >= P) Ft[seq_len(P), , drop = FALSE]
        else               Ft[((seq_len(P) - 1L) %% nrow(Ft)) + 1L, , drop = FALSE]
  dimnames(Ft) <- NULL
  abs(Ft)
} else NULL

factor_rows  <- list()
summary_rows <- list()

for (s in SEEDS) {
  d <- generate_multicohort_data(
    C = C, n_per = N_PER, p = P,
    K_shared = K_shared, K_specific = K_specific,
    F_templates = F_TEMPLATES,
    a_shared = A_SHARED, a_specific = A_SPECIFIC,
    specific_prognostic = FALSE,
    target_censoring = TARGET_CENS, seed = s
  )
  cohort_dummy <- as.numeric(d$cohort_id) - 1   # 0/1 for cohort 2 membership

  # NO cohort_id fit -- reproduces Analysis B's Condition A, K_init=10, this seed.
  fit_base <- suppressMessages(fit_cox_on_yf(
    d$Y, d$time, d$status, K = K_INIT, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA, verbose = FALSE
  ))
  cls_base <- classify_factors(fit_base, d$Y, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
  active_idx <- which(cls_base$category == "survival_active")

  cat(sprintf("--- seed=%d: %d survival_active factors (true=1) ---\n", s, length(active_idx)))

  for (k in active_idx) {
    Lk <- fit_base$EL[, k]
    cor_shared <- suppressWarnings(cor(Lk, d$L_true[, 1]))
    # true specific L columns: 2:5 (2 for cohort 1, 2 for cohort 2)
    cor_specific <- sapply(2:ncol(d$L_true), function(j) suppressWarnings(cor(Lk, d$L_true[, j])))
    best_specific_cor <- max(abs(cor_specific), na.rm = TRUE)
    cor_cohort <- suppressWarnings(cor(Lk, cohort_dummy))

    cat(sprintf("    factor %2d: |beta|=%.4f  cor_vs_true_shared=%.3f  best|cor_vs_specific|=%.3f  cor_vs_cohort=%.3f\n",
                k, abs(fit_base$EBeta[k]), cor_shared, best_specific_cor, cor_cohort))

    factor_rows[[length(factor_rows) + 1]] <- data.frame(
      seed              = s,
      factor            = k,
      abs_beta          = round(abs(fit_base$EBeta[k]), 4),
      pve               = round(cls_base$PVE[k], 4),
      cor_vs_true_shared    = round(cor_shared, 4),
      best_cor_vs_specific  = round(best_specific_cor, 4),
      cor_vs_cohort_dummy   = round(cor_cohort, 4),
      stringsAsFactors  = FALSE
    )
  }

  # WITH cohort_id fit -- does explicit cohort correction change the count?
  fit_cohort <- suppressMessages(fit_cox_on_yf(
    d$Y, d$time, d$status, K = K_INIT, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA, verbose = FALSE, cohort_id = d$cohort_id
  ))
  cls_cohort <- classify_factors(fit_cohort, d$Y, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
  n_active_cohort <- sum(cls_cohort$category == "survival_active")
  cat(sprintf("    [with cohort_id] survival_active count = %d (vs %d without)\n\n",
              n_active_cohort, length(active_idx)))

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    seed             = s,
    n_active_base    = length(active_idx),
    n_active_cohort  = n_active_cohort,
    stringsAsFactors = FALSE
  )
}

factor_results  <- do.call(rbind, factor_rows)
summary_results <- do.call(rbind, summary_rows)

out_csv_factors <- file.path(OUT_DIR, "k_recovery_diagnostic_factors.csv")
out_csv_summary <- file.path(OUT_DIR, "k_recovery_diagnostic_summary.csv")
write.csv(factor_results,  out_csv_factors, row.names = FALSE)
write.csv(summary_results, out_csv_summary, row.names = FALSE)
cat(sprintf("\nResults saved: %s, %s\n", out_csv_factors, out_csv_summary))
