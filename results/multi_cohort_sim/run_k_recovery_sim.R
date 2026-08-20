# ============================================================
# Script:  results/multi_cohort_sim/run_k_recovery_sim.R
# Purpose: Analysis B (docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md).
#          Validates ARD's ability to recover a KNOWN number of survival-active
#          (shared) and genomics-only (study-specific, non-prognostic) factors
#          when CAVI is started at K_init >> K_true.
#
#          Three conditions (via generate_multicohort_data(), C=2 cohorts,
#          specific_prognostic=FALSE so specific factors are genomics-only by
#          construction):
#            A: K_shared=1, K_specific=[2,2]  (K_true_total=5)
#            B: K_shared=2, K_specific=[2,2]  (K_true_total=6)
#            C: K_shared=2, K_specific=[3,3]  (K_true_total=8)
#
#          For each condition, fits YFB at K_init = K_true_total + {5,10,15},
#          5 seeds, and classify_factors() each fit. Metrics per fit:
#            K_eff_survival  vs. K_true_shared (ground truth)
#            K_eff_genomics  vs. K_true_specific (ground truth, summed over cohorts)
#            factor recovery: max-cor of estimated F vs. true F, shared vs.
#              specific programs scored separately (match_factors/sim_scoring.R)
#            beta RMSE for the shared (survival-active) factors, |beta_hat| vs
#              |beta_true|, aligned via the same factor matching
#
#   Output: results/multi_cohort_sim/outputs/k_recovery_sim_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-19
# Usage:   Rscript results/multi_cohort_sim/run_k_recovery_sim.R [--quick]
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

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
source("code/preprocess_desurv.R")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/select_K.R")   # classify_factors()

source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/build_ebmf_templates.R")
source("results/multi_cohort_sim/sim_scoring.R")   # match_factors(), beta_recovery()

`%+%` <- function(a, b) paste0(a, b)

# --------------------------------------------------------------------------
# 1. Configuration
# --------------------------------------------------------------------------
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
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
SEEDS       <- if (QUICK_MODE) unlist(mc$seeds)[1:2] else unlist(mc$seeds)
K_INIT_OFFSETS <- if (QUICK_MODE) c(5L) else c(5L, 10L, 15L)

# Three conditions per the plan (K_shared, K_specific per cohort).
CONDITIONS <- list(
  A = list(K_shared = 1L, K_specific = c(2L, 2L)),
  B = list(K_shared = 2L, K_specific = c(2L, 2L)),
  C = list(K_shared = 2L, K_specific = c(3L, 3L))
)

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" Analysis B — ARD K-Recovery Simulation\n")
cat(sprintf(" conditions: %s | K_init offsets: %s | seeds: %s\n",
            paste(names(CONDITIONS), collapse = ", "),
            paste(K_INIT_OFFSETS, collapse = ", "),
            paste(SEEDS, collapse = ", ")))
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# 2. EBMF templates (real gene-program shapes; NULL -> synthetic fallback)
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
# 3. Main loop: condition x K_init offset x seed
# --------------------------------------------------------------------------
rows <- list()

for (cond_name in names(CONDITIONS)) {
  cond       <- CONDITIONS[[cond_name]]
  K_shared   <- cond$K_shared
  K_specific <- cond$K_specific
  K_true_total    <- K_shared + sum(K_specific)
  K_true_specific <- sum(K_specific)

  cat(sprintf("--- Condition %s: K_shared=%d, K_specific=%s (K_true_total=%d) ---\n",
              cond_name, K_shared, paste(K_specific, collapse = ","), K_true_total))

  for (offset in K_INIT_OFFSETS) {
    K_init <- K_true_total + offset

    for (s in SEEDS) {
      d <- generate_multicohort_data(
        C = C, n_per = N_PER, p = P,
        K_shared = K_shared, K_specific = K_specific,
        F_templates = F_TEMPLATES,
        a_shared = A_SHARED, a_specific = A_SPECIFIC,
        specific_prognostic = FALSE,
        target_censoring = TARGET_CENS, seed = s
      )

      fit <- tryCatch(
        suppressMessages(fit_cox_on_yf(
          d$Y, d$time, d$status, K = K_init, max_iter = MAX_ITER, alpha = ALPHA,
          prior_beta = PRIOR_BETA, verbose = FALSE
        )),
        error = function(e) { message(sprintf("  [%s, K_init=%d, seed=%d] fit failed: %s",
                                              cond_name, K_init, s, e$message)); NULL }
      )
      if (is.null(fit)) next

      cls <- classify_factors(fit, d$Y, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
      K_eff_survival <- sum(cls$category == "survival_active")
      K_eff_genomics <- sum(cls$category == "genomics_only")
      K_eff_dead     <- sum(cls$category == "dead")

      mf     <- match_factors(fit$EF, d$F_true)
      is_sh  <- d$factor_labels == "shared"
      rec_shared   <- mean(mf$best_cor[is_sh])
      rec_specific <- mean(mf$best_cor[!is_sh])

      # beta recovery for the shared (true survival-active) factors: |beta_hat|
      # aligned to true factor order via the same match, vs. |beta_true|. Not
      # sim_scoring.R's beta_recovery() -- that compares |beta_hat| to a
      # threshold (TP/FP rate), not to the true beta magnitude (RMSE), which
      # this needs.
      bhat <- abs(fit$EBeta[mf$match])
      bhat[is.na(bhat)] <- 0
      beta_true_abs   <- abs(d$beta_true)
      beta_rmse_shared <- if (any(is_sh))
        sqrt(mean((bhat[is_sh] - beta_true_abs[is_sh])^2)) else NA_real_

      cat(sprintf("  K_init=%2d, seed=%d: K_eff_survival=%d (true=%d), K_eff_genomics=%d (true=%d), "
                  %+% "rec_shared=%.3f, rec_specific=%.3f, beta_rmse_shared=%.4f\n",
                  K_init, s, K_eff_survival, K_shared, K_eff_genomics, K_true_specific,
                  rec_shared, rec_specific, beta_rmse_shared))

      rows[[length(rows) + 1]] <- data.frame(
        condition        = cond_name,
        K_shared_true    = K_shared,
        K_specific_true  = K_true_specific,
        K_true_total     = K_true_total,
        K_init           = K_init,
        seed             = s,
        K_eff_survival   = K_eff_survival,
        K_eff_genomics   = K_eff_genomics,
        K_eff_dead       = K_eff_dead,
        rec_shared       = round(rec_shared, 4),
        rec_specific     = round(rec_specific, 4),
        beta_rmse_shared = round(beta_rmse_shared, 4),
        n_iter           = fit$history$n_iter,
        stringsAsFactors = FALSE
      )
    }
  }
  cat("\n")
}

results <- do.call(rbind, rows)

# --------------------------------------------------------------------------
# 4. Save + console summary
# --------------------------------------------------------------------------
out_csv <- file.path(OUT_DIR, "k_recovery_sim_results.csv")
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nResults saved: %s\n\n", out_csv))

cat("============================================================\n")
cat(" Mean over seeds, by condition x K_init\n")
cat("============================================================\n")
agg <- aggregate(cbind(K_eff_survival, K_eff_genomics, rec_shared, rec_specific, beta_rmse_shared) ~
                    condition + K_shared_true + K_specific_true + K_init,
                  data = results, FUN = function(x) mean(x, na.rm = TRUE))
agg <- agg[order(agg$condition, agg$K_init), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  %s (true shared=%d, specific=%d) K_init=%2d: K_eff_survival=%.1f K_eff_genomics=%.1f "
              %+% "rec_shared=%.3f rec_specific=%.3f beta_rmse=%.4f\n",
              agg$condition[i], agg$K_shared_true[i], agg$K_specific_true[i], agg$K_init[i],
              agg$K_eff_survival[i], agg$K_eff_genomics[i],
              agg$rec_shared[i], agg$rec_specific[i], agg$beta_rmse_shared[i]))
}
cat("============================================================\n")
