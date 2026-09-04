# ============================================================
# Script:  results/multi_cohort_sim/run_cohort_beta_recovery_sim.R
# Purpose: Stage 5 addition (DECISIONS.md 2026-09-04): ground-truth recovery
#          test for beta_cohort_id, motivated directly by the question "is
#          performance the right target, or should we be checking recovery
#          of cohort-specific factors?" External C-index cannot distinguish
#          "the model attributes signal to the right cohort" from "it
#          doesn't"; simulation with known ground truth can.
#
#          generate_multicohort_data() already supports specific_prognostic
#          = TRUE (built, never exercised by run_multicohort_sim.R): the
#          FIRST study-specific factor of each cohort gets a real, non-zero
#          survival coefficient (beta=0.8) that -- because that factor's
#          loadings are block-zeroed outside its owning cohort by
#          construction -- is a genuine, cohort-specific survival effect:
#          prognostic in one cohort, structurally silent in the other. This
#          is exactly the case beta_cohort_id was built to handle and a
#          shared beta cannot represent.
#
#          Two arms, same hybrid scenario (K_shared=2, K_specific=[2,2]),
#          same data per seed:
#            YFB_base        -- one beta_k per factor, shared across cohorts
#            YFB_beta_cohort -- beta_k^(c), cohort-specific (this session's
#                               Stage 2)
#
#          Metric (the point of this script): for the true cohort-specific
#          prognostic factor, does YFB_beta_cohort's per-cohort coefficient
#          correctly concentrate in the OWNING cohort (large |beta| there,
#          near-zero in the other), and does that give a materially
#          different (and more correct) picture than YFB_base's single
#          shared coefficient, which must average across a factor that is
#          prognostic in only half the data? classify_specificity() factor
#          recovery (genomics side) and external C-index are also reported,
#          but are secondary to the attribution question here.
#
#   Output: results/multi_cohort_sim/outputs/
#             cohort_beta_recovery_sim_results.csv  (one row per seed x arm)
#             cohort_beta_recovery_sim_attribution.csv (per-seed cohort attribution detail)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/multi_cohort_sim/run_cohort_beta_recovery_sim.R
#          Rscript results/multi_cohort_sim/run_cohort_beta_recovery_sim.R --quick
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({ library(yaml); library(survival) })
cfg <- yaml::read_yaml("config/globals.yml")

source("code/update_beta.R"); source("code/update_beta_cohort.R")
source("code/update_L.R"); source("code/update_F.R"); source("code/update_tau.R")
source("code/compute_elbo.R"); source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")
source("results/benchmark_sim/benchmark_helpers.R")
source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/build_ebmf_templates.R")
source("results/multi_cohort_sim/sim_scoring.R")

mc          <- cfg$synthetic_multicohort
C           <- mc$C
N_PER       <- unlist(mc$n_per)
P           <- mc$p
A_SHARED    <- mc$a_shared
A_SPECIFIC  <- mc$a_specific
TARGET_CENS <- mc$target_censoring
EBMF_KMAX   <- mc$ebmf_kmax
BETA_THRESH <- cfg$k_selection$beta_threshold
ALPHA       <- cfg$benchmark$alpha
PRIOR_BETA  <- "normal"
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
SEEDS       <- if (QUICK_MODE) 1:2 else 1:10
K_INIT_VALUES <- if (QUICK_MODE) c(6L) else c(6L, 12L)

# hybrid partition: 2 shared, 2 specific per cohort (matches
# synthetic_multicohort$scenarios$hybrid in globals.yml)
K_SHARED   <- 2L
K_SPECIFIC <- rep(2L, C)

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" Cohort-specific beta recovery simulation (hybrid, specific_prognostic=TRUE)\n")
cat(sprintf(" mode=%s | seeds=%s | K_init=%s\n",
            if (QUICK_MODE) "QUICK" else "FULL",
            paste(SEEDS, collapse = ","), paste(K_INIT_VALUES, collapse = ",")))
cat("============================================================\n\n")

templates <- tryCatch(build_ebmf_templates(Kmax = EBMF_KMAX),
                      error = function(e) { message("EBMF templates: ", e$message); NULL })
F_TEMPLATES <- if (!is.null(templates)) {
  Ft <- templates$F
  Ft <- if (nrow(Ft) >= P) Ft[seq_len(P), , drop = FALSE]
        else               Ft[((seq_len(P) - 1L) %% nrow(Ft)) + 1L, , drop = FALSE]
  dimnames(Ft) <- NULL
  abs(Ft)
} else NULL

results_rows <- list()
attribution_rows <- list()

for (s in SEEDS) {
  d <- generate_multicohort_data(
    C = C, n_per = N_PER, p = P, K_shared = K_SHARED, K_specific = K_SPECIFIC,
    F_templates = F_TEMPLATES, a_shared = A_SHARED, a_specific = A_SPECIFIC,
    specific_prognostic = TRUE,   # <-- the point of this script
    target_censoring = TARGET_CENS, seed = s
  )
  spl    <- stratified_split(d$status, test_frac = 0.25, seed = s)
  Ytr    <- d$Y[spl$train_idx, , drop = FALSE]
  ttr    <- d$time[spl$train_idx]
  str_   <- d$status[spl$train_idx]
  cid_tr <- droplevels(d$cohort_id[spl$train_idx])

  # The true cohort-specific prognostic factors: first specific factor per cohort.
  true_specific_prog_idx <- vapply(seq_len(C), function(cc) which(d$factor_owner == cc)[1], integer(1))

  for (K_here in K_INIT_VALUES) {
    cat(sprintf("--- seed %d, K_init=%d ---\n", s, K_here))

    fit_base <- tryCatch(
      suppressMessages(fit_cox_on_yf(Ytr, ttr, str_, K = K_here, max_iter = MAX_ITER,
                                       alpha = ALPHA, prior_beta = PRIOR_BETA, verbose = FALSE)),
      error = function(e) { message("YFB_base failed: ", e$message); NULL })
    fit_bc <- tryCatch(
      suppressMessages(fit_cox_on_yf(Ytr, ttr, str_, K = K_here, max_iter = MAX_ITER,
                                       alpha = ALPHA, prior_beta = PRIOR_BETA, verbose = FALSE,
                                       beta_cohort_id = cid_tr)),
      error = function(e) { message("YFB_beta_cohort failed: ", e$message); NULL })
    if (is.null(fit_base) || is.null(fit_bc)) next

    mf_base <- match_factors(fit_base$EF, d$F_true)
    mf_bc   <- match_factors(fit_bc$EF, d$F_true)
    est_base <- classify_specificity(fit_base$EL, cid_tr)
    est_bc   <- classify_specificity(fit_bc$EL, cid_tr)
    sa_base  <- specificity_accuracy(est_base, mf_base$match, d$factor_labels)
    sa_bc    <- specificity_accuracy(est_bc, mf_bc$match, d$factor_labels)

    # External (held-out) C-index for both arms.
    Yte <- d$Y[spl$test_idx, , drop = FALSE]; tte <- d$time[spl$test_idx]; ste <- d$status[spl$test_idx]
    pred_base <- predict_cox_on_yf(Yte, fit_base$EF, fit_base$EBeta, EF_norms = fit_base$EF_norms)
    c_base <- oriented_cindex(pred_base$risk_scores, tte, ste)
    cid_te <- droplevels(d$cohort_id[spl$test_idx])
    pred_bc <- predict_cox_on_yf(Yte, fit_bc$EF, fit_bc$EBeta, EF_norms = fit_bc$EF_norms,
                                  cohort_id_test = cid_te)
    c_bc <- oriented_cindex(pred_bc$risk_scores, tte, ste)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      seed = s, K_init = K_here,
      spec_acc_base = sa_base$accuracy, spec_acc_beta_cohort = sa_bc$accuracy,
      c_index_base = c_base, c_index_beta_cohort = c_bc,
      stringsAsFactors = FALSE
    )

    # --- The attribution question: does beta_cohort_id correctly localize
    # each true cohort-specific-prognostic factor's coefficient to its
    # OWNING cohort, vs. YFB_base's one shared (necessarily diluted) value?
    for (cc in seq_len(C)) {
      true_k <- true_specific_prog_idx[cc]
      if (is.na(true_k)) next
      matched_k_base <- mf_base$match[true_k]
      matched_k_bc   <- mf_bc$match[true_k]
      beta_base_val <- if (!is.na(matched_k_base)) abs(fit_base$EBeta[matched_k_base]) else NA_real_
      if (!is.na(matched_k_bc)) {
        owning_col <- as.character(cc)
        other_cols <- setdiff(colnames(fit_bc$EBeta), owning_col)
        beta_owning <- abs(fit_bc$EBeta[matched_k_bc, owning_col])
        beta_other  <- mean(abs(fit_bc$EBeta[matched_k_bc, other_cols]))
      } else {
        beta_owning <- NA_real_; beta_other <- NA_real_
      }
      attribution_rows[[length(attribution_rows) + 1]] <- data.frame(
        seed = s, K_init = K_here, true_cohort = cc, true_beta = 0.8,
        yfb_base_shared_beta = beta_base_val,
        yfb_beta_cohort_owning_col = beta_owning,
        yfb_beta_cohort_other_col_mean = beta_other,
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, results_rows)
attribution <- do.call(rbind, attribution_rows)
write.csv(results, file.path(OUT_DIR, "cohort_beta_recovery_sim_results.csv"), row.names = FALSE)
write.csv(attribution, file.path(OUT_DIR, "cohort_beta_recovery_sim_attribution.csv"), row.names = FALSE)

cat("\n=== Mean over seeds ===\n")
agg <- aggregate(cbind(spec_acc_base, spec_acc_beta_cohort, c_index_base, c_index_beta_cohort) ~ K_init,
                  data = results, FUN = mean, na.action = na.pass)
print(agg)

cat("\n=== Attribution: true cohort-specific-prognostic factor (true beta=0.8) ===\n")
attr_agg <- aggregate(cbind(yfb_base_shared_beta, yfb_beta_cohort_owning_col, yfb_beta_cohort_other_col_mean) ~ K_init,
                       data = attribution, FUN = function(x) mean(x, na.rm = TRUE))
print(attr_agg)
cat(sprintf("\nResults: %s\n", file.path(OUT_DIR, "cohort_beta_recovery_sim_results.csv")))
cat(sprintf("Attribution detail: %s\n", file.path(OUT_DIR, "cohort_beta_recovery_sim_attribution.csv")))
