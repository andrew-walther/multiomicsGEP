# ============================================================
# Script:  results/multi_cohort_sim/run_multicohort_sim.R
# Purpose: Multi-cohort simulation study.  For each of three scenarios
#          (all-shared, hybrid, nothing-shared), several seeds, and three
#          K_init settings, fit two arms and score how well each recovers the
#          shared vs. study-specific factor structure and the survival
#          coefficients:
#
#            YFB_base    — YFB model (η = (YF)β), no cohort indicator
#            EBMF        — unsupervised benchmark (flashier, survival-blind)
#
#          K_init settings mirror the real-data K-selection procedure (over-
#          specified K_init + ARD pruning) rather than fitting only at the
#          oracle true K:
#            oracle_k6 — K_init = 6 = the true total K in every scenario;
#                        retained as an internal reference only (not a deck
#                        figure), to compare ARD-based recovery against it
#            ard_k12   — K_init = 12 (2x the true K)
#            ard_k20   — K_init = 20 (= ebmf_kmax)
#          LB_base, LB_cohort, and YFB_cohort arms are dropped: LB is out of
#          the narrative entirely, and 6/18 already found no meaningful YFB
#          vs. YFB+cohort difference in simulation.
#
#          Metrics: shared/specific factor recovery (|cor| of gene programs),
#          specificity-classification accuracy, β recovery (TP/FP prognostic
#          rates, using the same |β|>0.001 ARD threshold as the real-data
#          procedure), and orientation-free held-out C-index.
#
#          Output: results/multi_cohort_sim/outputs/
#            multicohort_sim_results.csv   (one row per scenario × K_init × arm × seed)
#            multicohort_sim_example.rds   (data + fits, first seed, per K_init, for figures)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-14
# Updated: 2026-08-27 -- ARD-based K_init sweep, YFB_base + EBMF arms only
# Usage:   Rscript results/multi_cohort_sim/run_multicohort_sim.R [--quick]
# ============================================================

# --------------------------------------------------------------------------
# 0. Setup
# --------------------------------------------------------------------------
args         <- commandArgs(trailingOnly = TRUE)
QUICK_MODE   <- "--quick" %in% args
n_seeds_arg  <- args[grepl("^--n-seeds=", args)]

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

# core model + utility code (mirrors run_synthetic.R source block)
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

# preprocessing + helpers (needed for the real-data EBMF templates)
source("code/preprocess_desurv.R")
source("results/benchmark_sim/benchmark_helpers.R")   # defines PDAC_DATA_ROOT etc. (uses cfg)

# simulation-specific code
source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/build_ebmf_templates.R")
source("results/multi_cohort_sim/sim_scoring.R")

HAVE_FLASHIER <- requireNamespace("flashier", quietly = TRUE)

# null-coalescing helper (base R only added %||% in 4.4; define for portability)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --------------------------------------------------------------------------
# 1. Configuration (from globals.yml — never hardcoded)
# --------------------------------------------------------------------------
mc          <- cfg$synthetic_multicohort
C           <- mc$C
N_PER       <- unlist(mc$n_per)
P           <- mc$p
K_FIT_TRUE  <- mc$k_fit
A_SHARED    <- mc$a_shared
A_SPECIFIC  <- mc$a_specific
OFFSET_SD   <- mc$offset_sd
TARGET_CENS <- mc$target_censoring
EBMF_KMAX   <- mc$ebmf_kmax
BETA_THRESH <- cfg$k_selection$beta_threshold
ALPHA       <- cfg$benchmark$alpha
PRIOR_BETA  <- "normal"          # avoids β→0 collapse (run_synthetic.R rationale)
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
# --n-seeds=N overrides config/globals.yml's synthetic_multicohort$seeds (still
# [1,2,3,4,5]) for THIS run only, rather than editing the shared config -- six
# other run_*.R scripts in this directory read the same seeds list, and
# changing it globally would silently double their runtime too.
SEEDS <- if (QUICK_MODE) { unlist(mc$seeds)[1]
} else if (length(n_seeds_arg) > 0) { seq_len(as.integer(sub("^--n-seeds=", "", n_seeds_arg[1])))
} else { unlist(mc$seeds) }

SCENARIOS <- lapply(mc$scenarios, function(s)
  list(K_shared = s$K_shared, K_specific = unlist(s$K_specific)))

# Arms: YFB_base (the recommended joint model) vs. EBMF (unsupervised
# two-step) only. LB_base, LB_cohort, and YFB_cohort are dropped -- LB is out
# of the narrative entirely, and 6/18 already found no meaningful YFB vs.
# YFB+cohort difference in simulation (cohort indicator is a future-directions
# item now, not a tested arm).
ARMS <- c("YFB_base", "EBMF")

# K: over-specified K_init + ARD pruning, mirroring the real-data procedure,
# instead of the oracle K_FIT = true total K in every scenario. K_init=12 (2x
# the true K=6) and K_init=20 (= ebmf_kmax) show whether ARD is stable across
# starting K, the same question the real-data K_init sweep asks. The oracle
# K=6 fit is RETAINED as an internal reference only (to compare ARD-based
# recovery against the 6/18 oracle-K figures for the deck-placement
# decision) -- it is not itself a deck figure.
#
# under_k2/3/4 (added 2026-09-04, DECISIONS.md): UNDER-specified K_init,
# below K_true=6, answering the advisors' question of whether multiple true
# factors merge into one estimated factor when K_init is too small to
# represent them separately. match_factors()/classify_specificity() already
# support this without change -- the "under-specified" case is just a
# smaller K passed to the same fitting/scoring pipeline.
K_SETTINGS <- c(under_k2 = 2L, under_k3 = 3L, under_k4 = 4L,
                oracle_k6 = K_FIT_TRUE, ard_k12 = 12L, ard_k20 = EBMF_KMAX)

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" Multi-Cohort Simulation — Shared vs. Study-Specific Factors\n")
cat(sprintf(" mode=%s | C=%d | n_per=%s | p=%d | K_settings=%s | seeds=%s\n",
            if (QUICK_MODE) "QUICK" else "FULL", C, paste(N_PER, collapse=","),
            P, paste(names(K_SETTINGS), K_SETTINGS, sep = "=", collapse = ", "),
            paste(SEEDS, collapse=",")))
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# 2. EBMF gene-program templates (real-data ground truth; NULL -> synthetic)
# --------------------------------------------------------------------------
templates  <- tryCatch(build_ebmf_templates(Kmax = EBMF_KMAX),
                       error = function(e) { message("EBMF templates: ", e$message); NULL })
F_TEMPLATES <- if (!is.null(templates)) {
  cat(sprintf("Using EBMF templates: %d programs (p=%d).\n\n",
              templates$K_ebmf, nrow(templates$F)))
  # templates are on the real gene set; resize rows to P by recycling if needed
  Ft <- templates$F
  Ft <- if (nrow(Ft) >= P) Ft[seq_len(P), , drop = FALSE]
        else               Ft[((seq_len(P) - 1L) %% nrow(Ft)) + 1L, , drop = FALSE]
  dimnames(Ft) <- NULL   # drop gene-name row names (recycling would duplicate them)
  # The supervised models use a non-negative prior on L and F (point_exponential,
  # NMF-style) — this is the established model assumption: see DECISIONS 2026-05-06
  # ("Phase 2 initialization constraints"), where the SVD init was set to abs()
  # precisely because it is "equally valid for the non-negative point_exponential
  # prior."  flashier's ldf$F is signed (~48% negative); using it directly would
  # impose a prior/data mismatch that handicaps the supervised models on programs
  # they structurally cannot represent.  Take abs() so the ground-truth programs
  # match the model's non-negative assumption (still realistic, sparse,
  # real-data-grounded).  The EBMF benchmark arm is likewise run non-negative.
  abs(Ft)
} else {
  cat("Using synthetic sparse-F fallback (no EBMF templates).\n\n")
  NULL
}

# --------------------------------------------------------------------------
# 3. Per-arm fitting + prediction helpers
# --------------------------------------------------------------------------

#' Fit one arm on the training split, at a given K_init.  Returns a list with
#' EF (p×K), EL (n_tr×K), EBeta (or NULL for EBMF), EF_norms, EF_cohort
#' (always NULL now -- the cohort-indicator arms are dropped), n_iter.
fit_arm <- function(arm, d, spl, K) {
  Ytr    <- d$Y[spl$train_idx, , drop = FALSE]
  ttr    <- d$time[spl$train_idx]
  str_   <- d$status[spl$train_idx]

  if (arm == "YFB_base") {
    f <- suppressMessages(fit_cox_on_yf(
      Ytr, ttr, str_, K = K, max_iter = MAX_ITER, alpha = ALPHA,
      prior_beta = PRIOR_BETA, verbose = FALSE, cohort_id = NULL))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = f$EF_norms,
         EF_cohort = f$EF_cohort, n_iter = f$history$n_iter)

  } else if (arm == "EBMF") {
    # Unsupervised benchmark, extended to a TWO-STAGE baseline:
    #   stage 1 — survival-blind non-negative EBMF factorization of training Y;
    #   stage 2 — a Cox model fit on the EBMF factor scores.
    # This is the "factorize first, attach survival afterward" pipeline that the
    # joint SSBMF model is meant to outperform: EBMF chooses its factors WITHOUT
    # seeing the outcome, so it cannot preferentially sharpen the prognostic ones.
    # Non-negative priors (point_exponential on L and F) match the supervised arm,
    # so the only difference being tested is joint vs. two-stage supervision.
    ebnm_args <- if (requireNamespace("ebnm", quietly = TRUE))
      list(ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential))
    else list()
    f <- do.call(flashier::flash,
                 c(list(Ytr, var_type = 2, greedy_Kmax = K,
                        backfit = TRUE, verbose = 0), ebnm_args))
    l  <- flashier::ldf(f, type = "2")
    EF <- as.matrix(l$F)                 # p × K unit-norm gene programs
    EL <- as.matrix(l$L)                 # n_tr × K loadings (for structural scoring)
    # Stage 2: project training Y onto the EBMF programs and fit a downstream Cox.
    # Using the least-squares projection (rather than ldf's L) keeps the train and
    # test factor scores on the SAME scale, because the held-out scores in
    # predict_test() are obtained by the identical projection operator.
    Ltr_proj <- Ytr %*% EF %*% solve(crossprod(EF) + diag(1e-8, ncol(EF)))
    cox_coef <- tryCatch(
      as.numeric(coef(survival::coxph(survival::Surv(ttr, str_) ~ Ltr_proj))),
      error = function(e) rep(NA_real_, ncol(EF)))
    list(EF = EF, EL = EL, EBeta = NULL, cox_coef = cox_coef,
         EF_norms = NULL, EF_cohort = NULL, n_iter = f$n_factors)
  }
}

#' Held-out test predictions for an arm.  Supervised arms use the joint model's
#' β; the EBMF arm uses its TWO-STAGE downstream Cox (stage 2) on held-out scores
#' obtained by projecting test Y onto the EBMF programs.  Returns the test risk
#' scores together with their survival outcome and cohort so the report can
#' compute the C-index AND draw a held-out risk-stratified Kaplan-Meier curve.
predict_test <- function(arm, fit, d, spl) {
  Yte <- d$Y[spl$test_idx, , drop = FALSE]
  tte <- d$time[spl$test_idx]
  ste <- d$status[spl$test_idx]

  if (arm == "EBMF") {
    # Two-stage baseline: project held-out Y onto the EBMF programs (same
    # projection used to build the training scores in fit_arm), then score with
    # the downstream Cox coefficients.  NA coefficients (rank-deficient fit) -> 0.
    if (is.null(fit$cox_coef) || all(is.na(fit$cox_coef))) return(NULL)
    EF <- fit$EF
    cc <- fit$cox_coef; cc[is.na(cc)] <- 0
    Lte_proj <- Yte %*% EF %*% solve(crossprod(EF) + diag(1e-8, ncol(EF)))
    return(list(risk = as.vector(Lte_proj %*% cc), time = tte, status = ste,
                cohort = d$cohort_id[spl$test_idx]))
  }

  # arm == "YFB_base" is the only other case now (cohort-indicator arms dropped).
  pr <- predict_cox_on_yf(Yte, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)

  list(risk = pr$risk_scores, time = tte, status = ste,
       cohort = d$cohort_id[spl$test_idx])
}

# --------------------------------------------------------------------------
# 4. Main loop: scenario × seed × arm
# --------------------------------------------------------------------------
rows        <- list()
example_fit <- list()   # store data + fits for the first seed of each scenario

for (sc in names(SCENARIOS)) {
  part <- SCENARIOS[[sc]]
  cat(sprintf("--- Scenario: %-15s (K_shared=%d, K_specific=%s) ---\n",
              sc, part$K_shared, paste(part$K_specific, collapse = ",")))

  for (s in SEEDS) {
    d <- generate_multicohort_data(
      C = C, n_per = N_PER, p = P,
      K_shared = part$K_shared, K_specific = part$K_specific,
      F_templates = F_TEMPLATES, a_shared = A_SHARED, a_specific = A_SPECIFIC,
      offset_sd = OFFSET_SD, target_censoring = TARGET_CENS, seed = s)

    spl    <- stratified_split(d$status, test_frac = 0.25, seed = s)
    cid_tr <- droplevels(d$cohort_id[spl$train_idx])

    for (k_name in names(K_SETTINGS)) {
      K_here <- K_SETTINGS[[k_name]]

      for (arm in ARMS) {
        if (arm == "EBMF" && !HAVE_FLASHIER) next

        fit <- tryCatch(fit_arm(arm, d, spl, K = K_here),
                        error = function(e) { message(arm, " (K=", K_here, ") failed: ", e$message); NULL })
        if (is.null(fit)) next

        mf   <- match_factors(fit$EF, d$F_true)
        est  <- classify_specificity(fit$EL, cid_tr)
        sa   <- specificity_accuracy(est, mf$match, d$factor_labels)
        # beta_recovery()'s BETA_THRESH cutoff (0.001, the same ARD
        # survival-active threshold used on real data) is exactly the
        # ARD classification this K-sweep is meant to test: at K_here > true
        # K, does a survival-irrelevant extra factor stay below threshold
        # (true negative) or get falsely activated (false positive)?
        br   <- beta_recovery(fit$EBeta, mf$match, d$factor_labels, BETA_THRESH)
        pred <- predict_test(arm, fit, d, spl)
        ci   <- if (is.null(pred)) NA_real_
                else oriented_cindex(pred$risk, pred$time, pred$status)

        is_sh <- d$factor_labels == "shared"
        REC_CUT <- 0.7   # |cor| above which a true factor counts as "recovered"
        rows[[length(rows) + 1]] <- data.frame(
          scenario      = sc,
          k_setting     = k_name,
          K_init        = K_here,
          arm           = arm,
          seed          = s,
          rec_shared    = if (any(is_sh))  mean(mf$best_cor[is_sh])  else NA_real_,
          rec_specific  = if (any(!is_sh)) mean(mf$best_cor[!is_sh]) else NA_real_,
          frac_shared   = if (any(is_sh))  mean(mf$best_cor[is_sh]  > REC_CUT) else NA_real_,
          frac_specific = if (any(!is_sh)) mean(mf$best_cor[!is_sh] > REC_CUT) else NA_real_,
          spec_acc      = sa$accuracy,
          beta_tp_rate  = br$tp_rate,
          beta_fp_rate  = br$fp_rate,
          beta_shared   = br$mean_abs_shared,
          beta_specific = br$mean_abs_specific,
          c_index       = ci,
          n_iter        = fit$n_iter %||% NA_integer_,
          stringsAsFactors = FALSE
        )

        # keep one example per scenario/K-setting/arm (first seed) for report
        # figures. store the split (so EL rows align to cohorts for the
        # loading heatmap) and the held-out predictions (for the risk-
        # stratified Kaplan-Meier figure).
        if (s == SEEDS[1]) {
          if (is.null(example_fit[[k_name]])) example_fit[[k_name]] <- list()
          if (is.null(example_fit[[k_name]][[sc]]))
            example_fit[[k_name]][[sc]] <- list(data = d, split = spl)
          example_fit[[k_name]][[sc]][[arm]] <- list(EF = fit$EF, EL = fit$EL,
                                           EBeta = fit$EBeta, match = mf,
                                           est_labels = est, pred = pred)
        }
      }
      cat(sprintf("  seed %d, K_init=%d (%s) done.\n", s, K_here, k_name))
    }
  }
  cat("\n")
}

results <- do.call(rbind, rows)

# --------------------------------------------------------------------------
# 5. Save outputs
# --------------------------------------------------------------------------
write.csv(results, file.path(OUT_DIR, "multicohort_sim_results.csv"), row.names = FALSE)
saveRDS(example_fit, file.path(OUT_DIR, "multicohort_sim_example.rds"))

# --------------------------------------------------------------------------
# 6. Console summary (mean over seeds)
# --------------------------------------------------------------------------
cat("============================================================\n")
cat(" Mean over seeds (rec_shared / rec_specific / spec_acc / C / FP)\n")
cat("============================================================\n")
agg <- aggregate(cbind(rec_shared, rec_specific, spec_acc, c_index, beta_fp_rate) ~ scenario + k_setting + K_init + arm,
                 data = results, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
agg <- agg[order(agg$scenario, agg$K_init, agg$arm), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  %-15s K_init=%-2d(%-9s) %-9s recS=%.3f recSp=%s specAcc=%.3f C=%s FP=%s\n",
              agg$scenario[i], agg$K_init[i], agg$k_setting[i], agg$arm[i],
              agg$rec_shared[i],
              ifelse(is.nan(agg$rec_specific[i]), "  NA", sprintf("%.3f", agg$rec_specific[i])),
              agg$spec_acc[i],
              ifelse(is.na(agg$c_index[i]), " NA", sprintf("%.3f", agg$c_index[i])),
              ifelse(is.nan(agg$beta_fp_rate[i]), " NA", sprintf("%.2f", agg$beta_fp_rate[i]))))
}
cat(sprintf("\nResults: %s\n", file.path(OUT_DIR, "multicohort_sim_results.csv")))
cat("============================================================\n")
