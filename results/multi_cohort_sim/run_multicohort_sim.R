# ============================================================
# Script:  results/multi_cohort_sim/run_multicohort_sim.R
# Purpose: Multi-cohort simulation study.  For each of three scenarios
#          (all-shared, hybrid, nothing-shared) and several seeds, fit five
#          arms and score how well each recovers the shared vs. study-specific
#          factor structure and the survival coefficients:
#
#            YFB_base    — YFB model (η = (YF)β), no cohort indicator
#            YFB_cohort  — YFB model + cohort_id (corner-point dummy)
#            LB_base     — LB model  (η = Lβ), no cohort indicator
#            LB_cohort   — LB model  + cohort_id
#            EBMF        — unsupervised benchmark (flashier, survival-blind)
#
#          Metrics: shared/specific factor recovery (|cor| of gene programs),
#          specificity-classification accuracy, β recovery (TP/FP prognostic
#          rates), and orientation-free held-out C-index.
#
#          Output: results/multi_cohort_sim/outputs/
#            multicohort_sim_results.csv   (one row per scenario × arm × seed)
#            multicohort_sim_example.rds   (data + fits, first seed, for figures)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-14
# Usage:   Rscript results/multi_cohort_sim/run_multicohort_sim.R [--quick]
# ============================================================

# --------------------------------------------------------------------------
# 0. Setup
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
K_FIT       <- mc$k_fit
A_SHARED    <- mc$a_shared
A_SPECIFIC  <- mc$a_specific
OFFSET_SD   <- mc$offset_sd
TARGET_CENS <- mc$target_censoring
EBMF_KMAX   <- mc$ebmf_kmax
BETA_THRESH <- cfg$k_selection$beta_threshold
ALPHA       <- cfg$benchmark$alpha
PRIOR_BETA  <- "normal"          # avoids β→0 collapse (run_synthetic.R rationale)
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
SEEDS       <- if (QUICK_MODE) unlist(mc$seeds)[1] else unlist(mc$seeds)

SCENARIOS <- lapply(mc$scenarios, function(s)
  list(K_shared = s$K_shared, K_specific = unlist(s$K_specific)))

ARMS <- c("YFB_base", "YFB_cohort", "LB_base", "LB_cohort", "EBMF")

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" Multi-Cohort Simulation — Shared vs. Study-Specific Factors\n")
cat(sprintf(" mode=%s | C=%d | n_per=%s | p=%d | K_fit=%d | seeds=%s\n",
            if (QUICK_MODE) "QUICK" else "FULL", C, paste(N_PER, collapse=","),
            P, K_FIT, paste(SEEDS, collapse=",")))
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

#' Fit one arm on the training split.  Returns a list with EF (p×K), EL (n_tr×K),
#' EBeta (or NULL for EBMF), EF_norms, EF_cohort (or NULL), n_iter.
fit_arm <- function(arm, d, spl) {
  Ytr    <- d$Y[spl$train_idx, , drop = FALSE]
  ttr    <- d$time[spl$train_idx]
  str_   <- d$status[spl$train_idx]
  cid_tr <- droplevels(d$cohort_id[spl$train_idx])

  if (arm %in% c("YFB_base", "YFB_cohort")) {
    cid <- if (arm == "YFB_cohort") cid_tr else NULL
    f <- suppressMessages(fit_cox_on_yf(
      Ytr, ttr, str_, K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
      prior_beta = PRIOR_BETA, verbose = FALSE, cohort_id = cid))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = f$EF_norms,
         EF_cohort = f$EF_cohort, n_iter = f$history$n_iter)

  } else if (arm %in% c("LB_base", "LB_cohort")) {
    cid <- if (arm == "LB_cohort") cid_tr else NULL
    f <- suppressMessages(fit_supervised_mf_modular(
      Ytr, ttr, str_, K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
      prior_beta = PRIOR_BETA, verbose = FALSE, cohort_id = cid))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = NULL,
         EF_cohort = f$EF_cohort, n_iter = f$history$n_iter)

  } else if (arm == "EBMF") {
    # Unsupervised benchmark, extended to a TWO-STAGE baseline:
    #   stage 1 — survival-blind non-negative EBMF factorization of training Y;
    #   stage 2 — a Cox model fit on the EBMF factor scores.
    # This is the "factorize first, attach survival afterward" pipeline that the
    # joint SSBMF model is meant to outperform: EBMF chooses its factors WITHOUT
    # seeing the outcome, so it cannot preferentially sharpen the prognostic ones.
    # Non-negative priors (point_exponential on L and F) match the supervised arms,
    # so the only difference being tested is joint vs. two-stage supervision.
    ebnm_args <- if (requireNamespace("ebnm", quietly = TRUE))
      list(ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential))
    else list()
    f <- do.call(flashier::flash,
                 c(list(Ytr, var_type = 2, greedy_Kmax = K_FIT,
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

  # subtract estimated per-cohort offset for cohort arms (matches run_synthetic.R)
  if (!is.null(fit$EF_cohort) && ncol(fit$EF_cohort) > 0) {
    lev <- levels(droplevels(d$cohort_id[spl$train_idx]))
    cid_te <- factor(d$cohort_id[spl$test_idx], levels = lev)
    Lct <- model.matrix(~ cid_te)[, -1, drop = FALSE]
    Yte <- Yte - Lct %*% t(fit$EF_cohort)
  }

  pr <- if (grepl("^YFB", arm))
    predict_cox_on_yf(Yte, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)
  else
    predict_supervised_mf(Yte, fit$EF, fit$EBeta)

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

    for (arm in ARMS) {
      if (arm == "EBMF" && !HAVE_FLASHIER) next

      fit <- tryCatch(fit_arm(arm, d, spl),
                      error = function(e) { message(arm, " failed: ", e$message); NULL })
      if (is.null(fit)) next

      mf   <- match_factors(fit$EF, d$F_true)
      est  <- classify_specificity(fit$EL, cid_tr)
      sa   <- specificity_accuracy(est, mf$match, d$factor_labels)
      br   <- beta_recovery(fit$EBeta, mf$match, d$factor_labels, BETA_THRESH)
      pred <- predict_test(arm, fit, d, spl)
      ci   <- if (is.null(pred)) NA_real_
              else oriented_cindex(pred$risk, pred$time, pred$status)

      is_sh <- d$factor_labels == "shared"
      REC_CUT <- 0.7   # |cor| above which a true factor counts as "recovered"
      rows[[length(rows) + 1]] <- data.frame(
        scenario      = sc,
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

      # keep one example per scenario/arm (first seed) for report figures.
      # store the split (so EL rows align to cohorts for the loading heatmap) and
      # the held-out predictions (for the risk-stratified Kaplan-Meier figure).
      if (s == SEEDS[1]) {
        if (is.null(example_fit[[sc]])) example_fit[[sc]] <- list(data = d, split = spl)
        example_fit[[sc]][[arm]] <- list(EF = fit$EF, EL = fit$EL,
                                         EBeta = fit$EBeta, match = mf,
                                         est_labels = est, pred = pred)
      }
    }
    cat(sprintf("  seed %d done.\n", s))
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
agg <- aggregate(cbind(rec_shared, rec_specific, spec_acc, c_index, beta_fp_rate) ~ scenario + arm,
                 data = results, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
agg <- agg[order(agg$scenario, agg$arm), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  %-15s %-11s recS=%.3f recSp=%s specAcc=%.3f C=%s FP=%s\n",
              agg$scenario[i], agg$arm[i],
              agg$rec_shared[i],
              ifelse(is.nan(agg$rec_specific[i]), "  NA", sprintf("%.3f", agg$rec_specific[i])),
              agg$spec_acc[i],
              ifelse(is.na(agg$c_index[i]), " NA", sprintf("%.3f", agg$c_index[i])),
              ifelse(is.nan(agg$beta_fp_rate[i]), " NA", sprintf("%.2f", agg$beta_fp_rate[i]))))
}
cat(sprintf("\nResults: %s\n", file.path(OUT_DIR, "multicohort_sim_results.csv")))
cat("============================================================\n")
