# ============================================================
# Script:  results/multi_cohort_sim/run_survival_strength_sweep.R
# Purpose: Phase 2 of docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md
#          -- joint model (YFB) vs. two-step baselines (EBMF+Cox, PCA+Cox) as
#          the true prognostic effect size is scaled from 0 (survival
#          completely off) to large, holding the genomics DGP fixed.
#
#          Success criteria (see the plan):
#            1. Equivalence at strength=0: joint (tuned alpha) vs both 2-step
#               baselines agree in C-index and factor recovery.
#            2. Separation growing with strength: joint's C-index advantage
#               over the 2-step baselines widens as the true effect grows.
#            3. Internal control: YFB's F is alpha-invariant by construction
#               (F is pure-genomics regardless of alpha in the YFB
#               reformulation) -- verified directly by comparing EF at
#               alpha=0 vs. the tuned alpha, not just asserted.
#
#          Comprehensive extension (2026-07-12, same-day follow-up): the
#          initial 5-seed/1-scenario run showed the right qualitative
#          pattern but was too thin to be strong evidence (wide SEs, no test
#          of whether the pattern is specific to one DGP structure). This
#          version adds:
#            - 10 seeds/condition (was 5) for tighter SEs.
#            - 3 additional DGP scenarios beyond the default (real EBMF
#              templates): sparse synthetic F, low SNR, more shared factors
#              (config/globals.yml, synthetic_multicohort.survival_strength_
#              sweep.scenarios) -- asks whether the joint-vs-2-step ordering
#              is a general pattern or an artifact of one specific DGP.
#          Deferred (flagged, not silently dropped): n/p variation toward
#          realistic PDAC scale, and a cohort-heterogeneity variant. Both are
#          worth doing but were scoped out of this pass to bound compute.
#
#          Output: results/multi_cohort_sim/outputs/
#            survival_strength_sweep_results.csv        (now has a `scenario` column)
#            survival_strength_sweep_alpha_invariance.csv (now has a `scenario` column)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-12
# Usage:   Rscript results/multi_cohort_sim/run_survival_strength_sweep.R [--quick]
# ============================================================

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
source("results/multi_cohort_sim/fit_pca_cox.R")

HAVE_FLASHIER <- requireNamespace("flashier", quietly = TRUE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --------------------------------------------------------------------------
# 1. Configuration
# --------------------------------------------------------------------------
mc  <- cfg$synthetic_multicohort
sss <- mc$survival_strength_sweep

C           <- mc$C
N_PER       <- unlist(mc$n_per)
P           <- mc$p
A_SHARED    <- mc$a_shared
TARGET_CENS <- mc$target_censoring
EBMF_KMAX   <- mc$ebmf_kmax
BETA_THRESH <- cfg$k_selection$beta_threshold
ALPHA_TUNED <- cfg$benchmark$alpha   # 0.5 -- the production default, our "tuned" alpha
PRIOR_BETA  <- "normal"
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
N_SEEDS     <- sss$n_seeds %||% length(mc$seeds)
SEEDS       <- if (QUICK_MODE) 1:2 else seq_len(N_SEEDS)

K_SPECIFIC  <- unlist(sss$K_specific)
STRENGTHS   <- unlist(sss$strength_values)
BETA_BASE_DEFAULT <- c(1.5, -1.2, 0.8, -0.5)   # default shared-factor coefficients, scaled by strength

# --- DGP scenario definitions (config/globals.yml, ...survival_strength_sweep.scenarios) ---
# Each scenario overrides a subset of {K_shared, a_shared, f_templates, active_rate}
# on top of the defaults below; asks whether the joint-vs-2-step ordering is a
# general pattern or an artifact of the single default DGP structure.
SCENARIOS_CFG <- sss$scenarios %||% list(default = list())
if (QUICK_MODE) SCENARIOS_CFG <- SCENARIOS_CFG["default"]   # quick mode: default scenario only

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" Survival-Strength Sweep — Joint (YFB) vs. Two-Step Baselines\n")
cat(sprintf(" strength grid: %s\n", paste(STRENGTHS, collapse = ", ")))
cat(sprintf(" scenarios: %s\n", paste(names(SCENARIOS_CFG), collapse = ", ")))
cat(sprintf(" alpha_tuned=%.2f | mode=%s | seeds=%s\n",
            ALPHA_TUNED, if (QUICK_MODE) "QUICK" else "FULL",
            paste(SEEDS, collapse = ",")))
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# 2. EBMF gene-program templates (real-data ground truth; NULL -> synthetic)
#    Built once; individual scenarios opt out via f_templates: false to force
#    the synthetic sparse fallback instead (see SCENARIOS_CFG below).
# --------------------------------------------------------------------------
templates  <- tryCatch(build_ebmf_templates(Kmax = EBMF_KMAX),
                       error = function(e) { message("EBMF templates: ", e$message); NULL })
F_TEMPLATES_REAL <- if (!is.null(templates)) {
  cat(sprintf("Using EBMF templates: %d programs.\n\n", templates$K_ebmf))
  Ft <- templates$F
  Ft <- if (nrow(Ft) >= P) Ft[seq_len(P), , drop = FALSE]
        else               Ft[((seq_len(P) - 1L) %% nrow(Ft)) + 1L, , drop = FALSE]
  dimnames(Ft) <- NULL
  abs(Ft)
} else {
  cat("No EBMF templates available -- 'default' scenario will use the synthetic fallback too.\n\n")
  NULL
}

# --------------------------------------------------------------------------
# 3. Per-arm fitting + prediction helpers
# --------------------------------------------------------------------------

#' Fit one arm on the training split. Returns EF (p x K), EL (n_tr x K or NULL),
#' EBeta (or NULL), plus arm-specific extras needed for held-out prediction.
fit_arm <- function(arm, d, spl) {
  Ytr  <- d$Y[spl$train_idx, , drop = FALSE]
  ttr  <- d$time[spl$train_idx]
  str_ <- d$status[spl$train_idx]

  if (arm == "YFB_tuned") {
    f <- suppressMessages(fit_cox_on_yf(
      Ytr, ttr, str_, K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA_TUNED,
      prior_beta = PRIOR_BETA, verbose = FALSE))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = f$EF_norms, n_iter = f$history$n_iter)

  } else if (arm == "YFB_alpha0") {
    # Internal control: alpha=0 should not change YFB's F at all (F is pure
    # genomics regardless of alpha under the YFB reformulation) -- verified
    # directly below (section 4b), not merely asserted.
    f <- suppressMessages(fit_cox_on_yf(
      Ytr, ttr, str_, K = K_FIT, max_iter = MAX_ITER, alpha = 0,
      prior_beta = PRIOR_BETA, verbose = FALSE, sign_correction = FALSE))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = f$EF_norms, n_iter = f$history$n_iter)

  } else if (arm == "LB_tuned") {
    f <- suppressMessages(fit_supervised_mf_modular(
      Ytr, ttr, str_, K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA_TUNED,
      prior_beta = PRIOR_BETA, verbose = FALSE))
    list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = NULL, n_iter = f$history$n_iter)

  } else if (arm == "PCA_Cox") {
    f <- fit_pca_cox(Ytr, ttr, str_, K = K_FIT)
    list(EF = f$EF, EL = NULL, EBeta = NULL, cox_coef = f$cox_coef, center = f$center,
         EF_norms = NULL, n_iter = NA_integer_)

  } else if (arm == "EBMF") {
    ebnm_args <- if (requireNamespace("ebnm", quietly = TRUE))
      list(ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential))
    else list()
    f <- do.call(flashier::flash,
                 c(list(Ytr, var_type = 2, greedy_Kmax = K_FIT,
                        backfit = TRUE, verbose = 0), ebnm_args))
    l  <- flashier::ldf(f, type = "2")
    EF <- as.matrix(l$F)
    EL <- as.matrix(l$L)
    Ltr_proj <- Ytr %*% EF %*% solve(crossprod(EF) + diag(1e-8, ncol(EF)))
    cox_coef <- tryCatch(
      as.numeric(coef(coxph(Surv(ttr, str_) ~ Ltr_proj))),
      error = function(e) rep(NA_real_, ncol(EF)))
    list(EF = EF, EL = EL, EBeta = NULL, cox_coef = cox_coef, EF_norms = NULL,
         n_iter = f$n_factors)
  }
}

#' Held-out test predictions. Supervised arms use their own beta; the 2-step
#' arms (PCA_Cox, EBMF) project held-out Y with the SAME operator built on
#' training data, then score with the stage-2 Cox coefficients.
predict_test <- function(arm, fit, d, spl) {
  Yte <- d$Y[spl$test_idx, , drop = FALSE]
  tte <- d$time[spl$test_idx]
  ste <- d$status[spl$test_idx]

  if (arm == "PCA_Cox") {
    pr <- predict_pca_cox(list(EF = fit$EF, center = fit$center, cox_coef = fit$cox_coef), Yte)
    return(list(risk = pr$risk_scores, time = tte, status = ste))
  }

  if (arm == "EBMF") {
    if (is.null(fit$cox_coef) || all(is.na(fit$cox_coef))) return(NULL)
    cc <- fit$cox_coef; cc[is.na(cc)] <- 0
    Lte_proj <- Yte %*% fit$EF %*% solve(crossprod(fit$EF) + diag(1e-8, ncol(fit$EF)))
    return(list(risk = as.vector(Lte_proj %*% cc), time = tte, status = ste))
  }

  pr <- if (grepl("^YFB", arm))
    predict_cox_on_yf(Yte, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)
  else
    predict_supervised_mf(Yte, fit$EF, fit$EBeta)

  list(risk = pr$risk_scores, time = tte, status = ste)
}

ARMS <- c("YFB_tuned", "YFB_alpha0", "LB_tuned", "PCA_Cox", "EBMF")

# --------------------------------------------------------------------------
# 4. Main loop: scenario x strength x seed x arm
# --------------------------------------------------------------------------
rows       <- list()
alpha_rows <- list()   # section 4b: YFB alpha-invariance check

for (scen_name in names(SCENARIOS_CFG)) {
  ov <- SCENARIOS_CFG[[scen_name]]

  # --- resolve this scenario's DGP parameters against the defaults ---
  K_SHARED    <- ov$K_shared %||% sss$K_shared
  K_FIT       <- K_SHARED + sum(K_SPECIFIC)
  A_SHARED    <- ov$a_shared %||% mc$a_shared
  ACTIVE_RATE <- ov$active_rate %||% mc$active_rate
  USE_TEMPLATES <- ov$f_templates %||% TRUE
  F_TEMPLATES <- if (USE_TEMPLATES) F_TEMPLATES_REAL else NULL
  BETA_BASE   <- rep_len(BETA_BASE_DEFAULT, K_SHARED)   # recycled for high_K (K_shared=8)

  cat(sprintf("============ scenario: %s (K_shared=%d, a_shared=%.1f, templates=%s) ============\n",
              scen_name, K_SHARED, A_SHARED, USE_TEMPLATES))

for (strength in STRENGTHS) {
  cat(sprintf("--- strength=%.2f (beta_shared = %s) ---\n",
              strength, paste(round(strength * BETA_BASE, 3), collapse = ", ")))

  for (s in SEEDS) {
    d <- generate_multicohort_data(
      C = C, n_per = N_PER, p = P,
      K_shared = K_SHARED, K_specific = K_SPECIFIC,
      F_templates = F_TEMPLATES, a_shared = A_SHARED, a_specific = A_SHARED,
      active_rate = ACTIVE_RATE,
      beta_shared = strength * BETA_BASE,
      target_censoring = TARGET_CENS, seed = s)

    spl <- stratified_split(d$status, test_frac = 0.25, seed = s)

    fits <- list()
    for (arm in ARMS) {
      if (arm == "EBMF" && !HAVE_FLASHIER) next
      fit <- tryCatch(fit_arm(arm, d, spl),
                      error = function(e) { message(arm, " failed: ", e$message); NULL })
      if (is.null(fit)) next
      fits[[arm]] <- fit

      mf   <- match_factors(fit$EF, d$F_true)
      pred <- predict_test(arm, fit, d, spl)
      ci   <- if (is.null(pred)) NA_real_ else oriented_cindex(pred$risk, pred$time, pred$status)

      rows[[length(rows) + 1]] <- data.frame(
        scenario      = scen_name,
        strength      = strength,
        arm           = arm,
        seed          = s,
        rec_mean      = mean(mf$best_cor),      # K_specific=[0,0] -> all factors "shared"
        frac_recov    = mean(mf$best_cor > 0.7),
        c_index       = ci,
        n_iter        = fit$n_iter %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    }

    # 4b. Internal control: YFB's EF should be IDENTICAL at alpha=0 vs. the
    # tuned alpha, since F never depends on alpha under YFB (F is pure
    # genomics; only beta responds to alpha). Verified directly with a
    # DEDICATED, controlled pair of fits -- tol=-1 forces both to run the
    # exact same number of iterations, isolating the true alpha effect from
    # a confound where different alpha values give different elbo_full
    # trajectories and hence different natural convergence timing (which by
    # itself makes EF differ trivially, since it's then a snapshot at a
    # different point in an otherwise-identical trajectory, not a real
    # alpha-driven difference). The main sweep's YFB_tuned/YFB_alpha0 fits
    # above use natural convergence (appropriate for the C-index/recovery
    # comparison) and are NOT reused here for that reason.
    ctrl_iter <- 40L
    f_ctrl_tuned <- suppressMessages(fit_cox_on_yf(
      d$Y[spl$train_idx, , drop = FALSE], d$time[spl$train_idx], d$status[spl$train_idx],
      K = K_FIT, max_iter = ctrl_iter, tol = -1, alpha = ALPHA_TUNED,
      prior_beta = PRIOR_BETA, verbose = FALSE, sign_correction = FALSE))
    f_ctrl_a0 <- suppressMessages(fit_cox_on_yf(
      d$Y[spl$train_idx, , drop = FALSE], d$time[spl$train_idx], d$status[spl$train_idx],
      K = K_FIT, max_iter = ctrl_iter, tol = -1, alpha = 0,
      prior_beta = PRIOR_BETA, verbose = FALSE, sign_correction = FALSE))
    alpha_rows[[length(alpha_rows) + 1]] <- data.frame(
      scenario = scen_name, strength = strength, seed = s,
      max_abs_EF_diff = max(abs(f_ctrl_tuned$EF - f_ctrl_a0$EF)),
      stringsAsFactors = FALSE
    )

    cat(sprintf("  seed %d done.\n", s))
  }
  cat("\n")
}
}   # end scenario loop

results    <- do.call(rbind, rows)
alpha_inv  <- do.call(rbind, alpha_rows)

# --------------------------------------------------------------------------
# 5. Save
# --------------------------------------------------------------------------
write.csv(results, file.path(OUT_DIR, "survival_strength_sweep_results.csv"), row.names = FALSE)
write.csv(alpha_inv, file.path(OUT_DIR, "survival_strength_sweep_alpha_invariance.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# 6. Console summary
# --------------------------------------------------------------------------
cat("============================================================\n")
cat(" Mean over seeds: C-index and factor recovery by scenario x strength x arm\n")
cat("============================================================\n")
agg <- aggregate(cbind(c_index, rec_mean) ~ scenario + strength + arm, data = results,
                 FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
agg <- agg[order(agg$scenario, agg$strength, agg$arm), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  [%-16s] strength=%.2f  %-11s  C=%.3f  rec=%.3f\n",
              agg$scenario[i], agg$strength[i], agg$arm[i], agg$c_index[i], agg$rec_mean[i]))
}

cat("\n============================================================\n")
cat(" Internal control: max|EF(alpha=0.5) - EF(alpha=0)| for YFB\n")
cat("============================================================\n")
print(alpha_inv)

cat(sprintf("\nResults: %s\n", file.path(OUT_DIR, "survival_strength_sweep_results.csv")))
cat("============================================================\n")
