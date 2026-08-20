# ============================================================
# Script:  results/multi_cohort_sim/run_k7_signal_sweep.R
# Purpose: YFB vs. two-step (EBMF + Cox) as true survival signal strength is
#          scaled from 0 (off) to a realistic magnitude, using a DGP that
#          matches the RECOMMENDED real-data configuration's factor structure
#          exactly: K_init=7 (over-specified, ARD-pruned), 2 factors with a
#          real, non-zero survival coefficient and 2 more that are real
#          gene-expression programs with NO survival effect (genomics-only)
#          -- i.e. K_eff_total=4 (DECISIONS.md 2026-08-19).
#
#          This is a more faithful analog of the recommended model's
#          structure than the existing survival_strength_sweep (which uses
#          K_shared=4 with ALL 4 factors prognostic, no genomics-only
#          factors, and fits at K_FIT=K_true directly rather than
#          over-specifying to K=7). Requested to make this comparison
#          directly for the 8/21 meeting (was never discussed at the 8/3
#          meeting despite being requested).
#
#          Single simulated cohort (C=1) via generate_multicohort_data() --
#          with C=1, "study-specific" factors are never block-zeroed (there
#          is only one cohort to be "specific" to), so they behave exactly
#          like normal, real, whole-sample gene-expression programs; their
#          beta stays fixed at 0 (specific_prognostic=FALSE) regardless of
#          strength, which is exactly the "genomics-only, never prognostic"
#          role Programs other than 3/7 play in the recommended model.
#
#          Both arms fit at K_INIT=7 (matches the recommended real-data
#          procedure of over-specifying and letting ARD prune), AND at
#          K_TRUE=4 directly (matching how the July 15 sweep fit both arms,
#          with no spare capacity) -- the same simulated dataset is used for
#          both K settings per seed, isolating whether over-specification
#          itself dilutes YFB's advantage once real non-prognostic factors
#          are present, separately from the mere presence of those factors.
#
#   Output: results/multi_cohort_sim/outputs/k7_signal_sweep_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-20
# Usage:   Rscript results/multi_cohort_sim/run_k7_signal_sweep.R [--quick]
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
source("code/select_K.R")   # classify_factors()

source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/build_ebmf_templates.R")
source("results/multi_cohort_sim/sim_scoring.R")

HAVE_FLASHIER <- requireNamespace("flashier", quietly = TRUE)
if (!HAVE_FLASHIER) stop("flashier package required for the EBMF two-step baseline.")

mc <- cfg$synthetic_multicohort
BETA_THRESH <- cfg$k_selection$beta_threshold
PVE_THRESH  <- cfg$k_selection$pve_threshold
ALPHA       <- cfg$benchmark$alpha
PRIOR_BETA  <- "normal"
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
SEEDS       <- if (QUICK_MODE) 1:2 else 1:10

# --- DGP: matches the recommended model's K_eff_total=4 structure exactly ---
C           <- 1L                        # single cohort -- "specific" factors are never zeroed
N_TOTAL     <- 300L                      # single-cohort sample size, comparable to n_per*C elsewhere
P           <- mc$p
A_SHARED    <- mc$a_shared               # amplitude for the 2 survival-active factors
A_SPECIFIC  <- mc$a_shared               # same amplitude for the 2 genomics-only factors (fair comparison)
K_SHARED    <- 2L                        # survival-active factor count (matches K_eff_survival=2)
K_SPECIFIC  <- 2L                        # genomics-only factor count (matches K_eff_genomics=2 at K_init=7)
K_INIT      <- 7L                        # over-specify to K=7, exactly like the real recommended fit
K_TRUE      <- K_SHARED + K_SPECIFIC     # =4; fit directly at the true count, no spare capacity --
                                          # isolates whether over-specification itself (not just the
                                          # presence of non-prognostic factors) dilutes YFB's advantage
BETA_BASE   <- c(1.5, -1.2)              # scaled by strength; K_SPECIFIC factors stay at beta=0 always
STRENGTHS   <- c(0.0, 0.25, 0.5, 1.0, 2.0)
TARGET_CENS <- mc$target_censoring

OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat(" K=7 Signal-Strength Sweep -- YFB vs. Two-Step (EBMF+Cox)\n")
cat(sprintf(" DGP: K_shared=%d (survival-active), K_specific=%d (genomics-only), K_INIT=%d\n",
            K_SHARED, K_SPECIFIC, K_INIT))
cat(sprintf(" strength grid: %s | seeds: %s\n", paste(STRENGTHS, collapse=", "), paste(SEEDS, collapse=",")))
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# EBMF templates (real gene-program shapes; NULL -> synthetic fallback)
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
# Per-arm fit + prediction (YFB and EBMF+Cox only, both at K_INIT=7)
# --------------------------------------------------------------------------
fit_yfb <- function(Ytr, ttr, str_, K) {
  f <- suppressMessages(fit_cox_on_yf(
    Ytr, ttr, str_, K = K, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA, verbose = FALSE))
  list(EF = f$EF, EL = f$EL, EBeta = f$EBeta, EF_norms = f$EF_norms, n_iter = f$history$n_iter)
}

fit_ebmf <- function(Ytr, ttr, str_, K) {
  ebnm_args <- if (requireNamespace("ebnm", quietly = TRUE))
    list(ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential))
  else list()
  f <- do.call(flashier::flash,
               c(list(Ytr, var_type = 2, greedy_Kmax = K, backfit = TRUE, verbose = 0), ebnm_args))
  l  <- flashier::ldf(f, type = "2")
  EF <- as.matrix(l$F)
  Ltr_proj <- Ytr %*% EF %*% solve(crossprod(EF) + diag(1e-8, ncol(EF)))
  cox_coef <- tryCatch(as.numeric(coef(coxph(Surv(ttr, str_) ~ Ltr_proj))),
                       error = function(e) rep(NA_real_, ncol(EF)))
  list(EF = EF, cox_coef = cox_coef, n_iter = f$n_factors)
}

predict_yfb  <- function(fit, Yte) predict_cox_on_yf(Yte, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
predict_ebmf <- function(fit, Yte) {
  if (all(is.na(fit$cox_coef))) return(rep(NA_real_, nrow(Yte)))
  cc <- fit$cox_coef; cc[is.na(cc)] <- 0
  Lte_proj <- Yte %*% fit$EF %*% solve(crossprod(fit$EF) + diag(1e-8, ncol(fit$EF)))
  as.vector(Lte_proj %*% cc)
}

# --------------------------------------------------------------------------
# Main loop: strength x seed
# --------------------------------------------------------------------------
rows <- list()

for (strength in STRENGTHS) {
  cat(sprintf("--- strength=%.2f (beta_shared = %s) ---\n",
              strength, paste(round(strength * BETA_BASE, 3), collapse = ", ")))

  for (s in SEEDS) {
    d <- generate_multicohort_data(
      C = C, n_per = N_TOTAL, p = P,
      K_shared = K_SHARED, K_specific = K_SPECIFIC,
      F_templates = F_TEMPLATES, a_shared = A_SHARED, a_specific = A_SPECIFIC,
      beta_shared = strength * BETA_BASE, specific_prognostic = FALSE,
      target_censoring = TARGET_CENS, seed = s
    )
    spl <- stratified_split(d$status, test_frac = 0.25, seed = s)
    Ytr <- d$Y[spl$train_idx, , drop = FALSE]; ttr <- d$time[spl$train_idx]; str_tr <- d$status[spl$train_idx]
    Yte <- d$Y[spl$test_idx, , drop = FALSE];  tte <- d$time[spl$test_idx];  ste_te <- d$status[spl$test_idx]

    # K_INIT=7 (over-specified, matches the real recommended procedure) and
    # K_INIT=K_TRUE=4 (fit directly at the true factor count, matching how
    # the July 15 sweep fit both arms) -- same simulated dataset for both,
    # to isolate whether over-specification itself dilutes the joint
    # model's advantage once real non-prognostic factors are present.
    for (k_setting in list(list(K = K_INIT, tag = "K7"), list(K = K_TRUE, tag = "Ktrue"))) {
      K   <- k_setting$K
      tag <- k_setting$tag

      fit_y <- tryCatch(fit_yfb(Ytr, ttr, str_tr, K), error = function(e) { message("YFB failed: ", e$message); NULL })
      fit_e <- tryCatch(fit_ebmf(Ytr, ttr, str_tr, K), error = function(e) { message("EBMF failed: ", e$message); NULL })

      if (!is.null(fit_y)) {
        cls <- classify_factors(fit_y, Ytr, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
        c_yfb <- oriented_cindex(predict_yfb(fit_y, Yte), tte, ste_te)
        rows[[length(rows) + 1]] <- data.frame(
          strength = strength, seed = s, arm = paste0("YFB_", tag), K_fit = K, c_index = c_yfb,
          k_survival_active = sum(cls$category == "survival_active"),
          k_genomics_only   = sum(cls$category == "genomics_only"),
          stringsAsFactors = FALSE
        )
      }
      if (!is.null(fit_e)) {
        c_ebmf <- oriented_cindex(predict_ebmf(fit_e, Yte), tte, ste_te)
        rows[[length(rows) + 1]] <- data.frame(
          strength = strength, seed = s, arm = paste0("EBMF_", tag), K_fit = K, c_index = c_ebmf,
          k_survival_active = NA_integer_, k_genomics_only = NA_integer_,
          stringsAsFactors = FALSE
        )
      }
      cat(sprintf("  seed %d [K=%d/%s]: YFB C=%.3f (K_surv=%s) | EBMF+Cox C=%.3f\n", s, K, tag,
                  if (!is.null(fit_y)) c_yfb else NA,
                  if (!is.null(fit_y)) sum(cls$category == "survival_active") else NA,
                  if (!is.null(fit_e)) c_ebmf else NA))
    }
  }
  cat("\n")
}

results <- do.call(rbind, rows)
out_csv <- file.path(OUT_DIR, "k7_signal_sweep_results.csv")
write.csv(results, out_csv, row.names = FALSE)

cat("============================================================\n")
cat(" Mean over seeds: C-index by strength x arm\n")
cat("============================================================\n")
agg <- aggregate(c_index ~ strength + arm, data = results, FUN = mean, na.action = na.pass)
agg <- agg[order(agg$strength, agg$arm), ]
for (i in seq_len(nrow(agg)))
  cat(sprintf("  strength=%.2f  %-10s  C=%.3f\n", agg$strength[i], agg$arm[i], agg$c_index[i]))

cat(sprintf("\nResults: %s\n", out_csv))
