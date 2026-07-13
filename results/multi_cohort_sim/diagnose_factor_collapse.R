# ============================================================
# Script:  results/multi_cohort_sim/diagnose_factor_collapse.R
# Purpose: Root-cause the factor-collapse failure mode discovered in the
#          `sparse_synthetic` scenario of run_survival_strength_sweep.R
#          (see Result 4, docs/reports/joint_vs_twostep_sweep_07_12_2026.qmd,
#          and DECISIONS.md 2026-07-12). Runs targeted checks per seed:
#            1. YFB at alpha=0 (survival term removed) -- still collapses?
#            2. LB (same update_L.R/update_F.R machinery) -- collapses worse?
#            3. Random init vs. SVD init -- does a better starting point help?
#            4. Deflation init (YFB and LB) -- does the K-parsimony follow-up
#               plan's Step 2 fix (code/deflation_init.R, DECISIONS.md
#               2026-07-13) resolve this collapse mode? Added 2026-07-13.
#          A factor is "dead" if every element of its EF column is <1e-8 in
#          absolute value (sd(x)==0 or all(abs(x)<1e-8)).
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   Rscript results/multi_cohort_sim/diagnose_factor_collapse.R
# ============================================================

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
source("results/multi_cohort_sim/generate_multicohort_data.R")
source("results/multi_cohort_sim/sim_scoring.R")

mc <- cfg$synthetic_multicohort

n_dead <- function(EF) sum(apply(EF, 2, function(x) sd(x) == 0 || all(abs(x) < 1e-8)))

SEEDS <- 1:5
rows <- list()

for (seed in SEEDS) {
  d <- generate_multicohort_data(
    C = mc$C, n_per = unlist(mc$n_per), p = mc$p,
    K_shared = 4, K_specific = c(0, 0), F_templates = NULL,
    a_shared = 12, a_specific = 12, active_rate = 0.02,
    beta_shared = 4 * c(1.5, -1.2, 0.8, -0.5),
    target_censoring = mc$target_censoring, seed = seed)
  spl <- stratified_split(d$status, test_frac = 0.25, seed = seed)
  Ytr <- d$Y[spl$train_idx, , drop = FALSE]
  ttr <- d$time[spl$train_idx]
  str_ <- d$status[spl$train_idx]

  # Check 1: YFB tuned alpha vs. alpha=0 (survival term removed)
  f_svd_tuned <- suppressMessages(fit_cox_on_yf(Ytr, ttr, str_, K = 4, max_iter = 300,
                    alpha = 0.5, prior_beta = "normal", verbose = FALSE, init_method = "svd"))
  f_svd_a0 <- suppressMessages(fit_cox_on_yf(Ytr, ttr, str_, K = 4, max_iter = 300,
                    alpha = 0, prior_beta = "normal", verbose = FALSE, init_method = "svd",
                    sign_correction = FALSE))

  # Check 2: LB (shared update_L.R/update_F.R machinery)
  f_lb <- suppressMessages(fit_supervised_mf_modular(Ytr, ttr, str_, K = 4, max_iter = 300,
                    alpha = 0.5, prior_beta = "normal", verbose = FALSE))

  # Check 3: random init vs. SVD init (YFB, tuned alpha)
  f_rnd_tuned <- suppressMessages(fit_cox_on_yf(Ytr, ttr, str_, K = 4, max_iter = 300,
                    alpha = 0.5, prior_beta = "normal", verbose = FALSE, init_method = "random"))

  # Check 4: deflation init (YFB and LB, tuned alpha) -- the candidate fix
  f_yfb_defl <- suppressMessages(fit_cox_on_yf(Ytr, ttr, str_, K = 4, max_iter = 300,
                    alpha = 0.5, prior_beta = "normal", verbose = FALSE, init_method = "deflation"))
  f_lb_defl  <- suppressMessages(fit_supervised_mf_modular(Ytr, ttr, str_, K = 4, max_iter = 300,
                    alpha = 0.5, prior_beta = "normal", verbose = FALSE, init_method = "deflation"))

  rows[[length(rows) + 1]] <- data.frame(
    seed = seed,
    yfb_svd_tuned_dead = n_dead(f_svd_tuned$EF), yfb_svd_tuned_niter = f_svd_tuned$history$n_iter,
    yfb_svd_a0_dead     = n_dead(f_svd_a0$EF),    yfb_svd_a0_niter    = f_svd_a0$history$n_iter,
    lb_svd_tuned_dead   = n_dead(f_lb$EF),        lb_svd_tuned_niter  = f_lb$history$n_iter,
    yfb_rnd_tuned_dead  = n_dead(f_rnd_tuned$EF), yfb_rnd_tuned_niter = f_rnd_tuned$history$n_iter,
    yfb_defl_dead       = n_dead(f_yfb_defl$EF),  yfb_defl_niter      = f_yfb_defl$history$n_iter,
    lb_defl_dead        = n_dead(f_lb_defl$EF),   lb_defl_niter       = f_lb_defl$history$n_iter,
    stringsAsFactors = FALSE
  )
  cat(sprintf("seed %d done.\n", seed))
}

results <- do.call(rbind, rows)
OUT_DIR <- "results/multi_cohort_sim/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(OUT_DIR, "factor_collapse_diagnostic.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat(" Factor-collapse diagnostic: dead factors (of 4) per check, per seed\n")
cat("============================================================\n")
print(results)
cat(sprintf("\nResults: %s\n", file.path(OUT_DIR, "factor_collapse_diagnostic.csv")))
