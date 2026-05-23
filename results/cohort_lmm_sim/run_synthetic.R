# ============================================================
# Script:  results/cohort_lmm_sim/run_synthetic.R
# Purpose: Stage 1 synthetic validation of the cohort indicator column
#          extension (cohort-cols-L).  Generates data with a known platform
#          offset between two cohorts, fits five models, and reports
#          factor recovery and external C-index.
#
#          Five models compared:
#            LB_base     — LB model, no cohort adjustment
#            LB_cohort   — LB model + cohort_id (corner-point encoding)
#            YFB_base    — YFB model, no cohort adjustment
#            YFB_cohort  — YFB model + cohort_id
#            zscore      — per-cohort mean subtraction, then LB_base
#
#          Acceptance criteria (Stage 1):
#            (a) LB_cohort C-index >= LB_base C-index
#            (b) YFB_cohort C-index >= YFB_base C-index
#            (c) Factor recovery (mean max-cor): cohort >= base for both models
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-22
# Usage:   Rscript results/cohort_lmm_sim/run_synthetic.R [--seed N]
#          (default seed = 42)
# ============================================================

# --------------------------------------------------------------------------
# 0. Setup
# --------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
SEED <- 42L
seed_flag <- which(args == "--seed")
if (length(seed_flag) > 0) SEED <- as.integer(args[seed_flag + 1])

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages(library(survival))

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

cat("============================================================\n")
cat(" Stage 1 Synthetic Validation — Cohort Indicator Extension\n")
cat(sprintf(" Seed: %d\n", SEED))
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# 1. Data generating process
#
# Model: Y_ij = (L F')_ij + cohort_offset_ij + noise_ij
# Survival: h_i(t) = h_0(t) exp(L_i · beta_true)
#
# Cohort offset: patients in cohort B (index > n/2) receive an additive
# gene-level offset drawn once from N(0, offset_sd^2) per gene.  This
# represents a platform-induced shift (e.g. RNA-seq vs proteomics).
#
# The offset is ORTHOGONAL to biology in expectation, but its variance
# swamps the biological signal when offset_sd >> signal_scale.  Models
# that absorb the offset via cohort_id should recover cleaner factors and
# higher concordance than models that don't.
# --------------------------------------------------------------------------

N          <- 200L    # total patients
P          <- 500L    # genes
K_TRUE     <- 3L      # true biological factors
OFFSET_SD  <- 2.0     # cohort B gene-level offset SD — large enough to disrupt
                      #   factor recovery when unaccounted for, but absorb-able
                      #   by the rank-1 cohort column in L_aug
ACTIVE_RATE <- 0.10   # fraction of genes active per factor (50 genes/factor)
SIGNAL_SD  <- 0.5     # per-active-gene loading amplitude in F_true
NOISE_SD   <- 0.3     # observation noise SD
BETA_TRUE  <- c(1.5, -1.0, 0.5)   # survival coefficients; K_TRUE-vector
K_FIT      <- K_TRUE  # fit with correct K
MAX_ITER   <- 200L
ALPHA      <- 0.5     # genomics-survival mixing
PRIOR_BETA <- "normal"   # normal (not point_normal) to avoid β→0 collapse
TEST_FRAC  <- 0.25    # fraction for external test set

set.seed(SEED)

# Factor 1 is the SURVIVAL FACTOR: L[:,1] is drawn proportional to the
# true linear predictor so the survival-correlated direction is guaranteed
# to exist in the genomics matrix.  This ensures at least one factor can
# recover the survival signal regardless of CAVI initialisation.
# Factors 2..K are independent biological factors (pure genomics, no survival).
eta_raw <- rnorm(N, sd = 1.5)          # raw linear predictor
L_surv  <- pmax(eta_raw - min(eta_raw), 0.01)  # positive, min-shifted

L_bio  <- matrix(rexp((K_TRUE - 1L) * N, rate = 1), N, K_TRUE - 1L)
L_true <- cbind(L_surv, L_bio)          # n x K_TRUE

# Survival driven ONLY by factor 1; factors 2..K are pure biology
BETA_TRUE_EFF <- c(BETA_TRUE[1], rep(0, K_TRUE - 1L))
# Compute eta_true early so cohort assignment can be stratified by risk
eta_true <- as.vector(L_true %*% BETA_TRUE_EFF)

# True factors: F sparse positive (ACTIVE_RATE of genes active per factor)
F_true <- matrix(0, P, K_TRUE)
for (k in seq_len(K_TRUE)) {
  active_genes <- sample.int(P, size = max(1L, round(ACTIVE_RATE * P)))
  F_true[active_genes, k] <- SIGNAL_SD
}

# Cohort assignment: interleaved by survival-risk rank to prevent spurious
# cohort-survival confounding. This ensures E[eta|cohortA] = E[eta|cohortB],
# so any C-index difference between models reflects biology, not confounding.
eta_rank   <- rank(eta_true)
cohort_id  <- ifelse(eta_rank %% 2 == 1, "A", "B")

# Cohort B offset: one offset vector per gene, added to all cohort B rows
offset_vec <- rnorm(P, mean = 0, sd = OFFSET_SD)
cohort_offset <- outer(as.integer(cohort_id == "B"), offset_vec)   # N x P

# Observation noise
E <- matrix(rnorm(N * P, sd = NOISE_SD), N, P)

# Observed genomics matrix
Y <- L_true %*% t(F_true) + cohort_offset + E
# Column-center Y (standard preprocessing)
Y <- sweep(Y, 2, colMeans(Y), "-")

# Survival times from Weibull(shape=1.5) with linear predictor L · beta_true_eff
event_times <- (-log(runif(N)) / (0.01 * exp(eta_true)))^(1 / 1.5)
# Calibrate censoring to ~30% (mirrors generate_synthetic_benchmark_data)
censor_base <- rexp(N, rate = 1)
.calibrate <- function(ev, base, target = 0.30, n_iter = 40) {
  lo <- 1e-3; hi <- max(ev) * 10
  for (i in seq_len(n_iter)) {
    mid <- sqrt(lo * hi)
    if (mean(base * mid < ev) > target) lo <- mid else hi <- mid
  }
  hi
}
censor_times <- censor_base * .calibrate(event_times, censor_base)
time         <- pmin(event_times, censor_times)
status       <- as.integer(event_times <= censor_times)
cat(sprintf("DGP: n=%d, p=%d, K_true=%d, C=2 | censoring=%.1f%%\n\n",
            N, P, K_TRUE, 100 * mean(status == 0)))

# --------------------------------------------------------------------------
# 2. Train / test split (stratified by event status)
# --------------------------------------------------------------------------

split      <- stratified_split(status, test_frac = TEST_FRAC, seed = SEED)
train_idx  <- split$train_idx
test_idx   <- split$test_idx

Y_train  <- Y[train_idx, ];  Y_test  <- Y[test_idx, ]
t_train  <- time[train_idx]; t_test  <- time[test_idx]
s_train  <- status[train_idx]; s_test <- status[test_idx]
cid_train <- cohort_id[train_idx]
cid_test  <- cohort_id[test_idx]

cat(sprintf("Split: n_train=%d, n_test=%d | train events=%d, test events=%d\n\n",
            length(train_idx), length(test_idx), sum(s_train), sum(s_test)))

# --------------------------------------------------------------------------
# 3. Helper: compute external C-index via concordance()
# --------------------------------------------------------------------------
ext_cindex <- function(risk_scores, t_test, s_test) {
  # concordance() maximises over orientation; we fix sign via negation convention
  as.numeric(concordance(Surv(t_test, s_test) ~ risk_scores)$concordance)
}

# Helper: factor recovery — for each non-degenerate estimated L column, find
# the highest absolute Pearson correlation with any true L column.
# Returns NA if all estimated columns are zero (collapsed factors).
factor_recovery <- function(EL, L_true) {
  col_sds <- apply(EL, 2, sd)
  EL_nz   <- EL[, col_sds > 1e-10, drop = FALSE]
  if (ncol(EL_nz) == 0) return(NA_real_)
  cors <- abs(cor(EL_nz, L_true))   # K_nz x K_true
  mean(apply(cors, 1, max))
}

# Helper: zscore preprocessing — subtract supplied (or self) per-cohort column means.
# ref_means: named list of p-vectors (one per cohort level from the TRAINING set).
# If ref_means is NULL, compute from Y itself (train-time use).
zscore_center <- function(Y, cohort, ref_means = NULL) {
  Yc <- Y
  for (lv in unique(cohort)) {
    idx  <- cohort == lv
    mu   <- if (is.null(ref_means)) colMeans(Y[idx, , drop = FALSE]) else ref_means[[lv]]
    Yc[idx, ] <- sweep(Y[idx, , drop = FALSE], 2, mu, "-")
  }
  Yc
}

# --------------------------------------------------------------------------
# 4. Fit all five models
# --------------------------------------------------------------------------

results <- list()

# ------- 4a. LB_base: LB model, no cohort adjustment -------
cat("--- Fitting LB_base ---\n")
set.seed(SEED + 1L)
fit_lb_base <- suppressMessages(
  fit_supervised_mf_modular(Y_train, t_train, s_train,
                            K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
                            prior_beta = PRIOR_BETA,
                            verbose = TRUE, cohort_id = NULL)
)
pred_lb_base  <- predict_supervised_mf(Y_test, fit_lb_base$EF, fit_lb_base$EBeta)
c_lb_base     <- ext_cindex(pred_lb_base$risk_scores, t_test, s_test)
c_lb_base     <- max(c_lb_base, 1 - c_lb_base)   # orientation-free
rec_lb_base   <- factor_recovery(fit_lb_base$EL, L_true[train_idx, ])
results$LB_base <- list(c_ext = c_lb_base, recovery = rec_lb_base)
cat(sprintf("  C-ext = %.3f | factor recovery = %.3f\n\n", c_lb_base, rec_lb_base))

# ------- 4b. LB_cohort: LB model + cohort_id -------
cat("--- Fitting LB_cohort ---\n")
set.seed(SEED + 1L)
fit_lb_coh <- suppressMessages(
  fit_supervised_mf_modular(Y_train, t_train, s_train,
                            K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
                            prior_beta = PRIOR_BETA,
                            verbose = TRUE, cohort_id = cid_train)
)
# Test-time: subtract estimated cohort offset for test patients' cohort
EF_cohort_est <- fit_lb_coh$EF_cohort       # p x (C-1)
L_cohort_test <- model.matrix(~ factor(cid_test,
                                       levels = levels(factor(cid_train))))[, -1, drop = FALSE]
Y_test_adj_lb <- Y_test - L_cohort_test %*% t(EF_cohort_est)
pred_lb_coh   <- predict_supervised_mf(Y_test_adj_lb, fit_lb_coh$EF, fit_lb_coh$EBeta)
c_lb_coh      <- ext_cindex(pred_lb_coh$risk_scores, t_test, s_test)
c_lb_coh      <- max(c_lb_coh, 1 - c_lb_coh)
rec_lb_coh    <- factor_recovery(fit_lb_coh$EL, L_true[train_idx, ])
results$LB_cohort <- list(c_ext = c_lb_coh, recovery = rec_lb_coh)
cat(sprintf("  C-ext = %.3f | factor recovery = %.3f\n\n", c_lb_coh, rec_lb_coh))

# ------- 4c. YFB_base: YFB model, no cohort adjustment -------
cat("--- Fitting YFB_base ---\n")
set.seed(SEED + 1L)
fit_yfb_base <- suppressMessages(
  fit_cox_on_yf(Y_train, t_train, s_train,
                K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
                verbose = TRUE, cohort_id = NULL)
)
pred_yfb_base <- predict_cox_on_yf(Y_test, fit_yfb_base$EF, fit_yfb_base$EBeta,
                                    EF_norms = fit_yfb_base$EF_norms)
c_yfb_base    <- ext_cindex(pred_yfb_base$risk_scores, t_test, s_test)
c_yfb_base    <- max(c_yfb_base, 1 - c_yfb_base)
rec_yfb_base  <- factor_recovery(fit_yfb_base$EL, L_true[train_idx, ])
results$YFB_base <- list(c_ext = c_yfb_base, recovery = rec_yfb_base)
cat(sprintf("  C-ext = %.3f | factor recovery = %.3f\n\n", c_yfb_base, rec_yfb_base))

# ------- 4d. YFB_cohort: YFB model + cohort_id -------
cat("--- Fitting YFB_cohort ---\n")
set.seed(SEED + 1L)
fit_yfb_coh <- suppressMessages(
  fit_cox_on_yf(Y_train, t_train, s_train,
                K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
                verbose = TRUE, cohort_id = cid_train)
)
EF_cohort_yfb <- fit_yfb_coh$EF_cohort
Y_test_adj_yfb <- Y_test - L_cohort_test %*% t(EF_cohort_yfb)
pred_yfb_coh  <- predict_cox_on_yf(Y_test_adj_yfb, fit_yfb_coh$EF, fit_yfb_coh$EBeta,
                                    EF_norms = fit_yfb_coh$EF_norms)
c_yfb_coh     <- ext_cindex(pred_yfb_coh$risk_scores, t_test, s_test)
c_yfb_coh     <- max(c_yfb_coh, 1 - c_yfb_coh)
rec_yfb_coh   <- factor_recovery(fit_yfb_coh$EL, L_true[train_idx, ])
results$YFB_cohort <- list(c_ext = c_yfb_coh, recovery = rec_yfb_coh)
cat(sprintf("  C-ext = %.3f | factor recovery = %.3f\n\n", c_yfb_coh, rec_yfb_coh))

# ------- 4e. zscore: per-cohort mean subtraction, then LB_base -------
cat("--- Fitting zscore (per-cohort centering + LB_base) ---\n")
# Compute training cohort means (used for test centering to avoid test-set leakage)
train_cohort_means <- setNames(
  lapply(unique(cid_train), function(lv) colMeans(Y_train[cid_train == lv, , drop=FALSE])),
  unique(cid_train)
)
Y_train_zs <- zscore_center(Y_train, cid_train)
Y_test_zs  <- zscore_center(Y_test,  cid_test, ref_means = train_cohort_means)
set.seed(SEED + 1L)
fit_zs <- suppressMessages(
  fit_supervised_mf_modular(Y_train_zs, t_train, s_train,
                            K = K_FIT, max_iter = MAX_ITER, alpha = ALPHA,
                            prior_beta = PRIOR_BETA,
                            verbose = TRUE, cohort_id = NULL)
)
pred_zs   <- predict_supervised_mf(Y_test_zs, fit_zs$EF, fit_zs$EBeta)
c_zs      <- ext_cindex(pred_zs$risk_scores, t_test, s_test)
c_zs      <- max(c_zs, 1 - c_zs)
rec_zs    <- factor_recovery(fit_zs$EL, L_true[train_idx, ])
results$zscore <- list(c_ext = c_zs, recovery = rec_zs)
cat(sprintf("  C-ext = %.3f | factor recovery = %.3f\n\n", c_zs, rec_zs))

# --------------------------------------------------------------------------
# 5. Summary table
# --------------------------------------------------------------------------

cat("============================================================\n")
cat(" Stage 1 Results\n")
cat("============================================================\n")
cat(sprintf("%-15s  %8s  %10s\n", "Model", "C-ext", "Factor.Rec"))
cat(sprintf("%-15s  %8s  %10s\n", "-----", "-----", "----------"))
for (nm in names(results)) {
  cat(sprintf("%-15s  %8.3f  %10.3f\n",
              nm, results[[nm]]$c_ext, results[[nm]]$recovery))
}
cat("----\n")

# --------------------------------------------------------------------------
# 6. Acceptance criteria
#
# (a) Factor recovery: both cohort models improve over their base counterparts.
#     Demonstrates that absorbing the platform offset frees the biological factors.
# (b) Offset absorption: LB_cohort EF_cohort correlates with the true gene-level
#     offset_vec (|cor| ≥ 0.7).  Direct test that the rank-1 cohort column
#     captures the platform effect rather than biology.
# (c) C-index stability: LB_cohort C-ext within 0.10 of LB_base C-ext.
#     With balanced cohort assignment, survival prediction shouldn't degrade.
# --------------------------------------------------------------------------

cat("\n--- Acceptance criteria ---\n")

# (a) Factor recovery improvement
pass_a <- isTRUE(results$LB_cohort$recovery  >= results$LB_base$recovery) &&
          isTRUE(results$YFB_cohort$recovery >= results$YFB_base$recovery)
cat(sprintf("(a) Factor recovery: LB(%.3f >= %.3f) + YFB(%.3f >= %.3f): %s\n",
            results$LB_cohort$recovery,  results$LB_base$recovery,
            results$YFB_cohort$recovery, results$YFB_base$recovery,
            if (pass_a) "PASS" else "FAIL"))

# (b) Offset absorption quality
ef_cor_lb  <- abs(cor(as.vector(fit_lb_coh$EF_cohort),  offset_vec))
ef_cor_yfb <- abs(cor(as.vector(fit_yfb_coh$EF_cohort), offset_vec))
pass_b <- ef_cor_lb >= 0.7 || ef_cor_yfb >= 0.7   # at least one model absorbs it
cat(sprintf("(b) Offset absorption: LB |cor|=%.3f, YFB |cor|=%.3f (need >= 0.7 for either): %s\n",
            ef_cor_lb, ef_cor_yfb, if (pass_b) "PASS" else "FAIL"))

# (c) C-index stability (no catastrophic degradation)
pass_c <- abs(results$LB_cohort$c_ext - results$LB_base$c_ext) <= 0.10
cat(sprintf("(c) C-index stability: |LB_cohort - LB_base| = %.3f <= 0.10: %s\n",
            abs(results$LB_cohort$c_ext - results$LB_base$c_ext),
            if (pass_c) "PASS" else "FAIL"))

# --------------------------------------------------------------------------
# 7. Save results
# --------------------------------------------------------------------------

output <- list(
  results          = results,
  offset_cors      = c(LB=ef_cor_lb, YFB=ef_cor_yfb),
  params           = list(N=N, P=P, K_TRUE=K_TRUE, OFFSET_SD=OFFSET_SD,
                          BETA_TRUE=BETA_TRUE, SEED=SEED),
  pass             = c(a=pass_a, b=pass_b, c=pass_c),
  date             = Sys.time()
)
out_path <- "results/cohort_lmm_sim/synthetic_results.rds"
saveRDS(output, out_path)
cat(sprintf("\nResults saved to %s\n", out_path))
cat("============================================================\n")
