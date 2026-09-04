# ============================================================
# Script:  results/benchmark_sim/run_k_init_multistart_dip_check.R
# Purpose: Stage 6 of the 9/4 plan. The K_init=2..20 sweep
#          (run_k_init_sweep.R) shows external C-index dips sharply at
#          K_init=11 and K_init=13 (0.5994, 0.5390) relative to the
#          otherwise flat K_init=7..20 plateau (~0.627). A single SVD init
#          at a given K is not guaranteed to reach that K's true ELBO
#          optimum (CAVI is non-convex), so a dip could be a genuine
#          feature of these K's, or a single-init CAVI factor-collapse
#          artifact -- exactly the failure mode already documented for
#          K=5/6 in run_k_init_multistart_check.R (2026-08-19) and for
#          K=2/K=4 in the K-parsimony follow-up (DECISIONS.md 2026-07-13).
#
#          Falsifiable prediction stated in the 9/4 plan BEFORE running
#          this (§1e): the two dips (K_init=11, K_init=13) coincide with
#          the two worst ELBOs in the single-init sweep (-825096, -834897)
#          AND with K_survival_active jumping from 2 to 3 -- the signature
#          of a worse local optimum. If that is right, multistart should
#          remove both dips (best-of-15-restart ELBO should reach a
#          solution with external C back on the K=7..20 plateau).
#
#          Method: fit_cox_on_yf_multistart() (same as
#          run_k_init_multistart_check.R, which did this for K in
#          {5,6,7,8,9,10}) -- n_init=15 restarts (1 SVD + 14 random inits),
#          best-ELBO selection -- at K_init in 10:15 (the plan's requested
#          range, bracketing both dips with K=10/12/14/15 as non-dip
#          controls already on the plateau). A separate script and separate
#          output files from run_k_init_multistart_check.R, so that
#          script's 2026-08-19 outputs (referenced in DECISIONS.md) are
#          left untouched.
#
#          Same D4 preprocessing as run_k_init_sweep.R (per-platform z-std,
#          combined_rank gene selection, top-3000 per cohort, no cohort_id).
#          K_init values fit in parallel (mclapply); n_init=15 restarts run
#          serially within each K's worker (fit_cox_on_yf_multistart() has
#          no internal parallelism).
#
#   Output: results/benchmark_sim/outputs/k_init_sweep/
#             k_init_multistart_dip_results.csv   (one row per K, best-of-multistart)
#             k_init_multistart_dip_restarts.csv  (one row per K x restart, full diagnostics)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_k_init_multistart_dip_check.R
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_multistart_dip_check.R --quick
#          K_SWEEP_CORES env var controls worker count (default 5), same as run_k_init_sweep.R.
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival); library(parallel) })

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/fit_modular_multistart.R")   # fit_cox_on_yf_multistart()
source("code/preprocess_desurv.R")
source("code/select_K.R")                 # classify_factors()

YML_PATH <- "config/globals.yml"
cfg      <- yaml::read_yaml(YML_PATH)
b        <- cfg$benchmark
pcfg     <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
N_INIT       <- if (QUICK_MODE) 3L  else 15L
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
PVE_THRESH   <- cfg$k_selection$pve_threshold
TOP_N_DESURV <- pcfg$top_n_genes_desurv

K_VALUES <- 10:15
CORES    <- max(1L, as.integer(Sys.getenv("K_SWEEP_CORES", "5")))

OUT_DIR <- "results/benchmark_sim/outputs/k_init_sweep"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | N_INIT: %d | Cores: %d | K values: %s\n",
            QUICK_MODE, MAX_ITER, N_INIT, CORES, paste(K_VALUES, collapse = ", ")))

# --------------------------------------------------------------------------
# 1. Load + preprocess training data (identical to run_k_init_sweep.R)
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

cat("--- Preprocessing training data (D4 config) ---\n")
pp <- preprocess_merged_cohorts(
  cohort_raw_list          = train_raw,
  log_transform_flags      = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n                    = TOP_N_DESURV,
  rank_transform           = FALSE,
  per_platform_standardize = TRUE,
  normalize_method         = "none",
  selection_per_cohort     = TRUE,
  selection_method         = "combined_rank"
)
Y_train     <- pp$Y
train_genes <- pp$gene_names
cat(sprintf("  n=%d, p=%d\n\n", nrow(Y_train), ncol(Y_train)))

# --------------------------------------------------------------------------
# 2. Load external cohorts ONCE before forking (same rationale as
#    run_k_init_sweep.R: load_pdac_raw()'s tempfile/symlink dance is not
#    fork-safe if done concurrently in multiple mclapply workers).
# --------------------------------------------------------------------------

cat("--- Loading external cohorts (5) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

ext_data <- list()
for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE
  )
  common <- intersect(train_genes, pre_ext$gene_names)
  if (length(common) < 100) {
    cat(sprintf("    Skipping %s: only %d common genes\n", ext_cohort, length(common)))
    next
  }
  ext_data[[ext_cohort]] <- list(
    Y_ext     = pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE],
    train_idx = match(common, train_genes),
    time      = raw_ext$time,
    status    = raw_ext$status
  )
}
cat(sprintf("  %d/%d external cohorts usable\n\n", length(ext_data), length(EXTERNAL_COHORTS)))

oriented_cindex <- function(risk, time, status) {
  if (sd(risk) == 0) return(NA_real_)
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 3. Per-K worker: multistart fit, external scoring, factor classification.
# --------------------------------------------------------------------------

run_one_K <- function(K) {
  tryCatch({
    set.seed(42L)
    ms <- fit_cox_on_yf_multistart(
      Y_train, time_train, status_train,
      K = K, max_iter = MAX_ITER, alpha = ALPHA, prior_beta = PRIOR_BETA,
      n_init = N_INIT, init_seed_base = 42L
    )
    fit <- ms$best
    ms$restarts$K <- K

    cohort_c <- list()
    for (ext_cohort in names(ext_data)) {
      d      <- ext_data[[ext_cohort]]
      EF_sub <- fit$EF[d$train_idx, , drop = FALSE]
      pred   <- predict_cox_on_yf(d$Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
      cohort_c[[ext_cohort]] <- oriented_cindex(pred$risk_scores, d$time, d$status)
    }
    mean_c <- if (length(cohort_c) > 0) mean(unlist(cohort_c), na.rm = TRUE) else NA_real_

    cls <- classify_factors(fit, Y_train, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
    K_survival_active <- sum(cls$category == "survival_active")
    K_genomics_only   <- sum(cls$category == "genomics_only")
    K_dead            <- sum(cls$category == "dead")
    K_eff_total       <- K_survival_active + K_genomics_only

    row <- data.frame(
      K                 = K,
      best_idx          = ms$best_idx,
      n_init            = N_INIT,
      best_elbo_full    = round(tail(fit$history$elbo_full, 1), 4),
      K_survival_active = K_survival_active,
      K_genomics_only   = K_genomics_only,
      K_dead            = K_dead,
      K_eff_total       = K_eff_total,
      mean_external_c   = round(mean_c, 4),
      stringsAsFactors  = FALSE
    )
    for (ext_cohort in EXTERNAL_COHORTS) {
      v <- cohort_c[[ext_cohort]]
      row[[paste0("c_", ext_cohort)]] <- if (is.null(v)) NA_real_ else round(v, 4)
    }

    list(K = K, status = "ok", row = row, restarts = ms$restarts, fit = fit)
  }, error = function(e) list(K = K, status = "error", error_msg = conditionMessage(e)))
}

cat(sprintf("=== Multistart fitting (n_init=%d, best-ELBO) at K_init in {%s} ===\n\n",
            N_INIT, paste(K_VALUES, collapse = ", ")))

raw_results <- if (CORES > 1) {
  mclapply(K_VALUES, run_one_K, mc.cores = min(CORES, length(K_VALUES)), mc.preschedule = FALSE)
} else {
  lapply(K_VALUES, run_one_K)
}

# --------------------------------------------------------------------------
# 4. Collect, report, save.
# --------------------------------------------------------------------------

results_rows <- list(); restart_rows <- list(); best_fits <- list()
for (i in seq_along(K_VALUES)) {
  K   <- K_VALUES[i]
  res <- raw_results[[i]]
  if (is.null(res) || inherits(res, "try-error") || !is.list(res) || is.null(res$status) ||
      identical(res$status, "error")) {
    msg <- if (is.list(res) && !is.null(res$error_msg)) res$error_msg else "worker failed/killed"
    cat(sprintf("K=%2d: FAILED — %s\n", K, msg))
    next
  }
  results_rows[[length(results_rows) + 1]] <- res$row
  restart_rows[[length(restart_rows) + 1]] <- res$restarts
  best_fits[[as.character(K)]] <- res$fit
  r <- res$row
  cat(sprintf("K=%2d: best_idx=%d/%d | best_elbo=%.4f | K_survival_active=%d, K_genomics_only=%d, K_dead=%d, K_eff_total=%d | mean external C=%.4f\n",
              K, r$best_idx, r$n_init, r$best_elbo_full,
              r$K_survival_active, r$K_genomics_only, r$K_dead, r$K_eff_total, r$mean_external_c))
}

results      <- do.call(rbind, results_rows)
restarts_all <- do.call(rbind, restart_rows)

out_csv          <- file.path(OUT_DIR, "k_init_multistart_dip_results.csv")
out_restarts_csv <- file.path(OUT_DIR, "k_init_multistart_dip_restarts.csv")
write.csv(results, out_csv, row.names = FALSE)
write.csv(restarts_all, out_restarts_csv, row.names = FALSE)
saveRDS(best_fits, file.path(OUT_DIR, "k_init_multistart_dip_best_fits.rds"))

cat(sprintf("\n=== Results saved: %s ===\n", out_csv))
cat(sprintf("=== Restarts saved: %s ===\n", out_restarts_csv))

# --------------------------------------------------------------------------
# 5. Falsifiable-prediction check (stated in the 9/4 plan before running):
#    does multistart remove the K_init=11/13 dips seen in the single-init
#    sweep (run_k_init_sweep.R)?
# --------------------------------------------------------------------------

sweep_csv <- file.path(OUT_DIR, "k_init_sweep_results.csv")
if (file.exists(sweep_csv)) {
  sweep <- read.csv(sweep_csv, stringsAsFactors = FALSE)
  cat("\n=== Single-init sweep vs. best-of-multistart, K_init in {10..15} ===\n")
  for (K in K_VALUES) {
    si <- sweep[sweep$K_init == K, ]
    ms <- results[results$K == K, ]
    if (nrow(si) == 0 || nrow(ms) == 0) next
    cat(sprintf("K=%2d: single-init C=%.4f (K_surv_active=%d) -> best-of-%d-multistart C=%.4f (K_surv_active=%d)\n",
                K, si$mean_external_c[1], si$K_survival_active[1],
                N_INIT, ms$mean_external_c[1], ms$K_survival_active[1]))
  }
}
