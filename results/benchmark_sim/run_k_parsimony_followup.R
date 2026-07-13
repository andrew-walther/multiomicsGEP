# ============================================================
# Script:  results/benchmark_sim/run_k_parsimony_followup.R
# Purpose: Step 1 of the K-parsimony follow-up plan
#          (docs/plans/ssbmf_k_parsimony_followup_plan_07_13_2026.md).
#          Phase 3 (DECISIONS.md 2026-07-13) found K=7 is not free to shrink
#          in a single-seed comparison, but flagged K=2/K=4's suspiciously
#          fast convergence to near-zero beta as consistent with the CAVI
#          factor-collapse artifact documented in Phase 2. This script
#          re-checks K in {2,3,4,5} with two improved-optimization strategies
#          -- warm-start (seeded from the converged K=7 fit's top-K PVE-ranked
#          columns) and best-ELBO multistart -- alongside the original
#          fresh-SVD single-init fit, for a direct three-way comparison.
#
#          Uses the same D4 config as run_k_parsimony_curve.R: per-platform
#          z-std, DeSurv combined-rank gene selection, top-3000 per cohort,
#          no cohort_id.
#
#          Decision rule (mechanical, plan-specified): for each K < 7,
#          compute the best external mean C-index across {fresh, warmstart,
#          multistart}. If ANY K < 7 reaches within 1 SE of K=7's reference
#          (0.6267, margin 0.6068) -> OPTIMIZATION-LIMITED (-> Step 2). If
#          NO K < 7 reaches it -> CAPACITY-LIMITED (-> Step 3).
#
#          Output: results/benchmark_sim/outputs/k_parsimony_followup/
#            k_parsimony_followup_results.csv   (one row per K x method x cohort)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_k_parsimony_followup.R
#          caffeinate -i Rscript results/benchmark_sim/run_k_parsimony_followup.R --quick
# Requires: PDAC_DATA_ROOT set (real data not in git) -- see CLAUDE.md /
#           tests/test_real_data_loading.R for the local/Longleaf paths.
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"), error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/fit_modular_multistart.R")   # fit_cox_on_yf_multistart()
source("code/warmstart_from_fit.R")       # extract_top_k_by_pve()
source("code/preprocess_desurv.R")

b <- cfg$benchmark
p <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L  else cfg$cavi$max_iter
N_INIT       <- if (QUICK_MODE) 3L   else 15L
K_GRID_SMALL <- if (QUICK_MODE) c(2L, 5L) else c(2L, 3L, 4L, 5L)
K_REF        <- 7L
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
TOP_N_DESURV <- p$top_n_genes_desurv

# Reference numbers from the Phase 3 fresh-SVD run (DECISIONS.md 2026-07-13),
# used for the plan's mechanical decision rule. Refit here too (below) so the
# comparison is self-contained and reproducible from this one script.
K7_REF_MEAN_C <- 0.6267
K7_REF_SE_C   <- 0.0199
K7_REF_MARGIN <- K7_REF_MEAN_C - K7_REF_SE_C

OUT_DIR <- "results/benchmark_sim/outputs/k_parsimony_followup"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d | N_INIT: %d | K grid: %s (+ K=%d reference)\n",
            QUICK_MODE, MAX_ITER, N_INIT, paste(K_GRID_SMALL, collapse = ","), K_REF))

# --------------------------------------------------------------------------
# 1. Load and preprocess training data (D4 config, identical to
#    run_k_parsimony_curve.R / run_desurv_comparison.R)
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d, events=%d\n\n", length(time_train), sum(status_train)))

cat("--- Preprocessing (D4 config) ---\n")
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

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. K=7 reference fit (fresh SVD, same settings as run_k_parsimony_curve.R)
#    -- provides the PVE-ranked columns for warm-starting K<7 fits.
# --------------------------------------------------------------------------

cat("=== Fitting K=7 reference (fresh SVD) ===\n\n")
set.seed(42L)
fit_k7 <- suppressMessages(fit_cox_on_yf(
  Y_train, time_train, status_train,
  K = K_REF, max_iter = MAX_ITER, alpha = ALPHA,
  prior_beta = PRIOR_BETA, verbose = TRUE))
cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
            sum(abs(fit_k7$EBeta) > BETA_THRESH), max(abs(fit_k7$EBeta)),
            fit_k7$history$n_iter))

# --------------------------------------------------------------------------
# 3. For each K < 7: fresh-SVD, warm-start (from K=7), and multistart fits
# --------------------------------------------------------------------------

fits <- list("7" = list(fresh = fit_k7))

for (K in K_GRID_SMALL) {
  cat(sprintf("=== K=%d ===\n", K))

  cat("  -- fresh SVD --\n")
  set.seed(42L)
  fit_fresh <- suppressMessages(fit_cox_on_yf(
    Y_train, time_train, status_train,
    K = K, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA, verbose = FALSE))
  cat(sprintf("     K_eff=%d | beta_max=%.4f | iters=%d\n",
              sum(abs(fit_fresh$EBeta) > BETA_THRESH), max(abs(fit_fresh$EBeta)),
              fit_fresh$history$n_iter))

  cat("  -- warm-start (from K=7, PVE-ranked) --\n")
  ws <- extract_top_k_by_pve(fit_k7, K)
  set.seed(42L)
  fit_warmstart <- suppressMessages(fit_cox_on_yf(
    Y_train, time_train, status_train,
    K = K, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA,
    init_method = "custom", EL_init = ws$EL_init, EF_init = ws$EF_init,
    verbose = FALSE))
  cat(sprintf("     K_eff=%d | beta_max=%.4f | iters=%d\n",
              sum(abs(fit_warmstart$EBeta) > BETA_THRESH), max(abs(fit_warmstart$EBeta)),
              fit_warmstart$history$n_iter))

  cat(sprintf("  -- multistart (n_init=%d, best ELBO) --\n", N_INIT))
  ms <- fit_cox_on_yf_multistart(
    Y_train, time_train, status_train,
    K = K, max_iter = MAX_ITER, alpha = ALPHA,
    prior_beta = PRIOR_BETA,
    n_init = N_INIT, init_seed_base = 42L)
  fit_multistart <- ms$best
  cat(sprintf("     best_idx=%d/%d | K_eff=%d | beta_max=%.4f | iters=%d\n",
              ms$best_idx, N_INIT,
              sum(abs(fit_multistart$EBeta) > BETA_THRESH), max(abs(fit_multistart$EBeta)),
              fit_multistart$history$n_iter))

  fits[[as.character(K)]] <- list(
    fresh = fit_fresh, warmstart = fit_warmstart, multistart = fit_multistart
  )
  cat("\n")
}

# --------------------------------------------------------------------------
# 4. External validation on the 5 held-out cohorts, for every K x method
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) x K x method ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows     <- list()

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
  Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
  train_idx <- match(common, train_genes)

  for (K_key in names(fits)) {
    for (method in names(fits[[K_key]])) {
      fit    <- fits[[K_key]][[method]]
      EF_sub <- fit$EF[train_idx, , drop = FALSE]
      pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
      c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

      results_rows[[length(results_rows) + 1]] <- data.frame(
        K              = as.integer(K_key),
        method         = method,
        cohort         = ext_cohort,
        c_index        = round(c_val, 4),
        k_eff          = sum(abs(fit$EBeta) > BETA_THRESH),
        beta_max       = round(max(abs(fit$EBeta)), 4),
        n_iter         = fit$history$n_iter,
        n_common_genes = length(common),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(results_rows) == 0)
  stop("No external validation rows: all cohorts had < 100 common genes.")

results <- do.call(rbind, results_rows)

out_csv <- file.path(OUT_DIR, "k_parsimony_followup_results.csv")
write.csv(results, out_csv, row.names = FALSE)

# --------------------------------------------------------------------------
# 5. Aggregate and apply the plan's mechanical decision rule
# --------------------------------------------------------------------------

agg <- aggregate(c_index ~ K + method, data = results,
                 FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
agg <- do.call(data.frame, agg)
names(agg)[3:4] <- c("mean_c", "se_c")
agg$k_eff <- mapply(function(k, m) unique(results$k_eff[results$K == k & results$method == m]),
                     agg$K, agg$method)
agg <- agg[order(agg$K, agg$method), ]

cat("\n============================================================\n")
cat(" K-parsimony follow-up: K x method vs. mean external C-index\n")
cat("============================================================\n")
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  K=%d  %-11s  mean C=%.4f (SE=%.4f)  K_eff=%d\n",
              agg$K[i], agg$method[i], agg$mean_c[i], agg$se_c[i], agg$k_eff[i]))
}

cat(sprintf("\nK=7 reference: mean C=%.4f, SE=%.4f, 1-SE margin=%.4f (plan-specified, DECISIONS.md 2026-07-13)\n",
            K7_REF_MEAN_C, K7_REF_SE_C, K7_REF_MARGIN))
k7_row <- agg[agg$K == K_REF & agg$method == "fresh", ]
if (nrow(k7_row) == 1) {
  cat(sprintf("K=7 refit in this script: mean C=%.4f, SE=%.4f (reproducibility check)\n",
              k7_row$mean_c, k7_row$se_c))
  if (abs(k7_row$mean_c - K7_REF_MEAN_C) > 1e-3) {
    warning(sprintf(
      "K=7 refit (%.4f) diverges from the hardcoded reference (%.4f) by more than 1e-3 -- the
decision rule below is being applied against a STALE reference margin. Update K7_REF_MEAN_C/
K7_REF_SE_C at the top of this script before trusting the OUTCOME line.",
      k7_row$mean_c, K7_REF_MEAN_C))
  }
}

cat("\n--- Decision rule: best-of-{fresh,warmstart,multistart} per K vs. K=7 margin ---\n")
best_per_k <- aggregate(mean_c ~ K, data = agg[agg$K != K_REF, ], FUN = max)
names(best_per_k)[2] <- "best_mean_c"
best_per_k <- best_per_k[order(best_per_k$K), ]
for (i in seq_len(nrow(best_per_k))) {
  reaches <- best_per_k$best_mean_c[i] >= K7_REF_MARGIN
  cat(sprintf("  K=%d  best mean C=%.4f  %s margin (%.4f)\n",
              best_per_k$K[i], best_per_k$best_mean_c[i],
              if (reaches) "REACHES" else "does not reach", K7_REF_MARGIN))
}

any_reaches <- any(best_per_k$best_mean_c >= K7_REF_MARGIN)
cat(sprintf("\nOUTCOME: %s\n",
            if (any_reaches) "OPTIMIZATION-LIMITED -> proceed to Step 2 (deflation-init fix)"
            else "CAPACITY-LIMITED -> proceed to Step 3 (joint K/alpha/penalty Bayesian optimization)"))

cat(sprintf("\nResults: %s\n", out_csv))
cat("============================================================\n")
