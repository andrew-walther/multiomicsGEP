# ============================================================
# Script:  results/benchmark_sim/run_k_init_sweep.R
# Purpose: Analysis A (docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md).
#          Fits YFB on the real TCGA+CPTAC training data across K_init =
#          2..20, and classifies factors at each fit with classify_factors(),
#          to check whether K_eff_survival stays stable regardless of the
#          starting K — i.e. whether ARD pruning from an over-specified K is
#          a stable alternative to CV-selecting K directly from held-out
#          C-index.
#
#          Also records each fit's full ELBO (elbo_full, alpha-weighted
#          genomics + survival + KL), a joint genomics+survival
#          log-likelihood and BIC (code/compute_bic.R -- an ELBO-style
#          bound, not an exact marginal likelihood; df = K_init*(n+p+1),
#          i.e. charged by the CAVI starting K, not by ARD's K_eff), and
#          final RMSE, so K_init can be selected by a *consensus* of ELBO /
#          BIC / log-likelihood / external C-index rather than any single
#          criterion. These are not guaranteed to agree -- report whatever
#          the sweep actually shows.
#
#          As of 2026-09-04, also records two genuinely held-out criteria
#          (code/compute_cv_loglik.R), answering directly whether the
#          in-sample log-likelihood/BIC above agree with a real held-out
#          check: cv_survival_loglik() (held-out Cox partial log-likelihood,
#          leakage-free) and bicv_genomics_loglik() (bi-cross-validated
#          genomics log-likelihood). These require fresh per-fold/per-block
#          fits that cannot be derived from the single cached main fit, so
#          they have their own cache (k_init_sweep_cv_results.rds,
#          --reuse-cv-cache) independent of --reuse-cache for the main fits.
#
#          Preprocessing matches the D4 configuration in
#          run_desurv_comparison.R exactly: YFB + per-platform z-std +
#          combined_rank gene selection (top-3000 per cohort, before
#          normalization) + no cohort_id.
#
#   Output: results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv
#           results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_fits.rds
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-19
# Updated: 2026-08-27 -- extended K_INIT_VALUES to 2:20, added BIC/log-
#          likelihood columns, parallelized the fit loop with mclapply, and
#          added a --reuse-cache back-fill path.
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R --quick
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R --serial
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R --reuse-cache
#          caffeinate -i Rscript results/benchmark_sim/run_k_init_sweep.R --k-init=5,6,7,8,9,10,15,20
#
#          Parallelism: K_SWEEP_CORES env var controls the worker count
#          (default 5), e.g. `K_SWEEP_CORES=8 Rscript ...`. This R build
#          links vecLib, which multithreads BLAS/LAPACK via Grand Central
#          Dispatch and ignores OMP_NUM_THREADS; the only lever on
#          intra-process BLAS threading is VECLIB_MAXIMUM_THREADS, which
#          must be set BEFORE the R process starts (it cannot be changed at
#          runtime), e.g.:
#            VECLIB_MAXIMUM_THREADS=1 K_SWEEP_CORES=5 Rscript results/benchmark_sim/run_k_init_sweep.R
#          Leaving VECLIB_MAXIMUM_THREADS unset lets each of the K_SWEEP_CORES
#          forked workers multithread its own BLAS calls, which can
#          oversubscribe the machine; pin it to 1 when running the full
#          sweep with mc.cores > 1.
# ============================================================

args           <- commandArgs(trailingOnly = TRUE)
QUICK_MODE     <- "--quick" %in% args
REUSE_CACHE    <- "--reuse-cache" %in% args
REUSE_CV_CACHE <- "--reuse-cv-cache" %in% args
SERIAL_MODE    <- "--serial" %in% args
k_init_arg     <- args[grepl("^--k-init=", args)]

# Navigate to project root if invoked from a subdirectory
if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival); library(parallel) })

cfg <- yaml::read_yaml("config/globals.yml")  # benchmark_helpers.R must be sourced after cfg exists

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/compute_bic.R")   # compute_joint_ll_bic -- needs calc_cox_taylor_yf from fit_cox_on_yf.R
source("code/compute_cv_loglik.R")  # cv_survival_loglik, bicv_genomics_loglik
source("code/preprocess_desurv.R")
source("code/select_K.R")

YML_PATH <- "config/globals.yml"
cfg      <- yaml::read_yaml(YML_PATH)
b        <- cfg$benchmark
p        <- cfg$preprocessing

ALPHA        <- b$alpha
MAX_ITER     <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA   <- "normal"
BETA_THRESH  <- cfg$k_selection$beta_threshold
PVE_THRESH   <- cfg$k_selection$pve_threshold
TOP_N_DESURV <- p$top_n_genes_desurv  # 3000 — DeSurv-aligned, D4 config

OUT_DIR      <- "results/benchmark_sim/outputs/k_init_sweep"
FITS_RDS     <- file.path(OUT_DIR, "k_init_sweep_fits.rds")
CV_CACHE_RDS <- file.path(OUT_DIR, "k_init_sweep_cv_results.rds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CVL <- cfg$cv_loglik

# --------------------------------------------------------------------------
# 0. Resolve K_INIT_VALUES and (if requested) load the cached fits to reuse.
# --------------------------------------------------------------------------

cached_fits <- NULL
if (REUSE_CACHE) {
  if (!file.exists(FITS_RDS)) stop("--reuse-cache given but no cache found at ", FITS_RDS)
  cached_fits <- readRDS(FITS_RDS)
  cat(sprintf("--reuse-cache: loaded %d cached fits from %s (K = %s)\n",
              length(cached_fits), FITS_RDS,
              paste(sort(as.integer(names(cached_fits))), collapse = ", ")))
}

# Separate cache for the held-out CV criteria (code/compute_cv_loglik.R):
# these need fresh per-fold/per-block fits that the single main fit above
# cannot supply, so they cannot be back-filled from cached_fits.
cached_cv <- NULL
if (REUSE_CV_CACHE) {
  if (!file.exists(CV_CACHE_RDS)) stop("--reuse-cv-cache given but no cache found at ", CV_CACHE_RDS)
  cached_cv <- readRDS(CV_CACHE_RDS)
  cat(sprintf("--reuse-cv-cache: loaded %d cached CV results from %s (K = %s)\n",
              length(cached_cv), CV_CACHE_RDS,
              paste(sort(as.integer(names(cached_cv))), collapse = ", ")))
}

if (length(k_init_arg) > 0) {
  K_INIT_VALUES <- sort(as.integer(strsplit(sub("^--k-init=", "", k_init_arg[1]), ",")[[1]]))
} else if (REUSE_CACHE) {
  # Default reuse-cache scope: exactly the K values already cached (the
  # back-fill verification gate is defined on this set).
  K_INIT_VALUES <- sort(as.integer(names(cached_fits)))
} else {
  K_INIT_VALUES <- 2:20
}

CORES <- if (SERIAL_MODE) 1L else max(1L, as.integer(Sys.getenv("K_SWEEP_CORES", "5")))

cat(sprintf("Quick mode: %s | Reuse cache: %s | Serial: %s | Cores: %d | MAX_ITER: %d | K_init values: %s\n",
            QUICK_MODE, REUSE_CACHE, SERIAL_MODE, CORES, MAX_ITER,
            paste(K_INIT_VALUES, collapse = ", ")))

# --------------------------------------------------------------------------
# 1. Load + preprocess training data (TCGA_PAAD + CPTAC) — D4 config.
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
# 2. Load + preprocess the 5 external cohorts ONCE, in the parent, BEFORE
#    any forking. load_pdac_raw() does tempfile() -> dir.create() ->
#    file.symlink(); forked children inherit the parent RNG state, so
#    concurrent calls to it from mclapply workers can collide on the
#    symlink path. Loading once here and sharing read-only via fork
#    copy-on-write avoids that entirely. Per-K external scoring below stays
#    serial inside each worker — it's cheap (no I/O), just matrix algebra
#    over 5 already-loaded cohorts.
# --------------------------------------------------------------------------

cat("--- Loading external cohorts (5) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

ext_data <- list()
for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  # top_n=NULL: keep all external genes; intersection with train_genes
  # controls the final gene set — matches run_desurv_comparison.R's D4
  # external eval.
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
  if (sd(risk) == 0) return(NA_real_)  # constant risk score: undefined orientation
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 3. Per-K_init worker: fit-or-reuse, BIC/LL, classify_factors, external
#    scoring. tryCatch wraps the whole body so one K's failure (e.g. a
#    degenerate fit that errors instead of just converging badly) doesn't
#    abort the run — it emits a failed row instead.
# --------------------------------------------------------------------------

run_one_K <- function(K_init, verbose_fit) {
  tryCatch({
    t0 <- proc.time()[["elapsed"]]

    cached <- cached_fits[[as.character(K_init)]]
    if (!is.null(cached)) {
      fit        <- cached
      fit_source <- "cached"
    } else {
      set.seed(42L)  # overrides the forked child's inherited RNG stream
      fit <- fit_cox_on_yf(
        Y_train, time_train, status_train,
        K = K_init, max_iter = MAX_ITER, alpha = ALPHA,
        prior_beta = PRIOR_BETA, verbose = verbose_fit
      )
      fit_source <- "fresh"
    }
    fit_secs <- proc.time()[["elapsed"]] - t0

    n_iter <- fit$history$n_iter
    beta_max <- max(abs(fit$EBeta))

    bic_res <- compute_joint_ll_bic(fit, Y_train, time_train, status_train)

    cv_cached <- cached_cv[[as.character(K_init)]]
    if (!is.null(cv_cached)) {
      cv_surv <- cv_cached$cv_surv
      cv_bicv <- cv_cached$cv_bicv
      cv_source <- "cached"
    } else {
      set.seed(42L)
      cv_surv <- tryCatch(
        cv_survival_loglik(Y_train, time_train, status_train, K = K_init,
                            n_folds = CVL$n_folds, seed = CVL$seed,
                            max_iter = MAX_ITER, prior_beta = PRIOR_BETA, alpha = ALPHA),
        error = function(e) { cat(sprintf("  [K=%d] cv_survival_loglik failed: %s\n", K_init, conditionMessage(e))); NULL }
      )
      cv_bicv <- tryCatch(
        bicv_genomics_loglik(Y_train, status_train, K = K_init,
                              n_row_folds = CVL$n_row_folds, n_col_folds = CVL$n_col_folds,
                              seed = CVL$seed, max_iter = MAX_ITER),
        error = function(e) { cat(sprintf("  [K=%d] bicv_genomics_loglik failed: %s\n", K_init, conditionMessage(e))); NULL }
      )
      cv_source <- "fresh"
    }

    cls <- classify_factors(fit, Y_train, beta_thresh = BETA_THRESH, pve_thresh = PVE_THRESH)
    K_survival_active <- sum(cls$category == "survival_active")
    K_genomics_only   <- sum(cls$category == "genomics_only")
    K_dead            <- sum(cls$category == "dead")
    K_eff_total       <- K_survival_active + K_genomics_only

    cohort_c <- list()
    for (ext_cohort in names(ext_data)) {
      d      <- ext_data[[ext_cohort]]
      EF_sub <- fit$EF[d$train_idx, , drop = FALSE]
      pred   <- predict_cox_on_yf(d$Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
      cohort_c[[ext_cohort]] <- oriented_cindex(pred$risk_scores, d$time, d$status)
    }
    mean_c <- if (length(cohort_c) > 0) mean(unlist(cohort_c), na.rm = TRUE) else NA_real_

    row <- data.frame(
      K_init             = K_init,
      K_total             = K_init,
      K_survival_active  = K_survival_active,
      K_genomics_only    = K_genomics_only,
      K_dead             = K_dead,
      K_eff_total        = K_eff_total,
      elbo_full          = round(fit$history$elbo_full[n_iter], 4),
      loglik_genomics    = round(bic_res$loglik_genomics, 4),
      loglik_survival    = round(bic_res$loglik_survival, 4),
      loglik_joint       = round(bic_res$loglik_joint, 4),
      df                 = bic_res$df,
      bic                = round(bic_res$bic, 4),
      rmse               = round(fit$history$rmse[n_iter], 6),
      n_iter             = n_iter,
      converged          = fit$history$converged,
      beta_max           = round(beta_max, 6),
      fit_status         = "ok",
      fit_source         = fit_source,
      fit_secs           = round(fit_secs, 2),
      n_train            = nrow(Y_train),
      p_genes            = ncol(Y_train),
      mean_external_c    = round(mean_c, 4),
      cv_survival_total_logPL     = if (is.null(cv_surv)) NA_real_ else round(cv_surv$total_logPL, 4),
      cv_survival_logPL_per_event = if (is.null(cv_surv)) NA_real_ else round(cv_surv$mean_logPL_per_event, 6),
      cv_survival_se_logPL        = if (is.null(cv_surv)) NA_real_ else round(cv_surv$se_logPL, 4),
      bicv_genomics_total_loglik  = if (is.null(cv_bicv)) NA_real_ else round(cv_bicv$total_loglik, 2),
      cv_source                   = cv_source,
      stringsAsFactors   = FALSE
    )
    for (ext_cohort in EXTERNAL_COHORTS) {
      v <- cohort_c[[ext_cohort]]
      row[[paste0("c_", ext_cohort)]] <- if (is.null(v)) NA_real_ else round(v, 4)
    }

    list(K_init = K_init, fit_status = "ok", row = row, fit = fit,
         cv_result = list(cv_surv = cv_surv, cv_bicv = cv_bicv))
  }, error = function(e) {
    list(K_init = K_init, fit_status = "error", error_msg = conditionMessage(e))
  })
}

# --------------------------------------------------------------------------
# 4. Dispatch. verbose=FALSE in parallel mode (forked stdout interleaves);
#    a per-K summary is printed in the parent after collection instead.
#    mc.preschedule=FALSE: K=20 costs far more than K=2, so round-robin
#    scheduling would strand large-K jobs on one worker; it also isolates a
#    C-level crash in one fit to one K's slot instead of taking down a
#    whole pre-scheduled chunk.
# --------------------------------------------------------------------------

cat(sprintf("=== Fitting/scoring YFB at %d K_init value(s) (%s) ===\n\n",
            length(K_INIT_VALUES), if (CORES > 1) "parallel" else "serial"))

if (CORES > 1) {
  raw_results <- mclapply(K_INIT_VALUES, run_one_K, verbose_fit = FALSE,
                           mc.cores = CORES, mc.preschedule = FALSE)
} else {
  raw_results <- lapply(K_INIT_VALUES, run_one_K, verbose_fit = TRUE)
}

# --------------------------------------------------------------------------
# 5. Post-collection sweep: a killed forked worker returns NULL (not an R
#    condition), and mclapply can also hand back a "try-error" object for a
#    worker that errored outside our own tryCatch. Neither can be caught
#    from inside run_one_K, so handle them here.
# --------------------------------------------------------------------------

fits          <- list()
cv_results    <- list()
results_rows  <- list()
n_ok <- 0L; n_failed <- 0L

for (i in seq_along(K_INIT_VALUES)) {
  K_init <- K_INIT_VALUES[i]
  res    <- raw_results[[i]]

  if (is.null(res) || inherits(res, "try-error") || !is.list(res) || is.null(res$fit_status)) {
    n_failed <- n_failed + 1L
    msg <- if (inherits(res, "try-error")) conditionMessage(attr(res, "condition"))
           else if (is.null(res)) "worker returned NULL (killed or crashed)"
           else "malformed worker result"
    cat(sprintf("K_init=%2d: FAILED — %s\n", K_init, msg))
    results_rows[[length(results_rows) + 1]] <- data.frame(
      K_init = K_init, fit_status = "worker_error", stringsAsFactors = FALSE
    )
    next
  }

  if (identical(res$fit_status, "error")) {
    n_failed <- n_failed + 1L
    cat(sprintf("K_init=%2d: FAILED — %s\n", K_init, res$error_msg))
    results_rows[[length(results_rows) + 1]] <- data.frame(
      K_init = K_init, fit_status = "error", stringsAsFactors = FALSE
    )
    next
  }

  n_ok <- n_ok + 1L
  fits[[as.character(K_init)]] <- res$fit
  cv_results[[as.character(K_init)]] <- res$cv_result
  results_rows[[length(results_rows) + 1]] <- res$row
  r <- res$row
  cat(sprintf("K_init=%2d [%s]: K_survival_active=%d, K_genomics_only=%d, K_dead=%d, K_eff_total=%d | mean external C=%.4f | elbo_full=%.4f | bic=%.4f | cv_surv_logPL/event=%s | bicv_genomics=%s | rmse=%.6f | %.1fs\n",
              K_init, r$fit_source, r$K_survival_active, r$K_genomics_only, r$K_dead, r$K_eff_total,
              r$mean_external_c, r$elbo_full, r$bic,
              ifelse(is.na(r$cv_survival_logPL_per_event), "NA", sprintf("%.4f", r$cv_survival_logPL_per_event)),
              ifelse(is.na(r$bicv_genomics_total_loglik), "NA", sprintf("%.1f", r$bicv_genomics_total_loglik)),
              r$rmse, r$fit_secs))
}

cat(sprintf("\n%d/%d K_init values fit successfully (%d failed)\n",
            n_ok, length(K_INIT_VALUES), n_failed))

# Row-bind with fill: failed rows have only K_init/fit_status, so align
# columns via a full union before rbind.
all_cols <- unique(unlist(lapply(results_rows, names)))
results_rows <- lapply(results_rows, function(r) {
  missing <- setdiff(all_cols, names(r))
  for (m in missing) r[[m]] <- NA
  r[all_cols]
})
results <- do.call(rbind, results_rows)
results <- results[order(results$K_init), ]

# --------------------------------------------------------------------------
# 6. Console summary: surface any disagreement between criteria rather than
#    burying it. ELBO/BIC/log-likelihood only compare fits that actually
#    converged, on this training set only; external C-index is the only
#    generalization signal here.
# --------------------------------------------------------------------------

ok_results <- results[results$fit_status == "ok" & !is.na(results$elbo_full), ]
if (nrow(ok_results) > 0) {
  elbo_best    <- ok_results[which.max(ok_results$elbo_full), ]
  bic_best     <- ok_results[which.min(ok_results$bic), ]
  c_best       <- ok_results[which.max(ok_results$mean_external_c), ]
  cv_surv_ok   <- ok_results[!is.na(ok_results$cv_survival_logPL_per_event), ]
  bicv_ok      <- ok_results[!is.na(ok_results$bicv_genomics_total_loglik), ]
  cat("\n=== K_init preferred by each criterion ===\n")
  cat(sprintf("  ELBO-preferred            K_init = %2d (elbo_full=%.4f, K_eff_total=%d)\n",
              elbo_best$K_init, elbo_best$elbo_full, elbo_best$K_eff_total))
  cat(sprintf("  BIC-preferred (in-sample) K_init = %2d (bic=%.4f, K_eff_total=%d)\n",
              bic_best$K_init, bic_best$bic, bic_best$K_eff_total))
  cat(sprintf("  External C-preferred      K_init = %2d (mean_external_c=%.4f, K_eff_total=%d)\n",
              c_best$K_init, c_best$mean_external_c, c_best$K_eff_total))
  if (nrow(cv_surv_ok) > 0) {
    cv_surv_best <- cv_surv_ok[which.max(cv_surv_ok$cv_survival_logPL_per_event), ]
    cat(sprintf("  Held-out survival LL-pref K_init = %2d (cv_survival_logPL_per_event=%.4f)\n",
                cv_surv_best$K_init, cv_surv_best$cv_survival_logPL_per_event))
  }
  if (nrow(bicv_ok) > 0) {
    bicv_best <- bicv_ok[which.max(bicv_ok$bicv_genomics_total_loglik), ]
    cat(sprintf("  Bi-CV genomics LL-pref    K_init = %2d (bicv_genomics_total_loglik=%.1f)\n",
                bicv_best$K_init, bicv_best$bicv_genomics_total_loglik))
  }
  if (elbo_best$K_init != bic_best$K_init || elbo_best$K_init != c_best$K_init) {
    cat("  NOTE: criteria disagree on the preferred K_init — see DECISIONS.md before presenting a single number.\n")
  }
}

out_csv <- file.path(OUT_DIR, "k_init_sweep_results.csv")
write.csv(results, out_csv, row.names = FALSE)
if (length(fits) > 0) saveRDS(fits, FITS_RDS)
if (length(cv_results) > 0) saveRDS(cv_results, CV_CACHE_RDS)

cat(sprintf("\n=== Results saved: %s ===\n", out_csv))
