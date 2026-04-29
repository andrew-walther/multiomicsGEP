# ============================================================
# Script: run_cluster_a_external.R
# Purpose: External-cohort C-index comparison for Cluster A.
#          Trains TWO models on merged TCGA_PAAD + CPTAC v2 data:
#            * cluster_a: N_burnin=10, normalize_AB=TRUE (post-fix)
#            * baseline:  N_burnin=0,  normalize_AB=FALSE (pre-fix opt-ins
#                         disabled; inner-loop reorder is unconditional)
#          Projects each of the 5 external cohorts via predict_supervised_mf()
#          and reports Harrell's C-index per cohort for both fits.
#          Pass condition (§4.9): >= 1 external cohort with C-index improved
#          vs. baseline.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-29
# Inputs:  PDAC raw data; results/benchmark_sim/outputs/cluster_a_smoke/
#          tables/smoke_fit.rds (cluster A fit; reused if present)
# Outputs: results/benchmark_sim/outputs/cluster_a_external/
#            tables/external_cindex.csv
#            tables/baseline_fit.rds
# ============================================================

suppressPackageStartupMessages({
  library(survival)
})

# ── resolve repo root ──────────────────────────────────────────────────────────
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R")) {
  if (file.exists("../../code/fit_modular.R")) setwd("../..")
}

source("results/benchmark_sim/run_ssbmf_benchmark.R")  # loaders + EXTERNAL_COHORTS

# Use the single-cohort preprocessing (per-cohort log/top-N/rank) for external
# cohorts, matching how run_ssbmf_benchmark.R projects them.
PDAC_DATA_ROOT <- Sys.getenv(
  "PDAC_DATA_ROOT",
  unset = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Liu)/PDAC_data"
)
PLATFORM_LOG_TRANSFORM <- c(TCGA_PAAD = TRUE, CPTAC = FALSE,
                            Dijk = FALSE, Moffitt_GEO_array = FALSE,
                            PACA_AU_array = FALSE, PACA_AU_seq = TRUE,
                            Puleo_array = FALSE)
TOP_N        <- 2000
K            <- 20
MAX_ITER     <- 60
ALPHA        <- 0.5
LAMBDA       <- 1.0
PRIOR_LF     <- "point_normal"
PRIOR_BETA   <- "point_normal"

OUT_ROOT  <- "results/benchmark_sim/outputs/cluster_a_external"
TABLE_DIR <- file.path(OUT_ROOT, "tables")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Reload Cluster A fit (from Step 5 smoke run) + re-fit baseline
# ==============================================================================
cat("=== Cluster A external-cohort C-index ===\n")

smoke_rds <- "results/benchmark_sim/outputs/cluster_a_smoke/tables/smoke_fit.rds"
if (!file.exists(smoke_rds))
  stop("Smoke fit not found at ", smoke_rds, " — run run_cluster_a_smoke.R first.")
smoke    <- readRDS(smoke_rds)
fit_A    <- smoke$fit
Y_train  <- smoke$Y
time_train   <- smoke$time
status_train <- smoke$status
training_gene_names <- smoke$gene_names
cat(sprintf("  Loaded Cluster A fit: n_train=%d, p_train=%d\n",
            nrow(Y_train), ncol(Y_train)))
cat(sprintf("  Cluster A EBeta range: [%.4f, %.4f] | active: %d/%d\n",
            min(fit_A$EBeta), max(fit_A$EBeta),
            sum(abs(fit_A$EBeta) > 0.05), K))

cat("\n  Fitting BASELINE (N_burnin=0, normalize_AB=FALSE) on same data ...\n")
set.seed(42)
fit_B <- fit_supervised_mf_modular(
  Y = Y_train, time = time_train, status = status_train,
  K = K, max_iter = MAX_ITER,
  prior_LF = PRIOR_LF, prior_beta = PRIOR_BETA,
  alpha = ALPHA, lambda = LAMBDA,
  init_method = "svd",
  N_burnin = 0,
  normalize_AB = FALSE,
  verbose = FALSE
)
cat(sprintf("  Baseline EBeta range:  [%.4f, %.4f] | active: %d/%d\n",
            min(fit_B$EBeta), max(fit_B$EBeta),
            sum(abs(fit_B$EBeta) > 0.05), K))

saveRDS(list(fit = fit_B, Y = Y_train, time = time_train,
             status = status_train, gene_names = training_gene_names),
        file.path(TABLE_DIR, "baseline_fit.rds"))

# ==============================================================================
# 2. External cohort C-index for both fits
# ==============================================================================
EXTERNAL_COHORTS <- c("Dijk", "Moffitt_GEO_array", "PACA_AU_array",
                      "PACA_AU_seq", "Puleo_array")

eval_external <- function(ds, fit, label, training_gene_names) {
  res <- tryCatch({
    raw <- load_pdac_raw(ds, PDAC_DATA_ROOT)
    log_flag <- PLATFORM_LOG_TRANSFORM[[ds]]
    if (is.null(log_flag)) log_flag <- FALSE
    preproc <- preprocess_desurv_cohort(
      Y = raw$Y, gene_names = raw$gene_names,
      top_n = TOP_N, log_transform = log_flag,
      cohort_name = ds
    )
    common_idx <- match(training_gene_names, preproc$gene_names)
    Y_ext <- matrix(0.0, nrow = raw$n, ncol = length(training_gene_names))
    present <- !is.na(common_idx)
    Y_ext[, present] <- preproc$Y[, common_idx[present], drop = FALSE]

    pred <- predict_supervised_mf(Y_ext, fit$EF, fit$EBeta)
    if (sd(pred$risk_scores) < 1e-10) {
      # Constant risk score -> C-index undefined; report NA
      list(c_index = NA_real_, n = raw$n, events = sum(raw$status),
           p_intersect = sum(present), constant_risk = TRUE)
    } else {
      cidx <- as.numeric(concordance(
        Surv(raw$time, raw$status) ~ I(-pred$risk_scores)
      )$concordance)
      list(c_index = round(cidx, 4), n = raw$n, events = sum(raw$status),
           p_intersect = sum(present), constant_risk = FALSE)
    }
  }, error = function(e) {
    cat(sprintf("    [%s] %s: ERROR %s\n", label, ds, conditionMessage(e)))
    list(c_index = NA_real_, n = NA_integer_, events = NA_integer_,
         p_intersect = NA_integer_, constant_risk = NA)
  })
  res$dataset <- ds
  res$model   <- label
  res
}

rows <- list()
for (ds in EXTERNAL_COHORTS) {
  cat(sprintf("\n  === %s ===\n", ds))
  res_A <- eval_external(ds, fit_A, "cluster_a", training_gene_names)
  res_B <- eval_external(ds, fit_B, "baseline",  training_gene_names)
  cat(sprintf("    cluster_a C=%s   baseline C=%s   (n=%s, events=%s, p_int=%s)\n",
              format(res_A$c_index, nsmall=4),
              format(res_B$c_index, nsmall=4),
              res_A$n, res_A$events, res_A$p_intersect))
  rows[[length(rows) + 1L]] <- res_A
  rows[[length(rows) + 1L]] <- res_B
}

# ==============================================================================
# 3. Write summary table
# ==============================================================================
df_long <- do.call(rbind, lapply(rows, function(r) {
  data.frame(
    dataset = r$dataset, model = r$model,
    n = r$n, events = r$events, p_intersect = r$p_intersect,
    c_index = r$c_index, constant_risk = r$constant_risk,
    stringsAsFactors = FALSE
  )
}))

# Wide form: one row per cohort with both C-indices and Δ
ci_A <- df_long$c_index[df_long$model == "cluster_a"]
ci_B <- df_long$c_index[df_long$model == "baseline"]
df_wide <- data.frame(
  dataset      = df_long$dataset[df_long$model == "cluster_a"],
  n            = df_long$n[df_long$model == "cluster_a"],
  events       = df_long$events[df_long$model == "cluster_a"],
  p_intersect  = df_long$p_intersect[df_long$model == "cluster_a"],
  c_cluster_a  = ci_A,
  c_baseline   = ci_B,
  delta        = round(ci_A - ci_B, 4)
)

write.csv(df_long, file.path(TABLE_DIR, "external_cindex_long.csv"), row.names = FALSE)
write.csv(df_wide, file.path(TABLE_DIR, "external_cindex.csv"),      row.names = FALSE)

cat("\n=== External-cohort C-index summary ===\n")
print(df_wide, row.names = FALSE)

# Pass: at least 1 cohort improved vs baseline
improvements <- sum(!is.na(df_wide$delta) & df_wide$delta > 0)
cat(sprintf("\n  Cohorts improved vs baseline: %d / %d\n",
            improvements, nrow(df_wide)))
cat(sprintf("  PASS Step 5b: %s\n",
            if (improvements >= 1) "YES" else "NO"))

cat(sprintf("\nWrote: %s/{external_cindex.csv, external_cindex_long.csv, baseline_fit.rds}\n",
            TABLE_DIR))
