# ============================================================
# Script:       compute_ph_diagnostics.R
# Purpose:      Re-fit SSBMF at the saved optimal alpha for
#               TCGA-only / point_normal and run cox.zph() on
#               each external PDAC cohort.  Writes
#               ph_diagnostics_table.csv to the benchmark
#               output directory without re-running alpha CV.
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-04-24
# Dependencies: results/benchmark_sim/run_ssbmf_benchmark.R
#               (sources all helpers automatically)
# Run from repo root:
#   Rscript results/benchmark_sim/compute_ph_diagnostics.R
# ============================================================

# ---- 0. Locate repo root and load helpers --------------------
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
} else if (file.exists("../code/fit_modular.R")) {
  setwd("..")
}

# Source runner with sys.nframe() > 0 so the entry-point block is skipped
source("results/benchmark_sim/run_ssbmf_benchmark.R")

# ---- 1. Configuration ----------------------------------------
training_mode <- "tcga_only"
prior_beta    <- "point_normal"

summary_path  <- file.path(
  "results/benchmark_sim/outputs/real_data",
  training_mode, prior_beta, "tables", "realdata_benchmark_summary.csv"
)
if (!file.exists(summary_path))
  stop("Saved benchmark summary not found. Run run_ssbmf_benchmark.R first.")

saved <- read.csv(summary_path, stringsAsFactors = FALSE)
alpha_opt  <- saved$alpha_opt[1]
K_max      <- saved$K_max[1]
p_genes    <- saved$p_genes[1]
cat(sprintf("Using saved alpha_opt=%.2f, K_max=%d, p_genes=%d\n",
            alpha_opt, K_max, p_genes))

top_n    <- 2000
max_iter <- 300
tol      <- 1e-5

output_dir <- file.path(
  "results/benchmark_sim/outputs/real_data",
  training_mode, prior_beta, "tables"
)

# ---- 2. Load + preprocess training cohort --------------------
cat("Loading TCGA_PAAD ...\n")
raw_train <- load_pdac_raw("TCGA_PAAD", PDAC_DATA_ROOT)
preproc_train <- preprocess_desurv_cohort(
  Y           = raw_train$Y,
  gene_names  = raw_train$gene_names,
  top_n       = top_n,
  log_transform = PLATFORM_LOG_TRANSFORM[["TCGA_PAAD"]],
  cohort_name = "TCGA_PAAD"
)
training_gene_names <- preproc_train$gene_names
cat(sprintf("Training set: n=%d, p=%d\n", nrow(preproc_train$Y), ncol(preproc_train$Y)))

# ---- 3. Fit SSBMF at saved alpha_opt -------------------------
cat(sprintf("Fitting SSBMF (alpha=%.2f, K=%d, max_iter=%d) ...\n",
            alpha_opt, K_max, max_iter))
final_fit <- fit_supervised_mf_modular(
  preproc_train$Y, raw_train$time, raw_train$status,
  K          = K_max,
  alpha      = alpha_opt,
  max_iter   = max_iter,
  tol        = tol,
  prior_beta = prior_beta,
  verbose    = TRUE
)
cat(sprintf("Converged: %s  Iterations: %d\n",
            final_fit$history$converged, final_fit$history$n_iter))

# Save fitted model (enables future downstream scripts to skip re-fitting)
saveRDS(
  list(EF = final_fit$EF, EBeta = final_fit$EBeta,
       alpha_opt = alpha_opt, training_gene_names = training_gene_names,
       training_mode = training_mode, prior_beta = prior_beta),
  file.path(output_dir, "final_model.rds")
)

# ---- 4. Project + cox.zph() for each external cohort ---------
cat("Running external cohort projections + PH diagnostics ...\n")

ph_results <- lapply(EXTERNAL_COHORTS, function(ds) {
  cat(sprintf("  %s ...\n", ds))
  tryCatch({
    raw     <- load_pdac_raw(ds, PDAC_DATA_ROOT)
    preproc <- preprocess_desurv_cohort(
      Y           = raw$Y,
      gene_names  = raw$gene_names,
      top_n       = top_n,
      log_transform = PLATFORM_LOG_TRANSFORM[[ds]],
      cohort_name = ds
    )

    common_idx <- match(training_gene_names, preproc$gene_names)
    Y_ext      <- matrix(0.0, nrow = raw$n, ncol = length(training_gene_names))
    present    <- !is.na(common_idx)
    Y_ext[, present] <- preproc$Y[, common_idx[present], drop = FALSE]

    pred   <- predict_supervised_mf(Y_ext, final_fit$EF, final_fit$EBeta)
    lp_vec <- as.vector(pred$risk_scores)

    c_idx <- as.numeric(concordance(
      Surv(raw$time, raw$status) ~ I(-lp_vec)
    )$concordance)

    # Grambsch-Therneau PH test
    coxfit <- coxph(Surv(raw$time, raw$status) ~ lp_vec, ties = "efron")
    zph    <- cox.zph(coxfit)
    tbl    <- zph$table
    ph_p   <- round(tbl["lp_vec", "p"], 4)

    data.frame(
      Cohort   = ds,
      n        = raw$n,
      Events   = sum(raw$status),
      Platform = PLATFORM_MAP[[ds]],
      C_index  = round(c_idx, 4),
      PH_chisq = round(tbl["lp_vec", "chisq"], 3),
      PH_df    = as.integer(tbl["lp_vec", "df"]),
      PH_p     = ph_p,
      PH_flag  = if (ph_p < 0.05) "FLAG (p<0.05)" else "PASS",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning(sprintf("%s failed: %s", ds, conditionMessage(e)))
    data.frame(
      Cohort   = ds, n = NA_integer_, Events = NA_integer_,
      Platform = NA_character_, C_index = NA_real_,
      PH_chisq = NA_real_, PH_df = NA_integer_, PH_p = NA_real_,
      PH_flag  = "ERROR", stringsAsFactors = FALSE
    )
  })
})

ph_df <- do.call(rbind, ph_results)

cat("\n=== PH Diagnostic Results ===\n")
print(ph_df[, c("Cohort", "n", "Events", "C_index", "PH_chisq", "PH_p", "PH_flag")])

out_path <- file.path(output_dir, "ph_diagnostics_table.csv")
write.csv(ph_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_path))
