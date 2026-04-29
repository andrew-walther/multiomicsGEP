# ============================================================
# Script: run_lambda_sweep.R
# Purpose: Sensitivity analysis for the lambda survival-scaling parameter on
#          the merged TCGA_PAAD + CPTAC v2-preprocessed SSBMF fit.
#
#          Rationale: with p ≈ 2000 genes and n ≈ 230 samples, the genomics
#          reconstruction ELBO term dominates the survival term by a factor of
#          ~p/n ≈ 10 when lambda = 1. This sweep tests whether amplifying the
#          survival objective (lambda > 1) is sufficient to yield non-zero β
#          coefficients on the merged cross-platform cohort.
#
#          Lambda grid  : {1, 5, 10, 20}  (p/n ≈ 10 is the scale-equalized target)
#          Prior grid   : {point_normal, point_laplace, normal}
#          Preprocessing: v2 for all runs (intersect → log2 → QN → top-2000 → rank)
#
#          lambda = 1 results already exist at outputs/real_data/merged/v2_{prior}/.
#          Those are included in the summary table by reading the saved RDS files,
#          so only lambda ∈ {5, 10, 20} trigger new benchmark runs (9 fits total).
#
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-29
# Dependencies: yaml, survival (via run_ssbmf_benchmark.R)
# Inputs:  results/benchmark_sim/run_ssbmf_benchmark.R
#          results/benchmark_sim/run_phase1_diagnostics.R
# Outputs: results/benchmark_sim/outputs/real_data/merged/v2_lambda{X}_{prior}/
#            tables/final_model.rds, training_beta_summary.csv, ...
#          results/benchmark_sim/outputs/diagnostic_heatmaps/
#            merged_lambda{X}_{prior}_phase1_loading_heatmap.{pdf,png}
#          results/benchmark_sim/outputs/real_data/lambda_sweep_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
})

# ── resolve repo root ──────────────────────────────────────────────────────────
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R")) {
  if (file.exists("../../code/fit_modular.R")) setwd("../..")
}

# ── source dependencies ────────────────────────────────────────────────────────
# Sourcing run_ssbmf_benchmark.R loads all helper functions, cfg, and
# PDAC_DATA_ROOT.  The if (sys.nframe() == 0) guard prevents the entry-point
# block from executing on source().
source("results/benchmark_sim/run_ssbmf_benchmark.R")
source("results/benchmark_sim/run_phase1_diagnostics.R")

# ── sweep parameters ───────────────────────────────────────────────────────────
LAMBDA_GRID <- c(1, 5, 10, 20)
PRIOR_GRID  <- c("point_normal", "point_laplace", "normal")

REAL_DATA_ROOT <- file.path("results/benchmark_sim/outputs/real_data")

# ── helper: output directory for a given (lambda, prior) combination ──────────
#' Returns the output directory for a (lambda, prior) combination.
#' lambda = 1 uses the existing v2_{prior} directories; lambda > 1 uses
#' v2_lambda{X}_{prior} so existing results are never overwritten.
lambda_out_dir <- function(lambda, prior,
                           real_data_root = REAL_DATA_ROOT) {
  tag <- if (lambda == 1) paste0("v2_", prior) else
    sprintf("v2_lambda%g_%s", lambda, prior)
  file.path(real_data_root, "merged", tag)
}

# ── helper: read key scalars from a completed run directory ───────────────────
#' Extract summary scalars from a completed run directory.
#' Returns NULL (with a warning) if key files are missing.
read_run_summary <- function(out_dir, lambda, prior, beta_thresh = 0.05) {
  rds_path   <- file.path(out_dir, "tables", "final_model.rds")
  elbo_path  <- file.path(out_dir, "tables", "training_elbo_trace.csv")

  if (!file.exists(rds_path) || !file.exists(elbo_path)) {
    warning(sprintf("Missing output files in %s — skipping.", out_dir))
    return(NULL)
  }

  m    <- readRDS(rds_path)
  elbo <- read.csv(elbo_path, stringsAsFactors = FALSE)

  data.frame(
    lambda        = lambda,
    prior_beta    = prior,
    ELBO          = round(tail(elbo$ELBO_Full, 1), 2),
    n_active_beta = sum(abs(m$EBeta) > beta_thresh),
    max_abs_beta  = round(max(abs(m$EBeta)), 6),
    stringsAsFactors = FALSE
  )
}

# ── main sweep loop ────────────────────────────────────────────────────────────
cat(sprintf(
  "\n=== Lambda sweep: %d lambda × %d prior = %d combinations ===\n",
  length(LAMBDA_GRID), length(PRIOR_GRID),
  length(LAMBDA_GRID) * length(PRIOR_GRID)
))
cat(sprintf("  Lambda grid : %s\n", paste(LAMBDA_GRID, collapse = ", ")))
cat(sprintf("  Prior grid  : %s\n", paste(PRIOR_GRID,  collapse = ", ")))
cat(sprintf("  lambda = 1  : reading from existing v2 run directories (no re-fit)\n"))
cat(sprintf("  lambda > 1  : %d new fits\n\n",
            sum(LAMBDA_GRID > 1) * length(PRIOR_GRID)))

summary_rows <- list()

for (lam in LAMBDA_GRID) {
  for (prior in PRIOR_GRID) {

    out_dir    <- lambda_out_dir(lam, prior)
    sweep_label <- sprintf("merged_lambda%g_%s", lam, prior)

    # ── lambda = 1: existing results, skip re-fitting ─────────────────────────
    if (lam == 1) {
      cat(sprintf("  [lambda=%g, prior=%-14s] Using existing results at %s\n",
                  lam, prior, out_dir))

      # Generate heatmap if it hasn't been produced yet (idempotent)
      heatmap_png <- file.path(out_dir, "figures", "phase1_loading_heatmap.png")
      if (!file.exists(heatmap_png) && file.exists(file.path(out_dir, "tables", "final_model.rds"))) {
        run_phase1_diagnostic(
          training_mode = "merged",
          prior_beta    = prior,
          base_dir      = out_dir,
          label         = sweep_label
        )
      } else if (file.exists(heatmap_png)) {
        # Heatmap exists — still ensure a correctly labelled copy is in diagnostic_heatmaps/
        compare_dir  <- normalizePath(
          file.path(REAL_DATA_ROOT, "diagnostic_heatmaps"), mustWork = FALSE)
        dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
        compare_stub <- file.path(compare_dir,
                                  paste0(sweep_label, "_phase1_loading_heatmap"))
        out_stub     <- file.path(out_dir, "figures", "phase1_loading_heatmap")
        file.copy(paste0(out_stub, ".pdf"), paste0(compare_stub, ".pdf"), overwrite = TRUE)
        file.copy(paste0(out_stub, ".png"), paste0(compare_stub, ".png"), overwrite = TRUE)
      }

      summary_rows[[sweep_label]] <- read_run_summary(out_dir, lam, prior)
      next
    }

    # ── lambda > 1: run new benchmark fit ─────────────────────────────────────
    cat(sprintf("\n===== lambda=%g | prior=%s =====\n", lam, prior))

    run_real_data_benchmark(
      training_mode         = "merged",
      prior_beta            = prior,
      preprocessing_version = "v2",
      output_root           = out_dir,
      lambda                = lam,
      max_iter              = cfg$cavi$max_iter,
      tol                   = cfg$cavi$tol,
      alpha_grid            = cfg$cavi$alpha_grid,
      K_max                 = cfg$cavi$k_max
    )

    # Phase 1 loading heatmap with lambda-tagged comparison filename
    run_phase1_diagnostic(
      training_mode = "merged",
      prior_beta    = prior,
      base_dir      = out_dir,
      label         = sweep_label
    )

    summary_rows[[sweep_label]] <- read_run_summary(out_dir, lam, prior)
  }
}

# ── summary table ──────────────────────────────────────────────────────────────
summary_df <- do.call(rbind, Filter(Negate(is.null), summary_rows))

# Order rows by lambda then prior for readability
prior_order  <- match(summary_df$prior_beta, PRIOR_GRID)
lambda_order <- match(summary_df$lambda,      LAMBDA_GRID)
summary_df   <- summary_df[order(lambda_order, prior_order), ]
rownames(summary_df) <- NULL

out_csv <- file.path(REAL_DATA_ROOT, "lambda_sweep_summary.csv")
write.csv(summary_df, out_csv, row.names = FALSE)

cat("\n=== Lambda sweep complete ===\n")
cat(sprintf("  Summary table: %s\n\n", out_csv))
print(summary_df)
