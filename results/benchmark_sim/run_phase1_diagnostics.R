# ============================================================
# Script: run_phase1_diagnostics.R
# Purpose: Phase 1 diagnostic visualizations for the merged-cohort
#          SSBMF fit. Produces a cohort-stratified L-loading heatmap
#          to identify factors that encode batch/study signal vs.
#          biological/prognostic signal.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-28
# Dependencies: pheatmap (CRAN)
# Inputs:  results/benchmark_sim/outputs/real_data/{mode}/{prior}/tables/final_model.rds
# Outputs: results/benchmark_sim/outputs/real_data/{mode}/{prior}/figures/
#            phase1_loading_heatmap.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
})

# Resolve repo root the same way run_ssbmf_benchmark.R does
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

# ------------------------------------------------------------
# Helper: build the cohort-stratified loading heatmap
# ------------------------------------------------------------

#' Cohort-stratified L-loading heatmap with β annotation strip.
#'
#' Rows are samples (sorted and grouped by cohort), columns are GEP factors.
#' Cell color = factor loading intensity. Left annotation = cohort label;
#' top annotation = |EBeta| per factor (highlighting prognostic vs. batch
#' factors). Mirrors the reference plot from the 4/27/26 meeting notes
#' (combined_microarray_snn_binary_K20.png).
#'
#' Cluster controls are intentionally OFF: row order is locked to the
#' cohort grouping (TCGA block, then CPTAC block). Column order preserves
#' GEP1..GEP_K_max as fit, so β annotation lines up with column positions.
#'
#' @param EL            n x K matrix of subject loadings (final_fit$EL)
#' @param EBeta         numeric vector of length K (final_fit$EBeta)
#' @param cohort_labels factor of length n, levels = unique cohort names
#' @param out_stub      file path prefix; '.pdf' and '.png' are appended
#' @param beta_thresh   |β| threshold separating "active" (red) from "shrunk" (gray)
#' @return invisibly: list with EL_sorted, cohort_order, beta_active_mask
plot_cohort_loading_heatmap <- function(EL, EBeta, cohort_labels,
                                        out_stub,
                                        beta_thresh = 0.05) {
  stopifnot(nrow(EL) == length(cohort_labels))
  stopifnot(ncol(EL) == length(EBeta))

  # Sort samples by cohort so the heatmap shows clean horizontal blocks
  cohort_order <- order(cohort_labels)
  EL_sorted    <- EL[cohort_order, , drop = FALSE]
  cohorts_ord  <- cohort_labels[cohort_order]

  K <- ncol(EL_sorted)
  factor_labels <- paste0("GEP", seq_len(K))
  rownames(EL_sorted) <- paste0("S", seq_len(nrow(EL_sorted)))
  colnames(EL_sorted) <- factor_labels

  # Row annotation: cohort label (one color per study)
  row_anno <- data.frame(Cohort = cohorts_ord)
  rownames(row_anno) <- rownames(EL_sorted)

  # Column annotation: |β| (continuous) + active/inactive (categorical strip)
  beta_active <- ifelse(abs(EBeta) > beta_thresh, "Active", "Shrunk")
  col_anno <- data.frame(
    Beta_Status = factor(beta_active, levels = c("Active", "Shrunk")),
    Abs_Beta    = abs(EBeta)
  )
  rownames(col_anno) <- factor_labels

  cohort_levels <- levels(cohorts_ord)
  cohort_palette <- setNames(
    c("#2166AC", "#D6604D", "#5AAE61", "#9970AB", "#FDB863",
      "#998EC3", "#80CDC1")[seq_along(cohort_levels)],
    cohort_levels
  )
  anno_colors <- list(
    Cohort      = cohort_palette,
    Beta_Status = c(Active = "#D62728", Shrunk = "#BDBDBD"),
    Abs_Beta    = c("#FFFFFF", "#08306B")
  )

  # Cell color palette: white → red, mirroring meeting reference image
  heatmap_palette <- colorRampPalette(c("#FFFFFF", "#FCBBA1", "#FB6A4A", "#A50F15"))(100)

  main_title <- sprintf(
    "L-loading heatmap (samples by cohort) | n=%d, K=%d, |b|>%.2f marks 'Active'",
    nrow(EL_sorted), K, beta_thresh
  )

  draw_fun <- function() {
    pheatmap(
      EL_sorted,
      cluster_rows        = FALSE,
      cluster_cols        = FALSE,
      show_rownames       = FALSE,
      show_colnames       = TRUE,
      color               = heatmap_palette,
      annotation_row      = row_anno,
      annotation_col      = col_anno,
      annotation_colors   = anno_colors,
      main                = main_title,
      fontsize_col        = 9,
      fontsize            = 9,
      gaps_row            = which(diff(as.integer(cohorts_ord)) != 0),
      border_color        = NA,
      legend              = TRUE,
      annotation_legend   = TRUE,
      silent              = FALSE
    )
  }

  pdf(paste0(out_stub, ".pdf"), width = 9, height = 7)
  draw_fun(); dev.off()

  png(paste0(out_stub, ".png"), width = 1100, height = 850, res = 140)
  draw_fun(); dev.off()

  invisible(list(
    EL_sorted        = EL_sorted,
    cohort_order     = cohort_order,
    beta_active_mask = abs(EBeta) > beta_thresh
  ))
}

# ------------------------------------------------------------
# Driver: load model, generate plot, write summary
# ------------------------------------------------------------

#' Run Phase 1 diagnostic for one (training_mode, prior_beta) combination.
#'
#' Reads final_model.rds saved by run_real_data_benchmark(), produces the
#' cohort-stratified loading heatmap, and writes a small summary CSV with
#' β values and active/shrunk status.
#'
#' @param training_mode  "merged", "tcga_only", or "cptac_only"
#' @param prior_beta     "point_normal", "point_laplace", or "normal"
#' @param output_root    benchmark output root (matches run_real_data_benchmark)
#' @param base_dir       optional override for the full directory containing
#'   tables/ and figures/. When supplied, output_root/training_mode/prior_beta
#'   is ignored. Useful for v2 preprocessing runs whose directory name differs
#'   from the standard {mode}/{prior} pattern.
#' @param beta_thresh    |β| threshold for "Active" vs "Shrunk" classification
#' @return invisibly: list returned by plot_cohort_loading_heatmap()
run_phase1_diagnostic <- function(training_mode = "merged",
                                  prior_beta    = "point_normal",
                                  output_root   = "results/benchmark_sim/outputs/real_data",
                                  base_dir      = NULL,
                                  beta_thresh   = 0.05) {

  if (is.null(base_dir))
    base_dir <- file.path(output_root, training_mode, prior_beta)
  rds_path <- file.path(base_dir, "tables", "final_model.rds")
  fig_dir  <- file.path(base_dir, "figures")

  if (!file.exists(rds_path))
    stop(sprintf("final_model.rds not found at %s — run run_real_data_benchmark() first.",
                 rds_path))
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  m <- readRDS(rds_path)

  required <- c("EL", "EBeta", "cohort_labels")
  missing_fields <- setdiff(required, names(m))
  if (length(missing_fields) > 0)
    stop(sprintf("RDS missing required fields: %s. Re-run benchmark with the updated saveRDS call.",
                 paste(missing_fields, collapse = ", ")))

  cat(sprintf("=== Phase 1 Diagnostic [mode=%s, prior=%s] ===\n",
              training_mode, prior_beta))
  cat(sprintf("  Loaded EL: %d x %d\n", nrow(m$EL), ncol(m$EL)))
  cat(sprintf("  Cohorts: %s\n",
              paste(sprintf("%s=%d", names(table(m$cohort_labels)),
                            as.integer(table(m$cohort_labels))),
                    collapse = ", ")))
  cat(sprintf("  Active factors (|EBeta| > %.2f): %d / %d\n",
              beta_thresh, sum(abs(m$EBeta) > beta_thresh), length(m$EBeta)))

  out_stub <- file.path(fig_dir, "phase1_loading_heatmap")
  result <- plot_cohort_loading_heatmap(
    EL            = m$EL,
    EBeta         = m$EBeta,
    cohort_labels = m$cohort_labels,
    out_stub      = out_stub,
    beta_thresh   = beta_thresh
  )

  # Copy to shared comparison folder with {mode}_{prior} prefix for easy side-by-side review
  compare_dir <- file.path(output_root, "..", "diagnostic_heatmaps")
  compare_dir <- normalizePath(compare_dir, mustWork = FALSE)
  dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
  compare_stub <- file.path(compare_dir,
                            sprintf("%s_%s_phase1_loading_heatmap", training_mode, prior_beta))
  file.copy(paste0(out_stub, ".pdf"), paste0(compare_stub, ".pdf"), overwrite = TRUE)
  file.copy(paste0(out_stub, ".png"), paste0(compare_stub, ".png"), overwrite = TRUE)

  # Save a small summary CSV
  beta_summary <- data.frame(
    Factor       = paste0("GEP", seq_along(m$EBeta)),
    EBeta        = round(m$EBeta, 6),
    Abs_EBeta    = round(abs(m$EBeta), 6),
    Beta_Status  = ifelse(abs(m$EBeta) > beta_thresh, "Active", "Shrunk")
  )
  write.csv(beta_summary,
            file.path(base_dir, "tables", "phase1_beta_summary.csv"),
            row.names = FALSE)

  cat(sprintf("  Heatmap written: %s.{pdf,png}\n", out_stub))
  cat(sprintf("  Comparison copy: %s.{pdf,png}\n", compare_stub))
  cat(sprintf("  Beta summary:    %s/tables/phase1_beta_summary.csv\n", base_dir))
  invisible(result)
}

# ------------------------------------------------------------
# Entry point: run for the merged + point_normal fit by default
# ------------------------------------------------------------
if (sys.nframe() == 0) {
  run_phase1_diagnostic(
    training_mode = "merged",
    prior_beta    = "point_normal"
  )
}
