# ============================================================
# Script:  results/multi_cohort_sim/build_ebmf_templates.R
# Purpose: Build realistic gene-program templates (F columns) for the
#          multi-cohort simulation by fitting EBMF (flashier) to the merged
#          real training cohorts (TCGA_PAAD + CPTAC).  These templates make the
#          simulated factors resemble real PDAC expression programs.
#
#          GUARDED: if PDAC_DATA_ROOT is unset or flashier is unavailable, the
#          function returns NULL and the DGP falls back to synthetic sparse F.
#          This keeps the whole pipeline runnable anywhere.
#
#          NB: this is the TEMPLATE-BUILDING EBMF (fit to REAL data to make
#          ground truth).  It is distinct from the BENCHMARK EBMF in the runner,
#          which is fit to each SIMULATED Y and scored as a competing method.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-14
# Dependencies: flashier (optional), code/preprocess_desurv.R,
#               results/benchmark_sim/benchmark_helpers.R (cfg must exist)
# ============================================================

# build_ebmf_templates ----
#' Fit EBMF to the merged real training cohorts and return unit-norm gene
#' programs to use as simulation templates.
#'
#' @param Kmax    greedy_Kmax for flashier (upper bound on factors).
#' @param top_n   genes selected per cohort before merging (D4-aligned).
#' @param cache   path to cache the result as .rds (skipped if NULL).
#' @return list(F = p×K matrix of unit-norm programs, gene_names, K_ebmf), or
#'   NULL when real data / flashier are unavailable (caller uses synthetic F).
#' @family multicohort-sim
build_ebmf_templates <- function(Kmax  = 20L,
                                  top_n = 1000L,
                                  cache = "results/multi_cohort_sim/outputs/ebmf_templates.rds") {

  # return cached templates if present
  if (!is.null(cache) && file.exists(cache)) {
    message("build_ebmf_templates: using cached templates at ", cache)
    return(readRDS(cache))
  }

  root <- Sys.getenv("PDAC_DATA_ROOT", unset = "")
  if (root == "") root <- tryCatch(path.expand(cfg$pdac$data_root_default),
                                   error = function(e) "")

  have_flashier <- requireNamespace("flashier", quietly = TRUE)
  have_data     <- root != "" && dir.exists(root)

  if (!have_flashier || !have_data) {
    message(sprintf("build_ebmf_templates: %s%s -> synthetic fallback.",
                    if (!have_flashier) "flashier missing " else "",
                    if (!have_data) "PDAC data not found" else ""))
    return(NULL)
  }

  message("build_ebmf_templates: fitting EBMF to merged training cohorts ...")
  train_cohorts <- cfg$pdac$training_cohorts
  raw <- lapply(setNames(train_cohorts, train_cohorts),
                function(d) load_pdac_raw(d, root))

  # D4-aligned preprocessing: per-platform z-std + combined-rank top-N per cohort.
  pp <- preprocess_merged_cohorts(
    cohort_raw_list          = raw,
    log_transform_flags      = PLATFORM_LOG_TRANSFORM[train_cohorts],
    top_n                    = top_n,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE,
    normalize_method         = "none",
    selection_per_cohort     = TRUE,
    selection_method         = "combined_rank"
  )

  fit  <- flashier::flash(pp$Y, var_type = 2, greedy_Kmax = Kmax,
                          backfit = TRUE, verbose = 0)
  ldfr <- flashier::ldf(fit, type = "2")     # $F is p × K_ebmf (unit-norm columns)

  out <- list(F = as.matrix(ldfr$F),
              gene_names = pp$gene_names,
              K_ebmf = fit$n_factors)

  if (!is.null(cache)) {
    dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
    saveRDS(out, cache)
    message("build_ebmf_templates: cached ", out$K_ebmf, " programs to ", cache)
  }
  out
}
