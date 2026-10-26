# ============================================================
# Script: prepare_real_single.R
# Purpose: Build the expression-only training and external data exactly as the
#          frozen single-modality benchmark (D4) uses them, align the frozen D4
#          risk scores and PACA-AU donor IDs, and save one git-ignored object.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: yaml; sources the benchmark loaders in code/ and results/benchmark_sim/
#               (analysis side only; the ssbmf package never sources them)
# Usage: PDAC_DATA_ROOT=... Rscript results/ssbmf/prepare_real_single.R  (from the repository root)
# Preprocessing: the current benchmark pipeline, unchanged (per-cohort z-scoring,
#   DeSurv combined-rank selection of the top 3,000 genes per training cohort,
#   intersected to 2,064 genes). TCGA, Dijk and PACA-AU RNA-seq are log2(x + 1)
#   transformed twice by this pipeline (recorded 2026-10-16); kept so that the
#   comparison with D4 is on identical inputs.
# ============================================================

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")

main_checkout <- "/Users/ajwalther/GithubProjects/multiomicsGEP"
snap <- file.path(main_checkout, "results/benchmark_snapshot_2026-10-14/single_modality")
out_dir <- "results/ssbmf/outputs/real_single"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
root <- normalizePath(PDAC_DATA_ROOT)

# Training: TCGA + CPTAC ----
TR <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TR, TR), load_pdac_raw, pdac_root = root)
pp <- preprocess_merged_cohorts(train_raw, PLATFORM_LOG_TRANSFORM[TR],
                                top_n = cfg$preprocessing$top_n_genes_desurv,
                                rank_transform = FALSE, per_platform_standardize = TRUE,
                                normalize_method = "none", selection_per_cohort = TRUE,
                                selection_method = "combined_rank")
genes <- pp$gene_names
frozen_genes <- readRDS(file.path(snap, "d4_gene_names.rds"))
stopifnot(identical(genes, frozen_genes))
Y <- pp$Y; colnames(Y) <- genes
time <- unlist(lapply(TR, function(d) train_raw[[d]]$time))
status <- unlist(lapply(TR, function(d) train_raw[[d]]$status))
cohort <- rep(TR, vapply(TR, function(d) train_raw[[d]]$n, 1L))
ids <- unlist(lapply(TR, function(d) {
  s <- train_raw[[d]]$sampID; if (is.null(s)) paste0(d, "_", seq_len(train_raw[[d]]$n)) else paste0(d, "_", s)
}))
rownames(Y) <- make.unique(ids)
stopifnot(nrow(Y) == 273, ncol(Y) == 2064)

# External cohorts ----
frozen_risk <- readRDS(file.path(snap, "desurv_comparison_riskscores.rds"))[["D4"]]
donor_map <- function(co, samp) {
  if (!co %in% c("PACA_AU_array", "PACA_AU_seq")) return(paste0(co, "_", seq_along(samp)))
  si <- readRDS(file.path(root, "original", paste0(co, ".rds")))$sampInfo
  s <- sub("_seq$", "", samp)
  d <- si$icgc_donor_id[match(s, si$icgc_sample_id)]
  if (anyNA(d)) stop("Unmatched PACA-AU samples in ", co)
  paste0("PACAAU_", d)
}
external <- list()
for (co in cfg$pdac$external_cohorts) {
  raw <- load_pdac_raw(co, root)
  pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
                                  log_transform = PLATFORM_LOG_TRANSFORM[[co]], cohort_name = co,
                                  rank_transform = FALSE, per_platform_standardize = TRUE)
  common <- intersect(genes, pre$gene_names)
  Ye <- pre$Y[, match(common, pre$gene_names), drop = FALSE]
  colnames(Ye) <- common
  fr <- frozen_risk[[co]]
  stopifnot(length(fr$risk) == nrow(Ye), all.equal(fr$time, raw$time), all.equal(fr$status, raw$status))
  donors <- donor_map(co, if (is.null(raw$sampID)) seq_len(nrow(Ye)) else raw$sampID)
  rownames(Ye) <- make.unique(paste0(co, "_", seq_len(nrow(Ye))))
  external[[co]] <- list(Y = Ye, time = raw$time, status = raw$status, donor = donors,
                         genes_present = length(common), frozen_d4_risk = -fr$risk)  # D4 eta is a good-prognosis score
}

# Check the frozen D4 orientation reproduces the reported mean external C (0.627)
cfun <- function(r, t, s) survival::concordance(survival::Surv(t, s) ~ r, reverse = TRUE)$concordance
c_d4 <- vapply(external, function(e) cfun(e$frozen_d4_risk, e$time, e$status), 0)
print(round(c_d4, 3)); cat(sprintf("Mean external C of frozen D4: %.4f\n", mean(c_d4)))
stopifnot(abs(mean(c_d4) - 0.627) < 0.001)

saveRDS(list(train = list(Y = Y, time = time, status = status, cohort = cohort),
             external = external, genes = genes, d4_external_C = c_d4),
        file.path(out_dir, "data_expression.rds"))
cat("Saved", file.path(out_dir, "data_expression.rds"), "\n")
