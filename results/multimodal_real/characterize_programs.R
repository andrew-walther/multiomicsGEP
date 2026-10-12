# ============================================================
# Script: characterize_programs.R
# Purpose: Characterize the programs of a cached TCGA multimodal fit: top
#          genes, how much of each modality's fitted signal each program
#          carries, survival classification, marginal survival association in
#          TCGA and ICGC, and correspondence with the expression-only fit.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: survival; code/multimodal_yfb_helpers.R, preprocess_multimodal_yfb.R
# Run: Rscript results/multimodal_real/characterize_programs.R [K_init]   (default 7)
# ============================================================
#
# Definitions
# - Fitted signal of program k in modality m: ||L_k||^2 ||F_mk||^2, the squared
#   Frobenius norm of its rank-one reconstruction L_k F_mk^T. "Methylation
#   share" is program k's fraction of the total fitted methylation signal
#   (summed over programs); likewise for expression. This avoids counting
#   "nonzero" CpGs, which depends on an arbitrary threshold under the
#   point-Laplace prior (posterior means are almost never exactly zero).
# - Marginal hazard ratio: univariable Cox model of survival on the program's
#   projection Z_k = sum_m Y_m F_mk, standardized within each cohort (per SD of
#   that cohort's projection). 7 programs x 2 cohorts = 14 tests, unadjusted:
#   descriptive only.
# - Correspondence: |Pearson r| between the expression loadings of the joint
#   fit and the expression-only fit (same 3,000 genes), best match per program.

suppressPackageStartupMessages(library(survival))
source("code/multimodal_yfb_helpers.R"); source("code/preprocess_multimodal_yfb.R")

args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args)) as.integer(args[1]) else 7L
fits_dir <- "results/multimodal_real/outputs/fits"
out_dir <- "results/multimodal_real/outputs"
d <- readRDS("data/processed_multiomics_tcga_icgc.rds")
J <- readRDS(file.path(fits_dir, sprintf("fits_K%02d_point_laplace-point_laplace_pruned.rds", K)))$joint
expr_only_path <- file.path(fits_dir, sprintf("fits_K%02d_point_laplace-point_laplace_pruned_expression-only.rds", K))
E <- if (file.exists(expr_only_path)) readRDS(expr_only_path)$joint else NULL
genes <- J$training_spec$feature_names$expression

projection <- function(fit, Y) Reduce(`+`, Map(function(y, f) y %*% f, Y[names(fit$EF)], fit$EF))
marginal_hr <- function(Z, cohort) {
  t(apply(Z, 2, function(z) {
    s <- summary(coxph(Surv(cohort$time, cohort$event) ~ as.vector(scale(z))))$coefficients
    c(HR = s[1, "exp(coef)"], p = s[1, "Pr(>|z|)"])
  }))
}
hr_tcga <- marginal_hr(projection(J, d$training$Y), d$training)
hr_icgc <- marginal_hr(projection(J, d$validation_primary$Y), d$validation_primary)

signal <- sapply(names(J$EF), function(m) colSums(J$EL^2) * colSums(J$EF[[m]]^2))
share <- sweep(signal, 2, colSums(signal), "/")

best_match <- if (!is.null(E)) {
  r <- abs(suppressWarnings(stats::cor(J$EF$expression, E$EF$expression)))
  r[!is.finite(r)] <- NA
  apply(r, 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
} else rep(NA, ncol(J$EL))

top_genes <- function(f, n, positive = TRUE) {
  o <- order(f, decreasing = positive)
  keep <- o[seq_len(n)]
  keep <- keep[if (positive) f[keep] > 0 else f[keep] < 0]
  paste(genes[keep], collapse = ", ")
}
tab <- data.frame(
  program = seq_len(ncol(J$EL)),
  class = J$diagnostics$factor_class,
  pve = round(J$diagnostics$factor_pve, 3),
  expression_share = round(share[, "expression"], 3),
  methylation_share = round(share[, "methylation"], 3),
  beta = round(J$EBeta, 3), beta_z = round(J$diagnostics$prognostic_z, 2),
  hr_tcga = round(hr_tcga[, "HR"], 2), p_tcga = signif(hr_tcga[, "p"], 2),
  hr_icgc = round(hr_icgc[, "HR"], 2), p_icgc = signif(hr_icgc[, "p"], 2),
  best_match_expression_only = round(best_match, 2),
  up_genes = apply(J$EF$expression, 2, top_genes, n = 8, positive = TRUE),
  down_genes = apply(J$EF$expression, 2, top_genes, n = 6, positive = FALSE)
)
print(tab[, 1:12], row.names = FALSE)
write.csv(tab, file.path(out_dir, sprintf("real_program_table_K%02d.csv", K)), row.names = FALSE)
cat("Wrote", file.path(out_dir, sprintf("real_program_table_K%02d.csv", K)), "\n")
