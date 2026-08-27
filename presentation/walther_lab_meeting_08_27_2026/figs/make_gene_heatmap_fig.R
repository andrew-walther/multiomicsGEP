# ============================================================
# Script:  presentation/walther_lab_meeting_08_27_2026/figs/make_gene_heatmap_fig.R
# Purpose: Build a presentation-sized gene-weight heatmap for the 4 ARD-kept
#          programs (K_init=7 D4 fit), TRANSPOSED relative to the pathway-
#          enrichment report's figure (programs x genes, not genes x
#          programs) so it fits a 16:9 slide: 4 rows, ~80 columns, instead of
#          ~80 rows and 4 columns. The report's figure (docs/progress_book/
#          figs/2026-08-21_k7_kept_factors_heatmap.png) is tall and narrow
#          and reads fine on a printed page; it does not fit a slide at
#          readable size, which is why this is a separate figure rather than
#          reusing that PNG.
#
#   Inputs:
#     results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds (D4 EF)
#     results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds
#     results/benchmark_sim/outputs/pathway_enrichment/K7_kept_factors_top_genes.csv
#       (defines the top-20-per-kept-program gene set to plot)
#   Output:
#     presentation/walther_lab_meeting_08_27_2026/assets/gene_program_heatmap.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-27
# Dependencies: pheatmap
# Usage:   Rscript presentation/walther_lab_meeting_08_27_2026/figs/make_gene_heatmap_fig.R
# ============================================================

suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer) })

if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

FITS_RDS   <- file.path(ROOT, "results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds")
GENES_RDS  <- file.path(ROOT, "results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")
TOPGENES_CSV <- file.path(ROOT, "results/benchmark_sim/outputs/pathway_enrichment/K7_kept_factors_top_genes.csv")
OUT_PNG    <- file.path(ROOT, "presentation/walther_lab_meeting_08_27_2026/assets/gene_program_heatmap.png")

for (f in c(FITS_RDS, GENES_RDS, TOPGENES_CSV)) if (!file.exists(f)) stop("Missing input: ", f)

d4         <- readRDS(FITS_RDS)$D4
gene_names <- readRDS(GENES_RDS)
top_genes  <- read.csv(TOPGENES_CSV, stringsAsFactors = FALSE)

stopifnot(nrow(d4$EF) == length(gene_names))
rownames(d4$EF) <- gene_names

# Kept factors in the order they should read left-to-right visually grouped:
# survival-active first (3 protective, 7 adverse), then genomics-only (5, 6).
KEPT_FACTORS <- c(3, 5, 6, 7)
FACTOR_LABEL <- c(`3` = "Program 3\n(protective)", `5` = "Program 5\n(stroma)",
                   `6` = "Program 6\n(immune)",     `7` = "Program 7\n(adverse)")

# Union of each kept factor's own top-20 genes, ordered by factor then rank so
# each program's block of genes stays contiguous (visually block-diagonal).
gene_order <- unique(top_genes$gene[order(match(top_genes$factor, KEPT_FACTORS), top_genes$rank)])

mat <- d4$EF[gene_order, KEPT_FACTORS, drop = FALSE]   # EF has no colnames; index by factor number
mat <- t(mat)                                   # TRANSPOSE: programs = rows, genes = columns
rownames(mat) <- FACTOR_LABEL[as.character(KEPT_FACTORS)]

n_genes <- ncol(mat)
png(OUT_PNG, width = 2400, height = 620, res = 150)
pheatmap(
  mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_colnames = TRUE, show_rownames = TRUE,
  angle_col = 90,
  fontsize_col = 7.6, fontsize_row = 13,
  color = colorRampPalette(c("#FFFFFF", "#FCBBA1", "#FB6A4A", "#A50F15"))(100),
  border_color = "grey55",  # visible cell/frame borders so the heatmap doesn't blend into the white slide background
  main = sprintf("Top-20 genes x 4 kept programs (K_init=7, D4 fit) -- %d genes total", n_genes),
  legend = TRUE
)
dev.off()

cat(sprintf("Wrote %s (%d programs x %d genes)\n", normalizePath(OUT_PNG), nrow(mat), n_genes))
