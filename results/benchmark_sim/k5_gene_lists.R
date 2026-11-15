# ============================================================
# Script: k5_gene_lists.R
# Purpose: Gene lists for the five programs of the rank-5 joint SSBMF
#          fit: the n genes with the largest positive and the n with the most negative
#          posterior mean loading E[f_jk] in each program. Source for Chapter 4 Appendix C.15.
# Author: Andrew Walther
# Created: 2026-11-15
# Dependencies: base R
# Run from the repository root:
#   Rscript results/benchmark_sim/k5_gene_lists.R [n_per_pole]   (default 50)
# Writes results/benchmark_sim/outputs/k5_working_2026-11-12/prelim_ch4/gene_lists_k5.csv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
n_pole <- if (length(args)) as.integer(args[1]) else 50L

fit <- readRDS("results/benchmark_sim/outputs/k5_working_2026-11-12/s11/fits/size_D4_N2218_w1_f0.rds")
EF <- fit$state$EF            # p x K posterior mean loadings (canonical scale)
genes <- fit$genes
stopifnot(nrow(EF) == length(genes))

## Top genes at each pole of each active program ----
rows <- list()
for (k in seq_len(ncol(EF))) {
  pos <- order(-EF[, k])[seq_len(n_pole)]   # largest positive loadings
  neg <- order(EF[, k])[seq_len(n_pole)]    # most negative loadings
  rows[[length(rows) + 1]] <- data.frame(program = k, pole = "positive", rank = seq_len(n_pole),
                                         gene = genes[pos], loading = EF[pos, k])
  rows[[length(rows) + 1]] <- data.frame(program = k, pole = "negative", rank = seq_len(n_pole),
                                         gene = genes[neg], loading = EF[neg, k])
}
out <- do.call(rbind, rows)
stopifnot(all(out$loading[out$pole == "positive"] > 0), all(out$loading[out$pole == "negative"] < 0))
write.csv(out, "results/benchmark_sim/outputs/k5_working_2026-11-12/prelim_ch4/gene_lists_k5.csv",
          row.names = FALSE)
for (k in seq_len(ncol(EF))) for (pl in c("positive", "negative")) {
  cat(sprintf("P%d %s: %s\n", k, pl, paste(out$gene[out$program == k & out$pole == pl], collapse = ", ")))
}
