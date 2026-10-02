# ============================================================
# Script: plot_program_loadings.R
# Purpose: Heatmap of the top-loading genes and CpGs for each program of a
#          cached TCGA multimodal fit, showing which programs use expression,
#          methylation, or both.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: ggplot2
# Run: Rscript results/multimodal_real/plot_program_loadings.R [K_init]  (default 7)
# ============================================================
#
# Each panel shows the top n features of each program by |loading|, within
# that modality. Loadings are scaled by the largest |loading| of that program
# in that modality, so each program's strongest feature is +/-1. Programs
# whose share of the modality's fitted signal is below 2% are shown greyed
# (their scaled loadings are not meaningful). Signal share = ||L_k||^2 ||F_mk||^2
# over its sum across programs (see characterize_programs.R).

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args)) as.integer(args[1]) else 7L
out_dir <- "results/multimodal_real/outputs"
J <- readRDS(sprintf("results/multimodal_real/outputs/fits/fits_K%02d_point_laplace-point_laplace_pruned.rds", K))$joint
tab <- read.csv(file.path(out_dir, sprintf("real_program_table_K%02d.csv", K)))
labels <- c("1" = "1 Classical", "2" = "2 Immune", "3" = "3 Mixed", "4" = "4 Basal-like*",
            "5" = "5 Stroma", "6" = "6 Exocrine", "7" = "7 Methylation-only")
if (K != 7L) labels <- setNames(paste("Program", seq_len(ncol(J$EL))), seq_len(ncol(J$EL)))
n_top <- 6L

panel <- function(modality, share_col) {
  F <- J$EF[[modality]]
  rownames(F) <- J$training_spec$feature_names[[modality]]
  keep_prog <- which(tab[[share_col]] >= 0.02)
  feats <- unique(unlist(lapply(keep_prog, function(k) rownames(F)[order(-abs(F[, k]))[seq_len(n_top)]])))
  scaled <- sweep(F[feats, , drop = FALSE], 2, apply(abs(F), 2, max) + 1e-12, "/")
  d <- expand.grid(feature = feats, program = seq_len(ncol(F)), stringsAsFactors = FALSE)
  d$value <- scaled[cbind(match(d$feature, feats), d$program)]
  d$active <- d$program %in% keep_prog
  d$value[!d$active] <- NA
  # order features by the program they were selected for
  owner <- sapply(feats, function(f) keep_prog[which.max(abs(scaled[f, keep_prog]))])
  d$feature <- factor(d$feature, levels = rev(feats[order(owner)]))
  d$program <- factor(labels[as.character(d$program)], levels = labels)
  d$modality <- if (modality == "expression") "Genes (expression)" else "CpGs (methylation)"
  d
}
d <- rbind(panel("expression", "expression_share"), panel("methylation", "methylation_share"))
d$modality <- factor(d$modality, levels = c("Genes (expression)", "CpGs (methylation)"))
share_lab <- data.frame(
  program = factor(labels, levels = labels),
  expr = sprintf("%.0f%%", 100 * tab$expression_share), meth = sprintf("%.0f%%", 100 * tab$methylation_share))

# Diverging scale: two poles (blue / red) with a neutral midpoint; NA = program
# carries < 2% of this modality's signal
p <- ggplot(d, aes(program, feature, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_wrap(~modality, scales = "free_y", nrow = 1) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f2f1ec", high = "#e34948", midpoint = 0,
                       limits = c(-1, 1), na.value = "#d9d8d3",
                       name = "Loading\n(scaled to the program's\nlargest |loading|)") +
  labs(x = NULL, y = NULL,
       title = sprintf("Top-loading genes and CpGs by program (TCGA fit, K_init = %d)", K),
       subtitle = "Grey: program carries < 2% of that modality's fitted signal.  * survival-active") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
ggsave(file.path(out_dir, sprintf("real_program_loadings_K%02d.png", K)), p, width = 10, height = 9, dpi = 150)

# Modality share of each program's fitted signal (companion bar chart)
s <- rbind(data.frame(program = factor(labels, levels = labels), modality = "Expression", share = tab$expression_share),
           data.frame(program = factor(labels, levels = labels), modality = "Methylation", share = tab$methylation_share))
q <- ggplot(s, aes(program, share, fill = modality)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c(Expression = "#2a78d6", Methylation = "#eb6834"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(x = NULL, y = "Share of the modality's fitted signal",
       title = "How each program divides across expression and methylation") +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                                        legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, sprintf("real_program_modality_share_K%02d.png", K)), q, width = 8, height = 4, dpi = 150)
cat("Wrote heatmap and modality-share figures to", out_dir, "\n")
