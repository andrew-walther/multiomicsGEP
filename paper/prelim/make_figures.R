# ============================================================
# Script: paper/prelim/make_figures.R
# Purpose: Publication-style versions of the proposal chapter's figures:
#          descriptive labels, no embedded titles, stated sign conventions.
#          Inputs are the same saved results the progress-book figures use.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: ggplot2, dplyr, tidyr; code/pathway_enrichment.R (load_d4_weights)
# Run (from the repo root): Rscript paper/prelim/make_figures.R
# ============================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("code/pathway_enrichment.R")
out <- "paper/prelim/figures"
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                  plot.title = element_blank()))

# 1. Gene weights of the four retained single-modality programs ----
d4 <- load_d4_weights()
kept <- c("3" = "Program 3\n(protective)", "5" = "Stromal", "6" = "Immune /\nquiescent stroma",
          "7" = "Program 7\n(adverse)")
top_n <- 12L
genes <- unique(unlist(lapply(as.integer(names(kept)), function(k) {
  d4$gene_names[order(d4$EF[, k], decreasing = TRUE)[seq_len(top_n)]]
})))
w <- d4$EF[match(genes, d4$gene_names), as.integer(names(kept)), drop = FALSE]
w <- sweep(w, 2, apply(d4$EF[, as.integer(names(kept)), drop = FALSE], 2, max), "/")
owner <- apply(w, 1, which.max)
hm <- data.frame(gene = rep(genes, length(kept)), program = rep(kept, each = length(genes)),
                 weight = as.vector(w))
hm$gene <- factor(hm$gene, levels = rev(genes[order(owner)]))
hm$program <- factor(hm$program, levels = kept)
p1 <- ggplot(hm, aes(program, gene, fill = weight)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#f2f1ec", high = "#104281", limits = c(0, 1),
                      name = "Loading (scaled to the\nprogram's largest)") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 7), panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(out, "single_programs_heatmap.png"), p1, width = 6.2, height = 7, dpi = 200)

# 2. Held-out concordance gain of the supervised model over two-step, by K_init ----
d <- read.csv("results/multi_cohort_sim/outputs/multicohort_sim_results.csv", stringsAsFactors = FALSE)
delta <- d |>
  select(scenario, K_init, arm, seed, c_index) |>
  pivot_wider(names_from = arm, values_from = c_index) |>
  mutate(delta = YFB_base - EBMF) |>
  group_by(scenario, K_init) |>
  summarise(mean = mean(delta, na.rm = TRUE), se = sd(delta, na.rm = TRUE) / sqrt(sum(!is.na(delta))),
            .groups = "drop")
scen <- c(all_shared = "All programs shared across cohorts",
          hybrid = "Shared and cohort-specific programs",
          nothing_shared = "All programs cohort-specific")
delta$scenario <- factor(scen[delta$scenario], levels = scen)
delta$K_init <- factor(delta$K_init, levels = sort(unique(delta$K_init)))
p2 <- ggplot(delta, aes(K_init, mean, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c("#2a78d6", "#eb6834", "#1baf7a"), name = NULL) +
  labs(x = "Starting rank", y = "Held-out C, supervised minus two-step") +
  guides(fill = guide_legend(nrow = 1))
ggsave(file.path(out, "single_delta_c_kinit.png"), p2, width = 8, height = 4.2, dpi = 200)

# 3. Cohort-specific survival coefficients (positive = higher hazard) ----
fits <- readRDS("results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_fits.rds")
B <- fits$joint_yfb_beta_c$EBeta
colnames(B) <- c(CPTAC = "CPTAC", TCGA_PAAD = "TCGA")[fits$joint_yfb_beta_c$beta_cohort_levels]
pb <- data.frame(program = factor(rep(paste("Program", seq_len(nrow(B))), ncol(B)),
                                  levels = paste("Program", seq_len(nrow(B)))),
                 cohort = rep(colnames(B), each = nrow(B)), beta = as.vector(B))
p3 <- ggplot(pb, aes(program, beta, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(CPTAC = "#2a78d6", TCGA = "#eb6834"), name = "Training cohort") +
  labs(x = NULL, y = "Survival coefficient\n(positive = higher hazard)")
ggsave(file.path(out, "single_percohort_beta.png"), p3, width = 7, height = 3.6, dpi = 200)

# 4. Multimodal: programs kept against starting rank ----
r <- read.csv("results/benchmark_sim/outputs/multimodal_yfb_pruning/pruning_comparison.csv")
r <- r[is.na(r$error) & r$variant %in% c("original", "pruned"), ]
r$variant <- factor(c(original = "Without the parsimony framework",
                      pruned = "With the parsimony framework")[r$variant],
                    levels = c("Without the parsimony framework", "With the parsimony framework"))
r$scenario <- factor(c(both_informative = "Two prognostic programs",
                       adverse_protective = "Adverse and protective programs",
                       low_variance_prognostic_factor = "Low-variance prognostic program")[r$scenario],
                     levels = c("Two prognostic programs", "Adverse and protective programs",
                                "Low-variance prognostic program"))
r$noise <- factor(ifelse(r$noise_scale == 1, "Each true program < 1% of variance",
                         "Each true program about 10% of variance"),
                  levels = c("Each true program about 10% of variance", "Each true program < 1% of variance"))
k <- aggregate(K_final ~ noise + scenario + variant + K_init, r,
               function(v) c(m = mean(v), lo = min(v), hi = max(v)))
k <- cbind(k[1:4], as.data.frame(k$K_final))
p4 <- ggplot(k, aes(K_init, m, colour = variant, shape = variant)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
  geom_hline(yintercept = 3, linetype = 2, colour = "grey40") +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6)) +
  geom_line(position = position_dodge(width = 0.6)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
  facet_grid(noise ~ scenario, labeller = label_wrap_gen(24)) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834"), name = NULL) +
  scale_shape_manual(values = c(16, 15), name = NULL) +
  scale_x_continuous(breaks = c(3, 5, 7, 10, 15)) +
  labs(x = "Starting rank", y = "Programs kept (mean and range over seeds)")
ggsave(file.path(out, "multimodal_pruning_K_final.png"), p4, width = 9, height = 6, dpi = 200)
cat("Proposal figures written to", out, "\n")
