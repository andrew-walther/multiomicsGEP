# ============================================================
# Script: make_writeup_figures.R
# Purpose: Draw the writeup's survival-curve, modality-share and
#          single-modality simulation figures at their printed size, as
#          vector PDF (PDF output) and 300-dpi PNG (HTML output), from
#          saved scores, the saved frozen fit and saved aggregates. No model
#          is fitted.
# Author: Andrew Walther
# Created: 2026-10-14
# Dependencies: survival, ggplot2, yaml; results/benchmark_sim/benchmark_helpers.R,
#               code/preprocess_desurv.R, code/predict_cox_on_yf.R
# Run (from the repository root; the program-curve figure needs the local
#      held-out expression data to recompute frozen program projections):
#   PDAC_DATA_ROOT=data/PDAC_data Rscript paper/writeup/figures/make_writeup_figures.R
# ============================================================
#
# Figures written to paper/writeup/figures/ (each as .pdf and .png):
#   km_risk_external       risk-score tertiles in the five held-out cohorts
#   km_programs_external   Program 7 and Program 3 projection tertiles, pooled
#   multimodal_share       modality share of each multimodal program's signal
#   single_sim_delta       model minus two-step test C in the cohort simulation
#
# Orientation. The saved single-modality fit's eta = Z beta is a
# good-prognosis score, so risk = -eta for every held-out cohort, fixed
# before any validation outcome is used (see results/benchmark_sim/make_km_figures.R,
# which this script follows for the curves).

suppressPackageStartupMessages({ library(survival); library(ggplot2) })
cfg <- yaml::read_yaml("config/globals.yml")   # benchmark_helpers.R reads cfg
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R"); source("code/predict_cox_on_yf.R")

out_dir <- "paper/writeup/figures"
width_in <- 6.5                                 # PDF text width (letter, 1-in margins)

#' Save one figure as vector PDF and 300-dpi PNG at its printed size
#'
#' @param plot ggplot object.
#' @param name file stem written to out_dir.
#' @param height printed height in inches.
#' @return Invisibly, the two paths written.
save_both <- function(plot, name, height) {
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width_in, height = height, device = "pdf", useDingbats = FALSE)
  ggsave(png_path, plot, width = width_in, height = height, dpi = 300, bg = "white")
  cat("Wrote", pdf_path, "and", png_path, "\n")
  invisible(c(pdf_path, png_path))
}

# Plot text is 8-9 pt at printed size (body text is 11 pt)
theme_set(theme_minimal(base_size = 9) +
            theme(panel.grid.minor = element_blank(),
                  text = element_text(colour = "black"),
                  axis.text = element_text(colour = "black"),
                  strip.text = element_text(face = "bold", size = 9, hjust = 0),
                  plot.background = element_rect(fill = "white", colour = NA)))
fmt_p <- function(p) if (p < 1e-4) "p < 0.0001" else
  sprintf("p = %s", format(signif(p, 2), scientific = FALSE))

#' Kaplan-Meier step data with censoring counts, restarted at (0, 1)
#'
#' @param d data frame with time, status and the grouping column.
#' @param group name of the grouping column.
#' @param facet panel label stored in the output.
#' @return Data frame of step coordinates per group.
km_steps <- function(d, group, facet) {
  s <- summary(survfit(stats::as.formula(paste("Surv(time, status) ~", group)), data = d),
               censored = TRUE)
  out <- do.call(rbind, lapply(split(data.frame(time = s$time, surv = s$surv,
                                                n.censor = s$n.censor, strata = s$strata),
                                     s$strata), function(x) {
    rbind(data.frame(time = 0, surv = 1, n.censor = 0, strata = x$strata[1]), x)
  }))
  out$group <- factor(sub(".*=", "", as.character(out$strata)),
                      levels = c("Low", "Middle", "High"))
  out$facet <- facet
  out
}
tertile <- function(x) cut(x, stats::quantile(x, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                           labels = c("Low", "Middle", "High"))

# Frozen single-modality scores and projections ----
dc <- "results/benchmark_sim/outputs/desurv_comparison"
fit <- readRDS(file.path(dc, "desurv_comparison_fits.rds"))[["D4"]]
train_genes <- readRDS(file.path(dc, "d4_gene_names.rds"))
saved <- readRDS(file.path(dc, "desurv_comparison_riskscores.rds"))[["D4"]]
cohort_labels <- c(Dijk = "Dijk (RNA-seq)", Moffitt_GEO_array = "Moffitt (microarray)",
                   PACA_AU_array = "PACA-AU (microarray)", PACA_AU_seq = "PACA-AU (RNA-seq)",
                   Puleo_array = "Puleo (microarray)")
pdac_root <- normalizePath(Sys.getenv("PDAC_DATA_ROOT", "data/PDAC_data"))

rows <- list()
for (cohort in cfg$pdac$external_cohorts) {
  raw <- load_pdac_raw(cohort, pdac_root)
  pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
                                  log_transform = cfg$pdac$platform_log_transform[[cohort]],
                                  cohort_name = cohort, rank_transform = FALSE,
                                  per_platform_standardize = TRUE)
  common <- intersect(train_genes, pre$gene_names)
  Y <- pre$Y[, match(common, pre$gene_names), drop = FALSE]
  EF <- fit$EF[match(common, train_genes), , drop = FALSE]
  eta <- predict_cox_on_yf(Y, EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
  # Fail loud if the rebuilt scores or outcomes differ from the saved ones
  stopifnot(length(eta) == saved[[cohort]]$n,
            max(abs(eta - saved[[cohort]]$risk)) < 1e-8,
            identical(as.numeric(raw$time), as.numeric(saved[[cohort]]$time)),
            identical(as.integer(raw$status), as.integer(saved[[cohort]]$status)))
  Z <- sweep(Y %*% EF, 2, fit$EF_norms, "/")       # frozen normalized projections
  d <- data.frame(cohort = cohort_labels[[cohort]], time = raw$time, status = raw$status,
                  risk = -eta, program7 = Z[, 7], program3 = Z[, 3])
  d$risk_group <- tertile(d$risk); d$p7_group <- tertile(d$program7)
  d$p3_group <- tertile(d$program3)
  rows[[cohort]] <- d
}
km_data <- do.call(rbind, rows)
km_data$cohort <- factor(km_data$cohort, levels = unname(cohort_labels))
stopifnot(nrow(km_data) == 616, sum(km_data$status) == 414)

# Figure: risk tertiles per held-out cohort ----
cols3 <- c(Low = "#2a78d6", Middle = "#eda100", High = "#e34948")
lab <- do.call(rbind, lapply(split(km_data, km_data$cohort), function(d) {
  p <- stats::pchisq(survdiff(Surv(time, status) ~ risk_group, data = d)$chisq,
                     df = 2, lower.tail = FALSE)
  C <- concordance(Surv(time, status) ~ risk, data = d, reverse = TRUE)$concordance
  data.frame(facet = as.character(d$cohort[1]),
             label = sprintf("n = %d, C = %.3f\nlog-rank %s", nrow(d), C, fmt_p(p)))
}))
steps <- do.call(rbind, lapply(split(km_data, km_data$cohort), function(d)
  km_steps(d, "risk_group", as.character(d$cohort[1]))))
steps$facet <- factor(steps$facet, levels = unname(cohort_labels))
lab$facet <- factor(lab$facet, levels = unname(cohort_labels))
p1 <- ggplot(steps, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.5) +
  geom_point(data = steps[steps$n.censor > 0, ], shape = 3, size = 0.8, show.legend = FALSE) +
  geom_text(data = lab, aes(x = Inf, y = 1, label = label), inherit.aes = FALSE,
            hjust = 1.03, vjust = 1.1, size = 2.6, lineheight = 0.95) +
  facet_wrap(~facet, ncol = 3) +
  scale_colour_manual(values = cols3, name = "Risk tertile\n(within cohort)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Months", y = "Overall survival probability") +
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.22),
        legend.key.width = unit(1.2, "lines"))
save_both(p1, "km_risk_external", height = 4.6)

# Figure: Program 7 and Program 3 tertiles, pooled across cohorts ----
seq3 <- c(Low = "#9ec5f4", Middle = "#3987e5", High = "#104281")
panels <- c(p7_group = "Program 7 (adverse; basal-like)",
            p3_group = "Program 3 (protective; classical)")
prog <- do.call(rbind, lapply(names(panels), function(g) km_steps(km_data, g, panels[[g]])))
prog$facet <- factor(prog$facet, levels = panels)
plab <- do.call(rbind, lapply(names(panels), function(g) {
  # cohort-stratified log-rank: tertiles are formed within cohort
  lr <- survdiff(stats::as.formula(paste("Surv(time, status) ~", g, "+ strata(cohort)")),
                 data = km_data)
  data.frame(facet = factor(panels[[g]], levels = panels),
             label = sprintf("n = %d, five cohorts\nstratified log-rank %s", nrow(km_data),
                             fmt_p(stats::pchisq(lr$chisq, df = 2, lower.tail = FALSE))))
}))
p2 <- ggplot(prog, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.5) +
  geom_point(data = prog[prog$n.censor > 0, ], shape = 3, size = 0.8, show.legend = FALSE) +
  geom_text(data = plab, aes(x = Inf, y = 1, label = label), inherit.aes = FALSE,
            hjust = 1.03, vjust = 1.1, size = 2.6, lineheight = 0.95) +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = seq3, name = "Projection tertile (within cohort)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Months", y = "Overall survival probability") +
  theme(legend.position = "bottom")
save_both(p2, "km_programs_external", height = 2.9)
print(lab); print(plab)

# Figure: modality share of each multimodal program ----
tab <- read.csv("results/multimodal_real/outputs/real_program_table_K07.csv")
stopifnot(nrow(tab) == 7)
prog_labels <- c("1 Classical/gastric", "2 Immune", "3 Mixed", "4 Basal-like*",
                 "5 Stroma", "6 Exocrine", "7 Methylation-only")
s <- rbind(data.frame(program = factor(prog_labels, levels = prog_labels),
                      modality = "Expression", share = tab$expression_share),
           data.frame(program = factor(prog_labels, levels = prog_labels),
                      modality = "Methylation", share = tab$methylation_share))
p3 <- ggplot(s, aes(program, share, fill = modality)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c(Expression = "#2a78d6", Methylation = "#eb6834"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Program", y = "Share of the modality's\nfitted signal") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top",
        panel.grid.major.x = element_blank())
save_both(p3, "multimodal_share", height = 2.7)

# Figure: single-modality simulation, model minus two-step test C ----
sim <- read.csv("results/multi_cohort_sim/outputs/multicohort_sim_results.csv")
stopifnot(nrow(sim) == 540, !anyNA(sim$c_index))
wide <- merge(sim[sim$arm == "YFB_base", c("scenario", "K_init", "seed", "c_index")],
              sim[sim$arm == "EBMF", c("scenario", "K_init", "seed", "c_index")],
              by = c("scenario", "K_init", "seed"), suffixes = c("_model", "_two_step"))
stopifnot(nrow(wide) == 270)                     # 3 scenarios x 6 ranks x 15 seeds, paired
wide$delta <- wide$c_index_model - wide$c_index_two_step
ds <- aggregate(delta ~ scenario + K_init, data = wide,
                FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
ds <- data.frame(ds[, 1:2], ds$delta)
scen <- c(all_shared = "All shared (6 shared programs)",
          hybrid = "Hybrid (2 shared, 2 specific per cohort)",
          nothing_shared = "Nothing shared (3 specific per cohort)")
stopifnot(all(ds$scenario %in% names(scen)))
ds$scenario <- factor(scen[ds$scenario], levels = scen)
ds$K_init <- factor(ds$K_init, levels = sort(unique(ds$K_init)))
p4 <- ggplot(ds, aes(K_init, mean, fill = scenario)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.4) +
  scale_fill_manual(values = c("#2a78d6", "#eb6834", "#1baf7a"), name = NULL) +
  labs(x = "Starting rank", y = "Test C, model minus two-step") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "top", panel.grid.major.x = element_blank())
save_both(p4, "single_sim_delta", height = 3.2)
