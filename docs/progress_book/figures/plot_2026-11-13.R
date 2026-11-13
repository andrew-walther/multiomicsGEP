# ============================================================
# Script: docs/progress_book/figures/plot_2026-11-13.R
# Purpose: Build the November 13 progress-book figures from committed
#          aggregate CSVs (rank curve, training-CV forest, false-active
#          programs, methylation CV C, external C, multimodal weight path).
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: ggplot2, patchwork
# Run (from the repo root): Rscript docs/progress_book/figures/plot_2026-11-13.R
# ============================================================
# No model fitting and no patient-level data. Every input is an aggregate
# table already tracked in Git, so the book's CI does not need this script;
# the PNGs it writes are committed.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# Setup ----

bench <- "results/benchmark_sim/outputs"
out_dir <- "docs/progress_book/figures/2026-11-13"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Two-colour palette (validated for CVD separation and contrast on a light
# surface): blue for the reference or working arm, orange for the alternative.
blue <- "#2a78d6"
orange <- "#eb6834"
ink <- "#0b0b0b"
ink_muted <- "#52514e"
grid_col <- "#e6e5e1"

#' Read a committed CSV and stop if it is missing
#'
#' @param path Path relative to the repository root.
#' @return A data frame.
#' @examples
#' read_csv_checked("results/benchmark_sim/outputs/lf_prior_2026-11-09/sim_cells.csv")
#' @family nov13-figures
read_csv_checked <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

#' Shared theme: recessive grid, muted axes, no legend box
#'
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
#' @examples
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_nb()
#' @family nov13-figures
theme_nb <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
      axis.text = element_text(colour = ink_muted),
      axis.title = element_text(colour = ink),
      plot.title = element_text(colour = ink, size = base_size, face = "bold"),
      strip.text = element_text(colour = ink, face = "bold", hjust = 0),
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

#' Save a figure as a PNG at book resolution
#'
#' @param plot A ggplot or patchwork object.
#' @param name File name without extension.
#' @param width,height Size in inches.
#' @return The output path, invisibly.
#' @examples
#' save_png(ggplot(), "example", 4, 3)
#' @family nov13-figures
save_png <- function(plot, name, width, height) {
  path <- file.path(out_dir, paste0(name, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = 200, bg = "white")
  message("Wrote ", path)
  invisible(path)
}

# (a) Rank curve averaged over five fold partitions ----

#' Plot mean within-cohort CV C against rank K with the one-SE threshold
#'
#' Mean C is averaged over 5 partitions × 5 folds. Bars are ±1 mean fold SE.
#' The one-SE rule picks the smallest K with mean C ≥ max(C) − SE(max).
#' @return A ggplot object.
#' @family nov13-figures
plot_rank_curve <- function() {
  d <- read_csv_checked(file.path(bench, "training_cv_redo_2026-11-12/seeds_k_average.csv"))
  stopifnot(nrow(d) == 11L, sum(d$selected) == 1L)
  d$colour <- ifelse(d$selected, orange, blue)
  thr <- unique(d$threshold)
  ggplot(d, aes(K, mean_C)) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = ink_muted, linewidth = 0.4) +
    annotate("text", x = 0.6, y = thr + 0.003, label = "One-SE\nthreshold",
             hjust = 0, vjust = 0, size = 2.7, colour = ink_muted, lineheight = 0.9) +
    geom_line(colour = blue, linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_C - mean_fold_se, ymax = mean_C + mean_fold_se,
                      colour = colour), width = 0, linewidth = 0.5) +
    geom_point(aes(colour = colour), size = 2.2) +
    annotate("text", x = 5, y = d$mean_C[d$selected] + d$mean_fold_se[d$selected] + 0.004,
             label = "K = 5 selected", colour = ink, size = 3, vjust = 0) +
    geom_text(aes(y = 0.565, label = folds_converged), size = 2.7, colour = ink_muted) +
    annotate("text", x = 1.6, y = 0.565, label = "Converged\n(of 25)", size = 2.5,
             colour = ink_muted, hjust = 1, lineheight = 0.9) +
    scale_colour_identity() +
    scale_x_continuous(breaks = 2:12, limits = c(0.6, 12.3)) +
    scale_y_continuous(limits = c(0.555, 0.69)) +
    labs(x = "Rank K", y = "Mean within-cohort CV C (±1 fold SE)") +
    theme_nb()
}

# (b) Training-CV paired differences from the working model ----

#' Assemble the paired training-CV differences behind the forest plot
#'
#' Each row is (alternative − working model) in within-cohort C with the
#' paired SE over five folds. Comparisons 1–3 are at K = 7. The LB row is at
#' K = 5 and is LB − YFB = −(YFB − LB).
#' @return A data frame with columns group, label, diff, se.
#' @family nov13-figures
forest_data <- function() {
  lf <- read_csv_checked(file.path(bench, "lf_prior_2026-11-09/real_paired.csv"))
  a2 <- read_csv_checked(file.path(bench, "lf_prior_2026-11-09/a2/a2_real_paired.csv"))
  a2 <- a2[a2$reference == "PL", ]
  be <- read_csv_checked(file.path(bench, "lf_prior_2026-11-09/beta_real_paired.csv"))
  cv <- read_csv_checked(file.path(bench, "training_cv_redo_2026-11-12/paired.csv"))
  lb <- read_csv_checked(file.path(bench, "lb_yfb_training_2026-11-12/cv_summary.csv"))
  lb <- lb[lb$metric == "C_within", ]
  pick <- function(df, key, value) {
    r <- df[df[[key]] == value, ]
    stopifnot(nrow(r) == 1L)
    c(r$diff_within, r$se_within)
  }
  g1 <- "Point-exponential prior on L and F (K = 7)"
  g2 <- "Prior on β (K = 7)"
  g3 <- "Gene screen and scaling (K = 7)"
  g4 <- "Survival parameterization (K = 5)"
  rows <- list(
    c(g1, "X/sd", pick(lf, "arm", "PE")),
    c(g1, "X/sd with cohort strata", pick(a2, "setup", "PES")),
    c(g1, "Shifted z-scores", pick(a2, "setup", "SZ")),
    c(g1, "Within-patient ranks", pick(lf, "arm", "RANK")),
    c(g2, "Shared point-Laplace", pick(be, "arm", "B_PLAP")),
    c(g2, "Shared point-normal", pick(be, "arm", "B_PN")),
    c(g2, "Shared Normal", pick(be, "arm", "B_NSH")),
    c(g3, "Variance screen after quantile normalization", pick(cv, "input", "D2_Z")),
    c(g3, "X/sd, centered projection", pick(cv, "input", "D4_XSD")),
    c(g3, "Within-patient ranks then z-scores", pick(cv, "input", "D4_RANKZ")),
    c(g4, "LB with fold-in", c(-lb$YFB_minus_LB, lb$paired_SE))
  )
  d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(d) <- c("group", "label", "diff", "se")
  d$diff <- as.numeric(d$diff)
  d$se <- as.numeric(d$se)
  d$group <- factor(d$group, levels = c(g1, g2, g3, g4))
  d$label <- factor(d$label, levels = rev(d$label))
  d
}

#' Forest plot of paired training-CV differences from the working model
#'
#' Each comparison is its own panel with its title above it, stacked on a
#' shared x-axis range.
#' @return A patchwork object.
#' @family nov13-figures
plot_forest_stacked <- function() {
  d <- forest_data()
  xr <- range(c(d$diff - d$se, d$diff + d$se, 0))
  xr <- xr + c(-0.005, 0.005)
  groups <- levels(d$group)
  panels <- lapply(seq_along(groups), function(i) {
    di <- droplevels(d[d$group == groups[i], ])
    p <- ggplot(di, aes(diff, label)) +
      geom_vline(xintercept = 0, colour = ink_muted, linewidth = 0.4) +
      geom_errorbar(aes(xmin = diff - se, xmax = diff + se), width = 0, orientation = "y",
                     colour = blue, linewidth = 0.5) +
      geom_point(colour = blue, size = 2) +
      scale_x_continuous(limits = xr) +
      labs(title = groups[i], x = NULL, y = NULL) +
      theme_nb() +
      theme(plot.title = element_text(size = 9.5, face = "bold", colour = ink))
    if (i < length(groups)) p <- p + theme(axis.text.x = element_blank())
    else p <- p + labs(x = "Difference in within-cohort CV C from the working model (±1 paired SE)")
    p
  })
  heights <- vapply(groups, function(g) sum(d$group == g), numeric(1)) + 1.2
  wrap_plots(panels, ncol = 1, heights = heights)
}

# (c) False-active programs by simulation cell ----

cell_levels <- c("0.05%\nnear threshold", "0.25%", "1%", "5%")

#' Line panel of mean false-active programs per fit across planted shares
#'
#' @param d Data frame with columns cell (factor), arm, value.
#' @param arms Two arm names; the first is drawn blue, the second orange.
#' @param title Panel title.
#' @return A ggplot object.
#' @family nov13-figures
false_active_panel <- function(d, arms, title) {
  d$colour <- ifelse(d$arm == arms[1], blue, orange)
  last <- d[d$cell == cell_levels[4], ]
  ggplot(d, aes(cell, value, group = arm, colour = colour)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = ink_muted, linewidth = 0.4) +
    annotate("text", x = 0.6, y = 1.08, label = "Limit: 1 per fit", hjust = 0, vjust = 0,
             size = 2.8, colour = ink_muted) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    geom_text(data = last, aes(label = arm), hjust = -0.35, size = 3.2, colour = ink) +
    scale_colour_identity() +
    scale_x_discrete(expand = expansion(add = c(0.4, 0.9))) +
    scale_y_continuous(limits = c(0, 4.6), breaks = 0:4) +
    labs(title = title, x = "Planted variance share", y = "Mean false-active programs per fit") +
    theme_nb()
}

#' Two panels: L/F prior comparison (K = 7) and LB vs YFB (K = 5, w_S = 1)
#'
#' @return A patchwork object.
#' @family nov13-figures
plot_false_active <- function() {
  lf <- read_csv_checked(file.path(bench, "lf_prior_2026-11-09/sim_cells.csv"))
  map_lf <- setNames(cell_levels, c("s14_0.05pct", "s9_0.25pct", "s9_1pct", "s9_5pct"))
  lf <- lf[lf$cell %in% names(map_lf) & lf$arm %in% c("PL", "PE"), ]
  stopifnot(nrow(lf) == 8L)
  d1 <- data.frame(cell = factor(map_lf[lf$cell], levels = cell_levels),
                   arm = lf$arm, value = lf$mean_false_active)
  lb <- read_csv_checked(file.path(bench, "lb_yfb_simulation_2026-11-12/simulation_summary.csv"))
  lb <- lb[lb$design %in% c("S14_strong", "S9") & lb$w_S == 1, ]
  stopifnot(nrow(lb) == 8L)
  map_lb <- setNames(cell_levels, c("5e-04", "0.0025", "0.01", "0.05"))
  d2 <- data.frame(cell = factor(map_lb[format(lb$share, scientific = FALSE, drop0trailing = TRUE)],
                                 levels = cell_levels),
                   arm = lb$model, value = lb$false_active_mean)
  if (anyNA(d2$cell)) {
    d2$cell <- factor(map_lb[as.character(lb$share)], levels = cell_levels)
  }
  stopifnot(!anyNA(d2$cell))
  (false_active_panel(d1, c("PL", "PE"), "Priors on L and F (K = 7)") |
     false_active_panel(d2, c("YFB", "LB"), "LB and YFB (K = 5, w_S = 1)") +
     labs(y = NULL))
}

# (d) Methylation CV C against CpG count ----

#' Methylation-only CV C by CpG count with expression reference lines
#'
#' Only fits whose CV ran are plotted (point-exponential full fits at 10k–50k
#' CpGs did not converge and were not cross-validated).
#' @return A ggplot object.
#' @family nov13-figures
plot_methylation <- function() {
  d <- read_csv_checked(file.path("results/multimodal_real/outputs",
                                  "methylation_single_modality_2026-11-12/K05/cv_summary_iteration1000.csv"))
  ref <- unique(d[, c("arm", "expression_mean_C_1000")])
  stopifnot(nrow(ref) == 2L)
  pts <- d[is.finite(d$methylation_mean_C), ]
  pts$colour <- ifelse(pts$arm == "PL", blue, orange)
  ref$colour <- ifelse(ref$arm == "PL", blue, orange)
  ref$label <- sprintf("Expression, %s: %.3f", ifelse(ref$arm == "PL", "point-Laplace", "point-exponential"),
                       ref$expression_mean_C_1000)
  last <- do.call(rbind, lapply(split(pts, pts$arm), function(x) x[which.max(x$count), ]))
  last$label <- ifelse(last$arm == "PL", "Methylation, point-Laplace", "Methylation, point-exponential")
  ggplot(pts, aes(count, methylation_mean_C, colour = colour, group = arm)) +
    geom_hline(data = ref, aes(yintercept = expression_mean_C_1000, colour = colour),
               linetype = "dashed", linewidth = 0.4) +
    geom_text(data = ref, aes(x = 9000, y = expression_mean_C_1000 + 0.004, label = label),
              hjust = 0, vjust = 0, size = 2.9, colour = ink) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    geom_text(data = last, aes(label = label), hjust = 0, nudge_x = 0.05, size = 2.9, colour = ink) +
    scale_colour_identity() +
    scale_x_log10(breaks = c(1e4, 2.5e4, 5e4, 1e5, 182934),
                  labels = c("10k", "25k", "50k", "100k", "all\n182,934"),
                  limits = c(9000, 1.6e6)) +
    scale_y_continuous(limits = c(0.49, 0.64)) +
    labs(x = "CpG sites in the fit (log scale)", y = "Mean CV C (144 TCGA patients)") +
    theme_nb()
}

# (e) External C for the working model ----

#' External C forest for the K = 5 working model with the equal-series mean
#'
#' @return A ggplot object.
#' @family nov13-figures
plot_external <- function() {
  ext <- read_csv_checked(file.path(bench, "lb_yfb_external_2026-11-12/external_cindex.csv"))
  ext <- ext[ext$model == "YFB", ]
  stopifnot(nrow(ext) == 5L)
  pool <- read_csv_checked(file.path(bench, "lb_yfb_external_2026-11-12/external_pooled.csv"))
  pool <- pool[pool$metric == "equal_YFB", ]
  stopifnot(nrow(pool) == 1L)
  names_map <- c(Dijk = "Dijk", Moffitt_GEO_array = "Moffitt array",
                 PACA_AU_array = "PACA-AU array", PACA_AU_seq = "PACA-AU RNA-seq",
                 Puleo_array = "Puleo array")
  ext$label <- sprintf("%s\n%d records, %d events, %d genes missing",
                       names_map[ext$cohort], ext$n, ext$events, ext$missing_genes)
  d <- rbind(
    data.frame(label = ext$label, C = ext$C, lower = ext$lower, upper = ext$upper, pooled = FALSE),
    data.frame(label = "Equal-series mean", C = pool$estimate, lower = pool$lower,
               upper = pool$upper, pooled = TRUE)
  )
  d$label <- factor(d$label, levels = rev(d$label))
  d$colour <- ifelse(d$pooled, orange, blue)
  ggplot(d, aes(C, label, colour = colour)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = ink_muted, linewidth = 0.4) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, orientation = "y", linewidth = 0.5) +
    geom_point(aes(shape = pooled), size = 2.4) +
    geom_text(aes(x = upper, label = sprintf("%.3f", C)), hjust = -0.3, size = 2.9, colour = ink) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0.44, 0.86)) +
    labs(x = "External C (95% interval)", y = NULL) +
    theme_nb() +
    theme(axis.text.y = element_text(lineheight = 0.9))
}

# (f) Multimodal weight path ----

#' Two panels: methylation share of the prognostic program and empty blocks
#'
#' Left: mean share of the prognostic program's reconstruction precision that
#' comes from methylation. Right: empty methylation blocks out of 50
#' (10 datasets × 5 programs). One line per methylation active fraction.
#' @return A patchwork object.
#' @family nov13-figures
plot_weight_path <- function() {
  wp <- file.path(bench, "feature_group_simulation_2026-11-12_weight_path")
  m <- read_csv_checked(file.path(wp, "summary_matched.csv"))
  m <- m[m$block == "A" & m$truth_program == 1 & m$modality == "methylation", ]
  e <- read_csv_checked(file.path(wp, "summary_empty.csv"))
  e <- e[e$block == "A" & e$modality == "methylation", ]
  stopifnot(nrow(m) == 10L, nrow(e) == 10L, all(e$denominator == 50))
  lab <- function(f) sprintf("%g%% CpGs active", 100 * f)
  m$colour <- ifelse(m$methylation_active_fraction == 0.03, blue, orange)
  e$colour <- ifelse(e$methylation_active_fraction == 0.03, blue, orange)
  m$label <- lab(m$methylation_active_fraction)
  e$label <- lab(e$methylation_active_fraction)
  brk <- c(0.2, 0.3, 0.4472, 0.7, 1)
  xs <- scale_x_continuous(breaks = brk, labels = c("0.2", "0.3", "0.45", "0.7", "1"),
                           limits = c(0.18, 1.32))
  mlast <- m[m$weight == 1, ]
  elast <- e[e$weight == 1, ]
  p1 <- ggplot(m, aes(weight, 100 * mean_share, group = label, colour = colour)) +
    geom_line(linewidth = 0.5) + geom_point(size = 2) +
    geom_text(data = mlast, aes(label = label), hjust = -0.12, size = 2.9, colour = ink) +
    scale_colour_identity() + xs +
    scale_y_continuous(limits = c(0, 70)) +
    labs(title = "Methylation share of the prognostic program",
         x = "Methylation likelihood weight", y = "Share of precision (%)") +
    theme_nb()
  p2 <- ggplot(e, aes(weight, empty, group = label, colour = colour)) +
    geom_line(linewidth = 0.5) + geom_point(size = 2) +
    geom_text(data = elast, aes(label = label), hjust = -0.12, size = 2.9, colour = ink) +
    scale_colour_identity() + xs +
    scale_y_continuous(limits = c(0, 50)) +
    labs(title = "Empty methylation loading blocks",
         x = "Methylation likelihood weight", y = "Empty blocks (of 50)") +
    theme_nb()
  p1 | p2
}

# Build ----

save_png(plot_rank_curve(), "rank-curve", 6.5, 3.6)
save_png(plot_forest_stacked(), "training-cv-differences", 7.5, 5.6)
save_png(plot_false_active(), "false-active", 8, 3.6)
save_png(plot_methylation(), "methylation-cv", 7, 3.6)
save_png(plot_external(), "external-c", 7, 3.8)
save_png(plot_weight_path(), "weight-path", 8.5, 3.4)
