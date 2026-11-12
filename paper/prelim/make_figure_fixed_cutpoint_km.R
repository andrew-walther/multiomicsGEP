# ============================================================
# Script: make_figure_fixed_cutpoint_km.R
# Purpose: Draw external Kaplan–Meier curves for the fixed training cutpoint of the supervised score.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: ggplot2
# ============================================================
# Run from the repository root: Rscript paper/prelim/make_figure_fixed_cutpoint_km.R
# Inputs are the committed aggregate outputs of run_prompt7_real.R: monthly
# Kaplan–Meier estimates (km_monthly.csv), counts at risk every 12 months
# (at_risk.csv), the cutpoint (cutpoints.csv) and the group HRs (group_hr.csv).
# No patient-level data are read. Styling follows make_figures_supervised_external_km.R.

suppressPackageStartupMessages(library(ggplot2))

## cairo_pdf embeds fonts but needs the cairo shared library (it fails to load without XQuartz's X11 libs).
## capabilities("cairo") can report TRUE when the library cannot load, so test by opening a device.
CAIRO_OK <- tryCatch({ f <- tempfile(fileext = ".pdf"); suppressWarnings(cairo_pdf(f)); dev.off(); file.exists(f) },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (!CAIRO_OK) warning("cairo_pdf unavailable: figures are written with pdf() and their fonts are NOT embedded",
                       call. = FALSE)
PDF_DEVICE <- if (CAIRO_OK) cairo_pdf else function(filename, ...) grDevices::pdf(file = filename, ..., useDingbats = FALSE)

# EXHIBIT_K (results/benchmark_sim/exhibit_rank.R) reads the risk groups of another rank
# (cutpoint chosen on that rank's out-of-fold training scores) and writes a "_k<K>" file.
source("results/benchmark_sim/exhibit_rank.R")
out7 <- exhibit_out("results/benchmark_sim/outputs/prompt7_desurv_style", "prompt7_desurv_style")
figure_dir <- "paper/prelim/figures"

#' Read a saved aggregate CSV, failing loudly when it is missing.
#' @param name File name under the Prompt 7 output directory.
#' @return data.frame
#' @examples
#' # rd7("km_monthly.csv")
rd7 <- function(name) {
  p <- file.path(out7, name)
  if (!file.exists(p)) stop("Missing saved evidence: ", p)
  read.csv(p, stringsAsFactors = FALSE)
}

cohort_names <- c(Dijk = "Dijk", Moffitt_GEO_array = "Moffitt",
                  PACA_AU_array = "PACA-AU array", PACA_AU_seq = "PACA-AU RNA-seq",
                  Puleo_array = "Puleo")
groups <- c("Low", "High")
colors <- c(Low = "#006878", High = "#B34226")
patterns <- c(Low = "solid", High = "dotdash")

# Supervised score only; the two-step curves are in the writeup.
km <- subset(rd7("km_monthly.csv"), model == "w1")
ar <- subset(rd7("at_risk.csv"), model == "w1")
hr <- subset(rd7("group_hr.csv"), model == "w1")
stopifnot(nrow(km) > 0, nrow(ar) == 5L * 2L * 6L, all(names(cohort_names) %in% km$cohort))

km$group <- factor(groups[km$high + 1L], levels = groups)
ar$group <- factor(groups[ar$high + 1L], levels = groups)
km$series <- factor(cohort_names[km$cohort], levels = unname(cohort_names))
ar$series <- factor(cohort_names[ar$cohort], levels = unname(cohort_names))
ar$y <- -0.20 - 0.12 * (as.integer(ar$group) - 1L)

# Panel titles carry the series HR (high vs low) and its Wald interval.
hr <- hr[hr$cohort %in% names(cohort_names), ]
title <- sprintf("%s, HR %.2f (%.2f to %.2f)", cohort_names[hr$cohort], hr$hr, hr$lo, hr$hi)
levels(km$series) <- levels(ar$series) <- title[match(levels(km$series), cohort_names[hr$cohort])]

group_labels <- subset(ar, month == 0)
tick_labels <- expand.grid(month = seq(0, 60, 12), series = levels(km$series))

p <- ggplot(km, aes(month, survival, colour = group, linetype = group)) +
  geom_hline(yintercept = seq(0, 1, by = 0.25), colour = "grey90", linewidth = 0.25) +
  geom_step(linewidth = 0.42) +
  geom_text(data = ar, aes(month, y, label = n, colour = group),
            inherit.aes = FALSE, size = 8.5 / .pt, show.legend = FALSE) +
  geom_text(data = group_labels, aes(x = -3, y, label = substr(group, 1, 1), colour = group),
            hjust = 1, inherit.aes = FALSE, size = 8.5 / .pt, show.legend = FALSE) +
  geom_text(data = tick_labels, aes(month, y = -0.075, label = month),
            inherit.aes = FALSE, size = 8.5 / .pt) +
  facet_wrap(vars(series), ncol = 2) +
  scale_colour_manual(name = "Risk group", values = colors, breaks = groups) +
  scale_linetype_manual(name = "Risk group", values = patterns, breaks = groups) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_cartesian(xlim = c(-8, 62), ylim = c(-0.36, 1.03), expand = FALSE) +
  labs(x = "Months", y = "Survival probability") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 3)), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.y = element_text(size = 8), legend.text = element_text(size = 8.5),
        legend.title = element_text(size = 8.5),
        panel.spacing = grid::unit(9, "pt"))
# PDF_DEVICE is cairo_pdf (embeds fonts) when available; all figure text is at least 8 pt at the printed width (6.5 in)
fig_name <- paste0("sup_external_km_fixed", exhibit_suffix(), ".pdf")
ggsave(file.path(figure_dir, fig_name), p, width = 6.5, height = 6.6, device = PDF_DEVICE)
cat("Wrote", fig_name, "\n")
