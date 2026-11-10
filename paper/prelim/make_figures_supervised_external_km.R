# ============================================================
# Script: make_figures_supervised_external_km.R
# Purpose: Summarize the fixed S-12 scores in survival curves and per-program hazard ratios.
# Created: 2026-11-06
# Dependencies: survival, ggplot2, code/concordance_ci.R
# Run from the repository root with Rscript.
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})
source("code/concordance_ci.R")

## cairo_pdf embeds fonts but needs the cairo shared library (it fails to load without XQuartz's X11 libs).
## capabilities("cairo") can report TRUE when the library cannot load, so test by opening a device.
CAIRO_OK <- tryCatch({ f <- tempfile(fileext = ".pdf"); suppressWarnings(cairo_pdf(f)); dev.off(); file.exists(f) },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (!CAIRO_OK) warning("cairo_pdf unavailable: figures are written with pdf() and their fonts are NOT embedded",
                       call. = FALSE)
PDF_DEVICE <- if (CAIRO_OK) cairo_pdf else function(filename, ...) grDevices::pdf(file = filename, ..., useDingbats = FALSE)

# Saved evidence and score validation ----
s11_dir <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04"
s12_dir <- "results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05"
out_dir <- "results/benchmark_sim/outputs/prelim_ch4_2026-11-06"
figure_dir <- "paper/prelim/figures"
table_dir <- "paper/prelim/tables"
score_file <- "data/yfb_external_patient_scores_ws1.rds"
if (!file.exists(score_file)) {
  stop("Missing ignored score cache. Run results/benchmark_sim/save_yfb_external_patient_scores.R.")
}

cohort_names <- c(Dijk = "Dijk", Moffitt_GEO_array = "Moffitt",
                  PACA_AU_array = "PACA-AU array", PACA_AU_seq = "PACA-AU RNA-seq",
                  Puleo_array = "Puleo")
scores <- readRDS(score_file)
scores_by_cohort <- split(scores, scores$cohort)
ref <- read.csv(file.path(s12_dir, "external_cindex.csv"))
programs <- read.csv(file.path(s11_dir, "programs_D4.csv"))
programs <- subset(programs, arm == "w1" & survival_active)
fit <- readRDS(file.path(s11_dir, "fits", "size_D4_N2218_w1_f0.rds"))
beta <- fit$state$EBeta
Z <- as.matrix(scores[paste0("Z", seq_along(beta))])

# Score direction is fixed by the fitted model. Check the projection against the
# fitted coefficients, and check both risk scores against the registered results.
stopifnot(setequal(names(scores_by_cohort), names(cohort_names)),
          all(is.finite(Z)), all(is.finite(scores$risk)), all(is.finite(scores$risk_twostep)),
          all(is.finite(scores$time) & scores$time >= 0), all(scores$event %in% c(0, 1)),
          all(abs(drop(Z %*% beta) - scores$risk) < 1e-10),
          all(abs(beta[programs$program] - programs$beta) < 1e-10))
for (i in seq_len(nrow(ref))) {
  d <- scores_by_cohort[[ref$cohort[i]]]
  risk <- if (ref$model[i] == "w1") d$risk else d$risk_twostep
  stopifnot(nrow(d) == ref$n[i], sum(d$event) == ref$events[i],
            abs(frozen_reverse_cindex(risk, d$time, d$event) - ref$C[i]) < 1e-10)
}

# Kaplan–Meier figures ----
months <- seq(0, 60, by = 12)
groups <- c("Low", "Middle", "High")
colors <- c(Low = "#006878", Middle = "#555555", High = "#B34226")
patterns <- c(Low = "solid", Middle = "dashed", High = "dotdash")

#' Draw within-series Kaplan–Meier curves and risk counts for one frozen score.
#' @param score_name Column name in the ignored external score cache.
#' @param direction Sign from training, applied before defining risk tertiles.
#' @param filename PDF basename in the Chapter 4 figures directory.
#' @return Invisibly returns the ggplot, after writing a vector PDF.
#' @examples
#' # plot_tertiles("risk", 1, "sup_external_km_risk.pdf")
plot_tertiles <- function(score_name, direction, filename) {
  curves <- list()
  risk_counts <- list()
  for (cohort in names(cohort_names)) {
    d <- scores_by_cohort[[cohort]]
    signed_score <- direction * d[[score_name]]
    cutoffs <- quantile(signed_score, c(1/3, 2/3), names = FALSE)
    stopifnot(cutoffs[1] < cutoffs[2])
    d$tertile <- cut(signed_score, c(-Inf, cutoffs, Inf), labels = groups, include.lowest = TRUE)
    stopifnot(all(groups %in% as.character(d$tertile)))

    for (group in groups) {
      g <- d[d$tertile == group, , drop = FALSE]
      km <- survfit(Surv(time, event) ~ 1, data = g)
      key <- paste(cohort, group)
      last_month <- max(months)
      in_window <- km$time <= last_month
      plot_time <- c(0, km$time[in_window])
      plot_surv <- c(1, km$surv[in_window])
      if (max(km$time) > last_month && tail(plot_time, 1) < last_month) {
        plot_time <- c(plot_time, last_month)
        plot_surv <- c(plot_surv, tail(plot_surv, 1))
      }
      curves[[key]] <- data.frame(cohort, group, time = plot_time, surv = plot_surv)
      counts <- vapply(months, function(t) sum(g$time >= t), integer(1))
      risk_counts[[key]] <- data.frame(cohort, group, time = months, n = counts,
                                      y = -0.20 - 0.12 * (match(group, groups) - 1))
    }
  }

  curves <- do.call(rbind, curves)
  risk_counts <- do.call(rbind, risk_counts)
  curves$series <- factor(cohort_names[curves$cohort], levels = unname(cohort_names))
  risk_counts$series <- factor(cohort_names[risk_counts$cohort], levels = unname(cohort_names))
  group_labels <- subset(risk_counts, time == 0)
  tick_labels <- expand.grid(time = months, series = factor(cohort_names, levels = cohort_names))

  p <- ggplot(curves, aes(time, surv, colour = group, linetype = group)) +
    geom_hline(yintercept = seq(0, 1, by = 0.25), colour = "grey90", linewidth = 0.25) +
    geom_step(linewidth = 0.42) +
    geom_text(data = risk_counts, aes(time, y, label = n, colour = group),
              inherit.aes = FALSE, size = 8.5 / .pt, show.legend = FALSE) +
    geom_text(data = group_labels, aes(x = -3, y, label = substr(group, 1, 1), colour = group),
              hjust = 1, inherit.aes = FALSE, size = 8.5 / .pt, show.legend = FALSE) +
    geom_text(data = tick_labels, aes(time, y = -0.075, label = time),
              inherit.aes = FALSE, size = 8.5 / .pt) +
    facet_wrap(vars(series), ncol = 2) +
    scale_colour_manual(name = "Risk tertile", values = colors, breaks = groups) +
    scale_linetype_manual(name = "Risk tertile", values = patterns, breaks = groups) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    coord_cartesian(xlim = c(-8, 62), ylim = c(-0.48, 1.03), expand = FALSE) +
    labs(x = "Months", y = "Survival probability") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.title.x = element_text(margin = margin(t = 3)), legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 9),
          axis.text.y = element_text(size = 8), legend.text = element_text(size = 8.5),
          legend.title = element_text(size = 8.5),
          panel.spacing = grid::unit(9, "pt"))
  # PDF_DEVICE is cairo_pdf (embeds fonts) when available; all figure text is at least 8 pt at the printed width (6.5 in)
  ggsave(file.path(figure_dir, filename), p, width = 6.5, height = 6.6, device = PDF_DEVICE)
  invisible(p)
}

plot_tertiles("risk", 1, "sup_external_km_risk.pdf")
for (k in programs$program) {
  plot_tertiles(paste0("Z", k), sign(beta[k]), sprintf("sup_external_km_program_%d.pdf", k))
}

# Cohort-specific program direction ----
rows <- list()
for (cohort in names(cohort_names)) {
  d <- scores_by_cohort[[cohort]]
  for (k in programs$program) {
    # Per-SD associations use the program's original orientation, so their
    # signs can be compared directly to the training coefficient.
    d$score_sd <- as.numeric(scale(d[[paste0("Z", k)]]))
    cox <- summary(coxph(Surv(time, event) ~ score_sd, data = d))
    ci <- cox$conf.int[1, ]
    stopifnot(all(is.finite(ci)), cox$n == nrow(d), cox$nevent == sum(d$event))
    rows[[paste(cohort, k)]] <- data.frame(cohort, program = k, n = cox$n, events = cox$nevent,
      training_sign = sign(beta[k]), log_hr = unname(cox$coefficients[1, "coef"]),
      hr = unname(ci["exp(coef)"]), lo = unname(ci["lower .95"]), hi = unname(ci["upper .95"]))
  }
}
direction <- do.call(rbind, rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(direction, file.path(out_dir, "external_program_direction.csv"), row.names = FALSE)

header <- paste(sprintf("Program %d", programs$program), collapse = " & ")
table_rows <- character()
for (cohort in names(cohort_names)) {
  d <- direction[direction$cohort == cohort, ]
  values <- sprintf("%.2f (%.2f to %.2f)", d$hr, d$lo, d$hi)
  table_rows <- c(table_rows, paste(c(cohort_names[[cohort]], values), collapse = " & "))
}
sign_words <- ifelse(beta[programs$program] > 0, "Adverse", "Protective")
training_row <- paste(c("Training direction", sign_words), collapse = " & ")
table_path <- file.path(table_dir, "tab_external_program_direction.tex")
writeLines(c("\\begin{table}[!htb]", "\\centering",
  "\\caption{External associations of the survival-active programs. Values are hazard ratios and 95\\% Wald intervals per cohort-specific standard deviation of each program score from separate Cox fits. The training direction comes from the fitted coefficient. These are descriptive, unadjusted associations.}",
  "\\label{tab:external-program-direction}", "\\begin{tabular}{lrrr}", "\\toprule",
  paste("Series &", header, "\\\\"), "\\midrule", paste(table_rows, "\\\\"),
  "\\midrule", paste(training_row, "\\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  table_path)

print(with(direction, table(program, agrees_with_training = sign(log_hr) == training_sign)))
