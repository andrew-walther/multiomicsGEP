# ============================================================
# Script: paper/prelim/make_figures_supervised_sim.R
# Purpose: Figures and LaTeX tables for the prelim chapter's survival-supervision results
#          (simulation studies S-9, S-13, S-14 and the real-data training/external results),
#          read from saved aggregate CSVs. Nothing is refitted.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: ggplot2, patchwork (figures); base R for tables
# Run (from the repository root): Rscript paper/prelim/make_figures_supervised_sim.R
# Outputs: paper/prelim/figures/sup_*.pdf and paper/prelim/tables/tab_*.tex
# ============================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# Paths and constants ----
OUT <- "results/benchmark_sim/outputs"
S9_DIR  <- file.path(OUT, "yfb_survival_weight_2026-11-04/grid")
S13_DIR <- file.path(OUT, "yfb_supervision_strength_2026-11-05")
S14_DIR <- file.path(OUT, "yfb_near_threshold_2026-11-05")

## Real-data inputs: the only section that changes if the real-data model is refitted (S-15).
## To switch, change these two folder variables (file names are assumed unchanged).
REAL_DIR     <- file.path(OUT, "yfb_supervised_real_2026-11-04")
EXTERNAL_DIR <- file.path(OUT, "yfb_supervised_external_2026-11-05")

FIG_DIR <- "paper/prelim/figures"
TAB_DIR <- "paper/prelim/tables"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

COL_SUP   <- "#0072B2"   # supervised (w_S = 1 or A1)
COL_TWO   <- "#E69F00"   # two-step
COL_OTHER <- "#8C8C8C"   # other weights / arms
COL_CELL  <- c("0.25%" = "#9E9AC8", "1%" = "#6A51A3", "5%" = "#3F007D")  # planted variance
W_ZERO_X  <- 0.3         # x position used for w_S = 0 on the log axis
FIG_W     <- 6.5         # printed width in inches

# Theme and helpers ----
#' Common theme: theme_bw with all text at least 8 pt at printed size.
#' @return A ggplot2 theme object.
theme_prelim <- function() {
  theme_bw(base_size = 9) +
    theme(axis.text = element_text(size = 8, colour = "grey15"),
          axis.title = element_text(size = 9),
          strip.text = element_text(size = 8.5),
          strip.background = element_rect(fill = "grey93", colour = "grey60"),
          legend.text = element_text(size = 8.5), legend.title = element_text(size = 8.5),
          legend.key.size = unit(0.9, "lines"),
          panel.grid.minor = element_blank(), plot.margin = margin(3, 5, 3, 3))
}

## cairo_pdf needs the cairo shared library, which fails to load on some installs (missing X11 libs)
CAIRO_OK <- tryCatch({ f <- tempfile(fileext = ".pdf"); suppressWarnings(cairo_pdf(f)); dev.off(); file.exists(f) },
                     error = function(e) FALSE, warning = function(w) FALSE)
#' Save a figure as vector PDF at printed size.
#' @param p ggplot/patchwork object. @param name file name. @param h height in inches.
save_fig <- function(p, name, h) {
  # cairo_pdf where the cairo library loads; plain pdf() (also vector, base-14 fonts) otherwise
  dev <- if (CAIRO_OK) cairo_pdf else grDevices::pdf
  ggsave(file.path(FIG_DIR, name), p, width = FIG_W, height = h, device = dev)
  cat("Wrote", file.path(FIG_DIR, name), sprintf("(%.1f x %.1f in)\n", FIG_W, h))
}

#' Compare a computed value with the number quoted in a findings file and print the verdict.
#' Expected numbers are typed only as checks; every plotted/tabulated value comes from CSVs.
#' @param label description. @param got value from CSV. @param expect value in findings.
#' @param tol absolute tolerance (default half a unit in the last quoted digit).
MISMATCHES <- character(0)
check <- function(label, got, expect, tol = 0.0051) {
  ok <- isTRUE(abs(got - expect) <= tol)
  if (!ok) MISMATCHES <<- c(MISMATCHES, sprintf("%s: csv %.4f vs findings %.4f", label, got, expect))
  cat(sprintf("  [%s] %-58s csv %9.4f | findings %9.4f\n", if (ok) "ok" else "MISMATCH", label, got, expect))
}

# F1: S-9 against the survival weight (p = 2,000) ----
#' Plot the fixed-weight S-9 arms at p = 2,000.
#' Reads saved metrics and calibration slopes recomputed with exact event times.
#' @return Invisibly, the result of writing a vector PDF.
#' @examples
#' \dontrun{make_f1()}
make_f1 <- function() {
  m  <- read.csv(file.path(S9_DIR, "metrics.csv"))
  cs <- read.csv(file.path(S9_DIR, "cal_slope_exact.csv"))
  fit_key <- c("seed", "p", "variance", "arm")
  stopifnot(!anyDuplicated(m[fit_key]), !anyDuplicated(cs[fit_key]))
  stopifnot(nrow(read.csv(file.path(S9_DIR, "design.csv"))) == 60)
  m <- merge(m, cs[, c(fit_key, "cal_slope_exact")], by = fit_key, all.x = TRUE)
  # The adaptive share-target arm has a different weight per dataset.
  m <- m[m$p == 2000 & m$arm != "share0.2", ]
  stopifnot(nrow(m) == 210, all(is.finite(m$cal_slope_exact[m$w_S > 0])),
            all(is.finite(m$c_validation[m$w_S > 0])))
  agg <- do.call(rbind, lapply(split(m, list(m$variance, m$w_S), drop = TRUE), function(d) {
    data.frame(variance = d$variance[1], w_S = d$w_S[1], n = nrow(d),
               recovered = sum(d$recovered), fa = mean(d$false_active), conv = sum(d$converged),
               share = median(d$share_match),
               # At w_S = 0 the predictor is constant and these quantities are undefined.
               slope = median(d$cal_slope_exact), C = mean(d$c_validation))
  }))
  stopifnot(all(agg$n == 10), all(is.finite(agg$share)))
  agg$cell <- factor(c("0.0025" = "0.25%", "0.01" = "1%", "0.05" = "5%")[as.character(agg$variance)], levels = names(COL_CELL))
  agg$x <- ifelse(agg$w_S == 0, W_ZERO_X, agg$w_S)
  cat("\nF1 (S-9, p = 2,000): key values\n")
  print(agg[order(agg$cell, agg$w_S), c("cell", "w_S", "n", "recovered", "fa", "conv", "slope", "share", "C")],
        row.names = FALSE, digits = 3)
  g <- function(v, w, col) agg[agg$variance == v & agg$w_S == w, col]
  for (v in c(0.0025, 0.01, 0.05)) check(sprintf("recovered at w_S=1, v=%g (10/10)", v), g(v, 1, "recovered"), 10, 0)
  for (i in 1:3) {
    v <- c(0.0025, 0.01, 0.05)[i]
    check(sprintf("mean FA at w_S=1, v=%g", v), g(v, 1, "fa"), c(0.9, 0.1, 0)[i])
    check(sprintf("slope at w_S=1, v=%g", v), g(v, 1, "slope"), c(1.45, 2.16, 3.23)[i])
    check(sprintf("median matched survival share at w_S=1, v=%g", v), g(v, 1, "share"), c(0.202, 0.049, 0.010)[i], 0.0006)
    check(sprintf("mean validation C at w_S=1, v=%g", v), g(v, 1, "C"), c(0.944, 0.972, 0.986)[i], 0.0006)
    check(sprintf("slope at w_S=10, v=%g", v), g(v, 10, "slope"), c(0.26, 0.71, 1.74)[i])
    check(sprintf("slope at w_S=100, v=%g", v), g(v, 100, "slope"), c(0.12, 0.17, 0.16)[i])
    check(sprintf("slope at w_S=10000, v=%g", v), g(v, 10000, "slope"), c(0.07, 0.08, 0.10)[i])
    check(sprintf("recovered at w_S=40, v=%g", v), g(v, 40, "recovered"), c(3, 4, 10)[i], 0)
    check(sprintf("recovered at w_S=100, v=%g", v), g(v, 100, "recovered"), c(0, 0, 8)[i], 0)
    check(sprintf("FA at w_S=10, v=%g", v), g(v, 10, "fa"), c(2.2, 1.8, 0.8)[i])
    check(sprintf("converged at w_S>=40 is 0, v=%g", v), sum(agg$conv[agg$variance == v & agg$w_S >= 40]), 0, 0)
  }
  brk <- c(W_ZERO_X, 1, 10, 40, 100, 1000, 10000); lab <- c("0", "1", "10", "40", "100", "1000", "10000")
  #' Draw a panel with the same weight and variance encoding as its neighbors.
  #' @param ylab Axis title.
  #' @param ycol Summary column to plot.
  #' @return A ggplot object.
  #' @examples
  #' \dontrun{base("Mean validation C", "C")}
  base <- function(ylab, ycol) {
    ggplot(agg, aes(x, .data[[ycol]], colour = cell, shape = cell, group = cell)) +
      geom_vline(xintercept = 0.55, linetype = "dotted", colour = "grey55") +
      geom_line(data = agg[agg$w_S > 0, ], linewidth = 0.6) + geom_point(size = 1.9, na.rm = TRUE) +
      scale_x_log10(breaks = brk, labels = lab, limits = c(0.22, 14000)) +
      scale_colour_manual(values = COL_CELL, name = "Planted variance") +
      scale_shape_manual(values = c(16, 17, 15), name = "Planted variance") +
      labs(x = expression(paste("Survival weight ", w[S])), y = ylab) + theme_prelim()
  }
  pA <- base("Planted program recovered (of 10)", "recovered") + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2))
  pB <- base("Mean false-active programs", "fa")
  pC <- base("Median calibration slope", "slope") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    scale_y_log10(breaks = c(0.1, 0.3, 1, 3), limits = c(0.05, 4))
  pD <- base("Fits converged (of 10)", "conv") + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2))
  pE <- base("Median survival precision share", "share")
  pF <- base("Mean validation C", "C")
  p <- (pA + pB + pC + pD + pE + pF) + plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", plot.tag = element_text(size = 10, face = "bold"))
  save_fig(p, "sup_sim_weight_s9.pdf", 7.0)
}

# F12: S-13, survival-block weight (arm A1) ----
## Data: weight_curves.csv (per-weight summaries, exact-times slope), cv_selected.csv, twostep.csv.
make_f12 <- function() {
  wc <- read.csv(file.path(S13_DIR, "weight_curves.csv")); wc <- wc[wc$arm == "A1", ]
  sel <- read.csv(file.path(S13_DIR, "cv_selected.csv")); sel <- sel[sel$arm == "A1" & !sel$null, ]
  ts <- read.csv(file.path(S13_DIR, "twostep.csv"))
  ws <- c(1, 2, 5, 10, 20, 50)
  cells <- c("v0.0025", "v0.01", "v0.05", "weak")
  cell_lab <- c("0.25% variance", "1% variance", "5% variance", "Weak survival")
  names(cell_lab) <- cells
  wc$cell <- factor(cell_lab[wc$cell], levels = cell_lab)
  ref <- do.call(rbind, lapply(cells, function(cc) { d <- ts[ts$cell == cc, ]
    data.frame(cell = factor(cell_lab[cc], levels = cell_lab), c_validation = mean(d$c_validation),
               slope_median = median(d$cal_slope), converged = sum(d$converged)) }))
  cnt <- as.data.frame(table(cell = factor(cell_lab[sel$cell], levels = cell_lab), w_S = factor(sel$w_S, levels = ws)))
  cnt$w_S <- as.numeric(as.character(cnt$w_S))
  cat("\nF12 (S-13, A1): key values\n")
  print(wc[, c("cell", "w_S", "recovered", "converged", "slope_median", "c_validation")], row.names = FALSE, digits = 3)
  print(ref, row.names = FALSE, digits = 3); print(cnt[cnt$Freq > 0, ], row.names = FALSE)
  # Checks against the selected-weight table in the findings
  sv <- read.csv(file.path(S13_DIR, "selected_vs_twostep.csv")); sv <- sv[sv$arm == "A1", ]
  for (i in seq_along(cells)) {
    cc <- cells[i]; r <- ref[i, ]
    check(sprintf("two-step C, %s", cc), r$c_validation, c(0.939, 0.972, 0.986, 0.668)[i])
    check(sprintf("two-step slope, %s", cc), r$slope_median, c(1.16, 1.10, 0.96, 1.04)[i])
    check(sprintf("two-step slope matches selected_vs_twostep, %s", cc), r$slope_median, sv$twostep_slope[sv$cell == cc], 1e-9)
    check(sprintf("two-step C matches selected_vs_twostep, %s", cc), r$c_validation, sv$twostep_c[sv$cell == cc], 1e-9)
  }
  g <- function(cc, w, col) wc[wc$cell == cell_lab[cc] & wc$w_S == w, col]
  check("A1 0.25% cell slope at w_S=1", g("v0.0025", 1, "slope_median"), 1.45)
  check("A1 5% cell slope at w_S=20 (1.03, not converged)", g("v0.05", 20, "slope_median"), 1.03)
  check("A1 5% cell converged at w_S=20", g("v0.05", 20, "converged"), 0, 0)
  check("selected w_S=10 in 5% cell (count)", cnt$Freq[cnt$cell == cell_lab["v0.05"] & cnt$w_S == 10], 10, 0)
  check("selected w_S=5 in 1% cell (count)", cnt$Freq[cnt$cell == cell_lab["v0.01"] & cnt$w_S == 5], 10, 0)
  check("selected w_S=2 in 0.25% cell (count)", cnt$Freq[cnt$cell == cell_lab["v0.0025"] & cnt$w_S == 2], 9, 0)
  check("weights 20/50 never selected", sum(cnt$Freq[cnt$w_S >= 20]), 0, 0)

  brk <- ws; xs <- scale_x_log10(breaks = brk, labels = brk)
  leg <- scale_colour_manual(values = c("Supervised, block-weighted fit" = COL_SUP, "Two-step (same data)" = COL_TWO), name = NULL)
  row <- function(ycol, refcol, ylab, top = FALSE, bottom = FALSE, logy = FALSE, ylim = NULL, h1 = FALSE) {
    p <- ggplot(wc, aes(w_S, .data[[ycol]])) +
      geom_hline(data = ref, aes(yintercept = .data[[refcol]], colour = "Two-step (same data)"), linewidth = 0.6) +
      (if (h1) geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.35)) +
      geom_line(aes(colour = "Supervised, block-weighted fit"), linewidth = 0.6) + geom_point(aes(colour = "Supervised, block-weighted fit"), size = 1.6) +
      facet_wrap(~cell, nrow = 1, scales = if (logy || !is.null(ylim)) "fixed" else "free_y") + xs + leg + labs(y = ylab) + theme_prelim()
    if (logy) p <- p + scale_y_log10(breaks = c(0.03, 0.1, 0.3, 1, 3), limits = c(0.012, 4))
    if (!is.null(ylim)) p <- p + scale_y_continuous(limits = ylim, breaks = seq(0, 10, 5))
    p + labs(x = if (bottom) expression(paste("Survival weight ", w[S])) else NULL) +
      theme(axis.text.x = if (bottom) element_text() else element_blank(),
            axis.ticks.x = if (bottom) element_line() else element_blank(),
            strip.text = if (top) element_text(size = 8.5) else element_blank(),
            strip.background = if (top) element_rect(fill = "grey93", colour = "grey60") else element_blank())
  }
  p1 <- row("c_validation", "c_validation", "Validation C", top = TRUE)
  p2 <- row("slope_median", "slope_median", "Median calib. slope", logy = TRUE, h1 = TRUE)
  p3 <- row("converged", "converged", "Converged (of 10)", ylim = c(0, 10.5))
  p4 <- ggplot(cnt, aes(w_S, Freq)) +
    geom_col(fill = COL_SUP, alpha = 0.55, width = 0.12, colour = NA) + facet_wrap(~cell, nrow = 1) +
    scale_x_log10(breaks = brk, labels = brk) + scale_y_continuous(limits = c(0, 10.5), breaks = seq(0, 10, 5)) +
    labs(x = expression(paste("Survival weight ", w[S])), y = "Chosen by CV (of 10)") + theme_prelim() +
    theme(strip.text = element_blank(), strip.background = element_blank())
  p <- (p1 / p2 / p3 / p4) + plot_layout(guides = "collect", heights = c(1, 1, 0.8, 0.8)) &
    theme(legend.position = "bottom")
  save_fig(p, "sup_sim_weight_s13.pdf", 7.0)
}

# F13: S-14, near-threshold planted programs ----
#' Plot S-14 recovery counts and paired validation C differences.
#' Reads saved metrics and pairs A1 w_S = 1 with two-step by dataset.
#' @return Invisibly, the result of writing a vector PDF.
#' @examples
#' \dontrun{make_f13()}
make_f13 <- function() {
  m <- read.csv(file.path(S14_DIR, "metrics.csv"))
  arms <- c(twostep = "Two-step", A1_w1 = "Supervised, w_S = 1", A1_w5 = "Supervised, w_S = 5")
  m <- m[m$arm %in% names(arms), ]
  m$share_lab <- factor(paste0(m$share * 100, "%"), levels = c("0.05%", "0.1%", "0.2%"))
  m$strength_lab <- factor(c(strong = "Strong survival", realistic = "Realistic survival")[m$strength],
                           levels = c("Strong survival", "Realistic survival"))
  m$armf <- factor(arms[m$arm], levels = arms)
  agg <- do.call(rbind, lapply(split(m, list(m$share_lab, m$strength_lab, m$armf), drop = TRUE), function(d)
    data.frame(share_lab = d$share_lab[1], strength_lab = d$strength_lab[1], armf = d$armf[1],
               recovered = sum(d$recovered), C = mean(d$c_validation), oracle = mean(d$oracle_c), n = nrow(d))))
  stopifnot(all(agg$n == 10))
  # Paired C difference (A1 w_S = 1 minus two-step), computed per dataset
  wd <- merge(m[m$arm == "A1_w1", c("seed", "share", "strength", "c_validation")],
              m[m$arm == "twostep", c("seed", "share", "strength", "c_validation")], by = c("seed", "share", "strength"))
  stopifnot(nrow(wd) == 60, !anyDuplicated(wd[c("seed", "share", "strength")]),
            all(is.finite(wd$c_validation.x)), all(is.finite(wd$c_validation.y)))
  wd$dC <- wd$c_validation.x - wd$c_validation.y
  dC <- aggregate(dC ~ share + strength, wd, mean)
  cat("\nF13 (S-14): key values\n"); print(agg[order(agg$strength_lab, agg$share_lab, agg$armf), ], row.names = FALSE, digits = 3)
  cat("Paired C difference, A1 w_S=1 minus two-step\n"); print(dC, row.names = FALSE, digits = 3)
  gg <- function(s, sh, a, col) agg[agg$strength_lab == s & agg$share_lab == sh & agg$armf == a, col]
  S <- "Strong survival"; R <- "Realistic survival"
  check("strong 0.05% recovered two-step", gg(S, "0.05%", "Two-step", "recovered"), 0, 0)
  check("strong 0.05% recovered A1 w1", gg(S, "0.05%", "Supervised, w_S = 1", "recovered"), 5, 0)
  check("strong 0.05% recovered A1 w5", gg(S, "0.05%", "Supervised, w_S = 5", "recovered"), 9, 0)
  check("strong 0.1% recovered two-step", gg(S, "0.1%", "Two-step", "recovered"), 8, 0)
  check("strong 0.1% recovered A1 w1", gg(S, "0.1%", "Supervised, w_S = 1", "recovered"), 10, 0)
  check("realistic 0.1% recovered A1 w5", gg(R, "0.1%", "Supervised, w_S = 5", "recovered"), 0, 0)
  check("realistic 0.1% recovered A1 w1", gg(R, "0.1%", "Supervised, w_S = 1", "recovered"), 9, 0)
  check("realistic 0.2% recovered A1 w5", gg(R, "0.2%", "Supervised, w_S = 5", "recovered"), 7, 0)
  check("strong 0.05% C two-step", gg(S, "0.05%", "Two-step", "C"), 0.644, 0.0006)
  check("strong 0.05% C A1 w1", gg(S, "0.05%", "Supervised, w_S = 1", "C"), 0.800, 0.0006)
  check("strong 0.05% oracle C", gg(S, "0.05%", "Two-step", "oracle"), 0.936, 0.0006)
  check("realistic 0.2% C A1 w1", gg(R, "0.2%", "Supervised, w_S = 1", "C"), 0.728, 0.0006)
  expect_dC <- list(strong = c(0.156, 0.041, 0.008), realistic = c(0.038, 0.008, 0.002))
  for (s in names(expect_dC)) for (i in 1:3) {
    sh <- c(5e-4, 1e-3, 2e-3)[i]
    check(sprintf("paired dC %s share %g", s, sh * 100), dC$dC[dC$strength == s & abs(dC$share - sh) < 1e-9], expect_dC[[s]][i], 0.0006)
  }
  cols <- c(COL_TWO, COL_SUP, COL_OTHER); names(cols) <- arms
  pa <- ggplot(agg, aes(share_lab, recovered, fill = armf)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = recovered), position = position_dodge(0.8), vjust = -0.3, size = 2.9) +
    facet_wrap(~strength_lab) + scale_fill_manual(values = cols, name = NULL) +
    scale_y_continuous(limits = c(0, 11.5), breaks = seq(0, 10, 2)) +
    labs(x = NULL, y = "Planted program recovered (of 10)") + theme_prelim()
  pb <- ggplot(wd, aes(factor(share * 100, labels = levels(m$share_lab)), dC)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1),
               size = 1.6, colour = COL_SUP, alpha = 0.65) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.8) +
    facet_wrap(~factor(strength, levels = c("strong", "realistic"), labels = levels(m$strength_lab))) +
    labs(x = "Planted share of variance", y = "Paired validation C difference") + theme_prelim()
  p <- (pa / pb) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", plot.tag = element_text(size = 10, face = "bold"))
  save_fig(p, "sup_sim_near_threshold.pdf", 5.4)
}

# Real-data section (separate; depends on REAL_DIR / EXTERNAL_DIR) ----
## F3: training CV C against the number of genes.
make_f3 <- function() {
  s <- read.csv(file.path(REAL_DIR, "summary.csv")); s <- s[s$part == "size" & s$arm %in% c("w1", "w0cox"), ]
  pref <- read.csv(file.path(REAL_DIR, "size_preferred.csv"))
  d <- data.frame(genes = s$n_genes, C = s$mean_C, se = s$se_C,
                  arm = factor(c(w1 = "Supervised, w_S = 1", w0cox = "Two-step")[s$arm], levels = c("Supervised, w_S = 1", "Two-step")))
  d <- d[order(d$arm, d$genes), ]
  cat("\nF3 (real-data size curve): key values\n"); print(d, row.names = FALSE, digits = 3)
  cat("Preferred sizes:\n"); print(pref[, c("arm", "preferred_vs_all_genes")], row.names = FALSE)
  g <- function(n, a, col) d[d$genes == n & d$arm == a, col]
  S1 <- "Supervised, w_S = 1"; T2 <- "Two-step"
  exp1 <- c(`2218` = 0.663, `3000` = 0.663, `5000` = 0.657, `8001` = 0.647, `12001` = 0.643, `17900` = 0.648)
  exp2 <- c(`2218` = 0.652, `3000` = 0.653, `5000` = 0.652, `8001` = 0.650, `12001` = 0.653, `17900` = 0.652)
  se1  <- c(0.024, 0.026, 0.030, 0.030, 0.037, 0.038); se2 <- c(0.026, 0.027, 0.031, 0.033, 0.037, 0.037)
  for (i in seq_along(exp1)) {
    n <- as.numeric(names(exp1)[i])
    check(sprintf("C w_S=1, %d genes", n), g(n, S1, "C"), exp1[i], 0.0006); check(sprintf("SE w_S=1, %d genes", n), g(n, S1, "se"), se1[i], 0.0006)
    check(sprintf("C two-step, %d genes", n), g(n, T2, "C"), exp2[i], 0.0006); check(sprintf("SE two-step, %d genes", n), g(n, T2, "se"), se2[i], 0.0006)
  }
  p <- ggplot(d, aes(genes, C, colour = arm, shape = arm)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey55") +
    geom_line(linewidth = 0.6, position = position_dodge(0.03)) +
    geom_errorbar(aes(ymin = C - se, ymax = C + se), width = 0.04, position = position_dodge(0.03), linewidth = 0.4) +
    geom_point(size = 2, position = position_dodge(0.03)) +
    scale_x_log10(breaks = c(2218, 3000, 5000, 8001, 12001, 17900), labels = c("2,218", "3,000", "5,000", "8,001", "12,001", "17,900")) +
    scale_colour_manual(values = c(COL_SUP, COL_TWO), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
    scale_y_continuous(limits = c(0.5, 0.72), breaks = seq(0.5, 0.7, 0.05)) +
    labs(x = "Number of genes (log scale)", y = "Training CV concordance") + theme_prelim() +
    theme(legend.position = "bottom")
  save_fig(p, "sup_real_size_curve.pdf", 3.3)
  invisible(d)
}

## F4 + F9: external C and calibration slope by cohort.
COHORT_LAB <- c(Dijk = "Dijk", Moffitt_GEO_array = "Moffitt (array)", PACA_AU_array = "PACA-AU (array)",
                PACA_AU_seq = "PACA-AU (RNA-seq)", Puleo_array = "Puleo (array)")
#' Plot concordance and calibration for the two fixed fits in each external series.
#' Reads bootstrapped C bounds and constructs Wald intervals for calibration.
#' @return Invisibly, the result of writing a vector PDF.
#' @examples
#' \dontrun{make_ext()}
make_ext <- function() {
  e <- read.csv(file.path(EXTERNAL_DIR, "external_cindex.csv"))
  e$cohort_lab <- factor(COHORT_LAB[e$cohort], levels = rev(COHORT_LAB))
  e$model_lab <- factor(c(w1 = "Supervised, w_S = 1", twostep = "Two-step")[e$model], levels = c("Supervised, w_S = 1", "Two-step"))
  cat("\nF4/F9 (external): key values\n"); print(e[, c("cohort", "model", "n", "events", "C", "lo", "hi", "cal_slope", "cal_se")], row.names = FALSE, digits = 3)
  ex <- list(Dijk = c(0.584, 0.567, 0.55, 0.36, 0.35), Moffitt_GEO_array = c(0.537, 0.531, 0.36, 0.36, 0.28),
             PACA_AU_array = c(0.656, 0.656, 0.78, 0.33, 0.73), PACA_AU_seq = c(0.681, 0.671, 1.57, 0.43, 1.26),
             Puleo_array = c(0.644, 0.643, 0.99, 0.20, 0.80))
  for (k in names(ex)) {
    a <- e[e$cohort == k & e$model == "w1", ]; b <- e[e$cohort == k & e$model == "twostep", ]
    check(paste(k, "C w1"), a$C, ex[[k]][1], 0.0006); check(paste(k, "C two-step"), b$C, ex[[k]][2], 0.0006)
    check(paste(k, "slope w1"), a$cal_slope, ex[[k]][3]); check(paste(k, "slope SE w1"), a$cal_se, ex[[k]][4])
    check(paste(k, "slope two-step"), b$cal_slope, ex[[k]][5])
  }
  labels <- e[c("cohort_lab", "model_lab")]
  cc <- rbind(data.frame(labels, panel = "Concordance C", x = e$C, lo = e$lo, hi = e$hi),
              data.frame(labels, panel = "Calibration slope", x = e$cal_slope, lo = e$cal_slope - 1.96 * e$cal_se,
                         hi = e$cal_slope + 1.96 * e$cal_se))
  cc$panel <- factor(cc$panel, levels = c("Concordance C", "Calibration slope"))
  refl <- data.frame(panel = factor(c("Concordance C", "Calibration slope"), levels = levels(cc$panel)), v = c(0.5, 1))
  p <- ggplot(cc, aes(x, cohort_lab, colour = model_lab, shape = model_lab)) +
    geom_vline(data = refl, aes(xintercept = v), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.55, position = position_dodge(0.5), na.rm = TRUE) +
    geom_point(size = 2.2, position = position_dodge(0.5)) +
    facet_wrap(~panel, scales = "free_x") +
    scale_colour_manual(values = c(COL_SUP, COL_TWO), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
    labs(x = NULL, y = NULL) + theme_prelim() + theme(legend.position = "bottom")
  save_fig(p, "sup_external_c.pdf", 3.0)
}

# LaTeX table helpers ----
#' Format a number with fixed digits and a LaTeX minus sign.
#' @param x numeric. @param d digits. @param plus prefix positive nonzero values with "$+$".
fn <- function(x, d = 3, plus = FALSE) {
  s <- formatC(abs(x), format = "f", digits = d)
  neg <- x < 0 && as.numeric(s) != 0
  paste0(if (neg) "$-$" else if (plus && as.numeric(s) != 0) "$+$" else "", s)
}
#' "est (lo, hi)" string.
fci <- function(est, lo, hi, d = 3, plus = FALSE) sprintf("%s (%s, %s)", fn(est, d, plus), fn(lo, d), fn(hi, d))
#' Write a complete table environment; returns invisibly.
#' @param body tabular code. @param caption caption text. @param label label.
write_table <- function(file, body, caption, label) {
  txt <- paste0("\\begin{table}[!htb]\n\\centering\n\\small\n", body, "\n\\caption{", caption, "}\n\\label{", label, "}\n\\end{table}")
  writeLines(txt, file.path(TAB_DIR, file)); cat("Wrote", file.path(TAB_DIR, file), "\n")
}

# T1: simulation designs ----
#' Summarize the simulation registrations and the S-9 cross-validation choice.
#' Reads saved design, fit and choice rows and writes a LaTeX table.
#' @return Invisibly, the result of writing a table.
#' @examples
#' \dontrun{make_t1()}
make_t1 <- function() {
  des <- read.csv(file.path(S9_DIR, "design.csv"))
  ps <- sort(unique(des$p)); vs <- sort(unique(des$variance)); nseed <- length(unique(des$seed))
  stopifnot(identical(ps, c(50L, 2000L)), identical(vs, c(0.0025, 0.01, 0.05)), nseed == 10, nrow(des) == 60)
  m14 <- read.csv(file.path(S14_DIR, "metrics.csv"))
  # oracle C is a per-dataset value; report the range of the six cell means per strength
  cm <- aggregate(oracle_c ~ strength + share, m14, mean)
  oc <- aggregate(oracle_c ~ strength, cm, range)
  m9 <- read.csv(file.path(S9_DIR, "metrics.csv"))
  stopifnot(nrow(m9) == 480)
  cv9 <- read.csv(file.path(S9_DIR, "cv_chosen.csv"))
  stopifnot(nrow(cv9) == nrow(des), !anyDuplicated(cv9[c("seed", "p", "variance")]))
  cv_w1 <- sum(cv9$arm == "w1")
  check("S-9 datasets where CV selected w_S = 1", cv_w1, 55, 0)
  m13 <- read.csv(file.path(S13_DIR, "metrics.csv"))
  n13_sig <- sum(!m13$null); n13_null <- sum(m13$null)
  cat("\nT1 values: S-9 datasets", nrow(des), ", fits", nrow(m9), "; S-13 signal fits", n13_sig, "null fits", n13_null,
      "; S-14 fits", nrow(m14), "; oracle C ranges\n"); print(oc)
  check("S-9 datasets", nrow(des), 60, 0); check("S-9 fits", nrow(m9), 480, 0)
  check("S-13 signal fits (A1+A2, with CV folds)", n13_sig, 480, 0); check("S-13 null fits", n13_null, 360, 0)
  check("S-14 fits (excluding two-step: 3 arms x 6 cells x 10)", sum(m14$arm != "twostep"), 180, 0)
  # The S-14 findings say 240 fits: 4 arms (A1 w=1, A1 w=5, A2 w=1, two-step) x 6 cells x 10 seeds
  check("S-14 fits, all arms", nrow(m14), 240, 0)
  sc <- oc$oracle_c[oc$strength == "strong", ]; rc <- oc$oracle_c[oc$strength == "realistic", ]
  or_s <- sprintf("%.2f--%.2f", sc[1], sc[2]); or_r <- sprintf("%.2f--%.2f", rc[1], rc[2])
  body <- paste0(
"\\setlength{\\tabcolsep}{4pt}
\\begin{tabular}{>{\\raggedright\\arraybackslash}p{0.82in}>{\\raggedright\\arraybackslash}p{1.72in}>{\\raggedright\\arraybackslash}p{1.72in}>{\\raggedright\\arraybackslash}p{1.72in}}
\\toprule
 & S-9 & S-13 & S-14 \\\\
\\midrule
Purpose & Test whether survival supervision recovers a planted low-variance prognostic program, and how the weight $w_S$ affects it
 & Choose $w_S$ from the data by cross-validation and test calibration against the two-step fit
 & Test supervision when the planted program is near the detection threshold \\\\
\\addlinespace
Patients and genes & 400 training and 120 validation patients. $p = 50$ and $2{,}000$ genes.
 & 400 training and 120 validation patients. $p = 2{,}000$.
 & 400 training and 120 validation patients. $p = 2{,}000$. \\\\
\\addlinespace
Planted program & Share of centered variance 0.25\\%, 1\\%, 5\\%
 & Same three shares, plus a weak-survival cell (1\\%)
 & Share 0.05\\%, 0.1\\%, 0.2\\% \\\\
\\addlinespace
Survival & Strong, with about 30\\% censoring
 & Strong and weak (oracle $C \\approx 0.68$ for weak). Null fits permute training outcomes.
 & Strong (oracle $C$ ", or_s, ") and realistic (oracle $C$ ", or_r, ") \\\\
\\addlinespace
Weights $w_S$ & 0, 1, 10, 100, 1000, 10000, $p/50$ (40 at $p = 2{,}000$), and an arm targeting a share of 0.2
 & 1, 2, 5, 10, 20 or 50, chosen by five-fold CV partial log-likelihood in A1 and A2
 & 1 and 5 for A1 and 1 for A2 \\\\
\\addlinespace
Comparator & Unsupervised fit ($w_S = 0$) and earlier unsupervised fit on uncentered input
 & Two-step Cox using scores from the $w_S = 0$ loadings
 & Two-step \\\\
\\addlinespace
Datasets per cell & 10, for 60 datasets and 480 fits
 & 10, for ", n13_sig, " fits with CV and ", n13_null, " null fits
 & 10, for ", nrow(m14), " fits \\\\
\\addlinespace
CV weight & $w_S = 1$ in ", cv_w1, "/", nrow(cv9), " datasets
 & Selected per dataset (see Table~\\ref{tab:strength})
 & Fixed weights \\\\
\\bottomrule
\\end{tabular}")
  write_table("tab_sim_design.tex", body,
    "Designs of the three simulation studies of survival supervision. All use the same shared-loading generator and fit $K = 7$ programs. The planted program carries a small share of the variance and its coefficient sets the strength of the survival effect. A1 weights the Cox term and the KL divergence for $\\beta$ by $w_S$. A2 also uses a joint Gaussian posterior for $\\beta$.",
    "tab:sim-design")
}

# T2: training CV at D4 ----
make_t2 <- function() {
  s <- read.csv(file.path(REAL_DIR, "summary.csv")); s <- s[s$part == "size" & s$n_genes == 2218, ]
  f <- read.csv(file.path(REAL_DIR, "fits.csv")); f <- f[f$part == "size" & f$n_genes == 2218 & f$fold >= 1, ]
  a <- s[s$arm == "w1", ]; b <- s[s$arm == "w0cox", ]
  fw <- f[f$arm == "w1", ]; ft <- f[f$arm == "w0cox", ]
  stopifnot(identical(fw$fold, ft$fold))
  dd <- fw$C - ft$C; dm <- mean(dd); dse <- sd(dd) / sqrt(length(dd))
  cat("\nT2 values\n"); print(s[, c("arm", "mean_C", "se_C", "pooled_oof_C", "cv_partial_loglik", "folds_converged", "full_active")], digits = 4)
  cat(sprintf("paired fold difference %.4f (SE %.4f)\n", dm, dse))
  check("w1 mean C", a$mean_C, 0.663, 0.0006); check("w1 SE", a$se_C, 0.024, 0.0006); check("two-step mean C", b$mean_C, 0.652, 0.0006)
  check("two-step SE", b$se_C, 0.026, 0.0006); check("w1 pooled OOF C", a$pooled_oof_C, 0.663, 0.0006); check("two-step pooled OOF C", b$pooled_oof_C, 0.647, 0.0006)
  check("w1 CV partial loglik", a$cv_partial_loglik, -775, 0.5); check("two-step CV partial loglik", b$cv_partial_loglik, -777, 0.5)
  check("w1 folds converged", a$folds_converged, 3, 0); check("w1 survival-active", a$full_active, 3, 0); check("two-step survival-active", b$full_active, 4, 0)
  check("paired difference", dm, 0.011, 0.0006); check("paired difference SE", dse, 0.005, 0.0006)
  row <- function(lab, r) sprintf("%s & %s (%s) & %s & %s & %d/5 & %d \\\\", lab, fn(r$mean_C), fn(r$se_C), fn(r$pooled_oof_C),
                                   fn(r$cv_partial_loglik, 1), r$folds_converged, r$full_active)
  body <- paste0(
"\\begin{tabular}{lccccc}
\\toprule
Fit & \\shortstack{Mean fold $C$\\\\(SE)} & \\shortstack{Pooled\\\\out-of-fold $C$} & \\shortstack{CV partial\\\\log-likelihood} & \\shortstack{Folds\\\\converged} & \\shortstack{Survival-active\\\\programs} \\\\
\\midrule
", row("$w_S = 1$", a), "\n", row("Two-step", b), "\n",
"\\midrule
\\multicolumn{6}{l}{Paired fold difference ($w_S = 1$ minus two-step) \\quad ", fn(dm, 3, TRUE), " (SE ", fn(dse, 3), ")} \\\\
\\bottomrule
\\end{tabular}")
  write_table("tab_training_cv.tex", body,
    "Five-fold training cross-validation on the 2,218 selected genes (273 patients, 139 events). Convergence refers to the factorization in each arm. Activity is measured in the fit to all training patients using $|E[\\beta_k]| \\geq 1.96$ posterior SD for $w_S = 1$ and the Wald test for the two-step fit.",
    "tab:training-cv")
}

# T3: external cohorts ----
#' Tabulate external concordance, paired differences and calibration.
#' Reads fixed-fit cohort summaries and donor and record bootstrap intervals.
#' @return Invisibly, the result of writing a LaTeX table.
#' @examples
#' \dontrun{make_t3()}
make_t3 <- function() {
  e <- read.csv(file.path(EXTERNAL_DIR, "external_cindex.csv")); d <- read.csv(file.path(EXTERNAL_DIR, "external_paired_diff.csv"))
  donor <- read.csv(file.path(EXTERNAL_DIR, "external_paired_diff_donor.csv"))
  historical <- read.csv(file.path(OUT, "../../benchmark_snapshot_2026-10-14",
                                   "single_modality", "external_cindex_ci.csv"))
  historical <- historical[historical$model == "YFB (recommended)", ]
  stopifnot(setequal(historical$cohort, names(COHORT_LAB)), !anyDuplicated(historical$cohort), nrow(donor) == 1)
  rows <- character(0)
  calibration_rows <- character(0)
  for (k in names(COHORT_LAB)) {
    a <- e[e$cohort == k & e$model == "w1", ]; b <- e[e$cohort == k & e$model == "twostep", ]; p <- d[d$cohort == k, ]
    rows <- c(rows, sprintf("%s & %d (%d) & %s & %s & %s & %s \\\\", COHORT_LAB[k], a$n, a$events,
                            fci(a$C, a$lo, a$hi), fn(b$C), fci(p$diff, p$lo, p$hi, 3, TRUE),
                            fn(historical$c_index[historical$cohort == k], 4)))
    calibration_rows <- c(calibration_rows, sprintf("%s & %s (%s) & %s (%s) \\\\", COHORT_LAB[k],
      fn(a$cal_slope, 2), fn(a$cal_se, 2), fn(b$cal_slope, 2), fn(b$cal_se, 2)))
  }
  pp <- d[match(c("pooled (stratified, 5 cohorts)", "pooled without PACA_AU_seq"), d$cohort), ]
  stopifnot(!anyNA(pp$diff))
  #' Format a pooled comparison across the columns of the external table.
  #' @param lab Text to print as the row label.
  #' @param r One row of the saved paired difference table.
  #' @return One LaTeX table row.
  #' @examples
  #' \dontrun{prow("Pooled, donors", donor)}
  prow <- function(lab, r) sprintf("\\multicolumn{4}{l}{%s} & %s & \\\\", lab, fci(r$diff, r$lo, r$hi, 3, TRUE))
  cat("\nT3 values\n"); print(d, digits = 3)
  print(donor, digits = 3)
  check("pooled stratified diff", pp$diff[1], 0.007, 0.0006); check("pooled stratified lo", pp$lo[1], -0.002, 0.0006); check("pooled stratified hi", pp$hi[1], 0.017, 0.0006)
  check("pooled without seq diff", pp$diff[2], 0.006, 0.0006); check("pooled without seq lo", pp$lo[2], -0.003, 0.0006); check("pooled without seq hi", pp$hi[2], 0.016, 0.0006)
  check("pooled donor point estimate agrees with records", donor$diff, pp$diff[1], 1e-10)
  check("pooled donor lower limit", donor$lo, -0.003, 0.0006); check("pooled donor upper limit", donor$hi, 0.017, 0.0006)
  exd <- list(Dijk = c(0.018, 0, 0.036), Moffitt_GEO_array = c(0.007, -0.014, 0.026), PACA_AU_array = c(0, -0.026, 0.025),
              PACA_AU_seq = c(0.010, -0.019, 0.038), Puleo_array = c(0.001, -0.007, 0.010))
  for (k in names(exd)) { p <- d[d$cohort == k, ]
    check(paste(k, "paired diff"), p$diff, exd[[k]][1], 0.0006); check(paste(k, "paired lo"), p$lo, exd[[k]][2], 0.0006); check(paste(k, "paired hi"), p$hi, exd[[k]][3], 0.0006) }
  for (k in names(COHORT_LAB)) { a <- e[e$cohort == k & e$model == "w1", ]
    check(paste(k, "CI lo/hi"), a$lo + a$hi, c(Dijk = 0.513 + 0.658, Moffitt_GEO_array = 0.459 + 0.611, PACA_AU_array = 0.562 + 0.751,
                                                PACA_AU_seq = 0.579 + 0.780, Puleo_array = 0.599 + 0.686)[k], 0.0012) }
  body <- paste0(
"\\setlength{\\tabcolsep}{3pt}
\\begin{tabular}{lccccc}
\\multicolumn{6}{l}{\\textbf{A. External concordance}} \\\\
\\toprule
Cohort & $n$ (events) & $w_S = 1$ (95\\% CI) & Two-step & \\shortstack{Paired difference\\\\(95\\% CI)} & \\shortstack{Historical\\\\$C$} \\\\
\\midrule
", paste(rows, collapse = "\n"), "
\\midrule
", prow("Pooled, records", pp[1, ]), "\n", prow("Pooled, donors", donor), "\n",
   prow("Pooled without PACA-AU (RNA-seq)", pp[2, ]), "
\\bottomrule
\\end{tabular}

\\medskip
\\begin{tabular}{lcc}
\\multicolumn{3}{l}{\\textbf{B. Calibration slopes}} \\\\
\\toprule
Cohort & $w_S = 1$ (SE) & Two-step (SE) \\\\
\\midrule
", paste(calibration_rows, collapse = "\n"), "
\\bottomrule
\\end{tabular}")
  write_table("tab_external.tex", body,
    "External scoring of the fixed fits. The paired difference is $w_S = 1$ minus two-step. Pooled concordance uses pairs within series. The donor bootstrap resamples shared PACA-AU donors jointly. The last pooled row excludes PACA-AU RNA-seq. Calibration slopes come from a Cox regression of survival on the risk score. A slope of 1 indicates calibration. Historical $C$ is descriptive because the frozen benchmark used different preprocessing.",
    "tab:external")
}

# T6: supervision strength summary ----
make_t6 <- function() {
  sv <- read.csv(file.path(S13_DIR, "selected_vs_twostep.csv")); sv <- sv[sv$arm == "A1", ]
  m14 <- read.csv(file.path(S14_DIR, "metrics.csv"))
  cell_lab <- c(v0.0025 = "0.25\\% variance", v0.01 = "1\\% variance", v0.05 = "5\\% variance", weak = "Weak survival")
  wfmt <- function(s) { p <- strsplit(strsplit(s, " ")[[1]], ":"); paste(vapply(p, function(z) sprintf("%s (%s)", z[1], z[2]), ""), collapse = ", ") }
  ra <- vapply(sv$cell, function(cc) { r <- sv[sv$cell == cc, ]
    sprintf("%s & %s & %d / %d & %s / %s & %s / %s & %s / %s \\\\", cell_lab[cc], wfmt(r$weights), r$recovered, r$twostep_recovered,
            fn(r$mean_false_active, 1), fn(r$twostep_false_active, 1), fn(r$c_validation), fn(r$twostep_c), fn(r$slope_median, 2), fn(r$twostep_slope, 2)) }, "")
  cat("\nT6 panel A (S-13 A1 at CV-selected weight)\n"); print(sv[, c("cell", "weights", "recovered", "twostep_recovered", "mean_false_active", "twostep_false_active", "c_validation", "twostep_c", "slope_median", "twostep_slope")], digits = 3)
  # Checks against the findings table
  ec <- list(v0.0025 = c(10, 10, 0.9, 1.7, 0.946, 0.939, 1.24, 1.16), v0.01 = c(10, 10, 0, 0.9, 0.972, 0.972, 1.27, 1.10),
             v0.05 = c(10, 10, 0, 0.2, 0.986, 0.986, 1.45, 0.96), weak = c(10, 10, 1.2, 0.5, 0.671, 0.668, 1.05, 1.04))
  for (cc in names(ec)) { r <- sv[sv$cell == cc, ]; got <- unlist(r[, c("recovered", "twostep_recovered", "mean_false_active", "twostep_false_active", "c_validation", "twostep_c", "slope_median", "twostep_slope")])
    for (i in 1:8) check(sprintf("S-13 %s %s", cc, names(got)[i]), got[i], ec[[cc]][i], c(0, 0, 0.051, 0.051, 0.0006, 0.0006, 0.0051, 0.0051)[i]) }
  agg <- do.call(rbind, lapply(split(m14[m14$arm %in% c("twostep", "A1_w1"), ], list(m14$strength[m14$arm %in% c("twostep", "A1_w1")], m14$share[m14$arm %in% c("twostep", "A1_w1")], m14$arm[m14$arm %in% c("twostep", "A1_w1")]), drop = TRUE), function(d)
    data.frame(strength = d$strength[1], share = d$share[1], arm = d$arm[1], rec = sum(d$recovered), C = mean(d$c_validation), oracle = mean(d$oracle_c))))
  rb <- character(0)
  for (s in c("strong", "realistic")) for (sh in c(5e-4, 1e-3, 2e-3)) {
    a <- agg[agg$strength == s & agg$share == sh & agg$arm == "A1_w1", ]; b <- agg[agg$strength == s & agg$share == sh & agg$arm == "twostep", ]
    rb <- c(rb, sprintf("%s & %s\\%% & %d / %d & %s / %s & %s \\\\", if (s == "strong") "Strong" else "Realistic", formatC(sh * 100, format = "fg"),
                        a$rec, b$rec, fn(a$C), fn(b$C), fn(a$oracle)))
  }
  cat("T6 panel B (S-14)\n"); print(agg[order(agg$strength, agg$share, agg$arm), ], digits = 3, row.names = FALSE)
  eb <- list(`strong_5e-04` = c(5, 0, 0.800, 0.644), `strong_0.001` = c(10, 8, 0.901, 0.861), `strong_0.002` = c(10, 10, 0.936, 0.928),
             `realistic_5e-04` = c(0, 0, 0.630, 0.592), `realistic_0.001` = c(9, 8, 0.712, 0.705), `realistic_0.002` = c(10, 10, 0.728, 0.727))
  for (nm in names(eb)) { pr <- strsplit(nm, "_")[[1]]; a <- agg[agg$strength == pr[1] & abs(agg$share - as.numeric(pr[2])) < 1e-9 & agg$arm == "A1_w1", ]
    b <- agg[agg$strength == pr[1] & abs(agg$share - as.numeric(pr[2])) < 1e-9 & agg$arm == "twostep", ]
    check(paste("S-14", nm, "rec A1"), a$rec, eb[[nm]][1], 0); check(paste("S-14", nm, "rec two-step"), b$rec, eb[[nm]][2], 0)
    check(paste("S-14", nm, "C A1"), a$C, eb[[nm]][3], 0.0006); check(paste("S-14", nm, "C two-step"), b$C, eb[[nm]][4], 0.0006) }
  body <- paste0(
"\\setlength{\\tabcolsep}{4pt}
\\begin{tabular}{lccccc}
\\multicolumn{6}{l}{\\textbf{A. Weight chosen by cross-validation (S-13, A1)}} \\\\
\\toprule
Cell & \\shortstack{Chosen $w_S$\\\\(datasets)} & \\shortstack{Recovered\\\\(of 10)} & \\shortstack{False-active\\\\programs} & \\shortstack{Validation\\\\$C$} & \\shortstack{Calibration\\\\slope} \\\\
\\midrule
", paste(ra, collapse = "\n"), "
\\bottomrule
\\end{tabular}

\\medskip
\\begin{tabular}{lcccc}
\\multicolumn{5}{l}{\\textbf{B. Near the detection threshold (S-14, A1 at $w_S = 1$)}} \\\\
\\toprule
Survival & Planted share & \\shortstack{Recovered\\\\(of 10)} & Validation $C$ & Oracle $C$ \\\\
\\midrule
", paste(rb, collapse = "\n"), "
\\bottomrule
\\end{tabular}")
  write_table("tab_supervision_strength.tex", body,
    "Supervised fits against the matched two-step fit in simulation (10 datasets per cell, 2,000 genes). Entries separated by a slash give the supervised fit first and the two-step fit second. Panel A counts the datasets that selected each weight by five-fold cross-validated partial log-likelihood and reports the median calibration slope. Panel B reports oracle $C$, the validation $C$ of the true predictor.",
    "tab:strength")
}

# Run ----
cat("=== Simulation figures ===\n")
make_f1(); make_f12(); make_f13()
cat("\n=== Real-data figures and tables (REAL_DIR =", REAL_DIR, "; EXTERNAL_DIR =", EXTERNAL_DIR, ") ===\n")
make_f3(); make_ext()
cat("\n=== Tables ===\n")
make_t1(); make_t2(); make_t3(); make_t6()
cat("\n=== Check summary ===\n")
if (length(MISMATCHES)) {
  stop("Saved results differ from the checked findings:\n", paste0("  ", MISMATCHES, collapse = "\n"), call. = FALSE)
}
cat("All checked values match the findings files.\n")
