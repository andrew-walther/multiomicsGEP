# ============================================================
# Script: summarize_yfb_supervision_strength.R
# Purpose: Attempt S-13 tables: weight selected by CV partial log-likelihood, the registered
#          criteria per arm, per-weight curves, program roles, null false activity, two-step.
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: base R, grDevices
# Run: Rscript results/benchmark_sim/summarize_yfb_supervision_strength.R
# ============================================================
# Reads checkpoints only; refits nothing. Registration: docs/plans/2026-11-05-supervision-
# strength-plan.md. Gate S (S-9 definition): a cell qualifies with recovery >= 3/10 above the
# S-7 unsupervised reference, mean false-active <= 1 and >= 9/10 converged; the arm passes with
# >= 2 of the 3 S-9 cells.

OUT <- "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05"
S7  <- "results/benchmark_sim/outputs/yfb_supervision_2026-10-30/attempt7"

# Load ----
files <- list.files(file.path(OUT, "fits"), pattern = "\\.rds$", full.names = TRUE)
ck <- lapply(files, readRDS)
m <- do.call(rbind, lapply(ck, function(x) { r <- x$metrics
  for (nm in c("n_survival_only", "n_prognostic")) if (is.null(r[[nm]])) r[[nm]] <- NA
  r }))
m$cell <- as.character(m$cell)
cat(sprintf("%d checkpoints (%d signal, %d null); %d with errors\n", nrow(m), sum(!m$null),
            sum(m$null), sum(!is.na(m$error) & m$error != "")))
sig <- m[!m$null, ]; nul <- m[m$null, ]
expected <- 10 * 4 * 2 * 6
if (nrow(sig) < expected) warning(expected - nrow(sig), " registered signal tasks are missing")
write.csv(m, file.path(OUT, "metrics.csv"), row.names = FALSE)

# Weight selection by CV partial log-likelihood ----
#' Per dataset and arm: the weight with the largest summed CV partial log-likelihood among
#' weights whose five folds all fit; ties to the smaller weight.
pick <- do.call(rbind, lapply(split(sig, list(sig$cell, sig$seed, sig$arm), drop = TRUE), function(x) {
  x <- x[x$cv_folds_ok == 5 & is.finite(x$cv_pl), ]
  if (!nrow(x)) return(NULL)
  x <- x[order(-x$cv_pl, x$w_S), ]; x[1, ]
}))
write.csv(pick, file.path(OUT, "cv_selected.csv"), row.names = FALSE)

# Criteria per arm ----
s7 <- read.csv(file.path(S7, "metrics.csv"))
ref <- subset(s7, arm == "unsupervised" & p == 2000)[, c("seed", "variance", "recovered")]
cells <- do.call(rbind, lapply(split(pick, list(pick$arm, pick$cell), drop = TRUE), function(x) {
  v <- if (x$cell[1] == "weak") NA else as.numeric(sub("v", "", x$cell[1]))
  rr <- if (is.na(v)) NA else sum(ref$recovered[abs(ref$variance - v) < 1e-12])
  data.frame(arm = x$arm[1], cell = x$cell[1], n = nrow(x),
             weights = paste(names(table(x$w_S)), table(x$w_S), sep = ":", collapse = " "),
             recovered = sum(x$recovered), reference_recovered = rr,
             mean_false_active = mean(x$false_active), converged = sum(x$converged),
             slope_median = median(x$cal_slope, na.rm = TRUE),
             c_validation = mean(x$c_validation, na.rm = TRUE), cv_c = mean(x$cv_c, na.rm = TRUE),
             match_survival_only = sum(x$match_role == "survival_only_direction", na.rm = TRUE),
             mean_survival_only = mean(x$n_survival_only, na.rm = TRUE))
}))
cells$gate_qualifies <- with(cells, cell != "weak" & recovered - reference_recovered >= 3 &
                               mean_false_active <= 1 & converged >= .9 * n)
cells$slope_ok <- cells$slope_median >= .8 & cells$slope_median <= 1.25
write.csv(cells, file.path(OUT, "selected_cells.csv"), row.names = FALSE)
print(cells, digits = 3)

# Null false activity at the modal selected weight ----
modal <- tapply(pick$w_S[pick$cell != "weak"], pick$arm[pick$cell != "weak"],
                function(w) as.numeric(names(which.max(table(w)))))
nullc <- if (nrow(nul)) do.call(rbind, lapply(split(nul, list(nul$arm, nul$cell, nul$w_S), drop = TRUE),
  function(x) data.frame(arm = x$arm[1], cell = x$cell[1], w_S = x$w_S[1], n = nrow(x),
                         mean_false_active = mean(x$false_active), converged = sum(x$converged)))) else NULL
if (!is.null(nullc)) {
  nullc$at_modal_weight <- nullc$w_S == modal[nullc$arm]
  write.csv(nullc, file.path(OUT, "null_cells.csv"), row.names = FALSE)
}

verdict <- do.call(rbind, lapply(split(cells, cells$arm), function(x) {
  s9 <- x[x$cell != "weak", ]
  nm <- if (!is.null(nullc)) subset(nullc, arm == x$arm[1] & at_modal_weight) else NULL
  data.frame(arm = x$arm[1], modal_weight = modal[[x$arm[1]]],
             gate_cells = sum(s9$gate_qualifies), gate_pass = sum(s9$gate_qualifies) >= 2,
             slope_cells = sum(s9$slope_ok), slope_pass = sum(s9$slope_ok) >= 2,
             convergence_pass = all(x$converged >= .9 * x$n),
             null_max_false_active = if (!is.null(nm) && nrow(nm)) max(nm$mean_false_active) else NA,
             null_pass = if (!is.null(nm) && nrow(nm)) all(nm$mean_false_active <= 1) else NA)
}))
write.csv(verdict, file.path(OUT, "verdict.csv"), row.names = FALSE)
print(verdict)

# Per-weight curves ----
curve <- do.call(rbind, lapply(split(sig, list(sig$arm, sig$cell, sig$w_S), drop = TRUE), function(x)
  data.frame(arm = x$arm[1], cell = x$cell[1], w_S = x$w_S[1], recovered = sum(x$recovered),
             mean_false_active = mean(x$false_active), converged = sum(x$converged),
             slope_median = median(x$cal_slope, na.rm = TRUE),
             c_validation = mean(x$c_validation, na.rm = TRUE), cv_pl = mean(x$cv_pl, na.rm = TRUE),
             cv_c = mean(x$cv_c, na.rm = TRUE), share_match = median(x$share_match, na.rm = TRUE),
             match_survival_only = sum(x$match_role == "survival_only_direction", na.rm = TRUE),
             mean_survival_only = mean(x$n_survival_only, na.rm = TRUE),
             sweeps = median(x$sweeps, na.rm = TRUE))))
curve <- curve[order(curve$arm, curve$cell, curve$w_S), ]
write.csv(curve, file.path(OUT, "weight_curves.csv"), row.names = FALSE)

png(file.path(OUT, "weight_curves.png"), width = 1600, height = 1000, res = 120)
op <- par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
mets <- c(recovered = "Recovered (of 10)", mean_false_active = "Mean false-active",
          slope_median = "Calibration slope (median)", converged = "Converged (of 10)",
          c_validation = "Validation C", cv_c = "Training CV C",
          share_match = "Survival share, matched", mean_survival_only = "Survival-only directions")
cols <- c(weak = "#d95f02", v0.0025 = "#1b9e77", v0.01 = "#7570b3", v0.05 = "#e7298a")
for (mt in names(mets)) {
  plot(NA, xlim = c(1, 50), ylim = range(curve[[mt]], na.rm = TRUE), log = "x", xlab = "w_S",
       ylab = mets[[mt]], main = mets[[mt]])
  if (mt == "slope_median") abline(h = c(.8, 1, 1.25), lty = c(3, 2, 3), col = "grey50")
  for (a in c("A1", "A2")) for (cl in names(cols)) {
    z <- curve[curve$arm == a & curve$cell == cl, ]
    lines(z$w_S, z[[mt]], col = cols[[cl]], lty = if (a == "A1") 1 else 2, type = "b", pch = 19, cex = .6)
  }
}
legend("topright", c(names(cols), "A1 solid, A2 dashed"), col = c(cols, "black"),
       lty = 1, bty = "n", cex = .8)
par(op); dev.off()

# Two-step comparison ----
ts <- file.path(OUT, "twostep.csv")
if (file.exists(ts)) {
  t2 <- read.csv(ts)
  tcell <- do.call(rbind, lapply(split(t2, t2$cell), function(x)
    data.frame(cell = x$cell[1], twostep_recovered = sum(x$recovered),
               twostep_false_active = mean(x$false_active),
               twostep_c = mean(x$c_validation), twostep_slope = median(x$cal_slope))))
  write.csv(merge(cells, tcell), file.path(OUT, "selected_vs_twostep.csv"), row.names = FALSE)
  print(merge(cells, tcell)[, c("arm", "cell", "recovered", "twostep_recovered", "c_validation",
                                "twostep_c", "slope_median", "twostep_slope")], digits = 3)
}
