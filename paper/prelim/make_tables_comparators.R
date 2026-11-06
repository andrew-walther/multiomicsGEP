# ============================================================
# Script: make_tables_comparators.R
# Purpose: Build the Chapter 4 comparator and program-rule tables from saved aggregate CSVs.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: base R only
# ============================================================
# Run from the repository root: Rscript paper/prelim/make_tables_comparators.R
# Inputs are the committed aggregate outputs of run_prompt6_comparators.R,
# run_prompt7_real.R and run_prompt7_ebmf.R. No patient-level data are read.

OUT <- "results/benchmark_sim/outputs"
TAB <- "paper/prelim/tables"

#' Read a saved aggregate CSV, failing loudly when it is missing.
#' @param ... Path components below the outputs directory.
#' @return data.frame
#' @examples
#' # rd("prompt6_comparators", "paired.csv")
rd <- function(...) {
  p <- file.path(OUT, ...)
  if (!file.exists(p)) stop("Missing saved evidence: ", p)
  read.csv(p, stringsAsFactors = FALSE)
}

#' Format a number with a LaTeX minus sign.
#' @param x Numeric vector. @param d Decimal places.
#' @return Character vector.
#' @examples
#' # f(-0.0123, 3)
f <- function(x, d = 3) {
  s <- formatC(x, format = "f", digits = d)
  ifelse(x < 0, paste0("$-$", sub("^-", "", s)), s)
}

#' Format a signed difference with an explicit plus sign.
#' @param x Numeric vector. @param d Decimal places.
#' @return Character vector.
#' @examples
#' # fs(0.007)
fs <- function(x, d = 3) ifelse(x > 0, paste0("$+$", f(x, d)), f(x, d))

# Comparator table ----
e <- rd("prompt6_comparators", "external.csv")
pr <- rd("prompt6_comparators", "paired.csv")
hr <- rd("prompt6_comparators", "pooled_hr.csv")
rk <- rd("prompt7_desurv_style", "rank_paired.csv")
rke <- rd("prompt7_desurv_style", "rank_external.csv")
rkh <- rd("prompt7_desurv_style", "rank_pooled_hr.csv")
lab <- c(w1 = "Supervised ($w_S = 1$)", twostep = "Two-step (rank-7 EBMF + Cox)",
         DeSurv = "DeSurv", elastic_net = "Elastic-net Cox (D4 genes)", PurIST = "PurIST",
         EBMF_15 = "EBMF rank 15 + Cox", EBMF_40 = "EBMF rank 40 + Cox")
allE <- rbind(e[, c("model", "C", "cal_slope")], rke[grepl("EBMF_(15|40)", rke$model), c("model", "C", "cal_slope")])
allH <- rbind(hr[, c("model", "hr", "donor_lo", "donor_hi")],
              rkh[grepl("EBMF_(15|40)", rkh$model), c("model", "hr", "donor_lo", "donor_hi")])
allP <- rbind(pr[, c("model", "diff", "lo", "hi")], rk[grepl("EBMF_(15|40)", rk$model), c("model", "diff", "lo", "hi")])
rows <- vapply(names(lab), function(m) {
  c_mean <- mean(allE$C[allE$model == m])
  slopes <- range(allE$cal_slope[allE$model == m])
  h <- allH[allH$model == m, ]
  d <- allP[allP$model == m, ]
  stopifnot(nrow(h) == 1L, sum(allE$model == m) == 5L)
  paste(lab[[m]], f(c_mean), sprintf("%s to %s", f(slopes[1], 2), f(slopes[2], 2)),
        sprintf("%s (%s, %s)", f(h$hr, 2), f(h$donor_lo, 2), f(h$donor_hi, 2)),
        if (nrow(d)) sprintf("%s (%s, %s)", fs(d$diff), f(d$lo), f(d$hi)) else "",
        sep = " & ")
}, character(1))
tex <- c("\\begin{table}[!htb]", "\\centering", "\\small", "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  "Model & Mean $C$ & \\shortstack{Calibration\\\\slope range} & \\shortstack{HR per SD\\\\(95\\% CI)} & \\shortstack{Supervised minus\\\\model (95\\% CI)} \\\\",
  "\\midrule", paste(rows, "\\\\"), "\\bottomrule", "\\end{tabular}",
  sprintf("\\caption{Fixed supervised fit and comparators trained on the same %d patients, scored in the five validation series (%d records from %d donors). Mean $C$ is the equal-weight mean over series. The hazard ratio is per pooled SD of each score from a Cox model stratified by series, with a robust interval clustered on donor. The difference in $C$ is averaged over series with a donor-bootstrap interval (%d replicates).}",
    273L, pr$n[1], pr$donors[1], pr$B[1]),
  "\\label{tab:comparators}", "\\end{table}")
writeLines(tex, file.path(TAB, "tab_comparators.tex"))

# Program rule table ----
r <- rd("prompt7_desurv_style", "program_rule.csv")
r <- r[r$active | r$ratio > .1, ]
rows <- sprintf("%d & %s & %s & %s & %d & %s \\\\", r$program, ifelse(r$active, "Yes", "No"),
  f(r$ratio, 2), f(r$variance_share), r$active_folds, ifelse(r$main, "Yes", "No"))
tex <- c("\\begin{table}[!htb]", "\\centering", "\\small", "\\begin{tabular}{lccccc}", "\\toprule",
  "Program & Active & \\shortstack{Coefficient\\\\ratio} & \\shortstack{Signed share of\\\\risk variance} & \\shortstack{Active\\\\folds (of 5)} & \\shortstack{Main\\\\program} \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
  "\\caption{Program classification on the supervised training fit. The coefficient ratio is $|E\\beta_k|\\,\\mathrm{SD}(YE[F_k])$ relative to the largest program. The signed share is $\\mathrm{Cov}(\\eta_k,\\eta)/\\mathrm{Var}(\\eta)$ and sums to one over programs. Active folds count training folds in which a loading matched with absolute correlation at least 0.7 was survival-active. A main program needs a ratio of at least 0.25, a share of at least 5\\% and at least four active folds. Programs with ratio below 0.1 and no activity are omitted.}",
  "\\label{tab:program-rule}", "\\end{table}")
writeLines(tex, file.path(TAB, "tab_program_rule.tex"))
cat("Wrote tab_comparators.tex and tab_program_rule.tex\n")
