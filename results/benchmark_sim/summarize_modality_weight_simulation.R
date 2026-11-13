# ============================================================
# Script: summarize_modality_weight_simulation.R
# Purpose: Summarize saved modality-weight simulation aggregates.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R
# ============================================================

#' Validate and summarize the registered modality-weight simulation
#'
#' Reads the existing aggregate tables without refitting. Checks every registered
#' seed, arm, precision block and pattern row, including consistency of the
#' three-pattern success flag. Undefined shares and reset counts remain NA.
#' @param out Directory containing the saved metrics, patterns and precision CSVs.
#' @return Invisibly, a named list of four supplemental summary data frames.
#' @examples
#' # summarize_modality_weight_simulation()
summarize_modality_weight_simulation <- function(out = "results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weights") {
  metrics <- read.csv(file.path(out, "metrics.csv"), stringsAsFactors = FALSE, na.strings = "NA",
                      colClasses = c(error = "character"))
  patterns <- read.csv(file.path(out, "patterns.csv"), stringsAsFactors = FALSE, na.strings = "NA")
  precision <- read.csv(file.path(out, "precision.csv"), stringsAsFactors = FALSE, na.strings = "NA")
  keys <- c("block", "feature_ratio", "methylation_active_fraction", "seed", "arm")
  cell <- c("block", "feature_ratio", "methylation_active_fraction", "arm")
  #' Form unique keys from a fixed set of table columns
  #' @param x Data frame.
  #' @param cols Key column names.
  #' @return Character row keys.
  #' @examples
  #' # key(metrics, keys)
  key <- function(x, cols) do.call(paste, c(x[cols], sep = "|"))
  expected <- do.call(rbind, lapply(seq_len(5), function(i) {
    cells <- data.frame(block = c("A", "A", "A", "B", "B"),
      feature_ratio = c(5, 5, 1, 1, 5), methylation_active_fraction = c(.15, .03, .03, .15, .15))
    arms <- if (i <= 3) c("grouped_unweighted", "grouped_count", "grouped_eigenvalue", "expression_only") else
      c("grouped_count", "grouped_eigenvalue", "pooled_count", "grouped_unweighted")
    cbind(cells[rep(i, 40), ], seed = rep(1:10, each = 4), arm = rep(arms, 10))
  }))
  #' Require exact registered row identities without duplication
  #' @param actual,wanted Observed and expected data frames.
  #' @param cols Key columns.
  #' @param label Table name in errors.
  #' @return Invisible NULL or an error.
  #' @examples
  #' # check_keys(metrics, expected, keys, "Metrics")
  check_keys <- function(actual, wanted, cols, label) {
    a <- key(actual, cols); w <- key(wanted, cols)
    if (anyDuplicated(a) || anyDuplicated(w) || length(a) != length(w) || !setequal(a, w))
      stop(label, " rows do not match the registered design.")
  }
  check_keys(metrics, expected, keys, "Metrics")
  if (nrow(metrics) != 200L || any(nzchar(metrics$error)) || anyNA(metrics$error))
    stop("Expected 200 successful metrics rows.")
  b <- metrics[metrics$block == "B", ]
  wanted_patterns <- merge(b[keys], data.frame(truth_program = 1:3))
  check_keys(patterns, wanted_patterns, c(keys, "truth_program"), "Pattern")
  for (i in seq_len(nrow(b))) {
    p <- patterns[key(patterns, keys) == key(b[i, ], keys), ]
    if (anyNA(p$loading_recovered) || anyNA(p$pattern_correct) ||
        is.na(b$full_pattern_success[i]) ||
        !identical(all(p$loading_recovered & p$pattern_correct), b$full_pattern_success[i]))
      stop("Full-pattern success disagrees with the three pattern rows.")
  }
  wanted_precision <- do.call(rbind, lapply(seq_len(nrow(metrics)), function(i) {
    mods <- if (metrics$arm[i] == "expression_only") "expression" else c("expression", "methylation")
    cbind(metrics[rep(i, 5L * length(mods)), keys], expand.grid(program = 1:5, modality = mods))
  }))
  check_keys(precision, wanted_precision, c(keys, "program", "modality"), "Precision")
  if (anyNA(precision$empty) || anyNA(patterns[c("loading_recovered", "pattern_correct", "active", "correct_sign", "survival_recovered")]) ||
      anyNA(metrics[c("recovered", "recovered2")])) stop("Required outcome fields contain missing values.")
  #' Split rows using the registered grouping columns
  #' @param x Data frame.
  #' @param cols Grouping column names.
  #' @return Named list of data frames.
  #' @examples
  #' # groups(metrics, cell)
  groups <- function(x, cols) split(x, interaction(x[cols], drop = TRUE, lex.order = TRUE))
  #' Retain identifying columns from the first row of a group
  #' @param x Data frame containing one group.
  #' @param cols Identifying column names.
  #' @return One-row data frame.
  #' @examples
  #' # identity_row(metrics, cell)
  identity_row <- function(x, cols) x[1, cols, drop = FALSE]
  #' Sum recorded reset events while keeping entirely unobserved totals missing
  #' @param x Numeric reset counts with optional NA values.
  #' @return Numeric total or NA with its separately reported denominator.
  #' @examples
  #' # sum_defined(c(1L, NA_integer_))
  sum_defined <- function(x) if (all(is.na(x))) NA_integer_ else sum(x, na.rm = TRUE)
  program_summary <- do.call(rbind, lapply(groups(patterns, c(cell, "truth_program")), function(x) {
    data.frame(identity_row(x, c(cell, "truth_program")), seeds = nrow(x),
      loading_recovered = sum(x$loading_recovered), pattern_correct = sum(x$pattern_correct),
      loading_and_pattern = sum(x$loading_recovered & x$pattern_correct),
      active = sum(x$active), correct_sign = sum(x$correct_sign), survival_recovered = sum(x$survival_recovered))
  }))
  precision_summary <- do.call(rbind, lapply(groups(precision, c(cell, "modality", "program")), function(x) {
    s <- x$share[is.finite(x$share)]
    data.frame(identity_row(x, c(cell, "modality", "program")),
      mean_share = if (length(s)) mean(s) else NA_real_, min_share = if (length(s)) min(s) else NA_real_,
      max_share = if (length(s)) max(s) else NA_real_, share_defined = length(s),
      empty_count = sum(x$empty), empty_denominator = nrow(x),
      point_mass_resets = sum_defined(x$point_mass_resets), reset_defined = sum(!is.na(x$point_mass_resets)))
  }))
  empty_summary <- do.call(rbind, lapply(groups(precision, c(cell, "modality")), function(x) {
    data.frame(identity_row(x, c(cell, "modality")), final_empty_count = sum(x$empty),
      block_denominator = nrow(x), point_mass_resets = sum_defined(x$point_mass_resets),
      reset_defined = sum(!is.na(x$point_mass_resets)))
  }))
  secondary <- do.call(rbind, lapply(groups(metrics[metrics$block == "A", ], cell), function(x) {
    data.frame(identity_row(x, cell), seeds = nrow(x), primary_recovered = sum(x$recovered),
      secondary_recovered = sum(x$recovered2))
  }))
  result <- list(summary_B_programs = program_summary, summary_precision = precision_summary,
                 summary_empty_blocks = empty_summary, summary_A_secondary = secondary)
  for (name in names(result)) write.csv(result[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "NA")
  invisible(result)
}

if (sys.nframe() == 0L) summarize_modality_weight_simulation()
