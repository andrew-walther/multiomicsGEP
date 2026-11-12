# ============================================================
# Script: run_methylation_iteration_extension.R
# Purpose: Check whether doubling the iteration cap resolves the 50k PL fit.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, survival, ebnm, Rcpp; frozen Task C runner and inputs
# ============================================================
# Approved numerical follow-up, 2026-11-12: K = 7, 50,000 CpGs, PL,
# 600 iterations, with the same initialization and every other option fixed.
# Run only after the registered study finishes, with one CPU/numerical thread.
# Memory: several GB. On Longleaf request one CPU and enough memory for copies
# of the 144 x 305,352 prepared universe and the selected feature matrix.
# Run from the repository root with /usr/bin/time -l to measure peak RSS:
# Rscript -e 'source("results/multimodal_real/run_methylation_iteration_extension.R")'
# Separate private source/config and checkpoints preserve all 300-iteration results.

# Reproduce the frozen baseline with a separate private iteration cap ----
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("data/methylation_single_modality_2026-11-12/source/runner.R")
base_local <- C_LOCAL
extension_local <- file.path(base_local, "iteration600")
base_path <- file.path(base_local, "methylation_50000_PL_f0.rds")
baseline <- readRDS(base_path)
config_base <- yaml::read_yaml(file.path(base_local, "source/config/globals.yml"))
iteration_cap <- 2L * config_base$cavi$max_iter
stopifnot(iteration_cap == 600L, baseline$row$status == "capped",
          ncol(baseline$fit$EF) == config_base$benchmark$k_merged_yfb_desurv)
dir.create(extension_local, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(extension_local, "design.rds"))) {
  stopifnot(file.copy(file.path(base_local, "design.rds"), extension_local),
            file.copy(file.path(base_local, "source"), extension_local, recursive = TRUE))
  config_extended <- config_base
  config_extended$cavi$max_iter <- iteration_cap
  yaml::write_yaml(config_extended, file.path(extension_local, "source/config/globals.yml"))
}
config_check <- yaml::read_yaml(file.path(extension_local, "source/config/globals.yml"))
stopifnot(config_check$cavi$max_iter == iteration_cap)
config_check$cavi$max_iter <- config_base$cavi$max_iter
if (!identical(config_check, config_base)) stop("Extension configuration differs from baseline.")
if (!identical(unname(tools::md5sum(file.path(base_local, "design.rds"))),
               unname(tools::md5sum(file.path(extension_local, "design.rds")))))
  stop("Extension design differs from baseline.")
code_paths <- list.files(file.path(base_local, "source/code"), recursive = TRUE)
if (!identical(unname(tools::md5sum(file.path(base_local, "source/code", code_paths))),
               unname(tools::md5sum(file.path(extension_local, "source/code", code_paths)))))
  stop("Extension model sources differ from baseline.")
C_LOCAL <- extension_local
# c_one() resets seed 42 and reads the copied model configuration. This reruns
# from the original initialization, not an incomplete posterior warm start.
c_one("methylation", 50000L, "PL", 0L)
extended <- readRDS(file.path(extension_local, "methylation_50000_PL_f0.rds"))
if (is.null(extended$fit)) {
  write.csv(extended$row, file.path(C_OUT, "iteration_extension.csv"), row.names = FALSE)
  stop("The 600-iteration diagnostic failed; its private checkpoint was saved.")
}

# Verify the same trajectory before assessing the additional iterations ----
prefix_n <- baseline$fit$history$n_iter
old_elbo <- baseline$fit$history$elbo_full
new_elbo <- extended$fit$history$elbo_full
if (length(new_elbo) < prefix_n) stop("Extended fit ended before the baseline trajectory prefix.")
prefix_error <- max(abs(new_elbo[seq_len(prefix_n)] - old_elbo) / pmax(1, abs(old_elbo)))
if (!is.finite(prefix_error) || prefix_error > 1e-10)
  stop("Extended trajectory does not reproduce the 300-iteration baseline.")
matched <- c_match(baseline$scores, extended$scores)
correlations <- rep(NA_real_, length(matched))
program_defined <- logical(length(matched))
for (k in seq_along(matched)) {
  old_score <- baseline$scores[, k]
  new_score <- extended$scores[, matched[k]]
  program_defined[k] <- all(is.finite(c(old_score, new_score))) &&
    stats::sd(old_score) > 0 && stats::sd(new_score) > 0
  if (program_defined[k]) correlations[k] <- stats::cor(old_score, new_score)
}
risk_finite <- all(is.finite(c(baseline$risk, extended$risk)))
risk_RMS_defined <- risk_finite && stats::sd(baseline$risk) > 0
risk_correlation_defined <- risk_RMS_defined && stats::sd(extended$risk) > 0
risk_correlation <- if (risk_correlation_defined) stats::cor(baseline$risk, extended$risk) else NA_real_
relative_risk_change <- if (risk_RMS_defined)
  sqrt(mean((extended$risk - baseline$risk)^2)) / stats::sd(baseline$risk) else NA_real_
row <- extended$row
row$K <- ncol(extended$fit$EF)
row$iteration_cap <- iteration_cap
row$baseline_iterations <- prefix_n
row$baseline_status <- baseline$row$status
row$baseline_fit_seconds <- baseline$row$seconds
row$prefix_relative_elbo_error <- prefix_error
row$final_relative_elbo_change <- tail(extended$fit$history$delta_elbo_rel, 1)
row$warning_count <- length(extended$warnings)
row$risk_correlation_defined <- risk_correlation_defined
row$risk_RMS_change_defined <- risk_RMS_defined
row$risk_correlation_with_300 <- risk_correlation
row$risk_RMS_change_over_baseline_SD <- relative_risk_change
# The external timing utility writes its resource summary after this script
# exits. Re-source after that exit to rebuild metrics from the saved checkpoint
# and completed log, without refitting. Pending measurements stay explicit.
log_path <- file.path(extension_local, "process.log")
lines <- if (file.exists(log_path)) readLines(log_path) else character()
memory <- grep("maximum resident set size", lines, value = TRUE)
wall <- grep("^[[:space:]]*[0-9.]+[[:space:]]+real", lines, value = TRUE)
row$peak_rss_MiB <- if (length(memory)) as.numeric(trimws(sub("maximum.*", "", memory[1]))) / 1024^2 else NA_real_
row$process_seconds <- if (length(wall)) as.numeric(strsplit(trimws(wall[1]), "[[:space:]]+")[[1]][1]) else NA_real_
row$memory_measured <- is.finite(row$peak_rss_MiB)
write.csv(row, file.path(C_OUT, "iteration_extension.csv"), row.names = FALSE)
write.csv(data.frame(baseline_program = seq_along(matched), extended_program = matched,
  comparison_defined = program_defined, score_correlation = correlations, baseline_active = baseline$activity$active,
  extended_active = extended$activity$active[matched]),
  file.path(C_OUT, "iteration_extension_programs.csv"), row.names = FALSE)
paths <- c(file.path(extension_local, "source/config/globals.yml"),
           "results/multimodal_real/run_methylation_iteration_extension.R")
write.csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))),
          file.path(C_OUT, "iteration_extension_hashes.csv"), row.names = FALSE)
print(row)
