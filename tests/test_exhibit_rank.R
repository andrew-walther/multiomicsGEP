# ============================================================
# Script: test_exhibit_rank.R
# Purpose: Check the EXHIBIT_K rank override: unset keeps the historical K = 7 rank and
#          paths; a set value changes the rank and moves every output to its own folder.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml; results/benchmark_sim/exhibit_rank.R
# Run: Rscript tests/test_exhibit_rank.R
# ============================================================

source("results/benchmark_sim/exhibit_rank.R")
cfg <- yaml::read_yaml("config/globals.yml")
fails <- 0L
check <- function(ok, msg) {
  cat(if (isTRUE(ok)) "PASS" else "FAIL", msg, "\n")
  if (!isTRUE(ok)) fails <<- fails + 1L
}
old <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04"

Sys.unsetenv("EXHIBIT_K")
check(is.na(exhibit_k()), "unset EXHIBIT_K gives NA")
check(exhibit_fit_k(cfg) == 7L && cfg$benchmark$k_merged_yfb_desurv == 7L,
      "unset EXHIBIT_K fits the historical rank 7")
check(identical(exhibit_out(old, "s11"), old), "unset EXHIBIT_K keeps the historical output folder")
check(identical(exhibit_data("data/x.rds", "y.rds"), "data/x.rds"), "unset EXHIBIT_K keeps the data path")
check(identical(exhibit_suffix(), "") && identical(exhibit_root(), "results/benchmark_sim/outputs"),
      "unset EXHIBIT_K adds no file suffix and keeps the output root")

Sys.setenv(EXHIBIT_K = "5")
check(exhibit_k() == 5L && exhibit_fit_k(cfg) == 5L, "EXHIBIT_K = 5 fits rank 5")
check(identical(exhibit_out(old, "s11"),
                "results/benchmark_sim/outputs/k5_working_2026-11-12/s11"),
      "EXHIBIT_K = 5 writes to the K = 5 folder")
check(identical(exhibit_data("data/x.rds", "y.rds"), "data/k5_working_2026-11-12/y.rds"),
      "EXHIBIT_K = 5 keeps patient-level files under data/")
check(exhibit_fit_k(cfg) == cfg$benchmark$k_working_yfb, "EXHIBIT_K = 5 equals the working rank")
check(identical(exhibit_suffix(), "_k5") &&
        identical(exhibit_root(), "results/benchmark_sim/outputs/k5_working_2026-11-12"),
      "EXHIBIT_K = 5 suffixes figure names and moves the output root")

for (bad in c("0", "-1", "5.5", "five", "05")) {
  Sys.setenv(EXHIBIT_K = bad)
  check(inherits(tryCatch(exhibit_k(), error = function(e) e), "error"),
        sprintf("EXHIBIT_K = '%s' is rejected", bad))
}
Sys.unsetenv("EXHIBIT_K")
cat(sprintf("%d failures\n", fails))
if (fails) quit(status = 1L)
