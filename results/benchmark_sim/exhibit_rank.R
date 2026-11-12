# ============================================================
# Script: exhibit_rank.R
# Purpose: Opt-in rank override for rerunning the reported-fit exhibits at another K.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R
# ============================================================
# The reported-fit exhibits (S-11 training fits, S-12 external scores, S-15 weight refit,
# the comparator and risk-group analyses and the Chapter 4 program builders) were written
# for K = 7 with fixed paths. Setting the environment variable EXHIBIT_K to a positive
# integer reruns them at that rank, reading and writing under separate folders:
#   results/benchmark_sim/outputs/k<K>_working_2026-11-12/<part>/   (aggregates)
#   data/k<K>_working_2026-11-12/                                    (patient-level)
# With EXHIBIT_K unset every script keeps its historical rank and paths, so the K = 7
# outputs reproduce unchanged. An environment variable (not a command-line argument) is
# used because these scripts source one another and start PSOCK workers, which inherit it.

#' Rank requested through EXHIBIT_K
#' @return Integer K, or NA_integer_ when EXHIBIT_K is unset (historical behaviour).
#' @examples
#' # Sys.setenv(EXHIBIT_K = 5); exhibit_k()  # 5
#' @family exhibit_rank
exhibit_k <- function() {
  v <- Sys.getenv("EXHIBIT_K", "")
  if (!nzchar(v)) return(NA_integer_)
  k <- suppressWarnings(as.integer(v))
  if (is.na(k) || k < 1L || as.character(k) != v) stop("EXHIBIT_K must be a positive integer, got: ", v)
  k
}

#' Rank to fit: EXHIBIT_K when set, otherwise the historical configuration value
#' @param cfg Parsed config/globals.yml.
#' @return Integer K.
#' @examples
#' # exhibit_fit_k(yaml::read_yaml("config/globals.yml"))  # 7 when EXHIBIT_K is unset
#' @family exhibit_rank
exhibit_fit_k <- function(cfg) {
  k <- exhibit_k()
  if (is.na(k)) as.integer(cfg$benchmark$k_merged_yfb_desurv) else k
}

#' Aggregate output folder for one exhibit part
#' @param default Historical folder (used when EXHIBIT_K is unset).
#' @param part Subfolder name under the rank-specific root, e.g. "s11".
#' @return Folder path.
#' @examples
#' # exhibit_out("results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04", "s11")
#' @family exhibit_rank
exhibit_out <- function(default, part) {
  k <- exhibit_k()
  if (is.na(k)) default else
    file.path(sprintf("results/benchmark_sim/outputs/k%d_working_2026-11-12", k), part)
}

#' Patient-level file path (git-ignored data/ folder)
#' @param default Historical path (used when EXHIBIT_K is unset).
#' @param name File name under data/k<K>_working_2026-11-12/.
#' @return File path.
#' @examples
#' # exhibit_data("data/yfb_external_patient_scores_ws1.rds", "external_patient_scores_ws1.rds")
#' @family exhibit_rank
exhibit_data <- function(default, name) {
  k <- exhibit_k()
  if (is.na(k)) default else file.path(sprintf("data/k%d_working_2026-11-12", k), name)
}

#' Root aggregate folder for the rank (holds the per-part subfolders)
#' @param default Historical root (used when EXHIBIT_K is unset).
#' @return Folder path.
#' @examples
#' # exhibit_root("results/benchmark_sim/outputs")
#' @family exhibit_rank
exhibit_root <- function(default = "results/benchmark_sim/outputs") {
  k <- exhibit_k()
  if (is.na(k)) default else sprintf("results/benchmark_sim/outputs/k%d_working_2026-11-12", k)
}

#' File-name suffix for figures and tables written next to the historical ones
#' @return "" when EXHIBIT_K is unset, otherwise "_k<K>" (e.g. "_k5").
#' @examples
#' # paste0("sup_programs_heatmap", exhibit_suffix(), ".pdf")
#' @family exhibit_rank
exhibit_suffix <- function() {
  k <- exhibit_k()
  if (is.na(k)) "" else sprintf("_k%d", k)
}
