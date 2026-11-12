# ============================================================
# Script: test_methylation_k5_rerun.R
# Purpose: Verify rank isolation, frozen dependencies and rerun coverage.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, parallel; methylation study runner
# ============================================================
local({
source("results/multimodal_real/run_methylation_single_modality.R", local = TRUE)
root <- tempfile("c-k5-"); dir.create(root)
on.exit(unlink(root, recursive = TRUE))
C_BASE_LOCAL <- file.path(root, "private")
C_BASE_OUT <- file.path(root, "aggregate")
dir.create(C_BASE_LOCAL); dir.create(C_BASE_OUT)
saveRDS(list(fold = rep(1:5, 2)), file.path(C_BASE_LOCAL, "design.rds"))
# Rank selection leaves the historical directory intact and links the saved design.
c_set_rank(7)
stopifnot(identical(C_LOCAL, C_BASE_LOCAL), C_K == 7L)
c_set_rank(5)
stopifnot(basename(C_LOCAL) == "K05", basename(C_OUT) == "K05", C_K == 5L,
  identical(readRDS(file.path(C_LOCAL, "design.rds"))$fold, rep(1:5, 2)))
stopifnot(inherits(try(c_set_rank(5.5), silent = TRUE), "try-error"))

# Freeze exactly the historical dependencies, while capturing the new runner.
baseline <- file.path(C_BASE_LOCAL, "source")
dir.create(file.path(baseline, "code"), recursive = TRUE)
dir.create(file.path(baseline, "config"))
writeLines("# frozen fitter", file.path(baseline, "code/fit_ssbmf_yfb.R"))
writeLines("// frozen kernel", file.path(baseline, "code/multimodal_yfb_sweep.cpp"))
yaml::write_yaml(list(benchmark = list(k_merged_yfb_desurv = 7L)),
                 file.path(baseline, "config/globals.yml"))
writeLines("# historical runner", file.path(baseline, "runner.R"))
files <- list.files(baseline, recursive = TRUE)
write.csv(data.frame(path = files, md5 = unname(tools::md5sum(file.path(baseline, files)))),
          file.path(C_BASE_OUT, "source_hashes.csv"), row.names = FALSE)
write.csv(data.frame(revision = "historical-revision"),
          file.path(C_BASE_OUT, "source_revision.csv"), row.names = FALSE)
c_freeze()
stopifnot(read.csv(file.path(C_OUT, "source_revision.csv"))$revision == "historical-revision",
  yaml::read_yaml(file.path(C_LOCAL, "source/config/globals.yml"))$benchmark$k_merged_yfb_desurv == 7L,
  identical(readLines(file.path(C_LOCAL, "source/code/fit_ssbmf_yfb.R")), "# frozen fitter"))
c_freeze() # Manifest validation on resume.
copy_path <- file.path(C_LOCAL, "source/code/fit_ssbmf_yfb.R")
writeLines("# changed", copy_path)
stopifnot(inherits(try(c_freeze(), silent = TRUE), "try-error"))
writeLines("# frozen fitter", copy_path)

# Child rank propagation and resumed rank labels, without starting R processes.
#' Simulate a child that checkpoints a successful rank-five fit
#' @param command,args Subprocess command and argument vector.
#' @param ... Remaining subprocess arguments.
#' @return Integer zero status.
#' @examples
#' # system2('unused', '--K=5')
system2 <- function(command, args, ...) {
  stopifnot("--K=5" %in% args)
  saveRDS(list(row = data.frame(status = "converged"), fit = list(EF = matrix(1, 2, 5)),
    warnings = character()), file.path(C_LOCAL, "expression_2218_PL_f0.rds"))
  0L
}
launched <- c_launch("expression", 2218L, "PL", 0L)
stopifnot(launched$K == 5L,
  !file.exists(file.path(C_BASE_LOCAL, "expression_2218_PL_f0.rds")))
rm(system2)

# Complete the expression reference, then screen the full sweep by actual timings.
#' Skip repeated snapshot work in scheduling tests
#' @return NULL.
#' @examples
#' c_freeze()
c_freeze <- function() NULL
calls <- list(); coverage_case <- "success"
#' Capture every requested job and return distinct synthetic fold metrics
#' @param modality,count,arm,fold Registered job identifiers.
#' @return One aggregate metric row.
#' @examples
#' # c_launch('expression', 2218L, 'PL', 1L)
c_launch <- function(modality, count, arm, fold) {
  if (fold == 1L && modality == "expression" && arm == "PL") {
    if (coverage_case == "missing") return(NULL)
    if (coverage_case == "duplicate") fold <- 2L
  }
  calls[[length(calls) + 1L]] <<- data.frame(modality, count, arm, fold)
  data.frame(modality, count, arm, fold, K = 5L, status = "converged", converged = TRUE,
    seconds = if (count >= 50000L) 2000 else 60,
    process_seconds = if (count >= 50000L) 2000 else 60,
    C = .5 + fold / 100 + if (arm == "PE") .05 else 0,
    nonzero = 5, survival_active = if (arm == "PL") 2 else 0)
}
#' Map jobs sequentially to keep mock state in this process
#' @param X,FUN Input tasks and function.
#' @param ... Ignored parallel settings.
#' @return List of results.
#' @examples
#' sequential_map(1:3, identity)
sequential_map <- function(X, FUN, ...) lapply(X, FUN)
ns <- asNamespace("parallel"); original <- get("mclapply", ns)
unlockBinding("mclapply", ns); assign("mclapply", sequential_map, ns); lockBinding("mclapply", ns)
tryCatch({
  c_rerun()
  jobs <- do.call(rbind, calls)
  stopifnot(nrow(jobs) == 12L, !anyDuplicated(jobs), all(jobs$modality == "expression"))
  ex <- read.csv(file.path(C_OUT, "expression_summary.csv"))
  stopifnot(all(abs(ex$mean_C - c(.53, .58)) < 1e-12), all(ex$converged_folds == 5),
    ex$full_survival_active[ex$arm == "PL"] == 2)
  # All full fits run, but counts whose two-setup CV projection exceeds 3h do not.
  saveRDS(list(counts = c(10000L, 25000L, 50000L, 100000L, 182934L)),
          file.path(C_LOCAL, "design.rds"))
  calls <- list()
  c_run()
  jobs <- do.call(rbind, calls)
  stopifnot(sum(jobs$fold == 0L & jobs$modality == "methylation") == 10L,
    sum(jobs$fold > 0L & jobs$modality == "expression") == 10L,
    sum(jobs$fold > 0L & jobs$modality == "methylation") == 20L,
    !any(jobs$fold > 0L & jobs$count >= 50000L), !anyDuplicated(jobs))
  projected <- read.csv(file.path(C_OUT, "cv_eligibility.csv"))
  stopifnot(all(projected$projected_both_setups_CV_seconds[1:2] == 600),
    all(projected$projected_both_setups_CV_seconds[3:5] == 20000),
    identical(projected$timing_eligible, c(TRUE, TRUE, FALSE, FALSE, FALSE)))
  coverage_case <- "missing"
  stopifnot(grepl("Missing/failed CV worker result",
    tryCatch({ c_rerun(); "" }, error = conditionMessage), fixed = TRUE))
  coverage_case <- "duplicate"
  stopifnot(grepl("CV coverage differs",
    tryCatch({ c_rerun(); "" }, error = conditionMessage), fixed = TRUE))
}, finally = {
  unlockBinding("mclapply", ns); assign("mclapply", original, ns); lockBinding("mclapply", ns)
})
cat("K=5 methylation rerun checks passed.\n")
})
