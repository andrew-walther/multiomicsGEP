# ============================================================
# Script: audit_share_targeted_ws.R
# Purpose: Audit saved S-16 precision shares without fitting models.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: survival, ebnm, yaml, existing YFB simulation helpers
# ============================================================
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
started <- proc.time()[["elapsed"]]
source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
s9_load()
source("code/share_targeted_ws.R")
out <- "results/benchmark_sim/outputs/share_targeted_ws_2026-11-13"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

#' Summarize saved-state shares and their approximations
#' @param fit Converged canonical fit or reconstructed reduced state.
#' @param Y,time,status Training matrix and survival outcomes.
#' @param label Source label. @param selected Program indices to report.
#' @return Aggregate data frame without patient-level information.
#' @examples
#' # audit_one(fit, Y, time, status, "S9", 1L)
audit_one <- function(fit, Y, time, status, label, selected) {
  x <- share_targeted_ws(fit, Y, time, status, .5)
  a <- x$diagnostics[selected, , drop = FALSE]
  a$source <- label
  a$loading_second_moment_bound <- x$loading_second_moment_bound
  a$directional_fold_error <- pmax(a$directional_ratio / a$closed_directional_ratio,
                                    a$closed_directional_ratio / a$directional_ratio)
  exact_gene_ratio <- a$per_gene_median_share / (1 - a$per_gene_median_share)
  a$gene_fold_error <- pmax(exact_gene_ratio / a$closed_gene_ratio,
                             a$closed_gene_ratio / exact_gene_ratio)
  a
}

# S-9: seed 1 at each variance, with the matched truth program.
rows <- list()
for (v in c(.0025, .01, .05)) {
  path <- sprintf("results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid/p2000_v%s_s01_w1.rds", v)
  ck <- readRDS(path)
  if (!isTRUE(ck$fit$history$converged)) stop(sprintf("S9 variance %s did not converge", v))
  d <- yfb_decision_data(1, 2000L, v, truth_mode = "shared")
  Y <- s9_zscore(d$Y)
  r <- abs(cor(d$truth, ck$fit$EF)); r[!is.finite(r)] <- 0
  k <- which.max(r)
  rows[[length(rows) + 1L]] <- audit_one(ck$fit, Y, d$time, d$event,
                                         sprintf("S9_v%s_seed1_w1", v), k)
}

# K=5 S-11: the registered EL² bound, not the S-15 state.
z <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
Y <- z$Y[, match(z$sets$D4, z$genes), drop = FALSE]
ck <- readRDS("results/benchmark_sim/outputs/k5_working_2026-11-12/s11/fits/size_D4_N2218_w1_f0.rds")
if (!isTRUE(ck$metrics$converged)) stop("K5 S11 checkpoint did not converge")
fit <- ck$state
fit$history <- list(converged = ck$metrics$converged)
fit$canonicalize <- TRUE
fit$w_S <- 1
beta_sd <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
active <- which(abs(fit$EBeta) >= 1.96 * beta_sd & beta_sd > 0 &
                  colSums(fit$EF^2) > 1e-8)
rows[[length(rows) + 1L]] <- audit_one(fit, Y, z$time, z$status, "K5_S11_w1_bound", active)
write.csv(do.call(rbind, rows), file.path(out, "saved_share_diagnostics.csv"), row.names = FALSE)


#' Count retained checkpoint fields without extracting patient-level values
#' @param directory Checkpoint directory. @param source Source label.
#' @param full_fit Whether state is stored under fit rather than state.
#' @return One aggregate inventory row.
#' @examples
#' # inventory_checkpoints("path/to/fits", "S13", FALSE)
inventory_checkpoints <- function(directory, source, full_fit = FALSE) {
  files <- list.files(directory, pattern = "[.]rds$", full.names = TRUE)
  if (!length(files)) return(data.frame(source = source, files = 0L, state = 0L, EL = 0L, EL2 = 0L, EF = 0L, EBeta = 0L, EBeta2 = 0L, Tau = 0L, full_trace = 0L, final_shareA = 0L))
  counts <- lapply(files, function(path) {
    x <- readRDS(path)
    z <- if (full_fit) x$fit else x$state
    c(state = !is.null(z), EL = !is.null(z$EL), EL2 = !is.null(z$EL2),
      EF = !is.null(z$EF), EBeta = !is.null(z$EBeta), EBeta2 = !is.null(z$EBeta2),
      Tau = !is.null(z$Tau), full_trace = !is.null(z$history$trace$shareA_median),
      final_shareA = !is.null(z$shareA))
  })
  data.frame(source = source, files = length(files),
             as.list(colSums(do.call(rbind, counts))), check.names = FALSE)
}
dirs <- c(S9_grid = "results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid",
          S13 = "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05/fits",
          S15_K7 = "results/benchmark_sim/outputs/yfb_supervised_refit_2026-11-06/fits",
          S15_K5 = "results/benchmark_sim/outputs/k5_working_2026-11-12/s15/fits",
          S11_K5 = "results/benchmark_sim/outputs/k5_working_2026-11-12/s11/fits")
inventory <- do.call(rbind, lapply(names(dirs), function(nm)
  inventory_checkpoints(dirs[[nm]], nm, full_fit = nm == "S9_grid")))
inventory$limitation <- c("Final trace uses intra-sweep state",
  "Tau and trace absent; L_norm holds sum EL squared",
  "Tau absent; inner states absent", "Tau absent; inner states absent",
  "EL only for full-data fits; EL squared bound; separate fit")
write.csv(inventory, file.path(out, "checkpoint_inventory.csv"), row.names = FALSE)
write.csv(data.frame(elapsed_seconds = proc.time()[["elapsed"]] - started,
                     n_new_fits = 0L), file.path(out, "audit_run.csv"), row.names = FALSE)
