# ============================================================
# Script: real_helpers.R
# Purpose: Shared helpers for the ssbmf real-data runners: loading the
#          package, concordance with a fixed direction, and donor-level paired
#          bootstrap comparisons of two risk scores.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: survival
# ============================================================

#' Load ssbmf: the installed copy in SSBMF_LIB if set, else the source tree
load_ssbmf <- function() {
  lib <- Sys.getenv("SSBMF_LIB")
  if (nzchar(lib)) {
    .libPaths(c(lib, .libPaths()))
    library(ssbmf)
  } else {
    suppressMessages(pkgload::load_all("ssbmf", quiet = TRUE))
  }
}

#' Harrell's C with the direction fixed in advance: higher risk = shorter survival
cindex <- function(risk, time, status) {
  unname(survival::concordance(survival::Surv(time, status) ~ risk, reverse = TRUE)$concordance)
}

#' Donor-level paired bootstrap of C differences across cohorts
#'
#' Resamples donors with replacement within each cohort group; a donor with
#' records in several cohorts (PACA-AU array and RNA-seq) is drawn once and
#' contributes all of their records. Returns per-cohort and mean-over-cohort
#' differences C(new) - C(reference) with percentile intervals.
#'
#' @param cohorts named list; each element has `risk_new`, `risk_ref`, `time`,
#'   `status`, `donor` (character IDs; shared IDs across cohorts mark the same person)
#' @param B bootstrap replicates
#' @param seed random seed
#' @return list(`per_cohort` data.frame, `mean` one-row data.frame, `draws` matrix)
paired_donor_bootstrap <- function(cohorts, B = 2000, seed = 1) {
  set.seed(seed)
  donors <- unique(unlist(lapply(cohorts, `[[`, "donor")))
  # group donors so that each cohort keeps its size in expectation: a donor's
  # group is the set of cohorts they appear in
  grp <- vapply(donors, function(d) paste(names(cohorts)[vapply(cohorts, function(x) d %in% x$donor, TRUE)],
                                          collapse = "+"), "")
  idx_by_cohort <- lapply(cohorts, function(x) split(seq_along(x$donor), x$donor))
  stat <- function(dsample) {
    vapply(names(cohorts), function(co) {
      x <- cohorts[[co]]
      rows <- unlist(idx_by_cohort[[co]][dsample[dsample %in% names(idx_by_cohort[[co]])]], use.names = FALSE)
      cindex(x$risk_new[rows], x$time[rows], x$status[rows]) -
        cindex(x$risk_ref[rows], x$time[rows], x$status[rows])
    }, 0)
  }
  obs <- vapply(names(cohorts), function(co) {
    x <- cohorts[[co]]; cindex(x$risk_new, x$time, x$status) - cindex(x$risk_ref, x$time, x$status)
  }, 0)
  draws <- do.call(rbind, lapply(seq_len(B), function(b) {
    ds <- unlist(lapply(split(donors, grp), function(g) sample(g, length(g), replace = TRUE)), use.names = FALSE)
    stat(ds)
  }))
  colnames(draws) <- names(cohorts)
  ci <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
  per <- data.frame(cohort = names(cohorts), diff = unname(obs),
                    lower = apply(draws, 2, function(v) ci(v)[1]), upper = apply(draws, 2, function(v) ci(v)[2]))
  m <- rowMeans(draws)
  list(per_cohort = per,
       mean = data.frame(cohort = "mean over cohorts", diff = mean(obs), lower = ci(m)[1], upper = ci(m)[2]),
       draws = draws)
}

#' Score an external cohort when some training features are missing
#'
#' Missing features are filled with their training mean, so they contribute
#' nothing to the projection (the frozen benchmark likewise dropped them).
score_external <- function(fit, Y_ext, block = names(fit$features)[1]) {
  feats <- fit$features[[block]]
  Yf <- matrix(rep(fit$center[[block]], each = nrow(Y_ext)), nrow(Y_ext), length(feats),
               dimnames = list(rownames(Y_ext), feats))
  common <- intersect(feats, colnames(Y_ext))
  Yf[, common] <- Y_ext[, common]
  nd <- stats::setNames(list(Yf), block)
  predict(fit, newdata = nd, orientation = "samples_by_features")
}
