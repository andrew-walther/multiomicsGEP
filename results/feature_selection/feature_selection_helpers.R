# ============================================================
# Script: feature_selection_helpers.R
# Purpose: Evaluation and gene-screening helpers for the survival-screen experiment.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: survival
# ============================================================
# Interfaces follow Task E1 of docs/plans/2026-10-16-feature-selection-experiment-plan.md.
# They live beside the runner rather than in code/ so this experiment does not touch
# files that the current-method work edits in parallel; move them to code/ if the
# experiment is kept.

# Concordance ----

#' Harrell's C for a risk score, with no direction flipping
#'
#' Higher risk means shorter survival. C = P(risk_i > risk_j | T_i < T_j, δ_i = 1).
#' @param risk Numeric risk score.
#' @param time,status Survival outcome (status 1 = event).
#' @return Scalar C, or NA_real_ when the risk is constant.
#' @examples harrell_c(c(3, 2, 1), c(1, 2, 3), c(1, 1, 1))  # 1
#' @family feature_selection_experiment
harrell_c <- function(risk, time, status) {
  if (stats::sd(risk) == 0) return(NA_real_)
  as.numeric(survival::concordance(survival::Surv(time, status) ~ risk,
                                   reverse = TRUE)$concordance)
}

#' Risk direction of a score, decided on training patients only
#'
#' d = +1 if C(η) ≥ 0.5 on the training patients, else −1; held-out risk is d·η.
#' @param eta Linear predictor on training patients.
#' @param time,status Training survival outcome.
#' @return 1L, -1L, or NA_integer_ when η is constant.
#' @examples training_direction(c(1, 2, 3), c(1, 2, 3), c(1, 1, 1))  # -1
#' @family feature_selection_experiment
training_direction <- function(eta, time, status) {
  C <- harrell_c(eta, time, status)
  if (is.na(C)) return(NA_integer_)
  if (C >= 0.5) 1L else -1L
}

# Gene screens ----

#' Rank genes by a cohort-stratified univariate Cox test
#'
#' For each column y_j fits coxph(Surv(time, status) ~ y_j + strata(cohort)) and
#' ranks by the Wald p-value, ascending. A gene that is constant or whose fit fails
#' or does not converge gets p = 1 and is counted in attr "n_failed".
#' @param Y Patients × genes matrix (fold-training rows only).
#' @param time,status Survival outcome for the same rows.
#' @param strata Cohort label per row (or NULL for an unstratified test).
#' @param n_keep Number of genes to return; must not exceed ncol(Y).
#' @return Integer column indices (length n_keep), with attributes "p_values"
#'   (length ncol(Y)) and "n_failed".
#' @examples # survival_screen_genes(Y, time, status, cohort, 500)
#' @family feature_selection_experiment
survival_screen_genes <- function(Y, time, status, strata, n_keep) {
  if (!is.matrix(Y) || n_keep < 1L || n_keep != as.integer(n_keep))
    stop("Y must be a matrix and n_keep a positive integer.")
  if (n_keep > ncol(Y))
    stop(sprintf("n_keep = %d exceeds the %d candidate genes.", n_keep, ncol(Y)))
  p <- rep(1, ncol(Y)); failed <- 0L
  for (j in seq_len(ncol(Y))) {
    y <- Y[, j]
    ok <- if (!is.null(strata)) all(tapply(y, strata, stats::sd) > 0) else stats::sd(y) > 0
    if (!isTRUE(ok)) { failed <- failed + 1L; next }
    # A convergence warning (e.g. infinite coefficient) counts as a failed fit.
    fit <- tryCatch(
      if (is.null(strata)) survival::coxph(survival::Surv(time, status) ~ y)
      else survival::coxph(survival::Surv(time, status) ~ y + strata(strata)),
      warning = function(w) NULL, error = function(e) NULL)
    if (is.null(fit) || !is.finite(stats::coef(fit)[1])) { failed <- failed + 1L; next }
    z <- stats::coef(fit)[1] / sqrt(stats::vcov(fit)[1, 1])
    p[j] <- 2 * stats::pnorm(-abs(z))
  }
  idx <- order(p, seq_along(p))[seq_len(n_keep)]
  attr(idx, "p_values") <- p
  attr(idx, "n_failed") <- failed
  idx
}

#' Smallest nested unsupervised gene set with at least n_shared genes
#'
#' Uses the entry ranks of feature_count_manifest() (a gene enters the per-cohort
#' intersection at the largest of its combined mean/variance rank positions). The
#' per-cohort screen size q is the n_shared-th smallest entry rank, so the set is
#' {genes : entry ≤ q}; ties can make it slightly larger than n_shared.
#' @param entry_rank Named integer vector of entry ranks (names are genes).
#' @param n_shared Requested number of shared genes.
#' @return List: genes, top_n_per_cohort (q), n_shared_actual.
#' @examples unsupervised_gene_set(c(a = 1L, b = 3L, c = 2L), 2)  # a, c
#' @family feature_selection_experiment
unsupervised_gene_set <- function(entry_rank, n_shared) {
  if (n_shared > length(entry_rank)) stop("n_shared exceeds the available genes.")
  q <- sort(entry_rank)[n_shared]
  genes <- names(entry_rank)[entry_rank <= q]
  list(genes = genes, top_n_per_cohort = unname(q), n_shared_actual = length(genes))
}
