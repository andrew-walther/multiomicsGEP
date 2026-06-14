# ============================================================
# Script:  results/multi_cohort_sim/sim_scoring.R
# Purpose: Scoring helpers for the multi-cohort simulation.  Quantify how well
#          a fitted model recovers (a) the gene programs, (b) the shared vs.
#          study-specific structure, and (c) the survival coefficients.
#
#          All three metrics are computed against the known ground truth
#          returned by generate_multicohort_data().
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-14
# Dependencies: base R (stats::cor)
# ============================================================

# match_factors ----
#' Greedy best-|correlation| matching of estimated to true gene programs.
#'
#' Both F_hat and F_true are p × K matrices (gene space), so the metric is
#' comparable across LB / YFB / EBMF.  For each TRUE factor (column of F_true)
#' we assign the unused estimated factor with the highest |Pearson cor|.
#'
#' @param F_hat   p × K_hat estimated gene programs.
#' @param F_true  p × K_true true gene programs.
#' @return list(match = K_true int vector of assigned estimated-column indices,
#'   best_cor = K_true numeric vector of achieved |cor|).
#' @family multicohort-sim
match_factors <- function(F_hat, F_true) {
  Cmat   <- abs(suppressWarnings(cor(F_hat, F_true)))   # K_hat × K_true
  Cmat[is.na(Cmat)] <- 0                                # constant columns -> 0 cor
  Ktrue  <- ncol(F_true)
  avail  <- rep(TRUE, nrow(Cmat))
  m      <- rep(NA_integer_, Ktrue)   # NA = no estimated factor left to match (K_hat < K_true)
  best   <- numeric(Ktrue)            # 0 for unmatched true factors
  # assign true factors in descending order of their best achievable correlation
  for (j in order(apply(Cmat, 2, max), decreasing = TRUE)) {
    cand <- which(avail)
    if (length(cand) == 0) next       # estimated factors exhausted -> leave NA / 0
    i    <- cand[which.max(Cmat[cand, j])]
    m[j] <- i
    best[j] <- Cmat[i, j]
    avail[i] <- FALSE
  }
  list(match = m, best_cor = best)
}

# classify_specificity ----
#' Label each ESTIMATED factor shared/specific from its per-cohort loading energy.
#'
#' For estimated loadings L_hat (n × K_hat) and a cohort_id factor, compute each
#' factor's fraction of squared-loading energy in each cohort.  A factor whose
#' energy is concentrated (> 1 - tol) in one cohort is labelled specific to that
#' cohort; otherwise shared.
#'
#' @param L_hat      n × K_hat estimated loadings (rows aligned to cohort_id).
#' @param cohort_id  length-n factor of cohort membership.
#' @param tol        concentration tolerance (default 0.15 -> >85% in one cohort).
#' @return character vector of length K_hat: "shared" or "specific_<c>".
#' @family multicohort-sim
classify_specificity <- function(L_hat, cohort_id, tol = 0.15) {
  cohort_id <- factor(cohort_id)
  vapply(seq_len(ncol(L_hat)), function(k) {
    e <- tapply(L_hat[, k]^2, cohort_id, sum)
    e <- e / (sum(e) + 1e-12)                 # fraction of energy per cohort
    if (max(e) > 1 - tol) paste0("specific_", levels(cohort_id)[which.max(e)])
    else                  "shared"
  }, character(1))
}

# specificity_accuracy ----
#' Accuracy of estimated specificity labels against the truth, via the factor
#' matching.  Each true factor j is compared to the estimated factor matched to
#' it (match[j]); "shared" vs "specific" agreement is scored (cohort identity of
#' specific factors is also required to match).
#'
#' @param est_labels    character vector (length K_hat) from classify_specificity.
#' @param match         K_true int vector from match_factors$match.
#' @param true_labels   character vector (length K_true) of ground-truth labels.
#' @return list(accuracy, shared_recall, specific_recall, per_factor data.frame).
#' @family multicohort-sim
specificity_accuracy <- function(est_labels, match, true_labels) {
  est_for_true <- est_labels[match]                    # estimated label aligned to truth (NA if unmatched)
  # coarse shared/specific agreement (ignoring which cohort); unmatched (NA) counts as wrong
  coarse_true <- ifelse(true_labels == "shared", "shared", "specific")
  coarse_est  <- ifelse(is.na(est_for_true), "unmatched",
                        ifelse(est_for_true == "shared", "shared", "specific"))
  agree <- coarse_true == coarse_est
  is_sh <- coarse_true == "shared"
  list(
    accuracy        = mean(agree),
    shared_recall   = if (any(is_sh))  mean(agree[is_sh])  else NA_real_,
    specific_recall = if (any(!is_sh)) mean(agree[!is_sh]) else NA_real_,
    per_factor = data.frame(true = true_labels, est = est_for_true,
                            correct = agree, stringsAsFactors = FALSE)
  )
}

# beta_recovery ----
#' Survival-coefficient recovery, grouped by the TRUE factor label.
#'
#' Using the factor matching, |β̂| is read off for the estimated factor matched
#' to each true factor.  Reports the magnitude on shared vs. specific factors and
#' the true-/false-positive prognostic rates against a threshold.
#'
#' @param EBeta         K_hat estimated coefficients (NULL for EBMF -> all NA).
#' @param match         K_true int vector from match_factors$match.
#' @param true_labels   character vector (length K_true) of ground-truth labels.
#' @param thresh        |β| threshold for "prognostic" (globals beta_threshold).
#' @return list(beta_shared, beta_specific, tp_rate, fp_rate,
#'   mean_abs_shared, mean_abs_specific).
#' @family multicohort-sim
beta_recovery <- function(EBeta, match, true_labels, thresh) {
  is_sh <- true_labels == "shared"
  if (is.null(EBeta)) {                                # EBMF: survival-blind
    return(list(beta_shared = NA_real_, beta_specific = NA_real_,
                tp_rate = NA_real_, fp_rate = NA_real_,
                mean_abs_shared = NA_real_, mean_abs_specific = NA_real_))
  }
  bhat <- abs(EBeta[match])                            # aligned to true factor order (NA if unmatched)
  bhat[is.na(bhat)] <- 0                               # unrecovered factor -> not prognostic
  list(
    beta_shared       = bhat[is_sh],
    beta_specific     = bhat[!is_sh],
    tp_rate           = if (any(is_sh))  mean(bhat[is_sh]  > thresh) else NA_real_,
    fp_rate           = if (any(!is_sh)) mean(bhat[!is_sh] > thresh) else NA_real_,
    mean_abs_shared   = if (any(is_sh))  mean(bhat[is_sh])  else NA_real_,
    mean_abs_specific = if (any(!is_sh)) mean(bhat[!is_sh]) else NA_real_
  )
}

# oriented_cindex ----
#' Orientation-free held-out concordance, max(C, 1 - C).
#'
#' Matches the convention in run_synthetic.R: the sign of the risk score is
#' resolved by taking the better of the two orientations, so a correctly-fit but
#' sign-flipped predictor is not penalised.
#'
#' @param risk   n_test vector of risk scores (higher = worse prognosis).
#' @param time   n_test survival/censoring times.
#' @param status n_test event indicators.
#' @return scalar C-index in [0.5, 1].
#' @family multicohort-sim
oriented_cindex <- function(risk, time, status) {
  if (all(is.na(risk)) || sd(risk) < 1e-12) return(0.5)
  cc <- as.numeric(survival::concordance(
    survival::Surv(time, status) ~ risk)$concordance)
  max(cc, 1 - cc)
}
