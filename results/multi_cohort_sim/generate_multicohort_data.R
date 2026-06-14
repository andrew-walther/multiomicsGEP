# ============================================================
# Script:  results/multi_cohort_sim/generate_multicohort_data.R
# Purpose: Data-generating process for the multi-cohort simulation study.
#          Generates genomics + survival data with a KNOWN partition of
#          latent factors into shared (present in all cohorts, prognostic)
#          and study-specific (private to one cohort, non-prognostic).
#
#          Model:
#            Y      = L Fᵀ + offset + E,   E_ij ~ N(0, 1/τ_j), τ_j ~ Gamma(2,2)
#            η_i    = Σ_{k∈shared} L_ik β_k          (specific factors: β = 0)
#            T_i    = (−log U / (scale0 · exp(η)))^(1/shape)   (Weibull-PH)
#
#          Specificity is encoded by BLOCK-ZERO loadings: a cohort-c-specific
#          factor has L_ik = 0 for every patient i not in cohort c.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-14
# Dependencies: base R only (flashier templates passed in optionally)
# ============================================================

# .calibrate_censor_mc ----
#' Bisection calibration of the censoring scale to a target censoring fraction.
#'
#' Finds a scalar multiplier `s` so that the fraction of patients with
#' `base * s < event_time` (i.e. censored) is approximately `target`.  Copied
#' from the established single-cohort DGP (`run_synthetic.R`) for consistency.
#'
#' @param ev      n-vector of event (failure) times.
#' @param base    n-vector of baseline censoring draws (Exp(1)).
#' @param target  target censoring fraction in (0, 1).
#' @param n_iter  bisection iterations.
#' @return scalar censoring scale `s`.
.calibrate_censor_mc <- function(ev, base, target = 0.30, n_iter = 40) {
  lo <- 1e-3; hi <- max(ev) * 10
  for (i in seq_len(n_iter)) {
    mid <- sqrt(lo * hi)                       # geometric-mean bisection (positive scale)
    if (mean(base * mid < ev) > target) lo <- mid else hi <- mid
  }
  hi
}

# generate_multicohort_data ----
#' Generate a multi-cohort genomics + survival dataset with known shared /
#' study-specific factor structure.
#'
#' Rows of Y and L are stacked by cohort: rows 1..n_per[1] are cohort 1, the
#' next n_per[2] rows are cohort 2, etc.  Factor columns are ordered
#' [shared..., specific_1..., specific_2..., ...].
#'
#' @param C            number of cohorts (default 2).
#' @param n_per        length-C integer vector of patients per cohort.
#' @param p            number of genes.
#' @param K_shared     number of shared factors (β ≠ 0).
#' @param K_specific   length-C integer vector: specific factors per cohort (β = 0).
#' @param F_templates  p × K_ebmf matrix of EBMF gene programs, or NULL to use a
#'                     synthetic sparse-positive fallback.  Columns are recycled
#'                     by index and renormalised to unit L2 norm.
#' @param a_shared     amplitude multiplier for shared F columns.
#' @param a_specific   amplitude multiplier for specific F columns.  The genomics
#'                     variance ratio scales as (a_shared / a_specific)^2.
#' @param active_rate  synthetic-fallback sparsity (fraction of genes active).
#' @param beta_shared  length-K_shared coefficient vector, or NULL to recycle
#'                     c(1.5, -1.2, 0.8, -0.5).
#' @param specific_prognostic  if TRUE, the FIRST specific factor of each cohort
#'                     also receives a non-zero β (within-cohort prognostic signal).
#'                     Default FALSE (faithful to the meeting notes: survival from
#'                     shared factors only).
#' @param offset_sd    SD of an explicit per-cohort platform offset (0 = off).
#' @param shape,scale0 Weibull-PH baseline parameters.
#' @param target_censoring  target censoring fraction.
#' @param seed         RNG seed.
#'
#' @return list with: Y (n×p), time, status, cohort_id (factor), L_true (n×K),
#'   F_true (p×K), beta_true (K), factor_labels (K char: "shared"/"specific_c"),
#'   factor_owner (K int: 0 = shared, c = owned by cohort c), eta (n), params.
#' @family multicohort-sim
generate_multicohort_data <- function(
    C                   = 2L,
    n_per               = rep(150L, C),
    p                   = 1000L,
    K_shared            = 2L,
    K_specific          = rep(2L, C),
    F_templates         = NULL,
    a_shared            = 1.0,
    a_specific          = 1.0,
    active_rate         = 0.05,
    beta_shared         = NULL,
    specific_prognostic = FALSE,
    offset_sd           = 0.0,
    shape               = 1.5,
    scale0              = 0.01,
    target_censoring    = 0.30,
    seed                = 1L) {

  # --- validation: fail loud on malformed factor partition ---
  if (length(n_per) != C)      stop("n_per must have length C.")
  if (length(K_specific) != C) stop("K_specific must have length C.")

  set.seed(seed)
  n  <- sum(n_per)
  K  <- K_shared + sum(K_specific)
  if (K < 1) stop("Total number of factors K must be >= 1.")
  cohort_id <- factor(rep(seq_len(C), times = n_per))

  # --- per-column label and owning cohort (owner 0 = shared) ---
  labels <- c(rep("shared", K_shared),
              unlist(lapply(seq_len(C),
                            function(cc) rep(paste0("specific_", cc), K_specific[cc]))))
  owner  <- c(rep(0L, K_shared),
              unlist(lapply(seq_len(C),
                            function(cc) rep(cc, K_specific[cc]))))

  # --- loadings L: Exponential(1), block-zeroed outside the owning cohort ---
  # E[L]=1, E[L^2]=2; matches the point_exponential prior used by the fitters.
  L <- matrix(rexp(n * K, rate = 1), n, K)
  for (k in seq_len(K)) {
    if (owner[k] != 0L) L[cohort_id != owner[k], k] <- 0   # study-specific: zero elsewhere
  }

  # --- gene programs F: EBMF templates (preferred) or synthetic sparse ---
  # Each column is normalised to unit L2 norm so the amplitude knob sets scale.
  Fm <- matrix(0, p, K)
  if (!is.null(F_templates)) {
    if (nrow(F_templates) != p)
      stop(sprintf("F_templates has %d rows but p = %d.", nrow(F_templates), p))
    idx <- ((seq_len(K) - 1L) %% ncol(F_templates)) + 1L   # recycle columns by index
    Fm  <- F_templates[, idx, drop = FALSE]
    Fm  <- sweep(Fm, 2, sqrt(colSums(Fm^2)) + 1e-10, "/")
  } else {
    for (k in seq_len(K)) {
      active <- sample.int(p, max(1L, round(active_rate * p)))
      Fm[active, k] <- 1 / sqrt(length(active))            # unit-norm sparse column
    }
  }
  amp <- ifelse(labels == "shared", a_shared, a_specific)  # signal-strength knob
  Fm  <- sweep(Fm, 2, amp, "*")

  # --- per-gene Gaussian noise: precision τ_j ~ Gamma(2,2), noise SD = 1/√τ_j ---
  tau <- rgamma(p, shape = 2, rate = 2)
  E   <- sweep(matrix(rnorm(n * p), n, p), 2, sqrt(tau), "/")

  # --- optional explicit per-cohort platform offset (off by default) ---
  offset <- matrix(0, n, p)
  if (offset_sd > 0) {
    for (cc in seq_len(C)) {
      ov <- rnorm(p, 0, offset_sd)                         # one offset vector per cohort
      rows <- which(cohort_id == cc)
      offset[rows, ] <- matrix(ov, length(rows), p, byrow = TRUE)
    }
  }

  # --- assemble genomics + column-center (standard preprocessing) ---
  Y <- L %*% t(Fm) + offset + E
  Y <- sweep(Y, 2, colMeans(Y), "-")
  dimnames(Y) <- NULL    # avoid duplicate gene names leaking from F_templates

  # --- survival coefficients: shared factors carry β; specific are 0 by default ---
  if (is.null(beta_shared)) {
    beta_shared <- if (K_shared > 0) rep_len(c(1.5, -1.2, 0.8, -0.5), K_shared) else numeric(0)
  } else if (length(beta_shared) != K_shared) {
    stop("beta_shared must have length K_shared.")
  }
  beta <- numeric(K)
  beta[labels == "shared"] <- beta_shared

  # optional within-cohort prognostic specific factor (first specific per cohort)
  if (specific_prognostic) {
    for (cc in seq_len(C)) {
      first_spec <- which(owner == cc)[1]
      if (!is.na(first_spec)) beta[first_spec] <- 0.8
    }
  }

  eta <- as.vector(L %*% beta)                             # == 0 in the nothing-shared scenario

  # --- Weibull-PH event times + bisection-calibrated censoring (~target) ---
  event       <- (-log(runif(n)) / (scale0 * exp(eta)))^(1 / shape)
  censor_base <- rexp(n, rate = 1)
  cs          <- .calibrate_censor_mc(event, censor_base, target = target_censoring)
  censor_time <- censor_base * cs
  time        <- pmin(event, censor_time)
  status      <- as.integer(event <= censor_time)

  list(
    Y = Y, time = time, status = status, cohort_id = cohort_id,
    L_true = L, F_true = Fm, beta_true = beta, eta = eta,
    factor_labels = labels, factor_owner = owner,
    params = list(C = C, n_per = n_per, p = p, K = K,
                  K_shared = K_shared, K_specific = K_specific,
                  a_shared = a_shared, a_specific = a_specific,
                  offset_sd = offset_sd, specific_prognostic = specific_prognostic,
                  censoring_rate = mean(status == 0), seed = seed)
  )
}
