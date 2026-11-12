# ============================================================
# Script: lb_yfb_study.R
# Purpose: Compare pooled LB fold-in and YFB projection on identical inputs.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm and project fitting/prediction modules
# ============================================================
# Serial orchestration only. Callers load fitters and predictors first.

#' Reject a comparison with different inputs, folds or initial posterior states
#' @param lb,yfb Lists with Y, fold and state entries.
#' @return Invisible TRUE when all entries are identical, otherwise an error.
#' @examples
#' # lb_yfb_study_assert_identical(design, design)
#' @family lb-yfb-study
lb_yfb_study_assert_identical <- function(lb, yfb) {
  for (nm in c("Y", "fold", "state")) {
    if (is.null(lb[[nm]]) || is.null(yfb[[nm]]) || !identical(lb[[nm]], yfb[[nm]]))
      stop("LB/YFB mismatch in ", nm, ".")
  }
  invisible(TRUE)
}

#' Scale one shared matrix using training-only cohort feature statistics
#' @param X Numeric patients by features matrix, without missing values.
#' @param cohort Cohort labels for all rows.
#' @param train Training row indices, default all patients.
#' @param scaling zscore subtracts the mean and divides by SD. sd only divides.
#' @return Matrix with the same dimensions and ordering as X.
#' @examples
#' # lb_yfb_study_scale(X, cohort, which(fold != 1), "sd")
#' @family lb-yfb-study
lb_yfb_study_scale <- function(X, cohort, train = seq_len(nrow(X)),
                               scaling = c("zscore", "sd")) {
  scaling <- match.arg(scaling)
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
      length(cohort) != nrow(X) || anyNA(cohort) || !length(train) ||
      anyNA(train) || anyDuplicated(train) || any(!train %in% seq_len(nrow(X))))
    stop("Finite matrix, cohort labels and unique valid training indices required.")
  Y <- X
  for (co in unique(cohort)) {
    idx <- which(cohort == co); tr <- intersect(idx, train)
    if (length(tr) < 2) stop("At least two training patients needed for cohort ", co)
    sdj <- apply(X[tr, , drop = FALSE], 2, stats::sd)
    if (any(!is.finite(sdj)) || any(sdj <= 0))
      stop("Zero or invalid training feature SD in cohort ", co)
    mu <- if (scaling == "zscore") colMeans(X[tr, , drop = FALSE]) else rep(0, ncol(X))
    Y[idx, ] <- sweep(sweep(X[idx, , drop = FALSE], 2, mu, "-"), 2, sdj, "/")
  }
  Y
}

#' Construct a shared canonical initial posterior with zero coefficient means
#' @param Y Shared training matrix. @param K Number of fitted programs.
#' @param prior_LF Molecular prior family.
#' @param center_projection Whether YFB centers its Cox projection.
#' @return Posterior-state list accepted by canonical YFB, also used by LB.
#' @examples
#' # lb_yfb_study_start(Y, 7, "point_laplace", FALSE)
#' @family lb-yfb-study
lb_yfb_study_start <- function(Y, K, prior_LF, center_projection) {
  if (length(K) != 1L || !is.finite(K) || K < 1 || K != as.integer(K) || K > min(dim(Y)))
    stop("K must be a positive integer at most min(dim(Y)).")
  s <- svd(Y, nu = K, nv = K); root <- sqrt(pmax(s$d[seq_len(K)], 0))
  L <- sweep(s$u, 2, root, "*"); F <- sweep(s$v, 2, root, "*")
  if (prior_LF == "point_exponential") { L <- abs(L); F <- abs(F) }
  Z <- if (center_projection) sweep(Y, 2, colMeans(Y)) else Y
  st <- canonicalize_yf_factors(Z, L, L^2, F, F^2, rep(0, K), rep(0, K), K)
  # Re-form squared moments after rescaling so both fitters receive exactly
  # the same zero-variance state, including floating-point rounding.
  list(EL = st$EL, EL2 = st$EL^2, EF = st$EF, EF2 = st$EF^2,
       EBeta = st$EBeta, EBeta2 = st$EBeta2,
       Tau = 1 / pmax(apply(Y, 2, stats::var), 1e-8),
       prior_F = vector("list", K), g_beta = NULL)
}

#' Score fixed Cox risks within cohorts and pooled across all patients
#' @param risk Fixed higher-is-worse Cox linear predictor. Never reoriented here.
#' @param time,status Outcomes. @param cohort Cohort labels, used only in scoring.
#' @return One-row data frame with C_within and C_pooled. Undefined C is NA.
#' @examples
#' # lb_yfb_study_concordance(risk, time, status, cohort)
#' @family lb-yfb-study
lb_yfb_study_concordance <- function(risk, time, status, cohort) {
  if (length(risk) != length(time) || length(status) != length(time) ||
      length(cohort) != length(time) || anyNA(cohort) ||
      any(!is.finite(c(risk, time, status))) || any(time <= 0) || any(!status %in% c(0, 1)))
    stop("Finite aligned risks, positive times, binary events and cohort labels required.")
  if (length(risk) < 2L || stats::sd(risk) == 0)
    return(data.frame(C_within = NA_real_, C_pooled = NA_real_))
  y <- survival::Surv(time, status)
  within <- survival::concordancefit(y, risk, strata = as.integer(factor(cohort)),
                                     reverse = TRUE, timefix = FALSE)$concordance
  pooled <- survival::concordancefit(y, risk, reverse = TRUE, timefix = FALSE)$concordance
  data.frame(C_within = if (is.finite(within)) within else NA_real_,
             C_pooled = if (is.finite(pooled)) pooled else NA_real_)
}

#' Correlate scores, leaving constant-program correlations undefined
#' @param a,b Numeric vectors for the same patients.
#' @return Pearson correlation or NA if either score is constant.
#' @examples
#' # lb_yfb_study_cor(1:4, 4:1)
#' @family lb-yfb-study
lb_yfb_study_cor <- function(a, b) {
  if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  stats::cor(a, b)
}

#' Fit a pooled LB/YFB pair under matched preprocessing, priors and starts
#'
#' LB uses G + w_S (S - KL_beta) - KL_L - KL_F, with no genomics tempering.
#' YFB uses canonicalize and survival_block_weight. Their historical update
#' orders remain model-specific. No pruning, cohort offsets or Cox strata.
#' @param Y Numeric shared training matrix.
#' @param time,status Survival outcomes. @param cohort Labels for scoring only.
#' @param K Fitted program count. @param prior_LF,prior_beta Prior families.
#' @param w_S Positive survival-block weight. @param max_iter,tol Fitting controls.
#' @param center_projection Center the YFB survival projection only.
#' @param fold Full split identity for the identical-design check, default zero.
#' @return Fits and initial state, predictive training risks, metrics, fitted
#'   linear-predictor correlation and within-fit correlations of L_k with (YF)_k.
#' @examples
#' # lb_yfb_study_pair(Y, time, status, cohort, 7, "point_laplace", "normal", 1, 300, 1e-5)
#' @family lb-yfb-study
#' @seealso predict_ssbmf_lb_foldin, predict_ssbmf_yfb
lb_yfb_study_pair <- function(Y, time, status, cohort, K, prior_LF, prior_beta,
                              w_S, max_iter, tol, center_projection = FALSE,
                              fold = rep(0L, nrow(Y))) {
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y)) ||
      length(time) != nrow(Y) || length(status) != nrow(Y) || length(cohort) != nrow(Y) ||
      anyNA(cohort) || any(!is.finite(c(time, status))) || any(time <= 0) ||
      any(!status %in% c(0, 1)) || length(fold) != nrow(Y) || anyNA(fold))
    stop("Invalid study matrix, outcomes, cohorts or fold labels.")
  st <- lb_yfb_study_start(Y, K, prior_LF, center_projection)
  # YFB's posterior-state interface bypasses recanonicalization and preserves
  # every supplied moment and precision. LB starts from the same means,
  # squared second moments, zero beta and sample-variance precisions.
  lb_design <- list(Y = Y, fold = fold, state = st)
  yf_design <- list(Y = Y, fold = fold, state = st)
  lb_yfb_study_assert_identical(lb_design, yf_design)
  lb <- fit_ssbmf_lb(lb_design$Y, time, status, K = K, max_iter = max_iter, tol = tol,
    prior_LF = prior_LF, prior_beta = prior_beta, w_S = w_S,
    EL_init = st$EL, EF_init = st$EF, EBeta_init = st$EBeta, verbose = FALSE)
  # Check the actual LB initialization as well as the requested arguments.
  lb_yfb_study_assert_identical(lb_design, list(Y = Y, fold = fold, state = lb$initial_state))
  yf <- fit_ssbmf_yfb(yf_design$Y, time, status, K = K, max_iter = max_iter, tol = tol,
    prior_LF = prior_LF, prior_beta = prior_beta, canonicalize = TRUE,
    posterior_init = yf_design$state, w_S = w_S, survival_block_weight = TRUE,
    center_projection = center_projection, verbose = FALSE)
  fits <- list(LB = lb, YFB = yf)
  fi <- predict_ssbmf_lb_foldin(Y, lb)
  py <- predict_ssbmf_yfb(Y, yf$EF, yf$EBeta, yf$EF_norms,
                         projection_center = yf$projection_center)
  predictions <- list(LB = fi$risk_scores, YFB = py$risk_scores)
  metrics <- do.call(rbind, lapply(names(fits), function(model) {
    f <- fits[[model]]
    cbind(data.frame(model = model, converged = isTRUE(f$history$converged),
                     n_iter = f$history$n_iter,
                     foldin_converged = if (model == "LB") isTRUE(fi$converged) else TRUE),
          lb_yfb_study_concordance(predictions[[model]], time, status, cohort))
  }))
  eta_lb <- drop(lb$EL %*% lb$EBeta)
  program <- do.call(rbind, lapply(names(fits), function(model) {
    f <- fits[[model]]; yf_score <- Y %*% f$EF
    data.frame(model = model, program = seq_len(K),
      cor_L_YF = vapply(seq_len(K), function(k) lb_yfb_study_cor(f$EL[, k], yf_score[, k]), numeric(1)))
  }))
  list(fits = fits, start = st, predictions = predictions, metrics = metrics,
       comparison = data.frame(eta_correlation = lb_yfb_study_cor(eta_lb, py$risk_scores)),
       program = program)
}

#' Five-fold pooled fitting with LB fold-in and YFB projection evaluation
#' @param X Raw shared patient-by-feature matrix.
#' @param time,status Outcomes. @param cohort Labels for folds and evaluation only.
#' @param fold Integer labels 1 to 5, with every patient assigned once.
#' @param scaling Input scaling, zscore or sd, learned within each training fold.
#' @param prior_LF,prior_beta Prior families. @param K Fitted program count.
#' @param w_S Positive survival-block weight. @param max_iter,tol Fitting controls.
#' @param checkpoint NULL or an existing ignored directory for per-fold checkpoints.
#' @return Fold metrics, paired summary with defined-count flags, and private
#'   paired fits and held-out scores. Save patient-level fields under data/ only.
#' @examples
#' # lb_yfb_study_cv(X, time, status, cohort, fold, "sd", "point_exponential", "normal", 7, 1, 300, 1e-5)
#' @family lb-yfb-study
lb_yfb_study_cv <- function(X, time, status, cohort, fold, scaling, prior_LF,
                            prior_beta, K, w_S, max_iter, tol, checkpoint = NULL) {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
      length(time) != nrow(X) || length(status) != nrow(X) ||
      length(cohort) != nrow(X) || anyNA(cohort) ||
      any(!is.finite(c(time, status))) || any(time <= 0) || any(!status %in% c(0, 1)))
    stop("CV requires a finite matrix and exactly aligned outcomes and cohort labels.")
  if (length(fold) != nrow(X) || anyNA(fold) || !setequal(fold, 1:5))
    stop("Five nonempty folds numbered 1 through 5 are required.")
  if (!is.null(checkpoint) && !dir.exists(checkpoint)) stop("Checkpoint directory must exist.")
  pairs <- vector("list", 5); rows <- vector("list", 5)
  for (j in seq_len(5)) {
    tr <- which(fold != j); te <- which(fold == j)
    result <- tryCatch({
    Y <- lb_yfb_study_scale(X, cohort, train = tr, scaling = scaling)
    z <- lb_yfb_study_pair(Y[tr, , drop = FALSE], time[tr], status[tr], cohort[tr],
      K, prior_LF, prior_beta, w_S, max_iter, tol, identical(scaling, "sd"), fold[tr])
    lb <- predict_ssbmf_lb_foldin(Y[te, , drop = FALSE], z$fits$LB)
    f <- z$fits$YFB
    yf <- predict_ssbmf_yfb(Y[te, , drop = FALSE], f$EF, f$EBeta, f$EF_norms,
                           projection_center = f$projection_center)
    z$heldout <- list(LB = lb$risk_scores, YFB = yf$risk_scores)
    z$test_indices <- te
    fold_metrics <- do.call(rbind, lapply(c("LB", "YFB"), function(model) {
      m <- z$metrics[z$metrics$model == model, c("model", "converged", "n_iter"), drop = FALSE]
      cbind(data.frame(fold = j), m,
        data.frame(foldin_converged = if (model == "LB") isTRUE(lb$converged) else TRUE,
                   error = ""),
        lb_yfb_study_concordance(z$heldout[[model]], time[te], status[te], cohort[te]))
    }))
    list(pair = z, metrics = fold_metrics)
    }, error = function(e) {
      msg <- conditionMessage(e)
      warning(sprintf("Study fold %d failed: %s", j, msg))
      list(pair = list(error = msg, test_indices = te),
           metrics = data.frame(fold = j, model = c("LB", "YFB"), converged = FALSE,
             n_iter = NA_integer_, foldin_converged = FALSE, error = msg,
             C_within = NA_real_, C_pooled = NA_real_))
    })
    pairs[[j]] <- result$pair; rows[[j]] <- result$metrics
    if (!is.null(checkpoint))
      saveRDS(result, file.path(checkpoint, sprintf("fold_%d.rds", j)))
  }
  metrics <- do.call(rbind, rows)
  # Do not silently remove undefined folds. Report the defined count, but the
  # requested five-fold mean and SE are NA unless all five estimates exist.
  summary <- do.call(rbind, lapply(c("C_within", "C_pooled"), function(metric) {
    a <- metrics[metrics$model == "LB", metric]; b <- metrics[metrics$model == "YFB", metric]
    delta <- b - a; ok <- is.finite(delta)
    data.frame(metric = metric, defined_LB = sum(is.finite(a)), defined_YFB = sum(is.finite(b)),
      defined_pairs = sum(ok), mean_LB = mean(a), mean_YFB = mean(b),
      YFB_minus_LB = mean(delta), paired_SE = stats::sd(delta) / sqrt(5))
  }))
  list(metrics = metrics, summary = summary, pairs = pairs,
       failed = any(nzchar(metrics$error)))
}
