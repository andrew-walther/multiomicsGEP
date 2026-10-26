# Fitting ----

#' Fit survival-supervised Bayesian matrix factorization
#'
#' Fits the model of the package documentation ([ssbmf-package]) by
#' coordinate ascent on one variational objective: the expected Gaussian
#' log-likelihood of each block, plus the Gaussian-moment bound on the
#' expected Cox log partial likelihood, minus KL divergences to empirical
#' Bayes priors. Every accepted step increases the objective; the number of
#' steps that lowered it by more than 1e-8 x |objective| is reported in
#' `step_decreases` and should be 0.
#'
#' @param data an [ssbmf_data()] object.
#' @param K number of programs (an upper bound; unused programs are switched
#'   off by the point-mass priors).
#' @param control settings from [ssbmf_control()].
#' @return An object of class `ssbmf_fit` with `L` (`mean`, `mean2`), `F`
#'   (by block, `mean`, `mean2`), `beta` (`mean`, `mean2`, `sd`), `mu`,
#'   `tau`, `center`, `priors`, `objective` (trace by sweep and final
#'   components), `step_decreases`, `halvings`, `rejections`, `converged`,
#'   `sweeps`, `roles`, `eta` (fitted predictor mean), `features`, `control`,
#'   `call`, `package_version`.
#' @family fit
#' @examples
#' sim <- ssbmf_simulate(n = 80, p = c(expression = 60), K = 2, seed = 1)
#' dd <- ssbmf_data(sim$X, sim$time, sim$status)
#' fit <- ssbmf_fit(dd, K = 3, control = ssbmf_control(max_sweeps = 30))
#' fit
#' @export
ssbmf_fit <- function(data, K, control = ssbmf_control()) {
  if (!inherits(data, "ssbmf_data")) stop("`data` must be created with ssbmf_data().")
  if (!inherits(control, "ssbmf_control")) stop("`control` must be created with ssbmf_control().")
  K <- as.integer(K)
  if (K < 1 || K >= data$n) stop("`K` must be between 1 and n - 1.")
  set.seed(control$seed)
  d <- .prep_data(data, control)
  st <- .init_state(d, K, control$init)
  p_eff <- NULL
  if (control$weight_rule == "effective_features") {
    p_eff <- vapply(seq_len(d$M), function(m) .effective_features(st$R[[m]] - outer(rep(1, d$n), colMeans(st$R[[m]])),
                                                                   d$n - K - 1), 0)
    d$w_block <- d$w_block * pmin(1, p_eff / d$p)
  }

  # Setup pass: give every column a fitted prior and KL so that the objective
  # is defined; afterwards every step is checked against it.
  st <- .update_mu(st, d); st <- .update_tau(st, d)
  for (k in seq_len(K)) st <- .update_L(st, d, k, safeguard = FALSE)
  for (m in seq_len(d$M)) for (k in seq_len(K)) st <- .update_F(st, d, m, k, safeguard = FALSE)
  st <- .objective(st, d)

  trace <- data.frame(sweep = 0L, total = st$obj$total)
  step_dec <- 0L; worst <- 0
  record <- function(old, new) {
    delta <- new$obj$total - old$obj$total
    if (delta < -1e-8 * abs(old$obj$total)) step_dec <<- step_dec + 1L
    worst <<- min(worst, delta / abs(old$obj$total))
    new
  }
  converged <- FALSE
  sweep_i <- 0L
  for (sweep_i in seq_len(control$max_sweeps)) {
    start <- st$obj$total
    st <- record(st, .update_mu(st, d))
    st <- record(st, .update_tau(st, d))
    for (k in seq_len(K)) st <- record(st, .update_L(st, d, k))
    for (m in seq_len(d$M)) for (k in seq_len(K)) st <- record(st, .update_F(st, d, m, k))
    st <- record(st, .update_beta(st, d))
    trace <- rbind(trace, data.frame(sweep = sweep_i, total = st$obj$total))
    rel <- abs(st$obj$total - start) / abs(st$obj$total)
    if (control$verbose > 0 && sweep_i %% control$verbose == 0) {
      message(sprintf("sweep %d  objective %.6f  rel change %.2e", sweep_i, st$obj$total, rel))
    }
    if (rel < control$tol) { converged <- TRUE; break }
  }
  fit <- .build_fit(st, d, data, control, trace, step_dec, worst, converged, sweep_i, match.call())
  fit$weights <- list(block = d$w_block, surv = d$w_surv, p_eff = p_eff)
  fit
}

#' Assemble the fitted object
#' @keywords internal
#' @noRd
.build_fit <- function(st, d, data, control, trace, step_dec, worst, converged, sweeps, call) {
  # recompute the projection scores exactly (the loop updates them incrementally)
  Zbar <- Reduce(`+`, lapply(seq_len(d$M), function(m) d$Yc[[m]] %*% st$EF[[m]]))
  K <- ncol(st$EL)
  prog <- paste0("P", seq_len(K))
  name_F <- function(A, m) { dimnames(A) <- list(data$features[[m]], prog); A }
  Fl <- lapply(seq_len(d$M), function(m) list(mean = name_F(st$EF[[m]], m), mean2 = name_F(st$EF2[[m]], m)))
  names(Fl) <- d$blocks
  st <- .objective(st, d)
  fit <- list(
    L = list(mean = `dimnames<-`(st$EL, list(data$samples, prog)),
             mean2 = `dimnames<-`(st$EL2, list(data$samples, prog))),
    F = Fl,
    beta = list(mean = stats::setNames(st$b, prog), mean2 = stats::setNames(st$b^2 + st$s2, prog),
                sd = stats::setNames(sqrt(st$s2), prog)),
    mu = stats::setNames(st$mu, d$blocks), tau = stats::setNames(st$tau, d$blocks),
    center = d$center,
    priors = list(L = st$gL, F = stats::setNames(st$gF, d$blocks), beta = st$gbeta),
    objective = list(trace = trace, final = c(total = st$obj$total,
                     stats::setNames(st$obj$rec, paste0("recon_", d$blocks)), surv = st$obj$S,
                     kl_L = sum(st$klL), kl_F = sum(unlist(st$klF)), kl_beta = st$klB)),
    step_decreases = step_dec, worst_relative_step = worst,
    halvings = st$counters$halvings, rejections = st$counters$rejections,
    converged = converged, sweeps = sweeps,
    eta = stats::setNames(drop(Zbar %*% st$b), data$samples),
    scores = `dimnames<-`(Zbar, list(data$samples, prog)),
    score_sd = stats::setNames(apply(Zbar, 2, stats::sd), prog),
    pve = .pve(st, d, prog),
    surv_share = .surv_share(st, d, prog),
    features = data$features, n = data$n, events = sum(data$status),
    strata_levels = if (is.null(data$strata)) NULL else levels(data$strata),
    control = control, call = call,
    package_version = as.character(utils::packageVersion("ssbmf")))
  class(fit) <- "ssbmf_fit"
  fit$roles <- ssbmf_classify(fit)
  fit
}

#' Proportion of each block's centered variance explained by each program
#' @keywords internal
#' @noRd
.pve <- function(st, d, prog) {
  out <- vapply(seq_len(d$M), function(m) {
    tot <- sum(d$Yc[[m]]^2)
    vapply(seq_along(prog), function(k) {
      lk <- st$EL[, k] - mean(st$EL[, k])
      sum(lk^2) * sum(st$EF[[m]][, k]^2) / tot
    }, 0)
  }, numeric(length(prog)))
  out <- matrix(out, length(prog), d$M, dimnames = list(prog, d$blocks))
  out
}

#' Survival share of the loading-update precision
#'
#' For each program and block, the median over features of
#' w_S A_jk / (a^R_jk + w_S A_jk): the fraction of the precision of
#' q(F_m,jk) that comes from the Cox term (derivation, Section 4.3).
#' @keywords internal
#' @noRd
.surv_share <- function(st, d, prog) {
  P <- .surv_state(st, d)$P
  out <- vapply(seq_len(d$M), function(m) {
    yP <- drop(crossprod(d$Yc2[[m]], P))
    vapply(seq_along(prog), function(k) {
      A <- d$w_surv * (st$b[k]^2 + st$s2[k]) * yP
      aR <- d$w_block[m] * st$tau[[m]] * sum(st$EL2[, k])
      stats::median(ifelse(A + aR > 0, A / (A + aR), 0))
    }, 0)
  }, numeric(length(prog)))
  matrix(out, length(prog), d$M, dimnames = list(prog, d$blocks))
}
