# Algorithm settings ----

#' Settings for ssbmf_fit()
#'
#' @param max_sweeps maximum number of full sweeps (default 1000).
#' @param tol relative change in the objective over one sweep below which the
#'   fit is declared converged (default 1e-6).
#' @param prior_L prior family for activities L (default "point_exponential").
#' @param prior_F prior family for loadings F, either one family for every
#'   block or a named vector by block (default "point_laplace"; also
#'   "point_exponential", "point_normal", "normal").
#' @param w_surv survival weight w_S multiplying the Cox term (default 1, the
#'   unweighted model; 0 gives unsupervised factorization with beta not fitted).
#' @param w_block optional named vector of block weights w_m (default 1 each).
#' @param init "svd" (default) or "survival" (last program started at the
#'   survival score direction).
#' @param seed random seed used for initialization.
#' @param max_halvings safeguard: number of step halvings before an update is
#'   rejected (default 8).
#' @param beta_inner maximum inner iterations of the beta update (default 25).
#' @param tau_cap upper bound on each precision, as a multiple of
#'   1 / Var(y_j) (default 1e6).
#' @param role_lfsr local false sign rate at or below which a program is
#'   survival-active (default 0.05).
#' @param role_pve proportion of a block's centered variance at or above which
#'   a program is reconstruction-active (default 1e-3).
#' @param verbose print progress every `verbose` sweeps (0 = silent).
#' @return A list of class `ssbmf_control`.
#' @family fit
#' @examples
#' ssbmf_control(max_sweeps = 200)
#' @export
ssbmf_control <- function(max_sweeps = 1000, tol = 1e-6, prior_L = "point_exponential",
                          prior_F = "point_laplace", w_surv = 1, w_block = NULL,
                          init = c("svd", "survival"), seed = 1, max_halvings = 8,
                          beta_inner = 25, tau_cap = 1e6, role_lfsr = 0.05, role_pve = 1e-3,
                          verbose = 0) {
  init <- match.arg(init)
  fams <- c("point_exponential", "point_laplace", "point_normal", "normal")
  if (!prior_L %in% c("point_exponential")) stop("`prior_L` must be 'point_exponential' (L is nonnegative).")
  if (!all(prior_F %in% fams)) stop("Unknown `prior_F` family.")
  if (w_surv < 0) stop("`w_surv` must be >= 0.")
  structure(as.list(environment()), class = "ssbmf_control")
}
