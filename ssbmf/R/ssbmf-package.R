#' ssbmf: survival-supervised Bayesian matrix factorization
#'
#' For data blocks \eqn{Y_m} (samples by features) the model is
#' \deqn{Y_m = 1\mu_m^T + L F_m^T + E_m,\quad (E_m)_{ij} \sim N(0, 1/\tau_{mj}),}
#' with patient activities \eqn{L} shared across blocks and block-specific
#' loadings \eqn{F_m}. Survival follows a Cox model with linear predictor
#' \eqn{\eta = \sum_m \tilde Y_m F_m \beta}, where \eqn{\tilde Y_m} is the
#' block centered at its training means. Priors on \eqn{L}, \eqn{F_m} and
#' \eqn{\beta} are fitted by empirical Bayes.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats coef predict sd
#' @importFrom survival Surv strata
## usethis namespace: end
NULL
