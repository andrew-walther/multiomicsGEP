# ============================================================
# Script: test_share_targeted_ws.R
# Purpose: Check directional precision weights against known quantities.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: code/share_targeted_ws.R
# ============================================================
local({
original_cox <- get0("calc_cox_taylor_yf", envir = .GlobalEnv, inherits = FALSE)
source("code/share_targeted_ws.R", local = TRUE)
# A fixed Cox working-weight stub isolates the algebra from risk-set calculation.
#' @param eta,time,status Unused training vectors.
#' @return A list of fixed working weights.
#' @examples
#' # calc_cox_taylor_yf(rep(0, 3), 1:3, c(1, 1, 0))
calc_cox_taylor_yf <- function(eta, time, status) list(w = c(1, 2, 3))
Y <- matrix(c(1, 2, 3, 2, 0, 1), 3, 2)
fit <- list(EL = matrix(c(1, 1, 1), 3, 1), EL2 = matrix(c(2, 2, 2), 3, 1),
            EF = matrix(c(1, 0), 2, 1), EBeta = 2, EBeta2 = 4.25,
            Tau = c(2, 4), canonicalize = TRUE, w_S = 1,
            history = list(converged = TRUE))
r <- share_targeted_ws(fit, Y, 1:3, c(1, 1, 0), 0.5)
expected <- 4.25 * (1 + 2 * 4 + 3 * 9) / (2 * 6)
stopifnot(isTRUE(all.equal(r$ratio, expected)), isTRUE(all.equal(r$w_S, 1 / expected)),
          r$selected_program == 1, !r$no_signal_guard, !r$capped,
          isTRUE(all.equal(r$diagnostics$directional_share, expected / (1 + expected))))
fit2 <- fit; fit2$w_S <- 2
r2 <- share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5, pilot_loading = r$pilot_loading)
stopifnot(isTRUE(all.equal(r2$ratio, 2 * expected)), isTRUE(all.equal(r2$w_S, r$w_S)))
fit2 <- fit; fit2$EBeta <- 0
guard <- share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5)
stopifnot(guard$no_signal_guard, guard$w_S == 1, is.na(guard$selected_program))
fit2 <- fit; fit2$EL2 <- NULL
stopifnot(share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5)$loading_second_moment_bound)
fit2 <- fit; fit2$Tau <- 1e10
cap <- share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5)
stopifnot(cap$capped, cap$w_S == 1e5)
# Factor rescaling preserves η, directional ratio and the selected weight.
fit2 <- fit; fit2$EF <- fit$EF * 2; fit2$EL <- fit$EL / 2
fit2$EL2 <- fit$EL2 / 4; fit2$EBeta <- fit$EBeta / 2; fit2$EBeta2 <- fit$EBeta2 / 4
r2 <- share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5)
stopifnot(isTRUE(all.equal(r2$ratio, r$ratio)), isTRUE(all.equal(r2$w_S, r$w_S)))
for (bad in list(0, 1, NA_real_)) {
  stopifnot(inherits(try(share_targeted_ws(fit, Y, 1:3, c(1, 1, 0), bad), silent = TRUE), "try-error"))
}
fit2 <- fit; fit2$history$converged <- FALSE
stopifnot(inherits(try(share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5), silent = TRUE), "try-error"))
# Exact per-gene shares use each feature's precision, then take their median.
gene_ratios <- 4.25 * c(36 / 12, 7 / 24)
stopifnot(isTRUE(all.equal(r$diagnostics$per_gene_median_share,
                           median(gene_ratios / (1 + gene_ratios)))))
fit2 <- fit; fit2$EL2 <- NULL
bound_result <- share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5)
stopifnot(bound_result$ratio > r$ratio,
          bound_result$diagnostics$per_gene_median_share > r$diagnostics$per_gene_median_share)
# Four programs: the second has the largest standardized coefficient, the
# third has zero beta SD, and the fourth has a constant score.
Y4 <- cbind(Y, rep(1, 3))
f4 <- fit
f4$EF <- cbind(c(1, 0, 0), c(0, 1, 0), c(1, 0, 0), c(0, 0, 1))
f4$EL <- matrix(1, 3, 4); f4$EL2 <- matrix(2, 3, 4)
f4$EBeta <- c(1.96, 5, 10, 20)
f4$EBeta2 <- f4$EBeta^2 + c(1, 1, 0, 1)
f4$Tau <- rep(2, 3)
r4 <- share_targeted_ws(f4, Y4, 1:3, c(1, 1, 0), 0.5)
stopifnot(identical(r4$diagnostics$active, c(TRUE, TRUE, FALSE, FALSE)),
          r4$selected_program == 2)
# Reordering and sign reversal preserve the pilot match even after activity is lost.
tracked <- f4
tracked$EF <- cbind(-f4$EF[, 2], f4$EF[, 1], f4$EF[, 3], f4$EF[, 4])
tracked$EBeta <- c(0, 9, 0, 0)
tracked$EBeta2 <- tracked$EBeta^2 + 1
tracked$w_S <- 2
rt <- share_targeted_ws(tracked, Y4, 1:3, c(1, 1, 0), .5,
                        pilot_loading = r4$pilot_loading)
stopifnot(rt$selected_program == 1, rt$lost_activity, !rt$no_signal_guard,
          isTRUE(all.equal(rt$matching_correlation, 1)),
          identical(rt$pilot_loading, r4$pilot_loading))
tracked$EBeta <- rep(0, 4); tracked$EBeta2 <- rep(1, 4)
stopifnot(!share_targeted_ws(tracked, Y4, 1:3, c(1, 1, 0), .5,
                             pilot_loading = r4$pilot_loading)$no_signal_guard)
stopifnot(inherits(try(share_targeted_ws(fit, Y, 1:3, c(1, 1, 0), .5,
                                       pilot_loading = c(1, 1)), silent = TRUE), "try-error"))
for (field in c("projection_center", "cohort_id", "beta_cohort_id", "modality_weight",
                "feature_group", "strata_id", "profile_cox_intercept")) {
  fit2 <- fit; fit2[[field]] <- TRUE
  stopifnot(inherits(try(share_targeted_ws(fit2, Y, 1:3, c(1, 1, 0), 0.5), silent = TRUE), "try-error"))
}
if (!is.null(original_cox)) {
  calc_cox_taylor_yf <- original_cox
  actual <- share_targeted_ws(fit, Y, 1:3, c(1, 1, 0), 0.5)
  stopifnot(is.finite(actual$ratio), actual$ratio > 0)
}
stopifnot(identical(original_cox, get0("calc_cox_taylor_yf", envir = .GlobalEnv, inherits = FALSE)))
cat("share_targeted_ws behavioral checks passed\n")
})
