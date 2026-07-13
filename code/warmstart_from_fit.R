# ============================================================
# Script:  warmstart_from_fit.R
# Purpose: PVE-ranked column extraction for warm-starting a smaller-K refit
#          from an already-converged, larger-K fit (fit_cox_on_yf() or
#          fit_supervised_mf_modular() output).
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Dependencies: none (operates on the returned fit list only)
# ============================================================

#' Extract Top-K Warm-Start Columns From a Converged Fit, Ranked by PVE
#'
#' Selects the `K_target` columns of a converged fit's EL/EF with the highest
#' final-iteration proportion of variance explained (PVE), for use as
#' `init_method = "custom"` warm-start (`EL_init`, `EF_init`) in a smaller-K
#' refit. Ranking by PVE rather than |beta| is deliberate: K serves both
#' reconstruction (EL/EF explain Y) and prediction (beta acts on ZF = Y·EF),
#' so a factor with a small beta but high reconstruction PVE can still be a
#' useful column to keep, and one with a large beta but negligible PVE is a
#' fragile warm-start signal.
#'
#' PVE for factor k is `sum(EL[,k]^2) * sum(EF[,k]^2) / ||Y||_F^2`, the same
#' quantity CAVI itself tracks each iteration in `fit$history$factor_pve`
#' (see code/fit_cox_on_yf.R). This function reads the LAST row of that
#' matrix (the converged fit's final ranking), not the first.
#'
#' @param fit      A fit_cox_on_yf() or fit_supervised_mf_modular() result
#'                 list; must contain `$EL` (n x K_source), `$EF`
#'                 (p x K_source), and `$history$factor_pve`
#'                 (iterations x K_source matrix).
#' @param K_target integer: number of columns to extract (1 <= K_target <=
#'                 K_source, the number of columns in `fit$EL`).
#'
#' @return list(EL_init = n x K_target matrix, EF_init = p x K_target matrix),
#'   columns ordered by descending final-iteration PVE (most-explanatory
#'   first). Pass directly as `fit_cox_on_yf(..., EL_init = out$EL_init,
#'   EF_init = out$EF_init)`.
extract_top_k_by_pve <- function(fit, K_target) {
  K_source <- ncol(fit$EL)
  if (!is.numeric(K_target) || length(K_target) != 1 ||
      K_target != round(K_target) || K_target < 1 || K_target > K_source) {
    stop(sprintf(
      "K_target must be a positive integer <= K_source (%d); got %s.",
      K_source, deparse(K_target)
    ))
  }

  pve_hist  <- fit$history$factor_pve
  pve_final <- pve_hist[nrow(pve_hist), ]
  top_idx   <- order(pve_final, decreasing = TRUE)[seq_len(K_target)]

  list(
    EL_init = fit$EL[, top_idx, drop = FALSE],
    EF_init = fit$EF[, top_idx, drop = FALSE]
  )
}
