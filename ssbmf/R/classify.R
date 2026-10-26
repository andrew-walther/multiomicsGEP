# Program roles ----

#' Classify every program of one fitted model by role
#'
#' Roles come from posterior quantities of the single fit, with nothing
#' deleted or refitted: *survival_active* if the local false sign rate of
#' beta is at most `control$role_lfsr`; otherwise *reconstruction_only* if
#' the program explains at least `control$role_pve` of some block's centered
#' variance; otherwise *inactive*.
#'
#' @param fit an `ssbmf_fit`.
#' @return data.frame with `program`, `role`, `pve_<block>` columns,
#'   `beta_mean`, `beta_sd` (canonical scale) and `lfsr`.
#' @family methods
#' @export
ssbmf_classify <- function(fit) {
  cf <- coef(fit)
  pve <- fit$pve
  active <- cf$lfsr <= fit$control$role_lfsr
  recon <- apply(pve, 1, max) >= fit$control$role_pve
  role <- ifelse(active, "survival_active", ifelse(recon, "reconstruction_only", "inactive"))
  out <- data.frame(program = cf$program, role = role)
  for (m in colnames(pve)) out[[paste0("pve_", m)]] <- unname(pve[, m])
  out$beta_mean <- cf$mean; out$beta_sd <- cf$sd; out$lfsr <- cf$lfsr
  out
}
