# ============================================================
# Script: multimodal_moment_audit.R
# Purpose: Summarize the existing F-update moment guard without modifying inputs.
# Created: 2026-10-15
# Dependencies: base R
# ============================================================

#' Describe which input violates the existing F-update moment guard
#' Negative h and EZ are valid; w must be nonnegative and VZ must be at least
#' -1e-10, exactly as in multimodal_yfb_update_F_mk(). Nonfinite values fail.
#' @param w Cox diagonal working curvature vector.
#' @param h_minus_k Cox working numerator excluding the current program.
#' @param EZ_k,VZ_k Projection mean and variance vectors for the current program.
#' @return Four aggregate rows with counts, finite ranges and guard status.
#' @examples
#' multimodal_moment_guard_summary(c(0,1),c(-1,2),c(-2,3),c(0,.1))
#' @family multimodal_moment_audit
#' @seealso multimodal_yfb_update_F_mk
multimodal_moment_guard_summary <- function(w,h_minus_k,EZ_k,VZ_k) {
  values<-list(w=w,h_minus_k=h_minus_k,EZ_k=EZ_k,VZ_k=VZ_k)
  rows<-lapply(names(values),function(nm) {
    x<-values[[nm]];finite<-x[is.finite(x)]
    lower<-switch(nm,w=0,VZ_k=-1e-10,-Inf)
    data.frame(field=nm,n=length(x),nonfinite=sum(!is.finite(x)),
      negative=sum(x<0,na.rm=TRUE),
      minimum=if(length(finite))min(finite) else NA_real_,
      maximum=if(length(finite))max(finite) else NA_real_,
      invalid=any(!is.finite(x)) || any(x<lower,na.rm=TRUE))
  })
  do.call(rbind,rows)
}
