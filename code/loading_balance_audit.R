# ============================================================
# Script: loading_balance_audit.R
# Purpose: Decompose reconstruction and survival information for a frozen loading coordinate.
# Created: 2026-10-15
# Dependencies: base R
# ============================================================

#' Decompose one program's loading likelihood at a frozen posterior state
#' A_Gj = τ_j sum_i E[L_ik²]; A_Sj = E[β_k²] sum_i w_i y_Sij².
#' B_Gj = τ_j sum_i E[L_ik] R_minus_k,ij;
#' B_Sj = E[β_k] sum_i y_Sij h_minus_k,i
#'   - E[β_k²] sum_i w_i y_Sij (E[Z_ik]-y_Sij E[F_jk]).
#' Each coordinate excludes its own current contribution while all other
#' moments stay frozen. These are not an executed sequential sweep or prior
#' contributions. Survival covariates may be profiled; reconstruction is not.
#' @param Y_gen,Y_surv Reconstruction and survival matrices, patients by features.
#' @param Tau Nonnegative feature precision vector.
#' @param EL,EL2 Program activity mean and second moment vectors.
#' @param R Reconstruction residual excluding this program, same shape as Y_gen.
#' @param EZ Current projection mean vector for this program, all modalities.
#' @param EF Current loading means for this modality/program.
#' @param beta,beta2 Coefficient mean and second moment.
#' @param w Frozen nonnegative Cox diagonal working weights.
#' @param h_minus Working numerator excluding the other programs' mean predictor.
#' @return Per-feature A/B contributions and survival shares. Shares with no
#'   information remain NA; signal share uses absolute terms before cancellation.
#' @examples
#' # loading_balance_terms(Y,Y,tau,L,L2,R,Z,F,b,b2,w,h_minus)
#' @family loading_balance_audit
#' @seealso multimodal_yfb_update_F_mk
loading_balance_terms <- function(Y_gen,Y_surv,Tau,EL,EL2,R,EZ,EF,
                                  beta,beta2,w,h_minus) {
  stopifnot(is.matrix(Y_gen),identical(dim(Y_gen),dim(Y_surv)),
    identical(dim(Y_gen),dim(R)),length(Tau)==ncol(Y_gen),
    length(EF)==ncol(Y_gen),length(EL)==nrow(Y_gen),length(EL2)==nrow(Y_gen),
    length(EZ)==nrow(Y_gen),length(w)==nrow(Y_gen),length(h_minus)==nrow(Y_gen),
    length(beta)==1L,length(beta2)==1L,
    all(is.finite(c(Y_gen,Y_surv,Tau,EL,EL2,R,EZ,EF,beta,beta2,w,h_minus))),
    all(Tau>=0),all(EL2>=0),beta2>=0,all(w>=0))
  wy2<-as.vector(crossprod(w,Y_surv^2))
  A_gen<-Tau*sum(EL2);A_surv<-beta2*wy2
  B_gen<-Tau*as.vector(crossprod(R,EL))
  B_surv<-beta*as.vector(crossprod(Y_surv,h_minus))-
    beta2*(as.vector(crossprod(Y_surv,w*EZ))-EF*wy2)
  precision<-A_gen+A_surv;signal<-abs(B_gen)+abs(B_surv)
  stopifnot(all(is.finite(c(A_gen,A_surv,B_gen,B_surv,precision,signal))))
  data.frame(A_gen=A_gen,A_surv=A_surv,B_gen=B_gen,B_surv=B_surv,
    precision_share=ifelse(precision>0,A_surv/precision,NA_real_),
    signal_share=ifelse(signal>0,abs(B_surv)/signal,NA_real_))
}
