# Proposed opt-in diagnostic: profile the working Cox intercept

Status: Andrew approved testing on2026-10-15, conditional on exact reproduction
of the unchanged failed-fold trace. Mathematical/software proposal only. No project source changed, no
new scientific fit, no scoring or adoption criterion changed.

The saved failed framework fold has beta=0 although loading program4 remains.
The diagonal surrogate uses A=Z'WZ+tr(W CovZ). Constant projection components
increase A despite the exact Cox partial likelihood's invariance to eta+c.
At this saved state, constant components account for 34–77% of each A (58%
for program4). Its pseudo-z is .9104; the ordinary exact-Cox null score z
on that same fixed loading is1.4164. Neither is claimed to be credible or
evidence of validation performance.

For frozen u,w,eta0 with 1'u=0, profile the arbitrary scalar c from
Q(eta+c)=constant+u'(eta+c-eta0)
 -1/2(eta+c-eta0)'W(eta+c-eta0).
The maximizing c removes the weighted mean of eta-eta0, giving
Hc=W-ww'/D, D=sum(w), hc=u+Hc eta0.
Thus Qp=constant+hc'eta-eta'Hc eta/2.
Hc is PSD and Hc 1=0. With w=sum_r d_r P_r and D=sum_r d_r,
Hc-Hexact=sum_r d_r(P_r-w/D)(P_r-w/D)' is PSD.
This profiles a nuisance shift in the existing local approximation; it does
not change observed matrices or the underlying Cox partial likelihood.

Use matrix-vector products Hc v=w*v-w*sum(w*v)/D, and
qj=yj'Hc yj=sum(w*yj^2)-sum(w*yj)^2/D.
No n-by-n covariance or feature-by-feature Gram matrix is needed.

Keep mean-field moments and existing priors. For factor k, define
r_k=hc-Hc eta_minus_k and CovZk=sum_j vF_j yj yj'.
The exact scalar updates in this profiled fixed quadratic are:
 A_beta=EZk'Hc EZk+sum_j vF_j qj;
 B_beta=EZk' r_k.
For loading j, including both modalities in EZk_minus_j:
 A_Fj=tau_j sum_i E[L_ik^2]+E[beta_k^2] qj;
 B_Fj=tau_j Lk'Rj+E[beta_k] yj' r_k
       -E[beta_k^2] yj'Hc EZk_minus_j.
The expected working quadratic uses
 E[eta'Hc eta]=mu_eta'Hc mu_eta+
 sum_k {E[beta_k^2]sum_j vF_j qj+
         Var(beta_k) EZk'Hc EZk}.
L/tau updates and raw prediction remain unchanged. Zero-event data must
retain the current zero-information behavior (D=0, u=0, Hc=0).
Default diagonal fitting remains unchanged; this would be an opt-in flag.

This is outside the handoff's listed R1–R3 prior remedies; Andrew explicitly
approved its prospective test. Preregistration, a commit, full tests and
unchanged-fold trace identity are still required before scientific execution. It does not resolve the separate
shared-prior/canonicalization question. A first test would use the historical
per-factor framework so that the shared-prior question is not combined with
this change. Preserve all original simulation, CV, no-collapse, cap and
benchmark rules. In particular, the existing collapsed reference fold stays
C=NA and cannot be silently replaced with .5 or a new reference.

Approved next step after the trace prerequisite: opt-in implementation with default-path
identity and independently verified scalar updates, then preregister and
commit one fixed simulation/TCGA diagnostic evaluation before fitting.
No threshold, test bed, default or external-selection rule may change.

Private software checks:8/8 model/math checks and6/6 new evaluation-protocol
checks pass (plus9/9 original R1 protocol checks). These are synthetic
software checks, not scientific candidate fits or the full project suite.
The expected full-suite total after applying/registering them is625; report
the actual count at that future checkpoint.
