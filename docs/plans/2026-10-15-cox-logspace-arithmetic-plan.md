# Stable arithmetic for Breslow working quantities

The25k CpG failure is a finite-precision error in the global-shift Cox
calculation, not a missing-value or projection-variance error. A stable
calculation preserves the same partial likelihood and diagonal working
approximation; it does not repair an unstable optimization trajectory.

For distinct event times t_g, define risk sums R_g=sum_{j:t_j>=t_g}exp(eta_j)
and deaths d_g. The existing working weight is
w_i=exp(eta_i) sum_{g:t_g<=t_i} d_g/R_g, with u_i=delta_i-w_i and h_i=u_i+w_i eta_i.
The same quantities can be evaluated as follows:

1. Sort patients by time. Compute each suffix log-risk-sum with log-add-exp.
   Equal-time patients use the suffix beginning at the group's first row.
2. Form log(d_g)-log(R_g), taking log(0)=-Inf for groups without deaths.
3. Compute cumulative log-add-exp of these increments, log H_i.
4. Set w_i=exp(eta_i+log H_i), with zero before the first event.
5. Evaluate logPL=sum over event patients of eta_i-log R(t_i). Selecting
   event rows avoids multiplying zero censoring indicators by undefined terms.

Every risk set retains its own finite log normalization even when globally
shifted exponentials round to zero. The mathematical total sum_i w_i is the
number of deaths. The algorithm uses sorting plus linear passes over patients,
not a dense feature or patient-risk-set matrix. Molecular fitting cost remains.
No clipping, rescaled survival likelihood, additional prior, tuning parameter,
new convergence rule or changed preprocessing is introduced.

The implementation is opt-in through control$cox_arithmetic="logspace".
Default/original fit objects must remain bit-for-bit equal to the protected
base-branch default. Original and log-space arithmetic agree to floating-point
tolerance on moderate predictors; output is expected to differ where original
arithmetic underflows. The same option is propagated to working evaluations,
monitoring and removal-score calculations so a fit does not mix arithmetic.

Task5-M3 first checks the captured failing point, then repeats just the25k
full-training case with all original scientific settings and caps. Passing
numerical checks does not mean convergence, credible recovery, validated
feature selection or benchmark adoption. Remaining failures are preserved.
