# Targeted YFB mechanism diagnostics

Approved in chat on October 30 following the four-step diagnostic plan.
Work remains on the existing supervision branch, without merging or external scoring.

1. Compare adaptive/fixed feature precisions with matched initialization.
   Three shared-loading synthetic cases: seed 1, p=2000, raw variance shares
   .0025/.01/.05, the saved lowest-share arm's γ and ρ, K/cap/tolerance from globals.
   Initialize τ in both arms from the paired saved unsupervised fit. Retain the
   original saved fit separately: it starts τ=1 and is not the matched control.
2. Repeat those pairs with a whole-sweep acceptance guard. Twelve new fits total.
   The guard uses a fixed diagonal upper curvature D_i=max(10^-6, .5 times the
   number of event risk sets containing patient i), calculated within strata.
   For fixed projection scaling, E[log PL(η)] >= log PL(Eη)-.5 ΣD_i Var(η_i).
   Define the normalized moment criterion by evaluating this expression with
   s_k=1/||E F_k|| at every candidate, plus expected Gaussian log likelihood
   and the matching q/g KL terms. This is an explicit normalized diagnostic
   criterion, not a proof that a q-dependent normalizer defines a fixed-model ELBO.
   Working Cox proposals use D instead of the local diagonal curvature. Try the
   original F damping and six successive halvings. Restore all moments, priors,
   variances and KL terms before retries. Acceptance begins after the first
   completed feasible sweep; initial placeholder moments have no valid KL.
   Reject terminally with backtrack_stalled if no candidate improves the same
   criterion within relative tolerance 10^-8. Unchanged rejected states cannot
   establish convergence. No interpolation of posterior moments or stale KL.
3. Isolate coefficient uncertainty with fixed standardized n=400,K=7 designs:
   orthogonal versus equicorrelated (.9 population correlation). Gaussian working
   observations with unit noise; β=0 or β=(1,0,...,0). Fifty independent response
   seeds, Cox-weight analogues 1 and the first saved γ, fixed point-normal prior
   π=1/K, slab SD=1. Compare independent-coordinate posterior updates (cap1000,
   tolerance10^-8) to exact enumeration of all inclusion patterns and Gaussian
   conditional covariances. This is a working-likelihood diagnostic, not full
   Cox coverage or an additional factor-recovery gate. Retain nonconvergence.
4. Save every fit and report reconstruction, τ, loading norms, survival activity,
   recovery, stopping, rejected trials and synthetic C. Screening success requires
   a finite converged fit, no guard stall, reconstruction ratios <=1 (better
   than zero prediction), recovered planted factor and false-active <=1.
   Guarded accepted objectives must not decrease beyond the registered tolerance.
   Compare objectives only within the same criterion. These fixed cases cannot
   demonstrate Gate S; its original complete-grid requirements remain unchanged.

All changes are opt-in. Verify default numerical parity against a pre-change
synthetic fixture, frozen τ equality, complete rollback and explicit stalled
stopping, the Cox curvature inequality, and exact joint-posterior formulas.
Run the full suite after implementation. Use PSOCK workers with single-threaded
BLAS; preserve all earlier outputs and the frozen benchmark. Do not choose an
arm from synthetic validation C. Update documentation with both failures and gains.
