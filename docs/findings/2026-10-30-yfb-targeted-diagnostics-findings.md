# Targeted YFB diagnostics: noise feedback, guarded updates and coefficient uncertainty

Written October 30, 2026. The user approved this small study in chat. Its
registration precedes execution in the attempts log and targeted diagnostic
plan. The fitting/reference source checkpoint is **0bf3740**; subsequent tests
strengthen independent checks without changing the fitted source.

## Result

**Partial improvement, not a solved supervised model.** All twelve registered
fits were saved without numerical exceptions. Holding feature precisions fixed
prevented the reconstruction explosion in these three cases. Guarded updates
also kept reconstruction bounded and, in the 1% case with adaptive noise,
recovered the planted program and converged. That fit still had two unmatched
survival-active factors. No case passes all screening criteria. Of twelve fits,
one converges, eight hit the sweep cap and three explicitly stall in backtracking.
No external cohort was scored; no full-grid Gate S claim is made or changed.

## Matched noise intervention

Each pair starts from the same training-only unsupervised precision vector.
Both use the same synthetic data, hybrid start, saved γ/ρ and existing priors.
The original saved S-7 fit starts at τ=1 and is not used as the matched control.
The ratio is expected residual energy divided by total expression energy;
it is not centered variance explained.

| variance | adaptive_expected_error | fixed_expected_error | error_reduction_factor | fixed_tau_difference |
|---|---|---|---|---|
| 0.0025 | 0.2057 | 0.07925 | 2.596 | 0 |
| 0.01 | 884.3 | 0.07861 | 1.125e+04 | 0 |
| 0.05 | 4.483e+14 | 0.07937 | 5.648e+15 | 0 |

The fixed precision vector is unchanged exactly. The intervention demonstrates
that allowing feature-noise adaptation enables the observed runaway feedback
in these cases. It does not prove that fixed noise is the best model or that
a weak Gamma prior would suffice. A finite Gamma prior alone still allows
precision to approach zero as residual energy diverges. No Gamma prior or
new noise default is adopted.

## Full twelve-fit inventory

Recovery here is the raw loading/sign/posterior-activity match; a recovered
loading with a capped or stalled fit is not screening success. All guard trials,
including rejected ones, are retained in each fit. No accepted guarded step
decreases its criterion beyond the registered tolerance. The numerical posterior
checks remain finite, and the full regression suite passes **749/749** with a
clean exit; the saved optimizer-fixture warning is intentional.

| variance | arm | expected_reconstruction_ratio | recovered | false_active | converged | flags | sweeps | rejected_trials | success |
|---|---|---|---|---|---|---|---|---|---|
| 0.0025 | adaptive | 0.2057 | FALSE | 6 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.0025 | fixed | 0.07925 | FALSE | 6 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.0025 | guard_adaptive | 0.07734 | FALSE | 4 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.0025 | guard_fixed | 0.07735 | FALSE | 4 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.01 | adaptive | 884.3 | FALSE | 4 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.01 | fixed | 0.07861 | FALSE | 4 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.01 | guard_adaptive | 0.07748 | TRUE | 2 | TRUE |  | 110 | 13 | FALSE |
| 0.01 | guard_fixed | 0.07753 | TRUE | 2 | FALSE | backtrack_stalled | 98 | 17 | FALSE |
| 0.05 | adaptive | 4.483e+14 | TRUE | 4 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.05 | fixed | 0.07937 | FALSE | 4 | FALSE | max_iter | 300 | 0 | FALSE |
| 0.05 | guard_adaptive | 0.07707 | TRUE | 3 | FALSE | backtrack_stalled | 6 | 7 | FALSE |
| 0.05 | guard_fixed | 0.07707 | TRUE | 3 | FALSE | backtrack_stalled | 6 | 7 | FALSE |

The guard is an update package: it changes working curvature to the global
risk-set bound as well as checking acceptance. It does not isolate acceptance
from that curvature change. Accepted-state trajectories are plotted separately
from rejected candidates; short lines can represent stalled fits.

## Coefficient experiment: information weighting is a major cause

This separate experiment uses fixed n=400,K=7 standardized Gaussian designs,
unit observational noise, a fixed point-normal prior π=1/K, slab SD=1,
orthogonal/equicorrelated designs and 50 paired response seeds. Null and
planted-coefficient responses are both included. The exact calculation sums
all inclusion patterns; the mean-field calculation updates independent
coordinates at the same target. No empirical prior differences confound them.

| design | signal | weight | method | n | mean_false | se_false | recovered | mean_sd | nonconverged |
|---|---|---|---|---|---|---|---|---|---|
| correlated | null | 1 | joint | 50 | 0 | 0 | NA | 0.01789 | 0 |
| orthogonal | null | 1 | joint | 50 | 0 | 0 | NA | 0.01047 | 0 |
| correlated | planted | 1 | joint | 50 | 0 | 0 | 50 | 0.04531 | 0 |
| orthogonal | planted | 1 | joint | 50 | 0 | 0 | 50 | 0.01614 | 0 |
| correlated | null | 6.223e+04 | joint | 50 | 6.86 | 0.04957 | NA | 0.0006058 | 0 |
| orthogonal | null | 6.223e+04 | joint | 50 | 6.92 | 0.03876 | NA | 0.0002003 | 0 |
| correlated | planted | 6.223e+04 | joint | 50 | 5.88 | 0.04642 | 50 | 0.000604 | 0 |
| orthogonal | planted | 6.223e+04 | joint | 50 | 5.92 | 0.03876 | 50 | 0.0002003 | 0 |
| correlated | null | 1 | mean_field | 50 | 0 | 0 | NA | 0.007649 | 0 |
| orthogonal | null | 1 | mean_field | 50 | 0 | 0 | NA | 0.01047 | 0 |
| correlated | planted | 1 | mean_field | 50 | 0 | 0 | 50 | 0.01162 | 0 |
| orthogonal | planted | 1 | mean_field | 50 | 0 | 0 | 50 | 0.01614 | 0 |
| correlated | null | 6.223e+04 | mean_field | 50 | 6.86 | 0.04957 | NA | 0.0001965 | 0 |
| orthogonal | null | 6.223e+04 | mean_field | 50 | 6.92 | 0.03876 | NA | 0.0002003 | 0 |
| correlated | planted | 6.223e+04 | mean_field | 50 | 5.88 | 0.04642 | 50 | 0.000197 | 0 |
| orthogonal | planted | 6.223e+04 | mean_field | 50 | 5.92 | 0.03876 | 50 | 0.0002003 | 0 |

At weight one, both methods recover the planted coefficient in all fifty
responses per design and show no false-active coefficients in these draws.
At the saved large weight, even exact inference flags about 6.9 of seven null
coefficients on average. Correlation broadens exact marginal uncertainty
relative to mean-field uncertainty, but does not repair the high-temperature
false activity. An always-zero estimator cannot explain these results because
the planted positive controls recover correctly.

The weighted Gaussian likelihood assumes effective noise variance 1/γ,
although the responses were generated with variance one. Thus the weighted
posterior is mathematically correct for its target yet overconfident for the
generating experiment. This isolates the mechanism of increased information
weight; it is not a full-Cox coverage study, a proof of universal zero error
at γ=1, or a recommendation to change reported SDs after fitting.

## Teach-back and implementation limits

**Purpose and interfaces.** `tau_init` supplies a positive length-p vector;
`tau_fixed` holds it fixed. `guard` enables diagnostic whole-sweep acceptance.
The existing fitter defaults retain their old numerical results against a
pre-change synthetic fixture. New returned fields include the trial history,
rejection count and KL terms. No destructive pruning or external model choice
is introduced. The exact Gaussian reference returns K means, a K×K covariance,
variances, inclusion probabilities and pattern weights; it is not integrated
into the YFB fitting updates.

**Update logic.** Snapshot all moments, separate variances, priors and KL terms.
Construct a proposal using the fixed risk-set curvature and current normalized
projection. Update L/F/β and their priors coherently; canonicalize, then hold
or update τ. Recompute normalization and the same criterion. Retry with six
successive halvings of the F survival damping, restoring the original state
before each trial. The first completed feasible sweep establishes the initial
comparison state because initialization has placeholder moments without valid KL.
If all trials fail, restore the last accepted state and return backtrack_stalled.
An unchanged rejected state cannot qualify as convergence.

**Criterion and mathematical choice.** For each event, the log-sum-exp Hessian
is bounded by I_R/2. Summing over event risk sets gives D_i equal to half their
membership count, floored at 10^-6. At fixed projection normalization,

E[log PL(η)] ≥ log PL(Eη) − (1/2) Σ_i D_i Var(η_i).

The guard adds this term, weighted by γ, to expected Gaussian log likelihood
minus matching q/g KLs. It recomputes s_k=1/||E F_k|| at each evaluated state.
This is one explicitly defined normalized moment criterion. The q-dependent
normalizer still needs a fixed-model interpretation; no fixed-generative-model
ELBO guarantee is asserted. Objectives under native and guarded modes are not
compared numerically. Finite, accepted and converged remain distinct concepts.

**Why some guarded fits stall.** F-only damping retries cannot always produce
an improving whole-state proposal: coefficient/prior updates and normalization
also change. Rejection is a useful safeguard, not a replacement for deriving
a better normalized update direction. No gradient stationarity test is claimed.
The remaining task is a consistent normalized F/β update or trust-region/MM
proposal, with coherent prior learning and uncertainty, rather than accepting
a lower objective or interpolating means while retaining stale KL values.

**Independent verification.** Tests cover default parity, unchanged fixed τ,
the exact Cox Hessian bound, full rejected-state rollback (including priors
and KL), an independently calculated objective and its rescaling invariance,
orthogonal spike/slab marginal formulas and correlated Gaussian covariance.
The exact reference uses covariance about the mixture mean to avoid subtraction
of large second moments. Mean-field capped runs are retained explicitly.

**R and computational choices.** Base R provides Cholesky factorizations and
finite-pattern enumeration; ashr constructs the existing prior class. No new
package dependency is added. The small exact reference is limited to K≤10
because its cost is exponential. Three PSOCK workers and one BLAS thread per
worker avoid macOS forking and oversubscription. All inputs and RDS checkpoints
are synthetic. The plotted comparisons use common axes and saved histories.

## Next test recommended

Keep YFB, nonnegative inputs and the sparse priors. Do not attempt to repair
false activity solely with a joint β covariance or a different slab family:
exact inference at the same large weight already exhibits the problem.
The next methodological choice is a single generalized-Bayes objective that
balances genomic information without treating the survival sample as tens of
thousands of repetitions, together with a normalized update that can make
accepted progress. Downweighting genomic information while keeping Cox weight
near one is a candidate to derive and register, not a change adopted here.
Frozen τ remains a useful controlled baseline for that next experiment.

## Evidence and reproduction

Outputs are under `results/benchmark_sim/outputs/yfb_targeted_2026-10-30/`:
summary, matched-noise comparison, all twelve checkpoints, Gaussian coefficient
metrics, accepted-step checks, hashes, plots and compressed logs. The README
explains their semantics. Run the registered scripts only after reading the
plan; existing checkpoints are reused and coefficient outputs refuse overwrite.
The original S-7 results, immutable benchmark and independent-review verdict
are untouched. The October 30 notebook includes this follow-up from saved CSVs.


## Review entry points

- Registration: `docs/plans/2026-10-30-yfb-targeted-diagnostics-plan.md`.
- Fitter/control changes: `code/fit_yfb_supervised.R`, `code/yfb_diagnostic_controls.R`.
- Exact/reference inference: `code/yfb_joint_beta.R`.
- Registered runners: `results/benchmark_sim/run_yfb_targeted_diagnostics.R` and `run_yfb_beta_uncertainty_screen.R` in the same directory.
- Behavioral checks: `tests/test_yfb_diagnostic_controls.R`, `tests/test_yfb_joint_beta.R`.
- All result semantics: the targeted output directory's `README.md`.
- The complete progress book was rendered successfully; figure resources match
  saved evidence. Default-parity and independent mathematical checks pass; the
  intentionally missing-feature red checks and both full-suite logs are retained.
