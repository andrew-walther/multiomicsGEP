# Fixed-posterior loading-prior consistency repair

Updated 2026-10-16. Andrew authorized proceeding with the recommended
prior-consistency implementation after the corrected observer. This is an
opt-in implementation and prospective evaluation specification. The earlier implementation block ran no scientific fit. Andrew subsequently
requested resuming the methodological plan after advisor preparation/audit;
Task2-9 prospectively specifies its single next case and commit-before-run gate.

## Exact change and derivation

For one modality-factor loading prior
$g(f)=(1-\pi)\delta_0+\pi(\lambda/2)\exp(-\lambda|f|)$,
keep the completed sequential coordinate posterior $q$ fixed. Define
$S=\sum_{j=1}^{p}\Pr_q(F_j\ne0)$ and $T=\sum_j E_q|F_j|$.
Its expected log prior is

$$
S\log\pi+(p-S)\log(1-\pi)+S\log(\lambda/2)-\lambda T.
$$

Differentiation gives $\pi_{new}=S/p$ and $\lambda_{new}=S/T$.
These are separate modality-factor parameters. There is no additional
survival weight or inverse-feature-count likelihood scaling.
For $S=0$, return delta-zero and rate NA (unidentified).
Positive slab mass requires positive absolute magnitude; a nonrepresentable
rate fails loudly rather than introducing a floor, cap or support truncation.
$\pi=1$ is allowed when every posterior has no spike mass.

The same posterior's KL must change atomically:

$$
KL(q\Vert g_{new})=KL(q\Vert g_{old})-
  E_q\log\{g_{new}(F)/g_{old}(F)\}.
$$

The old KL uses the actual incoming-prior normalizers and coordinate
pseudo-data from the sweep, not pseudo-data recomputed afterward.
$E|F|$ adds the positive/negative truncated-Normal pieces; $|EF|$ is wrong
for a symmetric posterior. Slab probabilities add those pieces directly,
avoiding cancellation in `1-spike`. No posterior/projection moment is
changed at this prior step; its next sweep consumes the new prior.

## Implementation and numerical limits

`code/loading_prior_fixed_q.R` implements the sufficient statistics, maximizer
and KL correction. `multimodal_yfb_update_F_mk()` selects it only with
`loading_prior_update="fixed_q_em"`. Single exposes the same trailing argument;
multimodal exposes it in `control`. Historical is the default and retains
its existing returned objects. Supported: point-Laplace F, per-factor Normal
beta, no pruning; multimodal also requires damping=1. Unsupported combinations
fail rather than return a KL for a different/damped posterior.

The posterior moment kernel remains unchanged. Its existing inverse-Mills
asymptotic branch is approximate: independent integration at A=1/8.01²,
B=0, pi=1, rate=1 gives E|F|=.9710125 versus approximately .9709770
from that kernel (absolute difference3.54e-5). At h=32 the error is6.35e-10.
This is an identified numerical limitation, not a new fit or an explanation
of a saved trajectory. Do not claim exact arithmetic or a globally improving
implemented ELBO. Tail-density/normalizer accuracy remains a prerequisite
for interpreting a future extreme-state diagnostic. No production tail
arithmetic is silently changed in this attempt.

## Fixed prospective arms: one existing simulation only

After renewed authorization to run a scientific fit, evaluate one candidate
against the saved profiled Single cell021, SHA256
`5b83bf02f24603e9be1ebc337b7cb3f9beba9bf53d7c2a0458885d429132b820`.
Use the same bed seed20261026, noise.25, n400/validation120, p50, K7,
point-Laplace L/F, per-program Normal beta, w_S1, profile TRUE,
no supervised start/pruning, original300 cap/ELBO convergence rule.
The exact single change is `loading_prior_update="fixed_q_em"`, including
its necessary KL correction. The reference is reused, not rerun. Do not
change initialization, tolerance, cap, working arithmetic, prior family,
monitor definition, benchmark or test bed. No real/external data loaded.

Before execution, prepare a fixed one-case driver with fresh output refusal,
source/reference hashes, warnings/errors retained, ignored full fit/telemetry,
and compact aggregates. Commit that driver and its final registration before
running. Observe actual fixed-q prior gains and normalizer/moment accuracy;
compare native convergence, predictor/reconstruction changes and planted
recovery. Different stopping time after corrected KL is a possible consequence,
not permission to redefine convergence. This is a mechanism check, not an
adoption panel; one seed cannot establish recovery rates.

## Original validation criteria and stopping

The future adoption panel remains the original80 simulations:20 seeds at
each noise level in each model, recovery correlation>=.7 with correct beta
sign in>=16/20 each, mean false prognostic programs<=original reference+.5,
and every fit converged under its original rule with no norm runaway.
Real training comparisons require>=1 credible program with positive posterior
SD and nonconstant risk, and five finite paired CV scores satisfying the
original mean-delta>=-SE(delta) rule. The saved multimodal reference has
C=NA in fold1; that comparison cannot pass unchanged. It is not replaced by.5.
Rank and final once-only external requirements remain in the original handoff;
Single external lower95% deltaC must exceed-.02. No external selection.
The present implementation does not authorize that panel or external scoring.

Required software checks: independent mixed-measure quadrature for E|F|/KL,
maximizer against alternative parameters, support/no-information boundaries,
scale invariance, missing-normalizer failure, unchanged entire posterior sweep,
actual caller forwarding, unsupported combination failures, protected-base
default identity, and the full suite. Mathematical maximization is local to
the prior block and does not certify L/beta reference priors or Cox refresh.

Stop after the one authorized case: any error, nonfinite moments/KL, incorrect
prior gain, hash mismatch or unresolved density accuracy blocks interpretation.
A cap or failed scientific gate is retained as failure. No retry, larger count,
damping, weight or second remedy follows automatically. Regroup on that result.

## Code walkthrough

- `loading_prior_validate_stats(prob, mean_abs)`: checks equal feature-vector
  lengths, support probabilities and finite nonnegative magnitude; impossible
  mass/magnitude combinations error rather than being clipped.
- `loading_prior_laplace_stats(A,B,prior)`: reconstructs the just-completed
  scalar posterior using its own incoming prior; returns two length-p vectors.
  It reuses the existing base-R truncated-Normal helper, including its tail limit.
- `loading_prior_laplace_em(prob,mean_abs)`: two sums produce the unique interior
  maximizer or the explicit delta-zero boundary; output is a prior list.
- `loading_prior_laplace_log(prob,mean_abs,prior)`: evaluates the mixed-measure
  expected log prior, retaining -Inf for incompatible support.
- `loading_prior_fixed_q(details,prior)`: validates the actual diagnostic vectors,
  computes the prior gain and subtracts it from old positive KL; returns the
  new prior, positive KL (up to numerical error) and prior gain.
- Existing F-update: identical R/C++ sweep first, opt-in prior/KL replacement
  second. Single/multimodal callers forward the selector and reject unsupported
  modes. Historical outputs gain no new metadata. No dependencies added.


### Additional software-only check

Running the focused fixtures with `options(multimodal_yfb.engine="R")`
exposed a pre-existing named-scalar lookup defect: A/B feature names changed
`pos` to `pos.gene` in the posterior weight vector. The named/unnamed
posterior regression failed before a two-line `unname(A/B)` correction and
passes afterward. Scalar values and the C++ kernel are unchanged. Both
engines pass the same focused checks. This is a separate metadata fix, not
an additional analytical remedy; the historical R fallback previously errored
on these named inputs. The protected-base usual-path identity remains tested.


## Task2-9 code walkthrough

The fixed case driver adds no configurable grid. Its controls helper changes
one field of the existing profiled arm. Its diagnostic integration uses
half-line Gaussian kernels in stable coordinates and does not change fitting.
Its metric helper reuses original recovery/sign/credibility and frozen scoring
rules. The executable block verifies the reference hash, refuses old output,
records source hashes, fits once, and saves compact comparison/local-prior/
four-selected-density checks. Full fit and update telemetry remain ignored.
Warnings/errors are retained; the wrapper is restored even after failure.
Four density checks cannot certify the full trajectory or a common objective.
No new package dependency, parallelism or real-data input is introduced.
