# Stronger supervision with a fixed prior on the survival coefficients

Date: 2026-11-14. Author: Andrew Walther. Status: registered before any fit.
Registration label: S-17.

## Why this study

S-16 ([plan](2026-11-13-share-targeted-survival-weight-plan.md),
[findings](../findings/2026-11-13-share-targeted-weight-findings.md),
[review notes](../reviews/2026-11-13-share-targeted-weight-handoff.md)) set the survival
weight w_S from a target share of the loading curvature along each program's direction. On
the real PDAC data that share is about 10⁻⁴ at w_S = 1. Weights for shares of 0.05 to 0.5
(w_S of about 430 to 8,200) never converged. In those fits the ELBO rose every sweep while
the largest |η| grew to about 95 and the largest |β| to about 26.

Canonicalization rescales each program score YF_k to unit SD, so β_k is the log hazard
ratio per SD of program k. The current β prior is a Normal whose variance is re-estimated
from the data every sweep. Nothing bounds β. With 2,218 genes and about 220 training
patients, a combination of genes can nearly order the survival times, so once survival
dominates the loading update the fit keeps increasing the risk score. This reading is
consistent with the traces but not proven by them.

DeSurv avoids this with an elastic-net penalty on β whose strength λ is fixed and chosen
by CV. The Bayesian counterpart here is a **fixed prior on the standardized coefficients**:
β_k ~ N(0, σ_b²) (ridge) or a Laplace prior with fixed scale (lasso). Because the scores
have unit SD after each canonicalization, σ_b limits how large a hazard ratio per SD the
model can claim at a fixed point. Within a sweep the loading update can still raise the
score SD, so whether a fixed point exists at large weights is what Step 1 tests. For
reference, the K = 5 real-data coefficients at w_S = 1 are 0.21 to 0.37 in absolute value,
and DeSurv's pooled external HR per SD (1.5) corresponds to 0.41.

**Question.** With β bounded, do large survival weights converge, do they change the
real-data programs, and does that change hold up in training CV and on a permuted-survival
null?

## Fixed choices

- Fitter: `fit_ssbmf_yfb(canonicalize = TRUE, survival_block_weight = TRUE)`, point-Laplace
  L and F, pooled training without strata, 500 sweeps maximum, as in S-16.
- **β prior:** `beta_prior_mode = "shared"`, `beta_prior_fix = TRUE`, with
  - Normal (primary): `prior_beta = "normal"`, `beta_prior_init = ashr::normalmix(1, 0, σ_b)`;
  - Laplace (secondary, Step 1 only): `prior_beta = "point_laplace"`,
    `beta_prior_init = ebnm::laplacemix(1, 0, σ_b / √2)`, which has the same variance.
    (`ebnm` rejects a Laplace `g_init` under `prior_family = "normal"`.) The S-16 runner
    hardcodes `prior_beta = "normal"`, so it needs an opt-in override.
  A fixed elastic-net prior (Normal × Laplace) is not available in `ebnm` and is not built
  here. The two arms are its two ends.
- σ_b ∈ {0.1, 0.25, 0.5, 1}. On the real data, with a coefficient standard error near
  0.1, these shrink an unpenalized estimate by factors of about 0.5, 0.86, 0.96 and 0.99.
  The larger values bound β only near separation, which is the case S-17 targets. σ_b = 0.1
  tests a prior that also shrinks ordinary estimates.
- **Caveat on activity counts.** The estimated per-program prior can set E[β_k²] near 0
  for an inactive program, which switches off its survival term in the F update. A fixed
  N(0, σ_b²) prior cannot, so at large weights every program receives survival curvature,
  while the weight is calibrated on one program. The 1.96-SD activity rule also becomes
  close to an unshrunk Wald test. Active and false-active counts and the program rule are
  therefore not directly comparable with the references. Report the per-program
  directional shares for every arm alongside the counts.
- Survival weight: the S-16 weight function `share_targeted_ws()` with the Amendment 1
  definitions, target directional shares s\* ∈ {0.05, 0.25, 0.5}, the pilot fitted at
  w_S = 1 **with the same fixed β prior**, up to 3 recalibrations, cap 10⁵, guard as in
  S-16. No ramp.
- Each σ_b also gets an arm at w_S = 1 with the fixed prior. It separates the effect of
  the prior from the effect of the weight.
- References: w_S = 1 with the current per-program estimated prior (the working model) and
  the two-step fit.
- Real data: the K = 5 working configuration, the S-11/S-15 folds and fold-specific
  training centers, as in S-16.
- Every quantity used to set a weight comes from training data only, inside CV from the
  fold's training data. σ_b is a fixed hyperparameter chosen by training CV. It is never
  estimated from the outcome inside a fit.
- Fitter defaults do not change. Any new control is opt-in. Work on the branch
  `s17-bounded-risk-supervision`, merged into `main` only after review.

## Step 1. Convergence screen

Fits, for each prior family (2) and σ_b (4):

- **Real full training data:** the fixed-prior w_S = 1 pilot, then each target s\* through
  the weight function with recalibration (up to 3 passes), every fit traced.
- **S-9 p = 2,000, 1% variance, seed 1:** fixed weights w_S ∈ {40, 100, 1,000}, which did not
  converge in S-16 (2 × 4 × 3 = 24 fits). The planted standardized coefficient there is
  about 6.7, so every σ_b ≤ 1 is badly misspecified. These fits are convergence traces only.
- **S-13 weak cell, seed 1:** the pilot and each target, as for the real data.

Record the S-16 per-sweep trace (`trace_s16 = TRUE`): ELBO and its change, max |η|, max
|E β|, ‖E F_k‖ and achieved directional share. Also record, per program and sweep, the
score SD c_k before canonicalization and β before and after it. Canonicalization runs only
at the end of a sweep, so the F update can raise c_k, and β ← c_k β follows. A fixed β prior
therefore does not guarantee a bounded η, and a fit can meet the relative-ELBO tolerance
while η is still moving, because the ELBO is dominated by w_S times the Cox term.

**Stable convergence** means the fitter's convergence flag plus, over the last 50 sweeps,
a relative change in η below 10⁻³ per sweep and a change in max |η| below 1%.

**Pass for a (family, σ_b) pair:** on the real data, the pilot and the final calibration pass
of every target converge stably within 500 sweeps, with max |E β| ≤ 4σ_b. The S-9 and weak
fits are reported but do not enter the pass rule. Only the Normal family goes on to
Steps 2 and 3. The Laplace results are reported for comparison. **If no Normal σ_b passes,**
report the traces and stop after writing the handoff.

Choosing which σ_b go forward uses the full training data, but only for convergence, not
for any outcome-based performance measure. The leakage into Step 2 CV is negligible.

## Step 2. Real data, training CV

For each passing σ_b, five folds: the w_S = 1 fixed-prior arm and each target s\*, plus
both references. Report CV partial log-likelihood, CV C (mean and SE), paired C differences
from the working model and from the two-step fit, folds converged, weights, achieved
shares and cap use.

**Eligibility:** at least 4 of 5 folds and the full fit converge stably (Step 1 definition).

**Selection rule.** Candidates are the working model, the fixed-prior w_S = 1 arms and the
target arms. The two-step fit is a reference and cannot be selected. Let the best value be
the highest mean CV partial log-likelihood among eligible candidates (the mean fold
contribution, as in S-16). Every eligible candidate within 0.5 of the best is tied. Among
tied candidates, choose the working model if it is tied, then a fixed-prior w_S = 1 arm,
then the smallest target. Remaining ties go to the larger σ_b, the weaker imposed prior.
S-16's arms differed by about 0.02 in this quantity, so the working model may well be
selected. In that case the result rests on the program-structure and CV C comparisons,
which are reported for every eligible arm.

Then fit every eligible arm on the full training data and run the S-16 Step 5 structure
analyses on it with `results/benchmark_sim/share_targeted_structure.R`: loading correlation
matrix against the two-step fit, per-program Δℓ and semi-partial R², survival-active
programs and the program rule, and signature annotation.

## Step 3. Simulation screen

Cells: the S-13 weak cell (directional share about 0.01 at w_S = 1, the closest to the
real data) and its permuted-survival null, seeds 1–10, at every σ_b that passed Step 1:
the fixed-prior w_S = 1 arm and each s\*, plus both references. In the S-14 strong 0.05%
cell (share about 0.56 at w_S = 1), fit only the fixed-prior w_S = 1 arms and the
references, to check that the fixed prior keeps the low-variance recovery that supervision
gives there. Use exact event times in every `coxph` call.

Report per arm: weight, achieved share, recovery, false-active count, the largest
standardized |β| among false actives, convergence, validation C, the exact-times
calibration slope and the matched program's loading correlation with the two-step fit.

**Null check (replaces the S-16 risk-variance gate).** Under permuted survival every
active program is false and carries all of the risk variance, so that fraction does not
measure harm. A (σ_b, s\*) arm passes if, over the ten null datasets, (a) the mean
false-active count is ≤ 1 and (b) the mean of the largest standardized |β| among false
actives (0 for a dataset with none) is ≤ 0.25. Report the same numbers for both
references. The weights come from each null dataset's own pilot, guard included.

- When several targets for one null dataset hit the cap, the larger targets are capped
  duplicates. They take the metrics of the fit they duplicate and count in the
  denominator.
- A fit that does not converge counts at its final state.
- The criterion (b) mean is meant to catch a single runaway fit, as in S-16 null seed 10
  (standardized |β| 28.5).
- The paired-weight null (each null dataset fitted at its signal dataset's weight) is run
  and reported as secondary.

## Adoption candidate

The Step 2 selected arm is a candidate for adoption only if its (σ_b, s\*) passes the Step 3
null check, its weak-cell recovery is at least that of the working model, and its fixed
prior keeps S-14 low-variance recovery at least that of the working model. Adoption is
Andrew's decision after review.

## Step 4. External scoring, once, after review

As in S-16 Step 6. It runs only after the reviewing session and Andrew approve. If the
working model is selected, nothing new is scored.

## Outcomes and how they will be read

- **Converges, programs change, CV C holds or improves, null passes:** candidate for
  adoption.
- **Converges and programs change, CV C falls:** strong supervision reorganizes the
  programs without improving prediction on these data.
- **Converges, programs do not change:** bounding β removes the instability, and the
  real-data loadings are insensitive to supervision at the tested shares.
- **Does not converge with β bounded:** the growth is not only in β. Report the traces.
  The next candidate is a bound on the loadings' fit to survival.

## Execution

GPT-6.1 Sol (Opus 5.5 in Claude Code) plans and orchestrates and makes every judgment
call. GPT-6.1 Terra (Sonnet 5.5) subagents do the bounded execution. The handoff prompt is
[`docs/prompts/2026-11-14-codex-bounded-risk-prompt.md`](../prompts/2026-11-14-codex-bounded-risk-prompt.md).

## Outputs

- Runners extending the S-16 runners under `results/benchmark_sim/`, and tests.
- Aggregate outputs under `results/benchmark_sim/outputs/bounded_risk_2026-11-14/`.
  Patient-level data and fit checkpoints under the ignored `data/` folder.
- `docs/findings/2026-11-14-bounded-risk-findings.md` and
  `docs/reviews/2026-11-14-bounded-risk-handoff.md`.
