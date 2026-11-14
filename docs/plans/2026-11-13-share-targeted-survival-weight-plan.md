# Survival weight set by a target share of loading precision

Date: 2026-11-13. Author: Andrew Walther. Status: registered before any fit.
Registration label: S-16.

## Why this study

The 2026-11-13 advisor meeting asked whether the survival weight w_S is reasonable and
whether it can give survival more influence on the fit. This plan answers both.

**What w_S is.** In the default fit the objective is

E[log p(Y | L, F, τ)] + w_S · E[log PL(YFβ)] − KL terms,

where PL is the Cox partial likelihood. Because w_S · log PL = log(PL^w_S), this targets a
power (generalized-Bayes) posterior in which the survival likelihood is raised to w_S and
the expression likelihood keeps power 1. At w_S = 1 it is the ordinary posterior.
Multiplying the likelihood itself by a constant would only add a constant to the
objective, so a power is the only reweighting that changes the posterior. The derivation
(`derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`, the
survival-weight box and Eq. AB-F) and the canonical code path (`code/fit_ssbmf_yfb.R`,
exact sequential sweep over genes) agree.

This study uses `survival_block_weight = TRUE` (the S-13 A1 arm), whose objective also
multiplies KL(q(β)‖g_β) by w_S. That keeps β at its untempered posterior while the
loadings see survival at weight w_S. It is not exactly the power posterior of the
derivation's box because the entropy of q(β) is scaled too. Its derivation is in
`derivations/SSBMF-YFB-model-derivations/components/survival_block_weight/`. With the
default (FALSE), the posterior SD of β shrinks by 1/√w_S, so at weights of 10³ to 10⁴ the
1.96-SD activity rule would mark programs active for that reason alone. At w_S = 1 the
two settings give the same fit.

**How much information survival gives the loadings.** In the sequential sweep, the
curvature of the objective in the loadings F_k of program k has two parts:
expression, τ̄ · Σ_i E[l_ik²] times the identity, and survival,
w_S · E[β_k²] · YᵀWY, where W holds the Cox working weights.

- *Per gene* (the share recorded in the fit trace), the survival-to-expression ratio is
  w_S · E[β_k²] · Σ_i W_ii y_ij² / (τ_j · Σ_i E[l_ik²]). Under canonical scaling and
  z-scored genes this is about w_S · w̄ · b_k² / (τ · (p · PVE_k)²), with b_k the
  standardized log hazard ratio per SD of the program score, w̄ the mean Cox working
  weight, PVE_k the program's share of expression variance and p the number of genes.
- *Along the program's own direction* f̂_k = E[F_k] / ‖E[F_k]‖, the ratio is

  r_k = w_S · E[β_k²] · f̂_kᵀ YᵀWY f̂_k / (τ̄ · Σ_i E[l_ik²]) ≈ w_S · w̄ · b_k² / (τ · p · PVE_k),

  which is p · PVE_k times the per-gene ratio. This is the quantity that decides whether
  survival can turn the program, so it is the one this plan targets. The directional
  share is s_k = r_k / (1 + r_k). At s_k = 0.5 survival and expression contribute equal
  curvature along the program axis.

Rough values from reported fits (to be checked in Step 1): in the S-9 p = 2,000 cells the
directional share at w_S = 1 is already near 0.5, which may explain why w_S = 1 to 10 works
there. On the K = 5 real-data fit the per-gene share is 4.9 × 10⁻⁷ to 1.2 × 10⁻⁶ (K = 5
exhibits findings), so the directional share is about 10⁻⁴. A directional share of 0.5
would need w_S of order 10⁴. The S-15 grid (w_S up to 50) could not reach it, and in S-9
and S-13 fits stopped converging from w_S = 20 to 40 at p = 2,000. Supervision strong
enough to change the real-data programs has not been tested.

**What this plan changes.** The weight is set from the data through a target directional
share s\*, which means the same thing in simulation and on real data. A raw w_S does not.

**Motivation from DeSurv.** Young et al. report that survival supervision and standard NMF
give different factor structures in PDAC at the same rank (k = 3). DeSurv concentrates
survival in one factor that is not the highest-variance factor. Its parsimony comes mainly
from the one-SE rule applied to a flattened C-index surface, so the evidence for
supervision is the comparison at the same rank plus external performance. Step 5 follows
its results sections.

## Fixed choices

- Fitter: `fit_ssbmf_yfb(canonicalize = TRUE, survival_block_weight = TRUE)`,
  point-Laplace L and F, per-program Normal β, pooled training without strata.
- Maximum sweeps: 500 for every S-16 fit, set by overriding `cfg$cavi$max_iter` in the
  runners. The w_S = 1 references are refitted at 500 so all arms match.
- Real data: the adopted K = 5 working configuration (TCGA + CPTAC, 273 patients, 139
  deaths, 2,218 genes by combined mean and variance rank, per-cohort z-scores), the same
  five training folds as S-11 and S-15, and fold-specific training centers.
- Target directional shares: s\* ∈ {0.01, 0.05, 0.1, 0.25, 0.5}. References: w_S = 1 and
  the two-step fit (w_S = 0 loadings, then an exact-times Cox fit on the live program
  scores).
- Every quantity used to set a weight comes from training data only. Inside CV, from the
  fold's training data only.
- Fitter defaults do not change. Every new control is opt-in.
- Work on the branch `s16-share-targeted-weight`. It is merged into `main` only after
  review.

## Step 1. Check the share formulas and write the weight function (no new fits)

1. From saved S-9, S-13, S-15 and K = 5 checkpoints, compute for the matched (simulation)
   or survival-active (real) programs: the per-gene median share, the directional ratio
   r_k and share s_k, and both closed forms from b_k, PVE_k, w̄ and τ. Compare the per-gene
   median share with the share recorded in the fit trace (`fit_ssbmf_yfb.R` around line
   1290). Where a checkpoint lacks E[L²], use E[L]² and say that the value is a bound.
2. Write `share_targeted_ws()` in `code/share_targeted_ws.R`, with roxygen and tests.
   Input: a converged fit on the training data and s\*. Target program: the
   survival-active program with the largest standardized |b_k|. Output: the w_S that gives
   directional share s\* for that program at the input fit's state. The ratio is linear in
   w_S at a fixed state, so w_S = w_current · [s\* / (1 − s\*)] / r_k(w_current).
   - **No-signal guard:** if no program is survival-active at the w_S = 1 pilot, return
     w_S = 1 and flag the fit. The method does not amplify supervision when the pilot
     shows no survival signal.
   - **Cap:** w_max = 10⁵. Record when it binds. If two targets for the same data both
     hit the cap, fit only the smaller target and report the others as capped.
3. **Recalibration.** b_k and PVE_k move as the weight grows. Starting from the w_S = 1
   pilot, compute the weight, fit (warm start from the previous fit), recompute the weight
   at the new converged state, and repeat, up to 3 recalibrations. Stop when the achieved
   directional share is within a factor of 2 of s\* on the ratio scale. Report the
   achieved share and the number of passes. Arms are labelled by target and reported with
   their achieved share. An arm that misses the tolerance after 3 passes is reported as
   such and stays eligible.

**Check:** the per-gene median share computed from checkpoints matches the trace to
rounding. The closed forms are within a factor of 2 of the exact values in every cell.
If a closed form is off by more than a factor of 4, report it and continue. The weight
function uses the exact quantities, not the closed forms.

## Step 2. Make large weights converge

1. Diagnose. Trace per-sweep ELBO, its change, max |η|, max |E β|, ‖E F_k‖, the change in
   the expansion point and the achieved directional share for: S-9 p = 2,000, 1% variance,
   seed 1, at w_S ∈ {10, 40, 100, 1,000}; and the real full training data at the weights
   Step 1 gives for s\* = 0.25 and 0.5. Name the failure: (a) the quadratic expansion of
   the Cox likelihood is not a lower bound, so large steps overshoot and the ELBO falls;
   (b) separation, with β and the loading scale growing without bound; or another
   mechanism that the traces show.
2. Fix with the smallest opt-in change that the diagnosis supports. Candidates: a
   geometric ramp of w_S from 1 to the target over the first sweeps, with convergence
   judged only at the target; step-halving of the survival part of the update when the
   ELBO falls; or a fixed-curvature quadratic bound on the Cox term. Write the change's
   derivation in the findings note, and run
   `derivations/SSBMF-YFB-model-derivations/verify_derivation.R` if an update changes.
3. Tests: at w_S = 1 with the new control off, results are identical to the current
   fitter. Add a regression test on a small case that diverges without the fix.
   `Rscript tests/run_tests.R` passes.

**Pass:** the diagnostic fits at every weight up to the real-data s\* = 0.5 weight
converge within 500 sweeps, with no ELBO decrease larger than the stopping tolerance after
the ramp. **If no fix passes,** report the diagnosis, and run Steps 3 to 5 only at targets
whose weights converge without a fix. Say which targets were dropped.

## Step 3. Simulation screen

Cells, seeds 1–10:

- **Weak survival** (the S-13 weak cell, oracle C about 0.68). Closest to the real-data
  regime.
- **Its permuted-survival null.**
- **Low-variance program:** the S-14 strong-survival 0.05% cell, where A1 at w_S = 5
  recovered 9 of 10 against 5 of 10 at w_S = 1 and 0 for the two-step fit.

For each dataset: the w_S = 1 pilot, then a fit at each s\* through the Step 1 weight
function (with recalibration), plus the two-step reference. No CV in the screen. Use exact
event times (`timefix = FALSE`) in every `coxph` call. A target below the pilot's own share
gives w_S < 1. Fit it and label it.

Report per arm: weight, achieved directional and per-gene shares, recovery and false
activity under the S-9 definitions, the false actives' share of risk-score variance and
largest standardized |β|, convergence, validation C, the exact-times calibration slope,
and the matched program's loading correlation with the two-step fit.

**Null check (registered).** In the null cell, weights come from each null dataset's own
pilot through the weight function, guard included, because that is how the method would
run on data without signal. A target passes if mean false-active is ≤ 1 and the false
actives carry under 1% of the risk-score variance. As a secondary check, also fit each
null dataset at the weight computed for the paired signal dataset (the S-13 approach) and
report it. A target that fails the primary null check is not carried to real data.

## Step 4. Real data, training CV

For each of the five folds and each target that passed Step 3: fit the w_S = 1 pilot on the
fold's training data, compute the weight with recalibration, fit, and score the held-out
fold. Report CV partial log-likelihood, CV C (mean and SE), folds converged, weights,
achieved shares, passes and cap use.

**Eligibility:** an arm needs at least 4 of 5 folds converged and a converged full fit.

**Selection rule:** among eligible arms (w_S = 1 and each s\*), the one with the highest
mean CV partial log-likelihood. Differences within 0.5 log-likelihood units go to the
weaker supervision. Also report paired CV C differences from w_S = 1 and from the two-step
fit, with SEs.

Then fit every eligible arm on the full training data for Step 5.

## Step 5. Structure of the programs (every eligible arm, full training data)

Against the K = 5 two-step (EBMF) fit at the same rank:

1. The K × K absolute loading correlation matrix with the EBMF loadings, and the matched
   correlation per program.
2. Per program, on the training data:
   - Δℓ_k = ℓ_full − ℓ_(−k), where ℓ_full is the maximized Cox partial log-likelihood
     (exact times) on all K program scores and ℓ_(−k) is the same model refitted without
     program k.
   - Semi-partial R²_k = (‖Y − Σ_(k'≠k) E[l_k'] E[f_k']ᵀ‖²_F − ‖Y − Σ_k' E[l_k'] E[f_k']ᵀ‖²_F)
     / ‖Y‖²_F, on the centered, scaled training input. This is the extra fraction of
     expression variance explained by program k given the others.
3. Survival-active programs, the program rule from the S-11 DeSurv-style analyses, and
   how many programs carry survival.
4. Annotation of changed programs against Classical, basal-like and CAF signatures with the
   existing functions in `code/pathway_enrichment.R` and the S-11 DeSurv-style scripts.

These are descriptive. They answer whether supervision at a given share changes which
programs the model finds.

**Stop here.** Write the findings note and the review handoff, and push the branch.

## Step 6. External scoring, once, after review

After the reviewing session has checked Steps 1 to 5 and Andrew has approved, score the
selected arm, w_S = 1 and the two-step fit on the external cohorts with
`frozen_reverse_cindex()` and the documented score direction: pooled C and pooled HR per
training SD, with all records from each PACA-AU donor kept together. If w_S = 1 is
selected, its external result is the existing one. External results for the other
eligible arms may be reported only in a separate table marked as post hoc sensitivity,
not used for any choice (as in DeSurv's k × α tables).

## Step 7 (if time remains before the handoff). Wider simulation and rank surface

Five-fold CV selection of s\* in the S-9 0.25%, 1% and 5% cells and the other S-14 cells.
On real data, the CV surface over K ∈ {3, 5, 7} × the eligible arms. Both are descriptive
and do not change the Step 4 selection.

## Outcomes and how they will be read

- **Supervision changes the programs, CV C holds or improves, the null check passes and
  fits converge:** a candidate for adoption, decided by Andrew after review.
- **Programs change but CV C falls:** reported as a finding. Strong supervision reorganizes
  the programs without improving prediction on this cohort.
- **Programs do not change up to s\* = 0.5:** the real-data posterior is insensitive to
  supervision at every share tested.
- **Large weights cannot be made to converge:** reported with the diagnosis. The share
  reached on real data is the limit of the current fitter.

## Execution

GPT-6.1 Sol (Opus 5.5 in Claude Code) plans and orchestrates and makes every judgment
call. GPT-6.1 Terra (Sonnet 5.5) subagents do the bounded execution: runners, tests, fit
grids, tables. The handoff prompt is
[`docs/prompts/2026-11-13-codex-share-targeted-weight-prompt.md`](../prompts/2026-11-13-codex-share-targeted-weight-prompt.md).

## Outputs

- Code: `code/share_targeted_ws.R`, any opt-in fitter control, tests, and runners under
  `results/benchmark_sim/`.
- Aggregate outputs under `results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/`.
  Patient-level inputs and outputs, and fit checkpoints, under the ignored `data/` folder.
- `docs/findings/2026-11-13-share-targeted-weight-findings.md`.
- `docs/reviews/2026-11-13-share-targeted-weight-handoff.md` for review before Step 6 and
  before any status file, progress-book chapter or Chapter 4 text is changed.
