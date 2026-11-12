# Independent review of tasks A and A2

Date: 2026-11-12. Author: Andrew Walther.

## Verdict

**Keep point-Laplace L and F on per-cohort z-scores.** Task A's point-exponential
arm passes recovery and simulation convergence, but fails false activity in every
simulation cell and fails the real-data criterion on both concordance and convergence.
**Neither A2 setup qualifies.** Shifted z-scores (SZ) and cohort-stratified X/sd (PES)
fail the same two criteria. Independent calculations reproduce the saved decisions.
All four requested fresh fits reproduce their saved metrics and loading matrices.

This supports retaining the comparator under the registered rule. It does not establish
that point-Laplace meets every acceptance criterion. Its near-threshold false-active
mean is 3.0, and only three of its five real-data fold fits converge. The registered
recovery, false-active and convergence gate still matters for scientific adoption.

Two qualifications should accompany the A2 result. Computing the SZ shift from all
273 patients can affect the fitted model through its reconstruction term, even though
a common shift cannot change concordance for fixed coefficients. PE and PES have nearly
identical loading directions, but substantially different column norms. They are not
identical fitted models.

## Scope and independence

The review followed [the review prompt](../prompts/2026-11-12-codex-task-a-review-prompt.md),
the reporting rules and task A section of the
[week plan](../plans/2026-11-09-week-plan.md), including the November 12 pooled-fit revision,
and the A and A2 registrations in the attempts log. The pooled-design decision was read
separately from the result entry. The findings notes and result-entry bodies were withheld
until a provisional verdict was written to a scratch file. Full blinding to the stated
choice was not achieved: the week-plan progress checklist and result headings encountered
while locating the registration already state that PL stays. The numerical verdict was
therefore established from the fit rows and code rather than from the findings narrative.

A separate reviewer audited A2's implementation, criteria and invariance argument without
reading the findings notes. The final review also received an independent statistical
check. Only this review and the reviews index are deliverables. No model, result CSV,
patient file or shared research-status document was changed by this review.

The substantive runners were read at repository commit `1bf8e6f6`. Other sessions advanced
main during the review. Source MD5 values recorded after the checks were:

| Source | MD5 |
|---|---|
| `results/benchmark_sim/run_lf_prior_comparison.R` | `5897c1a5d37635b40594344b8f3c4e45` |
| `results/benchmark_sim/run_lf_prior_followup.R` | `f2522a0474ee93132c87937cf577699b` |
| `code/fit_ssbmf_yfb.R` | `a6bff03fe8e783c07c9936149d304624` |
| `config/globals.yml` | `8981cbd3998476d45de0cb44c60fbbdd` |

## Task A: design, scoring and completeness

**Arms and inputs match the registration.** `lfp_fit()` calls the working model for PL,
with point-Laplace L/F, Normal β per program, canonicalization, w_S = 1, K = 7, no
pruning, 300 sweeps and tolerance 10⁻⁵. PE changes the L/F prior to point-exponential and
turns on `center_projection`. RANK does the same on the registered rank input and is
reported descriptively. No task A arm passes `strata_id`. The fitter centers the survival
design using the current training rows while retaining the original input for reconstruction.

The prepared inputs contain 273 patients and 2,218 genes. Recomputing X/sd within each
cohort reproduces PE, with nonnegative entries. Recomputing per-cohort z-scores reproduces
PL within 10⁻⁸. The prepared fold vector and outcomes are identical to S-11. Test-fold
sizes are 55, 55, 55, 55 and 53, with both cohorts represented in each fold. The three
S-9 cells and strong-survival S-14 cell use the specified generators and seeds 1 to 10.
The registered PL refits replace the original plan's reuse-only option.

**The primary concordance pairs patients within cohort.** `lfp_concordance()` uses
`Surv(time, status) ~ risk + strata(cohort)`, with `reverse = TRUE`. The risk is the fitted
Cox score. No held-out outcome chooses its sign. Independently counting comparable pairs,
including an event before a censor recorded at the same time, reproduces all 25 CV rows
across A and A2 within 1.2 × 10⁻¹⁶. Recomputing pooled and within-cohort held-out partial
log-likelihood differences from the saved scores gives maximum discrepancy zero.
PES pooled scores retain arbitrary cohort offsets and its pooled concordance is descriptive
only, as the A2 findings state.

**Recovery, false activity and calibration use identical definitions across arms.**
The best absolute loading correlation must reach 0.7, its coefficient must be active under
|E β| ≥ 1.96 SD with positive SD, and its sign must agree with the loading direction.
An unmatched program is not exempted from false activity. Calibration is the validation
Cox slope with `timefix = FALSE`. Independent rescoring of all 120 saved PL, PE and SZ
simulation states reproduces recovery and false-active counts, with maximum numerical
discrepancy zero for loading correlation, validation C, calibration slope and the reported
SZ maximum coefficient ratios. This is a re-scoring check, separate from the four fresh fits.

**Failures are retained in the present analysis.** There are exactly 80 A simulation rows,
18 A real-data rows, 40 A2 simulation rows and 12 A2 real-data rows, with unique design
keys and no missing fits. The saved fit records report no fit errors or warnings.
`s9_score()` returns unrecovered, nonconverged and seven false-active programs on a fit
error. The summaries expand the registered design and assign those same penalties to
missing simulation checkpoints. Missing real fits are marked nonconverged, and undefined
paired metrics cannot pass criterion 4. Nonconverged fits remain visible and fail the
convergence criterion. The separate cached-reference completeness caveat is listed below.

All A/A2 CSVs were read, including input, program, overlap, reuse and paired tables.
The task B CSVs sharing the parent directory were inventoried but are outside this review.

## Independently recomputed criteria

Cell order below is S-14 0.05%, S-9 0.25%, S-9 1% and S-9 5%. Counts and means come from
`sim_fits.csv` and `a2_sim_fits.csv`, not the cell summaries.

| Arm/setup | Recovery counts | Mean false-active counts | Converged counts |
|---|---|---|---|
| PL reference | 5, 10, 10, 10 | 3.0, 0.9, 0.1, 0.0 | 10, 10, 10, 10 |
| PE, also PES one-cohort reference | 8, 10, 10, 10 | 2.6, 3.0, 4.0, 4.3 | 10, 10, 10, 10 |
| SZ | 6, 10, 10, 10 | 2.8, 3.1, 4.1, 4.4 | 10, 10, 10, 10 |

The real-data calculations join each setup's five fold rows to task A's PL rows. Each SE
is the sample SD of the paired differences divided by √5, as registered.

| Arm/setup | Mean within-cohort C | Difference from PL | Paired SE | Converged folds |
|---|---|---|---|---|
| PL reference | 0.657841724 | 0 | | 3/5 |
| PE | 0.572979213 | −0.084862511 | 0.030769417 | 0/5 |
| RANK, descriptive | 0.614549718 | −0.043292006 | 0.016372667 | 0/5 |
| SZ | 0.635143498 | −0.022698227 | 0.007906220 | 1/5 |
| PES | 0.628271441 | −0.029570283 | 0.015307666 | 0/5 |

PE, SZ and PES each pass criteria 1 and 3 and fail criteria 2 and 4. These verdicts
agree exactly with `criteria.csv` and `a2_criteria.csv`. The paired within-cohort and pooled
differences and SEs agree with `real_paired.csv` and all four PL/PE-reference rows of
`a2_real_paired.csv` to less than 10⁻¹⁴. Descriptive cell medians also reproduce the raw
fit rows. The five-fold SE is the registered comparison
statistic, not an independent-sample uncertainty estimate, because training sets overlap.
None of these candidate comparisons is within one SE of PL, so the week-plan rule for
repeating close comparisons is not triggered here.

## Four fresh fits

The runners were sourced into a scratch verification session. `PRIV` and `PRIV2` were
redirected to `/private/tmp/task-a-review-20261112/A` and `/private/tmp/task-a-review-20261112/A2`.
`lfp_run()` ran one PE simulation and one PE fold. `a2_run()` then ran one SZ simulation
and one PES fold, each with two PSOCK workers and one BLAS/OpenMP thread per worker.
No saved checkpoint was eligible for reuse in those empty scratch folders. The initial
sandbox attempt failed before fitting because its local PSOCK socket could not open.
The same command succeeded with the required execution permission.

| Task | Fresh fit | Main reproduced metrics | Comparison with saved fit |
|---|---|---|---|
| A | PE, S-9 1%, seed 1 | recovered, 5 false-active, 19 sweeps, C 0.9649682, slope 2.184517, 7 nonzero programs | all metrics except elapsed time agree within 10⁻¹⁰, EF difference exactly 0 |
| A | PE, real fold 1 | within-C 0.6053269, pooled-C 0.5253837, 0 active, 7 nonzero programs, not converged at 300 sweeps | same |
| A2 | SZ, S-9 1%, seed 1 | recovered, 3 false-active, 38 sweeps, C 0.9649682, slope 1.523144, 5 nonzero programs | same |
| A2 | PES, real fold 1 | within-C 0.6416465, pooled-C 0.5430933, 3 active, 7 nonzero programs, not converged at 300 sweeps | same |

The scratch command was `Rscript /private/tmp/task-a-review-20261112/rerun.R`. Its log and
the aggregate-only verification logs remain in that scratch directory for local inspection.
Patient-level rerun states remain local and are not committed. No production code was
modified, so the full software suite was not rerun for this documentation-only review.

## Cohort and mean-expression explanations

Recomputing every A/A2 program's training correlation and every full-data top-100 gene
overlap reproduces the saved tables exactly. The full PE fit has seven nonzero programs,
none survival-active. Absolute cohort correlations range from 0.6400 to 0.9977 and
mean-expression correlations from 0.6696 to 0.8878. Across the five CV fits the reported
maximum cohort correlation is about 0.998. This supports substantial cohort tracking,
consistent with PE's retained gene means. The input's median gene means are 11.0 for
TCGA and 15.5 for CPTAC, with within-cohort gene SD 1.

The correlation evidence supports the proposed mechanism but does not isolate its causal
contribution to the concordance loss. Cohort stratification restores nonzero β while
leaving the reconstruction input unchanged. Simulation still has excessive false activity
without any cohort difference. The simulation also bundles the prior, input offset and
projection centering, so it does not isolate the prior alone as the cause.

PES and PE match uniquely by the same-index loading direction in all six real-data fits.
Their cosine similarities range from 0.9999274 to 1. Their column norms differ. In the
full fit, PE norms range from 0.00746 to 0.11312 and PES norms from 0.07333 to 0.14720.
The corresponding PES/PE norm ratios are 1.30 to 14.67, reaching 16.98 in a CV fit.
Canonicalization uses global score variation for PE and within-cohort score variation
for PES. Direction agreement therefore does not establish identical F, coefficient
contributions or predictions. All seven programs remain nonzero, so these saved PES fits
do not exhibit the complete-program zeroing described in the earlier stratified scale warning.

## PL reuse checks

`pl_reuse_check.csv` has identical recovery verdicts and validation C to three decimals
for all 40 fits. Its five false-active disagreements are S-14 seeds 1, 4, 5, 7 and 8.
The saved S-14 E[F] and β are exactly identical to the new PL checkpoint for all five.
The older S-14 scorer at `run_yfb_near_threshold.R:68` always subtracts the best-correlated
program from activity. The shared scorer at `run_yfb_survival_weight_p2000.R:169` subtracts
it only when |r| ≥ 0.7. All five disagreements are unmatched fits, and the new count is
one larger. This is a scoring-definition difference, not a fit difference.

`pl_s11_reuse_check.csv` reproduces pooled C, convergence and active/nonzero program counts
in all five CV folds. Its fold-0 task A fields are blank because the writer uses the
CV-only subset. Directly inspecting both full-data checkpoints gives the same 219 sweeps,
convergence and 3 active/7 nonzero programs. Maximum loading and coefficient differences
are 6.16 × 10⁻⁹ and 3.71 × 10⁻⁷, respectively. Thus the full fits agree numerically,
without claiming byte identity.

## A2 implementation and limitations

**SZ implements the registered transformation.** For each gene it adds
c_j = −min_i z_ij. The same constant is used in both cohorts. Synthetic training z-scores
use their own minima, and validation z-scores receive the training constants. Validation
observations can consequently fall below zero. They are used for projection, not for a
nonnegative refit, so this does not violate the training-input constraint.

**The global minimum can affect held-out metrics through fitting.** For fixed F and β,
adding c to every row changes each η by the same scalar c′Fβ. Pairwise concordance,
calibration with a Cox baseline, and Cox partial likelihood are invariant to that scalar.
Training projection centering also removes c from the Cox design. However, the fitter
still reconstructs Y + 1c′, initializes from that matrix, and estimates nonnegative L/F
under their priors. There is no translation-invariance guarantee for those estimates.
Held-out expression therefore contributes to the fitted representation when it determines
c. Outcome-free preprocessing does not make that contribution absent. The reported CV
is the registered globally preprocessed comparison, not a proof of training-fold-only
prediction performance. A causal or numerical claim about the size or direction of this
effect would require a separately authorized comparison.

The same general qualification applies to inherited full-training cohort standardization
and gene selection in the prepared S-11 inputs. Those were fixed by the registration.
This review does not refit preprocessing or change the registered study.

**PES simulation reuse is justified for the stated one-cohort setting.** In
`fit_ssbmf_yfb.R:596–607`, one stratum subtracts the same training mean as pooled
centering. In `fit_ssbmf_yfb.R:127–145`, one stratum calls the same Cox core on the same
patients. Other settings are identical. The mathematical fit is the same, subject to
floating-point summation order. Reusing PE for criteria 1 to 3 is valid, but offers no
simulation test of heterogeneous cohorts, separate baselines or stratified scaling.

**The false-active coefficient ratios are descriptive.** Independent rescoring reproduces
the ratios. The median of each fit's largest ratio is 0.280, 0.113, 0.201 and 0.226 in
the four cells. The largest values are 0.408, 0.181, 0.323 and 0.508. Only 6 of 10 S-14
fits have a matched coefficient, so four fits have no defined ratio. All 10 fits still
count in the false-active criterion. Canonical score scaling makes these coefficient
magnitudes interpretable on the fitted score scale, but correlated programs prevent them
from measuring independent contributions or replacing the registered activity criterion.

**Synthetic concordance retains the historical time correction.** Calibration explicitly
uses exact times, as required. `frozen_reverse_cindex()` at `code/concordance_ci.R:37`
uses survival's default time correction. With these extremely dispersed synthetic times,
that can change C slightly. For PE S-9 1%, seed 1, default C is 0.9649682 and exact-time
C is 0.9648145. The convention is identical across arms and simulation C is not one of
A's four adoption criteria. This is a reporting limitation, not a decision-changing defect.

## Comparison with the findings notes

[Task A findings](../findings/2026-11-12-lf-prior-comparison-findings.md): the decision,
cell counts, concordance differences, cohort tracking, gene overlaps and PL reuse explanation
agree with this review. Its nonconvergence and one-cohort limitations are supported.
Retention of PL should continue to be distinguished from passing an absolute acceptance gate.

[A2 findings](../findings/2026-11-12-lf-prior-followup-findings.md): the separate no-adoption
verdicts, criteria, coefficient-ratio summaries and PE/PES direction agreement are supported.
The discussion already acknowledges different norms. Its statements that the common shift
cannot affect held-out concordance need to be limited explicitly to a fixed fitted model.
The present wording does not establish that global preprocessing leaves refitted CV
performance unchanged. No observed numerical discrepancy reverses either decision.

## Proposed changes for the owning sessions

No changes below were applied to code, findings or shared status files.

1. **P2, qualify the SZ CV-invariance claim.**
   `docs/findings/2026-11-12-lf-prior-followup-findings.md:51–54` and `:177–178`
   should distinguish fixed-score invariance from changes to the fitted factorization.
   The all-patient transform is created at `run_lf_prior_followup.R:227` and passed to
   training at `:150–156`. State that the registered comparison uses held-out covariates
   in preprocessing. Do not silently replace it with an unregistered refit.
2. **P3, include the full-data PL fit in the reuse table.**
   `run_lf_prior_comparison.R:466–467` should use all PL rows for the S-11 comparison,
   rather than `cv`. Fold 0 currently looks missing despite its saved fit being present.
   The direct full-data check above resolves this omission for the review.
3. **P3, validate cached A reference completeness before A2 summary.**
   `run_lf_prior_followup.R:356–373` inner-joins cached cell summaries without asserting
   all four expected cells and unique reference keys. An incomplete reference CSV could
   remove a simulation cell from criteria 1 to 3. All current reference cells are present,
   and the generated checkpoint summaries penalize missing fits correctly. The real-data
   criterion already requires five converged paired rows, so one missing reference fold
   alone cannot pass it. Add explicit expected-cell/fold and uniqueness checks before
   future reuse rather than treating the reference file as self-validating.
4. Keep the numerical decisions unchanged. Carry forward the retained model's S-14
   false activity and real-data nonconvergence, the single-cohort simulation scope, and
   the preprocessing and fold-SE limitations when reporting these negative results.
   Label this evidence as the registered K = 7 comparison rather than a new evaluation
   of the later K = 5 working rank.
