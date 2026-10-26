# What still keeps the supervised factorization from working: ranked roadblocks

Written 2026-10-26 for the advisor meeting. Read-only analysis: no model was fitted and no code was
changed. Numbers come from the files cited. The rebuilt package and its results are on the
unmerged branch `rebuild/unified-ssbmf`. Paths prefixed `rebuild:` are read with
`git show rebuild/unified-ssbmf:<path>`.

## The question

The method factorizes molecular data as Y ≈ L Fᵀ and links survival through a Cox model on the
projection scores, η = (Ỹ F) β. It should (a) converge to a defined optimum, (b) let survival shape
the loadings F so that a prognostic program that explains little variance can be found,
(c) predict external survival at least as well as the frozen benchmark (mean external C 0.627),
(d) do this with expression alone or with expression and methylation, and (e) report calibrated
uncertainty for β. Andrew already tracks two problems: survival hardly moves F because the
reconstruction term carries far more precision than the Cox term, and the multimodal fits did not
converge. This document asks what else stands between the current state and a model that does
(a)–(e) reliably.

## Summary

| Rank | Roadblock | Type | Among the two known issues? | Fixed by the rebuild? |
|---|---|---|---|---|
| 1 | Survival barely shapes the loadings, and what recovery there was came from the starting values | Modeling | Yes (known issue 1), with a new part on initialization | No |
| 2 | No principled way to set the survival weight, and a weight chosen by CV breaks the β intervals | Modeling and inference | Partly: a consequence of issue 1 that is a decision in its own right | No |
| 3 | In the joint fit, methylation dominates and the β prior collapses | Modeling and data | **Changed**: convergence is fixed, but the joint model now converges to an uninformative risk score | Convergence yes, usefulness no |
| 4 | The target (a 0.25%-variance prognostic program) may be undetectable at these sample sizes, and the simulations disagree with each other about whether supervision helps | Design of the goal and of the simulations | **New** | No |
| 5 | The evaluation cannot tell the candidate models apart | Evaluation | **New** | No |
| 6 | Rank and program identity are not determined | Modeling and evaluation | **New** | No |
| 7 | β intervals undercover with two blocks | Inference | **New** | Partly (one block calibrated) |
| 8 | Preprocessing choices are open, especially for methylation | Data | **New** (audited 10/26, not decided) | Partly |
| 9 | Biology of the rebuilt programs is not characterized | Interpretation | **New** | No |
| 10 | Software and comparator gaps | Software and evaluation | **New** | Partly |

Roadblocks 1–5 decide whether the method's central claim can be made at all. Roadblocks 6–10 are
work that remains once a configuration is chosen.

Each entry separates **Verified** statements (a number in a cited file, or a calculation from
cited numbers) from **Hypotheses** (explanations not yet tested).

---

## 1. Survival barely shapes the loadings

**What it is.** In the update for a loading F_jk, the reconstruction contributes a precision of
order τ_j Σ_i E[L²_ik], summed over n samples for each of p features. The Cox term contributes
precision of order E[β²_k] Σ_i P_i ỹ²_ij, which is limited by the number of events. With every
matrix entry treated as an independent Gaussian observation, the loadings follow the
reconstruction.

**Verified.**
- Survival share of the loading-update precision in the rebuilt model at w_S = 1: median 4 × 10⁻⁴
  (one block), 1 × 10⁻⁵ (two blocks), 1 × 10⁻¹³ in the 0.25% cell (rebuild:
  `docs/findings/2026-10-27-rebuild-first-results.md` §4.1). On the real expression data it is
  4.9 × 10⁻⁷ to 5.1 × 10⁻⁷ for the two survival-active programs (rebuild:
  `results/ssbmf/outputs/real_single/frozen_fit_roles.csv`, column `surv_share`). In the old
  code it was 6.3 × 10⁻¹¹ to 3.0 × 10⁻⁷ (`docs/plans/2026-10-16-round2-return-note.md`).
- The rebuild computes the objective exactly, so this is a property of the model, not a coding
  error.
- O4 fails. At 1% of variance every arm recovers the program (15/15), so supervision adds nothing.
  At 0.25%, the rebuilt model recovers it in 2/15 and 3/15 datasets against 4/15 for the two-step
  method (rebuild findings §2, §4.1).
- The survival message to F_·k is proportional to the current b_k = E[β_k] (rebuild:
  `derivations/unified_ssbmf/ssbmf_derivation.qmd`, "Absorbing zero"). A program that does not
  already point at the prognostic direction gets no pull toward it. The rebuild implements a
  survival-score start (`init = "survival"`) but it is **off by default and was not run** in
  any panel (derivation §6, item 4, and results §5, item 7).
- The older single-modality code used a supervised start (one program replaced by Ỹᵀu₀, the
  null-model martingale-residual direction). In its calibrated simulation (n = 400, β₃ = 3.3) the
  supervised model recovered the planted program in 17/20 and 9/10 datasets, against 0/20 and
  0/10 for EBMF followed by Cox (`docs/findings/2026-10-15-m6-m7-summary.md`, section B).

**Hypothesis.** The older recoveries came from the supervised start, not from survival acting in
the loading updates, because the survival share in that code was 10⁻¹¹ to 10⁻⁷. If this holds,
the method's one demonstrated advantage over the two-step method is an initialization, which a
two-step method could also use. Test: fit the rebuilt model with `init = "survival"` at w_S = 1
and w_S = 0 on cells C and C3. If both recover the program equally often, the start does the work.

**Fixed by the rebuild?** No. The rebuild fixed the arithmetic around the problem and showed that
it remains.

**Candidate remedies and costs.**
- Survival-score start in the rebuild (cheap: one panel rerun). Tests the hypothesis above. Does
  not by itself make the loadings survival-informed.
- Supervise a smaller survival-screened feature set. The 10/26 screen kept both benchmark
  survival programs with 251 genes, while an unsupervised screen of the same size lost them
  (`docs/findings/2026-10-26-survival-screen-cv.md`). Cost: the supervision then lives in
  preprocessing, and the claim becomes "screen, then factorize". An outside review warned this
  weakens the method's thesis. That warning is a judgement, not a finding.
- Survival weight w_S (see roadblock 2).
- Change the likelihood so it does not overcount correlated features, for example a
  low-rank-plus-diagonal noise model or a feature-group likelihood. Cost: a new derivation, and the
  one tuning-free version tried (effective-feature weight) made recovery worse.

---

## 2. No principled way to set the survival weight, and a CV-chosen weight breaks calibration

**What it is.** Requirement 5 of the rebuild plan asks the model to set the balance itself, with no
tuning parameter (`docs/plans/2026-10-22-unified-method-rebuild-plan.md`, Requirements). The
evidence so far says it cannot.

**Verified.**
- Effective-feature block weight, the only tuning-free rule: w_m ≈ 0.23 (one block) and 0.06 (two
  blocks). Recovery fell to 3/15 (C), 0/15 (C2) and 0/15 (C3). Dropped (rebuild findings §4.2).
- Survival weight from {1, 10, 100} by training CV: in C3, the program became survival-active in
  12/15 datasets and validation C rose to 0.603 (two-step 0.584). The loadings still matched the
  planted program only partly (mean |r| 0.40, threshold 0.7), and 90% β intervals covered the
  truth in **27%** of fits. In cells C and C2, coverage fell to 67% and 63% (rebuild findings §4.2).
- Large weights cause many damped or rejected steps (median 22 rejected steps per fit in C3)
  because the column-wise (Jacobi) message step overshoots when correlated features move together.
- DeSurv, the lab's frequentist counterpart, does not set the balance from the model either. It
  weights the reconstruction and Cox losses by (1 − α)/(np) and α·2/n_event, then rescales the Cox
  gradient by a balancing factor ζ, the ratio of the squared gradient norms, clipped at 10⁶
  (DeSurv SI Appendix, Section 2.3, Eq. 4). α is chosen by Bayesian optimization of CV C
  (selected α = 0.334 in PDAC). Its validated model also truncates each program to its top
  n_top = 270 genes before projection. With all 1,970 genes, 5 of 35 (k, α) settings passed
  adjusted external validation, against 18 of 35 with n_top = 270 (DeSurv SI, Table S2).

**Hypothesis.** A power likelihood exp{w_S ℓ_PL} with w_S > 1 treats the events as w_S times more
informative than they are, so the posterior for β is too narrow by roughly √w_S. This would
account for the 27% coverage. It follows from the form of the tempered posterior but has not
been checked against coverage at fixed w_S.

**Fixed by the rebuild?** No. The rebuild made the question precise.

**Decision needed.** Either (a) accept a CV-chosen weight, as DeSurv does, and give up calibrated β
intervals, or report them from a separate step (for example a Cox refit on the frozen programs, or
a bootstrap). Or (b) keep w_S = 1 and accept that the loadings are unsupervised. Or (c) change the
model (roadblock 1, last remedy). Calibration methods for tempered posteriors (general Bayesian
updating, SafeBayes) are named in the derivation but not checked or implemented. Option (a)
conflicts with requirement 5 and with O7 as written.

---

## 3. In the joint model, methylation dominates and the β prior collapses

**What it is.** Known issue 2 was non-convergence. The rebuild converges, so that problem has
changed form: the joint fit now converges to programs that carry no survival information.

**Verified.**
- Both selected real-data fits converged (joint 147 sweeps, expression-only 343) with no objective
  decreases. CV fold fits converged in 17/20 (joint) and 19/20 (expression only)
  (rebuild: `results/ssbmf/outputs/real_multimodal/frozen_fit_diagnostics.csv`).
- Joint fit: all three programs reconstruction-only, β posterior means about 10⁻⁸, shared prior
  variance about 2 × 10⁻⁷, survival share 10⁻¹⁵. Program 1 explains 35% of methylation variance and
  8% of expression variance (`frozen_fit_roles_joint.csv`).
- ICGC primary C: joint 0.554, expression only 0.707. Paired donor-level difference −0.153
  (95% interval −0.266 to −0.051). Training C of the joint fit is 0.577 (`icgc_C.csv`,
  `frozen_fit_diagnostics.csv`).
- Inputs: 3,000 genes on the log scale against 10,000 unscaled CpG β-values, one Gaussian
  observation per entry (rebuild findings §4.4).
- In the rebuild, a β of zero is not absorbing: the prior variance is fitted from pseudo-data, so
  it becomes positive once any program's score exceeds its noise. The joint collapse therefore
  means no program carries survival signal at the converged loadings. It is not the absorbing-prior
  defect of the old code (derivation, "Absorbing zero").

**Hypotheses.**
- With 10,000 CpGs and the same entry counting, shared activities L follow methylation variance.
  The prognostic expression programs then never form, so β has nothing to fit. This matches the
  numbers above but is not proven to be the only cause.
- Part of the −0.153 gap is a K effect: the joint fit was selected at K_init = 3, the bottom of
  its grid, and the expression-only fit at 7. The comparison is not at matched K.

**Fixed by the rebuild?** Convergence, yes. A useful joint model, no.

**Candidate remedies and costs.**
- Per-block weights chosen by training CV. Cost: one more tuning parameter and the calibration
  problem of roadblock 2.
- Reduce the methylation set toward the size of the expression set, or screen it by survival.
  Cost: the CpG-count question (10k to 305k) stays open.
- Separate activities per block plus shared ones, or a block-specific β. Cost: new derivation.
- Methylation scale: β-values against M-values (roadblock 8).

---

## 4. The goal may be out of reach at these sample sizes, and the simulations disagree

**What it is.** O4 asks for recovery of a program at 0.25% of variance, with n = 200 training
patients. Whether that is possible, and whether it is the regime that matters for PDAC, has not
been examined.

**Verified.**
- The rebuild's simulation (`rebuild:results/ssbmf/sim_design.md`): n = 200, 30% censoring, β = 1
  per SD of the true score, programs on disjoint sets of 10% of features, nuisance noise in groups
  of 20 features with within-group correlation 0.3.
- The old simulation that favoured supervision used n = 400, β₃ = 3.3 and a supervised start
  (roadblock 1). The two panels differ in sample size, effect size and start, and they reach
  opposite conclusions about supervision.
- The rebuild panel has 5 datasets per cell, 15 per arm across K_init values. A Fisher exact test
  of 2/15 against 4/15 gives p = 0.65, and the 95% interval for 4/15 runs from 0.08 to 0.55.
  At this size the panel cannot show that one method recovers the program more often than another.
  The full 20-dataset panel and the semi-synthetic designs were planned but not run.
- The real survival-active programs are not low-variance in this sense. In the rebuilt
  expression-only fit they explain 10.1% and 4.1% of variance (`frozen_fit_roles.csv`). DeSurv's
  prognostic program D1 explains 7% (DeSurv paper). In the 1% cell every method already succeeds.
- DeSurv's own simulation used p = 3,000, n = 200, k = 3, 150 marker genes per program with
  Gamma(3, 0.8) loadings, and β = (2, 0, 0). Prognostic programs "explain low variance relative to
  outcome-neutral signals", but the share is not stated (DeSurv paper, Methods, Simulation Studies).

**Hypotheses.**
- Back-of-envelope, not computed from data: in the rebuild design each 20-feature nuisance group
  has a leading eigenvalue of about 1 + 19 × 0.3 = 6.7 noise variances, and there are 100 such
  groups. A program carrying 0.25% of the total variance of 2,000 features has an eigenvalue of
  roughly 0.0025 × 2,000 × (average feature variance), about 5 or more. It is therefore of the same
  size as 100 nuisance directions. Variance cannot pick it out, and survival information from
  about 140 events must do all the work over a 2,000-dimensional search.
- The method may add the most where the program is between 0.25% and 1%. That range has not been
  simulated.

**Fixed by the rebuild?** No. This is a question about the target, not the algorithm.

**Candidate remedies and costs.**
- Map recovery against program variance share (0.25%, 0.5%, 1%) and effect size (β = 0.5, 1, 2
  per SD) for supervised and two-step arms with 20 datasets per cell. Cost: one Longleaf panel.
  It turns O4 into a statement of where supervision helps.
- Calibrate a semi-synthetic simulation to the real data: real L and F from the frozen fit, real
  noise covariance, cohort structure, n = 273.
- Decide with the advisors whether O4 should be judged at 0.25%, at 1%, or on the prognostic score
  instead of the loadings. The rebuild findings already raise this (§7).

---

## 5. The evaluation cannot tell the candidate models apart

**Verified.**
- Training CV over K: mean C 0.640, 0.652, 0.649, 0.656 with SE about 0.04
  (`rebuild:results/ssbmf/outputs/real_single/cv_by_K.csv`). The joint CV means range from 0.548
  to 0.586 with SE 0.036–0.048 (`cv_by_K_joint.csv`).
- External expression-only differences are small against their intervals: rebuilt − benchmark
  +0.002 (−0.021 to 0.024), rebuilt − two-step −0.001 (−0.010 to 0.009) (rebuild findings §4.3).
- ICGC primary has 50 donors and 30 deaths. Intervals are about ±0.10 (writeup Limitations).
  Reasonable fits on these same 50 donors have given C of 0.42 (stacked EBMF then Cox), 0.59
  (old selected expression-only), 0.69–0.70 (old joint framework), 0.554 (rebuilt joint) and
  0.707 (rebuilt expression-only) (`docs/findings/2026-10-15-m6-m7-summary.md`, section C;
  rebuild findings §4.4; `paper/writeup/writeup.qmd`).
- A fit with no survival-active program has a constant risk score, and C computed from rounding
  noise came out at 0.70–0.75 in one CV fold (`docs/findings/2026-10-26-survival-screen-cv.md`).
  C is undefined there and must be handled explicitly.
- The benchmark's starting rank was informed by external concordance, which makes its external
  evaluation development-informed (writeup Limitations). O5's threshold of 0.617 is anchored to
  that number.
- PACA-AU contributes 616 records from 568 people. 47 donors appear on both platforms and one
  appears twice in RNA-seq. The ICGC multimodal validation set overlaps the RNA-seq series
  (writeup Limitations). The rebuild's pooled bootstraps resample shared donors together, but
  older record-level intervals do not.

**Hypothesis.** External C at this sample size cannot rank the configurations that are now in
question (survival weight, screening, block balance). A decision rule built on C alone will follow
noise.

**Fixed by the rebuild?** Partly. The rebuild uses donor-level paired bootstraps, fixes the risk
direction in training and handles constant scores. Sample size is unchanged.

**Candidate remedies and costs.**
- Choose configurations on simulation recovery and training-CV partial likelihood, and use
  external C only once, for a frozen final model. Cost: discipline.
- Report one pooled donor-level estimate across the five external series, stratified by cohort,
  as DeSurv does with its pooled HR per SD. It gives more power than five separate C values.
- Compare joint and expression-only fits at matched K.

---

## 6. Rank and program identity are not determined

**Verified.**
- The one-standard-error rule chose the bottom of the grid in both real-data analyses
  (K_init = 5 of {5, 7, 10, 15}, joint K_init = 3 of {3, 5, 7, 10}), because CV C is flat in K
  within its SE (`cv_by_K.csv`, `cv_by_K_joint.csv`).
- K is meant to be an upper bound, with point masses switching programs off. On the real data no
  program was switched off: all 5, all 3 and all 7 programs are survival-active or
  reconstruction-only (the three `frozen_fit_roles*.csv` files). The old multimodal fits never
  removed a factor either (writeup).
- Program loadings change with the feature set. At fixed K = 7, median loading correlation with
  the 2,064-gene programs was 0.30, 0.33 and 0.57 at 5,001, 10,000 and 17,900 genes
  (`docs/plans/2026-10-16-round2-return-note.md`). Between 10k and 25k CpGs, three of seven
  methylation programs had loading r of 0.12, 0.04 and −0.09 (same file, unconverged fits).
- The role thresholds (lfsr ≤ 0.05, 0.1% of variance) are provisional (rebuild findings §5).
- The derivation states that the objective is invariant to (L/c, cF, β/c). The rebuild reports
  on a canonical scale and never rescales inside the loop, so this scale freedom is handled.

**Hypotheses.** Rotations and splits among correlated programs are only weakly identified under
sparse priors, so program identity may not be stable across seeds or K_init. No seed or K_init
stability study has been run on the rebuilt real fits.

**Fixed by the rebuild?** Scale, yes. Rank and stability, no.

**Candidate remedies and costs.** A stability study across seeds and K_init on the real data
(cheap, no new code). Report programs that recur, with their matching correlations. Rank chosen by
CV partial likelihood or by recurrence rather than CV C.

---

## 7. β intervals undercover with two blocks

**Verified.** 90% intervals cover the truth in 82% of fits overall: 85–90% with one block, 70–75%
with two blocks (O7, rebuild findings §2). Under the CV-chosen weight, 27% to 67% (roadblock 2).

**Hypotheses.** (i) Mean field ignores the posterior correlation between F and β, so the
uncertainty in the program score does not reach β. (ii) The derivation notes that the Jensen bound
slightly overstates the precision of β when Z_k is centered. (iii) Overcounted reconstruction
precision makes q(F) too narrow, which shrinks V_ik and thus the β variance. Any of these grows
with more features, which fits the two-block result. None has been isolated.

**Fixed by the rebuild?** Partly. One-block coverage meets O7.

**Candidate remedies and costs.** Report β intervals from a Cox refit on frozen scores, or from a
bootstrap over patients, alongside the variational ones. Cost: separate from the model's own
posterior. A joint (F_k, β_k) update would address (i). It needs a derivation and was listed
untried in round 2.

---

## 8. Preprocessing is open, especially for methylation

**Verified** (`docs/findings/2026-10-26-preprocessing-audit.md`).
- The two analyses differ in twelve places, including the TCGA expression source (RSEM against
  TMM-CPM), scaling (per-cohort z-scores against raw log scale) and how validation cohorts are
  scaled. Eleven alignments are proposed and none is adopted.
- Each external cohort is z-scored on its own mean and SD, so scoring uses the validation cohort's
  distribution.
- Methylation is used as raw β-values with KNN imputation whose neighbours are limited to probes
  with missing values. Cross-reactive and SNP-probe filtering is undocumented. The CpG screen keeps
  10,000 of 305,352 eligible probes by variance.

**Hypothesis.** The methylation scale and probe set feed directly into roadblock 3. Unfiltered
cross-reactive or SNP probes and β-value heteroscedasticity could create high-variance programs
that carry no survival signal.

**Fixed by the rebuild?** Partly. Feature intercepts and signed loadings remove the sign floor
(below), and the double log is shown harmless. The methylation choices and the cross-analysis
alignment are not decided.

**Remedies and costs.** Decide the eleven alignments (decisions, not experiments). Check probe QC
(already listed as a backlog item). Refit the joint model on M-values. Each refit changes the
inputs, so the frozen benchmark numbers must not be quoted for it.

---

## 9. Biology of the rebuilt programs is not characterized

**Verified.** The published labels (Program 7 adverse and basal-like, Program 3 protective and
classical, Program 5 stromal) belong to the frozen benchmark fit. The writeup notes that the
enrichment analyses share inputs, include a correction designed after an earlier result, use a
within-training subtype check and report unadjusted tests (`paper/writeup/writeup.qmd`,
Limitations). The rebuilt programs have not been matched to these or characterized. No CpG has
been annotated to genes.

**Hypothesis.** The two rebuilt survival-active programs correspond to the basal-like and classical
axes. This is untested.

**Remedies and costs.** Match rebuilt loadings to the frozen programs (cheap). Rerun the enrichment
pipeline on the chosen configuration only, once it is frozen.

---

## 10. Software and comparator gaps

**Verified.**
- `ssbmf` passes `R CMD check` with 216 tests. Vignettes are not written, the branch is not merged
  or pushed, and the old repository suite was not run on it (rebuild findings §2, §6).
- The full 20-dataset simulation panel and semi-synthetic designs have not run (sim design file).
- No runtime is reported at 305,000 CpGs. The rebuild plan leaves compiled code to be decided after
  timings (plan, Open questions 3).
- The DeSurv head-to-head on the same data was deferred by Andrew. DeSurv is installed in
  `data/round2_R_library/`, and the remaining blocker is authorization, not software.
- The column-wise message step is the source of rejected steps at large w_S. A feature-wise
  (Gauss–Seidel) update is proposed and untried (rebuild findings §7).

**Fixed by the rebuild?** Partly: the package exists and passes check.

---

## What is not a roadblock any more

| Former problem | Status and evidence |
|---|---|
| Objective not monotone; working Cox quadratic not a bound (refreshed objective fell in 18 of 50 refreshes) | Resolved. The rebuild maximizes a Jensen bound on E_q log PL, checked below Monte Carlo in 24/24 states. 180/180 simulation fits and every real fit had no step decrease beyond 1e-8 × \|objective\|; worst relative step −9.9 × 10⁻¹³ (rebuild findings §2, §3). |
| Non-convergence, including the multimodal fits (known issue 2, in its original form) | Resolved for simulation (180/180, max 279 sweeps) and every selected real fit (297, 147 and 343 sweeps). Residual: one expression CV fold fit at K_init = 7 hit 1,000 sweeps, and 3 of 20 joint CV fold fits did not converge. |
| KL computed against stale priors; rescaling inside the loop that did not transform the β prior | Resolved by design: every prior refit returns its posterior and KL in one `ebnm` call, and nothing is rescaled inside the loop (derivation §4.5, §2). |
| β absorbed at zero by the shared prior | Resolved as a mechanism: the prior variance is fitted to score-based pseudo-data, and β leaves zero within 5 sweeps in tests. The joint-fit collapse is a different problem (roadblock 3). |
| Risk-set underflow at 25,000 CpGs | Resolved: log-space risk sets, finite at \|η\| = 10⁴ (rebuild findings §3). |
| Double log on TCGA, Dijk and PACA-AU RNA-seq | Negligible for the benchmark: mean external C 0.625 against 0.627, risk r ≥ 0.99 (`docs/findings/2026-10-26-single-log-check.md`). Should still be corrected in new work. |
| Nonnegative loadings on z-scored data (47% of entries negative; 69% of reconstruction error on them, a floor of 55% of total variation) | Resolved in the rebuild by feature intercepts and signed (point-Laplace) loadings (preprocessing audit §4; derivation §2). |
| Expression-only prediction | Met: 0.628 against 0.627 (O5). |
| Two separate fitters with about 30 option flags each | Resolved: one `ssbmf_fit()`, one block as the M = 1 case. |
| Scale non-identifiability between L, F and β | Handled: canonical reporting scale only, prediction unaffected (derivation §2). |

---

## Claims from the outside review (2026-10-17), checked

| Claim | Check | Verdict |
|---|---|---|
| DeSurv rescales the Cox gradient by a balancing factor ζ | DeSurv SI Section 2.3, Eq. 4, with ζ clipped at 10⁶ | Verified |
| DeSurv validation depends on gene focusing (n_top = 270); with all genes most (k, α) settings failed | DeSurv SI Table S2: 5/35 settings significant with all 1,970 genes, 18/35 with n_top = 270 | Verified |
| DeSurv rank-transforms within sample and standardizes Z with training mean and SD | DeSurv paper, Methods, and SI, Survival Analysis | Verified |
| The imbalance arises because β must shrink to O(1/P) | The rebuild measures the share directly. Its stated mechanism is entry counting in the reconstruction (rebuild findings §7). The review's β-scale argument was not tested. | Imbalance verified, mechanism unconfirmed |
| Weighting each block by P̄/P_m fixes methylation dominance | Not tested. The related effective-feature weight lowered recovery. | Untested |
| A point-normal β prior will prune non-survival programs | Not tested. In the joint fit no program carried signal, so a different β prior would not create one. | Untested, and not a remedy for roadblock 3 |
| A supervised start "guarantees" the fit locks onto low-variance programs | The old supervised start recovered 17/20 in its simulation. Not run in the rebuild. "Guarantees" is not supported. | Partly supported |
| Deduplicating PACA-AU donors removes the dependence caveat | Counts confirmed (568 people, 47 shared, 1 repeated). It removes the caveat for pooled analyses only. | Verified for pooled analyses |
| A survival screen undermines the thesis | A judgement, not a finding. The 10/26 screen result is relevant evidence. | Opinion |

---

## Suggested order, for discussion

1. Run the survival-score start in the rebuild at w_S = 1 and w_S = 0 (tests the main hypothesis
   of roadblock 1, one panel).
2. Run the recovery map across variance share and effect size, with 20 datasets per cell
   (roadblock 4). This shows where supervision can help before choosing a balance rule.
3. Decide the survival-weight question (roadblock 2) and how β uncertainty will be reported
   (roadblock 7).
4. Only then work on the joint model's block balance (roadblock 3), after the methylation
   preprocessing decisions (roadblock 8).
5. Keep external cohorts for one final scoring (roadblock 5).
