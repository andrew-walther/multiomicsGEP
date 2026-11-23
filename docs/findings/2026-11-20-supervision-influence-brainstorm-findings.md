# Making survival supervision influence the factorization: diagnosis and candidate solutions

Date: 2026-11-20. Author: Andrew Walther.

This note collects the brainstorm and literature review on why survival has almost no
influence on the SSBMF loadings at the working weight, why larger weights fail, and
which changes are most likely to help. It adopts nothing. The candidates are ranked for
discussion before any derivation or registration. An independent statistical review
corrected the first draft, and its points are incorporated below.

## Summary

- **At w_S = 1 the model behaves as specified.** Along each survival-active program,
  survival supplies at most about 1.5 × 10⁻⁴ of the loading precision on the real data
  (per gene, about 10⁻⁶). The joint posterior is therefore the unsupervised
  factorization. Changing the loadings requires a deliberate departure from the joint
  posterior, which is what DeSurv's weight α does.
- **Large weights overfit through low-variance expression directions.** With 2,218 genes
  and 273 patients, some gene combinations with little expression variance order the
  training survival times almost perfectly. Above a weight of roughly 100 the weighted
  likelihood itself prefers such a solution: moving from the w_S = 1 fit to the
  degenerate fit gains 403 nats of Cox log likelihood and loses about 44,600 nats of
  expression log likelihood. The per-gene noise precision then collapses on a few genes
  and amplifies the failure. The fitting algorithm also stops increasing its objective,
  so part of the problem is numerical.
- **Upweighting alone moves the loadings the wrong way.** In the real-data-like
  simulation, w_S = 10 converges and moves the loadings (correlation with the two-step fit
  about 0.77), but it loses a weak program that both w_S = 1 and the two-step fit find, and
  validation C falls by about 0.08. Larger weights do not converge.
- **DeSurv avoids this by construction:** nonnegative unit-norm programs that cannot form
  gene contrasts, a fixed common noise level, a survival push capped at a share of the
  reconstruction gradient, and steps that must lower the stated loss.
- **Leading candidates:** restrict or shrink the low-variance directions that survival can
  use; make the algorithm increase its objective at every step; guard the noise
  precision. Separate survival-driven programs are a structural alternative. Any candidate
  must improve recovery and held-out C in the real-data-like test bed, not only change the
  loadings.

## 1. Why the working weight has no influence

Each gene's loading update combines an expression estimate and a survival estimate,
weighted by their precisions (YFB derivation, F-combination equation). The survival
share is ρ/(1 + ρ) with ρ = w_S · A_surv / A_gen. Here A_gen = τ_j Σ_i E[l_ik²], which grows
with the number of patients and the program's strength, and A_surv = E[β_k²] Σ_i W_ii y_ij²,
which is bounded by the Cox working weights. On the real data at w_S = 1, the share along
each survival-active program's loading direction is at most about 1.5 × 10⁻⁴, and the
median per-gene share is about 10⁻⁶ (share-targeted weight findings, 2026-11-13). The
share falls roughly like 1/p² as genes are added.

The near-threshold simulation of 2026-11-19 shows the same from the outcome side. With a
real-data-like survival signal (true-risk C 0.69), the supervised and two-step fits give
programs with median loading correlation 0.98 to 0.999. Only with true-risk C of 0.86 or
more does supervision move the loadings and recover a weak program the two-step fit misses.

This is the correct Bayesian answer when 2,218 Gaussian gene likelihoods share parameters
with one Cox likelihood. The partial factor model (Hahn, Carvalho and Mukherjee, JASA
2013) describes the same imbalance: the joint likelihood is dominated by the
high-dimensional block even when prediction is the goal.

## 2. Why large weights fail: evidence from the saved fits

The fixed-β-prior study (2026-11-14) saved final fits and per-sweep traces on the real
training data (273 patients, 2,218 genes, K = 5). The checks use the Normal prior with
σ_b = 0.5. Its pilot (w_S = 1) and 5% target (w_S ≈ 390) converged. Its 25% and 50%
targets (w_S ≈ 2,500 and 7,500) did not. The traces record loading norms, score SDs, β and
the linear predictor, but not τ, and only final states were saved.

| Final fit | Genes with τ_j < 10⁻³ | Loading norm² on those genes | Loading norm in top 5 / top 50 PCs | Loading norm in null space of Y | Share of η variance from PCs 51–271 | Largest \|E[l_ik]\| | Largest correlation between columns of E[L] |
|---|---|---|---|---|---|---|---|
| w_S = 1 (pilot) | 0 | 0 | 94–100% / 100% | 0% | – | 925 | 0.03 |
| 5% target | 0 | 0 | 80–100% / 98–100% | 0–1% | – | 942 | 0.06 |
| 25% target | 6 | 85% | 0% / 2–3% | 47–86% | 21–46% | 6 × 10⁷ | 0.82 |
| 50% target | 42 | 94% | 0% / 2–3% | 54–83% | 30–63% | 2 × 10⁸ | 0.90 |

Ranges are over the five programs. The centered expression matrix has rank 271, and its
null space holds gene combinations with zero variance across patients. The η-variance
column was computed in the independent review.

**Sequence in the traces (25% target).** The program scores YF keep SD 1 throughout,
because canonical rescaling holds them there, so the score scale is not what runs away.
Separation comes first. By sweeps 15 to 30 the largest linear predictor is already 15 to
20 and β is 3 to 4 per SD. Reconstruction error rises steadily from the first sweep (RMSE
0.75 at sweep 1, 1.9 by sweep 25, 3.5 by sweep 50). The loading norms grow gradually, and
the largest jump comes between sweeps 100 and 200. The monitored objective rises until
about sweep 43 and then mostly declines: 75% of later steps decrease it at the 25% target
and all of them at the 50% target.

**The weighted likelihood prefers the degenerate fit.** From the pilot to the 25%-target
fit, the Cox partial log likelihood improves from −647 to −245 (403 nats). The expression
log likelihood with τ_j profiled falls from about −108,400 to −153,000 (about 44,600 nats).
The weighted sum favors the degenerate fit whenever w_S exceeds about 44,600 / 403 ≈ 110.
This ignores the KL terms and compares two points, so it is not a proof that the optimum
is degenerate. It does explain why weights near 100 and above failed in every study.
With τ_j estimated per gene, a poorly reconstructed gene costs only (n/2) log of its
residual sum of squares, so the penalty for abandoning a few genes is small.

**Interpretation.** With p > n, a large weight makes fitting training survival through
low-variance row-space directions worth more than describing expression. The noise
precision then collapses on the few genes those directions use, after which nothing but
the prior constrains them, and the loading mass in the null space of Y grows. The τ
collapse and the null-space mass look like amplifiers and consequences, not the trigger.
The objective's decline after its peak means the updates are no longer ascent steps.
Candidate causes are the quadratic surrogate for the Cox term, the canonical rescaling,
the empirical-Bayes prior refits and the mean-field β update. The October 30 diagnostics
found the mean-field β information weighting to be a major cause of instability.

**Checks that would settle the order of causes:** evaluate the full objective (with KL
terms) at the pilot and degenerate solutions under the same w_S; rerun the 25% target with
τ held at the pilot values and see whether η and β still blow up through low-variance PCs;
log τ_j, loading norms, the largest η and the Cox log likelihood at every sweep.

**What this rules out.** The working explanation recorded on 2026-11-14, that the score
scale grows through F, is not supported. A prior on β alone cannot help, because β is
already per SD and the separation runs through F. A geometric weight ramp does not change
which solution the weighted likelihood prefers.

## 3. How DeSurv handles supervision

Read from the package source (`src/functions.cpp`, github.com/rashidlab/DeSurv, version
1.0.1 installed locally) and the paper.

- **Loss.** (1 − α)[‖X − WH‖² / (2‖X‖²) + λ_H pen(H)] − α[2 ℓ_Cox / n_events − λ pen(β)]
  + λ_W pen(W). Both terms are O(1), and the β elastic net acts on standardized scores.
- **Stated weight.** Multiplying through by ‖X‖² / ((1 − α)σ²) turns the reconstruction
  term into a Gaussian log likelihood with common variance σ² and gives the Cox term a
  weight of [α / (1 − α)] · 2‖X‖² / (σ² n_events). ‖X‖² is the uncentered sum of squares of
  nonnegative rank data, which inflates the ratio, so the stated weight is plausibly
  10⁴ or more for the selected α = 0.334. This is far above where our fits break.
- **But the stated weight does not set the steps.** The Cox gradient on W is rescaled so
  that its Frobenius norm is α times the reconstruction gradient's (the ratio is capped at
  10⁶). The direction is therefore not the gradient of the stated loss, and the comparable
  quantity is a gradient share: the survival push is about one third of the
  reconstruction push. The stated loss is enforced only through backtracking.
- **Bounded, monotone steps.** Multiplicative updates with the ratio clamped to [0.2, 5],
  damping θ = 0.5 that shrinks by 0.9 for up to 30 backtracks, and acceptance only when the
  full loss does not increase. If no step is accepted, W stays where it was. Convergence is
  declared when the cosine change of W and H falls below the tolerance, which a stalled W
  also satisfies. The β update backtracks the same way.
- **Fixed common noise level.** The reconstruction loss has no per-gene precision, so the
  τ collapse in Section 2 cannot happen.
- **Scale and sign constraints.** W is nonnegative with unit-norm columns, H and β are
  rescaled so η is unchanged, and the input is nonnegative (within-patient ranks).
  Nonnegative programs cannot form gene contrasts, which removes most low-variance
  directions that could separate survival.
- **Search.** 20 random initializations, warm starts across α, and Bayesian optimization of
  (k, α, λ) by CV C-index. The paper reports k = 3 and α = 0.334 for PDAC.

DeSurv's results may therefore owe as much to implicit regularization (bounded steps, warm
starts, nonnegativity, stalls counted as convergence) as to the weight α.

## 4. Other literature

- **Prediction-constrained training** (Hughes et al., AISTATS 2018). Raising the response
  likelihood to a power in a joint model ("label replication") can fail at any power,
  because it does not control how latent variables are inferred from the data alone. YFB
  already has the prediction-constrained form: η depends only on the observed Y and the
  global F. LB (η = Lβ) is the replication form. They also found that continuation from
  small weights fails and that restarts are needed.
- **Generalized Bayes and modular inference** (Bissiri, Holmes and Walker 2016; Holmes and
  Walker 2017; Grünwald and van Ommen; Jacob et al. 2017; Carmona and Nicholls 2020). These
  calibrate a power for a misspecified module, and the calibrations down-weight. Here the
  clearly misspecified block is the Gaussian genomics likelihood with independent genes,
  which argues for tempering genomics rather than inflating the Cox term. In a mean-field
  fit the two differ in how the KL terms and empirical-Bayes priors respond.
- **Partial factor model** (Hahn, Carvalho and Mukherjee 2013). Lets the outcome load on
  low-variance residual directions through separate shrunken coefficients, so the factors
  that describe the predictors stay unsupervised.
- **Supervised factor models with a response row** (sJIVE, Palzer et al. 2021; SupSVD, Li et
  al. 2016). The outcome enters as an extra weighted row of the data matrix, with a
  unit-norm joint loading so the outcome weight cannot inflate the coefficients.
- **Supervised PCA** (Bair and Tibshirani 2004). Screen genes on univariate Cox scores with
  a threshold chosen by CV, then factorize the retained genes.
- **Firth-penalized Cox** (Heinze and Schemper 2001). Keeps β finite under separation, a
  symptom here.
- **CoxNMF** (Huang et al. 2020). An unnormalized weight up to 1,000 with no scale
  constraint, reporting the iterate with the best C. An example of what to avoid.

## 5. Candidate solutions, ranked

1. **Restrict or shrink the low-variance directions available to survival.** Let survival
   reweight genes only within the top r right singular vectors of the training expression
   (r < rank(Y), fixed in advance or chosen by training CV), or shrink the survival
   contribution along low-variance directions in proportion to their variance. Restricting
   to r < rank(Y) is what removes an exact fit to the training survival times when p > n.
   This targets the trigger in Section 2. Risk: a weak program that lies below the top r
   directions cannot be reached, which bears directly on the near-threshold case.
2. **Make every update increase the objective.** Damp the updates of q(F) and q(β) and
   accept a step only if the monitored objective does not decrease, as DeSurv does. This is
   a precondition for evaluating any weight, because the traces show the current updates
   stop being ascent steps. It also exposes whether the decline comes from the Cox
   surrogate, the rescaling, the prior refits or the β update. It does not change where the
   optimum is.
3. **Guard the noise precision.** Hold τ_j at the unsupervised fit's values, or put a
   shared prior on the τ_j that bounds them away from zero (the current `update_tau` caps τ
   from above at 10⁸ but has no lower floor). This removes the amplifier in Section 2. The
   October evidence is mixed: fixed τ bounded the reconstruction but none of three fits
   converged or recovered the program, while the guarded-update arm did both.
4. **Separate survival-driven programs.** Add one or two programs that enter the Cox model
   with strong shrinkage and leave the unsupervised programs alone, in the spirit of the
   partial factor model. Survival then has a channel that does not bend the expression
   programs. Risk: those programs may be pure noise unless shrinkage and held-out tuning
   control them.
5. **Share-targeted weight or a DeSurv-style gradient share, on top of 1 to 3.** The
   2026-11-13 study failed at targets of 5% and above without these guards. It becomes
   worth repeating once they hold.
6. **Temper the genomics block instead of inflating survival.** A power below 1 on the
   Gaussian genomics likelihood with explicitly weighted KL terms (an open ROADMAP item),
   the reading the generalized-Bayes literature supports. Risk: in the multimodal weight
   path a tempered likelihood made empirical-Bayes priors empty the loadings.
7. **A Gaussian survival response row.** Append a weighted working response (martingale or
   deviance residuals, refreshed each sweep) as an extra feature, as in sJIVE and SupSVD.
   Conjugate and bounded, but an approximation to the Cox likelihood with its own weight.
8. **Survival gene screening inside CV.** Restrict to genes with the strongest univariate
   Cox association, recomputed within each training fold. This raises the survival share
   mechanically, but brings outcome-based preprocessing that must stay inside the folds.

**Selection of the weight.** Whatever candidate is chosen, w_S (or r, or the share) should
be tuned by held-out training-CV partial likelihood, as the derivation already proposes.
That rule picks w_S between 1 and 10 in simulation and 1 on the real data under the
current model.

**Technical note.** With ρ much larger than 1, the per-gene mean-field factor for q(f_k)
ignores strong correlation among genes induced by the survival term and understates
variance. A structured q(f_k) with covariance in the n-dimensional row space (cheap via
the Woodbury identity) would be more accurate if large weights are pursued.

Not sufficient alone: a prior or Firth penalty on β (a reasonable companion to 1). Not
recommended: continuation from small weights without restarts (Hughes et al. found it
unreliable) and fixed large weights without guards (CoxNMF).

## 6. How to evaluate a candidate

The 2026-11-19 near-threshold simulation at K = 5 reproduces the real-data behavior. Its
ceilings limit what it can show: no arm recovers the 0.05% program in the real-data-like
setting, and the two-step fit already recovers the 0.1% program in 8 of 10 datasets.
An evaluation should therefore:

- use more seeds (at least 30 per cell), paired comparisons against the two-step fit and
  prespecified thresholds;
- count success as better recovery and better or equal held-out C and calibration slope;
  report loading correlation with the two-step fit only descriptively, since a candidate
  that recovers better need not diverge, and divergence alone can be overfitting
  (Section 7);
- include a positive control in the strong setting, where the current model already gains;
- check permuted survival for both false-active programs and held-out C;
- tune the weight by inner CV, as it would be on the real data;
- require a nondecreasing objective as part of convergence, with no gene's τ_j collapsing
  when τ is estimated;
- save the fits so τ and the loadings can be inspected afterwards;
- add a supervised PCA or principal-component Cox comparator to show how much any
  method can gain in that setting.

Only a candidate that passes there moves to training CV on the real data, with the
external cohorts scored once after the model is fixed.

### Measures of success beyond CV C

A large gain in CV C is not expected (DeSurv's CV C was 0.655, close to the current
model's), so these measures carry the evaluation alongside recovery. A factorization
that differs from EBMF counts only if it moves toward survival-relevant structure, since
the spike in Section 7 shows that divergence alone can be overfitting.

In simulation, where the truth is known:

- **Loading accuracy as a continuous measure.** The correlation of the matched program
  with the true prognostic loadings, not only the pass or fail recovery threshold.
- **Gene-level recovery.** Precision and recall of the true program's genes among a
  fitted program's top genes. DeSurv reported 0.50 against 0.07 for NMF.
- **Specificity.** Under permuted survival, no spurious survival-active programs and no
  gain in held-out C.
- **Rank needed.** The smallest K at which the prognostic program appears. Supervision
  should find it at a smaller K than the two-step fit.

On the real data, where there is no truth:

- **Parsimony.** Similar C with fewer programs. DeSurv's main argument was that it needed
  3 factors where NMF needed 7 for the same per-cohort C.
- **Concentration.** Survival signal concentrated in one or two programs, measured by each
  program's share of the risk score or its drop in partial log likelihood (DeSurv's
  leading factor had Δℓ = 66.4).
- **Low-variance prognostic programs.** A program that is clearly prognostic while
  explaining little expression variance. DeSurv's main factor explained 7% of the
  variance. The current Program 5 explains 3.6%.
- **Localized change.** Supervision changes the survival-relevant programs while the
  others stay close to EBMF. Every program changing a little would look like noise.
- **Stability.** The same program appears across CV folds and random starts.
- **External replication per program.** Each program's hazard ratio has the same direction
  and similar size in the five external cohorts, not only the combined risk score.
- **Added value over known subtypes.** The program stays prognostic after adjusting for
  PurIST and DeCAF. DeSurv reported a hazard ratio of 1.33 after that adjustment.
- **Biological coherence.** Enrichment for known PDAC biology (basal and classical,
  stromal subtypes, immune).
- **Held-out partial log likelihood.** More sensitive than C, so it can show small real
  gains that C hides.

Several of these already exist in the repository: the program rule with its risk-score
share and fold stability, the external per-program curves and the PurIST comparisons.

## 7. Weight-grid spike in the test bed

A spike on 2026-11-20 (script and per-fit metrics in
`results/benchmark_sim/outputs/near_threshold_weight_spike_2026-11-20/`) refit the supervised
model at w_S of 10, 100 and 1,000 in the real-data-like setting and the moderate setting
(named `realistic` in the CSV), for the 0.05% and 0.1% shares, with the same 10 seeds,
K = 5 and a 300-sweep cap. It used the committed data generator and scoring and compared
each fit with the saved two-step fit of the same dataset. All 120 fits ran without
errors. Means over 10 datasets, with w_S = 1 from the committed run:

| Setting, share | w_S | Recovered (of 10) | Converged (of 10) | C gain over two-step | Median loading r with two-step | Mean largest \|β\| per SD |
|---|---|---|---|---|---|---|
| Real-data-like, 0.05% | 1 / 10 / 100 / 1,000 | 0 / 0 / 0 / 0 | 10 / 10 / 0 / 0 | +0.017 / +0.026 / +0.025 / −0.001 | 0.98 / 0.77 / 0.72 / 0.63 | – / 16 / 21 / 21 |
| Real-data-like, 0.1% | 1 / 10 / 100 / 1,000 | 8 / 0 / 0 / 0 | 10 / 10 / 0 / 0 | +0.005 / −0.082 / −0.077 / −0.079 | 1.00 / 0.78 / 0.75 / 0.66 | – / 19 / 20 / 21 |
| Moderate, 0.05% | 1 / 10 / 100 / 1,000 | 0 / 0 / 0 / 0 | 10 / 10 / 0 / 0 | +0.030 / +0.036 / +0.029 / +0.007 | 0.98 / 0.76 / 0.65 / 0.68 | – / 17 / 21 / 22 |
| Moderate, 0.1% | 1 / 10 / 100 / 1,000 | 9 / 0 / 0 / 0 | 10 / 10 / 0 / 0 | +0.008 / −0.081 / −0.076 / −0.087 | 1.00 / 0.80 / 0.72 / 0.63 | – / 19 / 21 / 21 |

The two-step fit also recovers the 0.1% program in 8 of 10 datasets in both settings, so
w_S = 1 does not add it. Larger weights destroy a program the unsupervised factorization
already finds. A weight of 10 converges and moves the loadings, but away from the
prognostic program, and validation C falls about 0.08 below the two-step fit.
Coefficients of 16 to 21 per SD, against a true risk score SD of 0.8, show ordinary p ≫ n
overfitting of the predictive direction. Weights of 100 and above never converge, in line
with the break-even weight of about 110 in Section 2. The fits were not saved, so τ was not
inspected.

## 8. Decisions so far (Andrew, 2026-11-23)

- **What counts as success.** Supervision succeeds if the model identifies
  survival-associated programs on the real data and gives a survival-informed
  factorization that differs from the unsupervised EBMF fit. The simulation analog is
  recovering a low-variance prognostic program at a realistic survival signal, which the
  current model does not do (Section 7). A gain in CV C is not expected to be large and is
  not the measure of progress, though C should not degrade.
- **Updates that always increase the objective** (candidate 2) make sense and should be
  part of any approach.
- **Status.** Tabled for now. Next, discuss any further ideas, then test the candidates
  systematically, each on its own branch as needed, in the real-data-like test bed.

## 9. Open questions for Andrew

- Is a restriction of survival to the high-variance expression directions (candidate 1)
  acceptable, given that it limits which weak programs can be found?
- Should survival-driven programs (candidate 4) be separate from the expression programs,
  or must supervision change the expression programs themselves?
- Is holding τ at unsupervised values acceptable for the reported model, or should τ stay
  estimated with a shrinkage prior?
