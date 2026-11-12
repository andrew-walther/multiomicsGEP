# Matched LB and YFB comparison at K = 5

Date: 2026-11-12. Author: Andrew Walther.

## Result

All 104 paired simulation cells and the 12 training fits finished. No captured
fit error occurred. All training fits and all 160 S-9/S-14 fits converged.
Seven matched-generator fits reached the iteration cap. The simulation process
also reported at least 50 warnings whose messages were not persisted. This is
a material numerical-diagnostics limitation, separate from the saved convergence flags.

The registered relative rule excludes LB because its mean false activity is
above 1 in all three S-9 cells at w_S = 1, where YFB passes. Neither model
passes the near-threshold S-14 cell at either weight. YFB is not excluded by
this relative comparison, but the S-14 failures and unclassified warnings
remain limits on the scientific interpretation.

On the training patients, the fitted predictors correlate 0.9978 and corresponding
programs have identical top-20 gene sets. Mean within-cohort CV concordance is
0.6516 for LB fold-in and 0.6546 for YFB projection, with paired difference
0.0029 and SE 0.0044. The donor-bootstrap external interval for the pair-weighted
pooled difference includes zero.

## Scope and registered choices

This study implements Part 1 of the [follow-up prompt](../prompts/2026-11-12-codex-d-e-followup-prompt.md).
The fixed inputs are D4 genes, per-cohort training-fold z-scores, point-Laplace
L and F, per-program Normal β, K = 5 and pooled fits without Cox strata.
Training uses w_S = 1. Simulation uses w_S = 1 and 10. No pruning is used.
The maximum iteration count and tolerance come from `config/globals.yml`.
The saved S-11 folds in `data/lf_prior_real_2026-11-12.rds` are retained.

## Pre-fit software corrections

`code/lb_yfb_study.R` returns within-cohort, pooled and individual-cohort
concordance, with each comparable-pair count. `survival::concordancefit()` uses
`reverse = TRUE` and `timefix = TRUE`. Larger scores always mean higher hazard.
Counts sum concordant, discordant and tied-score pairs, including every stratum.
Constant scores remain undefined under this study's reporting convention, with
pair counts retained. Near-equal event times are treated as ties, matching task A's
current concordance calculation. No existing test required keeping `timefix = FALSE`.
The hand-computed unequal-cohort case has six concordant pairs in one cohort
and one discordant pair in the other, giving within-cohort C = 6/7.

The runner accepts a positive integer `--K` and a positive finite `--w_S`.
Rank defaults to the historical configuration value in every mode. Training
weight defaults to 1. The scientific commands explicitly set K = 5. Simulation
retains its registered weights 1 and 10. Actual settings are written once to
`settings.txt`. The smoke test checks the saved models as well as this file.

When a full posterior start is supplied, `code/fit_ssbmf_yfb.R` now saves its
actual moments, precisions and persistent priors immediately before the first
sweep. Ordinary fit objects are unchanged. The paired-study assertion compares
this snapshot with the common start using exact identity. A mutation test changes
the recorded β start and requires the assertion to reject it. Andrew explicitly
authorized this small addition to the fitter shared with Part 2 before E's fits.

Failed folds and cells retain aggregate rows with explicit error messages and
undefined metrics. An incomplete five-fold comparison has an undefined mean and
paired SE. Patient-level fits and scores remain under ignored `data/`.

## Verification and independent review

The two focused files pass 10/10 and 5/5 checks. Independent source review found
one mismatch with the specification: smoke mode still defaulted to rank 2.
That exception was removed so all modes use the configuration default unless
`--K` is supplied. Review found no other defect in concordance, pair counts,
initialization checks, weight propagation or failure schemas. The first full-suite run picked up Part 2's temporarily registered feature-group tests
and reported 821 passed and seven failed. All seven failures were the intentionally
unimplemented feature-group checks. Andrew authorized a temporary copy of the committed
manifest at `19f910bd`, with current Part 1 code and tests. That full suite passed
819 tests with zero failures. It emitted the existing bounded β optimizer
profile-likelihood fallback warning. These checks passed before scientific fitting started.

The external checks pass 3/3. The summary checks pass 4/4, including an actual
CSV round trip. That test caught R's conversion of an all-blank error column
to logical NA. The summarizer now reads error messages as character, so
successful seeds are not counted as failures. Independent source review
approved this correction. Neither the fits nor the gate thresholds changed.
The final combined run used the same approved committed manifest plus the
external and summary test files, and passed 826 tests with zero failures.
The existing bounded β optimizer fallback warning remained the only diagnostic
warning reported by the suite.

## Simulation

The [per-seed metrics](../../results/benchmark_sim/outputs/lb_yfb_simulation_2026-11-12/simulation_metrics.csv)
and [cell summaries](../../results/benchmark_sim/outputs/lb_yfb_simulation_2026-11-12/simulation_summary.csv)
retain all 104 paired cells, or 208 model fits. There were no captured hard
errors. Molecular loading recovery, coefficient activity and the sign check
follow the registered S-9/S-14 definition. These generators simulate survival
from a centered molecular projection under the planted loading. The matched
generator instead uses latent Lβ.

Each S-9/S-14 row below has all ten expected seeds. False-activity means use
ten defined values. Validation-C means use only the defined values, with that
count shown beside the mean. No undefined score is replaced by 0.5. The
binary gate always uses all ten expected seeds, including recovery failures.
The simulated validation set has one cohort, so its individual-cohort,
within-cohort and pooled concordance agree. Comparable-pair counts are saved
in the per-seed table.

| Design | Planted share | w_S | Model | Recovery | False-active mean | Converged | Validation C [defined n] | Gate |
|---|---|---|---|---|---|---|---|---|
| S9 | 0.25% | 1 | LB | 10/10 | 1.8 | 10/10 | 0.9376 [10] | Fail |
| S9 | 0.25% | 1 | YFB | 10/10 | 0.9 | 10/10 | 0.9435 [10] | Pass |
| S9 | 0.25% | 10 | LB | 9/10 | 0.2 | 10/10 | 0.9349 [10] | Pass |
| S9 | 0.25% | 10 | YFB | 9/10 | 0.9 | 10/10 | 0.9211 [10] | Pass |
| S9 | 1% | 1 | LB | 10/10 | 1.5 | 10/10 | 0.9671 [10] | Fail |
| S9 | 1% | 1 | YFB | 10/10 | 0.1 | 10/10 | 0.9724 [10] | Pass |
| S9 | 1% | 10 | LB | 10/10 | 0.0 | 10/10 | 0.9661 [10] | Pass |
| S9 | 1% | 10 | YFB | 10/10 | 0.0 | 10/10 | 0.9693 [10] | Pass |
| S9 | 5% | 1 | LB | 10/10 | 1.2 | 10/10 | 0.9793 [10] | Fail |
| S9 | 5% | 1 | YFB | 10/10 | 0.0 | 10/10 | 0.9856 [10] | Pass |
| S9 | 5% | 10 | LB | 10/10 | 0.0 | 10/10 | 0.9783 [10] | Pass |
| S9 | 5% | 10 | YFB | 10/10 | 0.0 | 10/10 | 0.9859 [10] | Pass |
| S14_strong | 0.05% | 1 | LB | 5/10 | 1.4 | 10/10 | 0.8022 [6] | Fail |
| S14_strong | 0.05% | 1 | YFB | 4/10 | 2.9 | 10/10 | 0.7488 [10] | Fail |
| S14_strong | 0.05% | 10 | LB | 7/10 | 0.1 | 10/10 | 0.8498 [8] | Fail |
| S14_strong | 0.05% | 10 | YFB | 3/10 | 1.7 | 10/10 | 0.8121 [10] | Fail |

At w_S = 1, LB fails the false-activity limit at shares 0.25%, 1% and 5%.
YFB passes those three cells. Both models pass all three S-9 cells at w_S = 10.
In S-14, LB recovers 5/10 programs at weight 1 and 7/10 at weight 10. YFB
recovers 4/10 and 3/10. Both therefore fail the recovery limit there, and YFB
also exceeds the false-activity limit at both weights. Applying only the
registered relative rule removes LB from carry-forward. The rule does not
remove YFB because no cell has LB passing while YFB fails. This does not turn
the S-14 failures into a successful acceptance result.

### Matched-generator diagnostics

No binary recovery gate was registered for these six-seed cells. The score
correlation is the mean absolute, one-to-one matched correlation for the three
true survival programs, evaluated on validation patients. LB uses fold-in
scores and YFB uses projection scores. All correlation means use six defined
seeds. Means include saved final states from unconverged fits, whose counts
are shown explicitly. These values are descriptive.

| Amplitude | w_S | Model | Matched score correlation | False-active mean | Converged | Validation C [defined n] |
|---|---|---|---|---|---|---|
| 0.5 | 1 | LB | 0.8006 | 0.83 | 6/6 | 0.7795 [6] |
| 0.5 | 1 | YFB | 0.8430 | 1.17 | 6/6 | 0.7830 [6] |
| 0.5 | 10 | LB | 0.8226 | 0.17 | 1/6 | 0.7318 [6] |
| 0.5 | 10 | YFB | 0.8491 | 1.00 | 5/6 | 0.7627 [6] |
| 1 | 1 | LB | 0.6992 | 1.17 | 6/6 | 0.7521 [6] |
| 1 | 1 | YFB | 0.7937 | 1.33 | 6/6 | 0.7867 [6] |
| 1 | 10 | LB | 0.6504 | 1.17 | 6/6 | 0.7910 [5] |
| 1 | 10 | YFB | 0.8120 | 1.50 | 5/6 | 0.7640 [6] |

### Numerical diagnostics and undefined scores

All 208 validation fold-ins or projections completed, and every LB fold-in
converged. Of the 208 model fits, 201 converged. The seven iteration-cap
failures were all at w_S = 10 in the matched generator:

- LB, amplitude 0.5: seeds 1, 2, 3, 4 and 6.
- YFB: seed 5 at amplitude 0.5 and seed 5 at amplitude 1.

The seven undefined validation concordances have positive comparable-pair
counts but constant risk scores. They are LB S-14 seeds 1, 4, 7 and 10 at
weight 1, LB S-14 seeds 1 and 8 at weight 10, and LB matched amplitude 1,
seed 5 at weight 10. Their undefined C values remain in the tables. In
particular, S-14 LB means summarize fewer seeds than the YFB means and
should not be interpreted as a paired prediction advantage.

The private run log ends with `There were 50 or more warnings (use warnings()
to see the first 50)`. The runner did not persist individual warning messages,
and the R process has exited. Their types and exact frequency cannot be
recovered from the saved output. The seven iteration caps above are established
from saved fit histories. They do not explain or classify the remaining
warnings. No claim is made that these warnings were benign, or that a recorded
CAVI convergence flag rules out an internal optimizer diagnostic. Warning
provenance remains a material limitation of this run. No scientific cell was
rerun to change its outcome or obtain a more favorable fit.

## Training cross-validation

All five folds converged for each model. Every held-out LB fold-in converged.
The fits are pooled without Cox strata. Within-cohort scoring counts only
TCGA–TCGA and CPTAC–CPTAC pairs. TCGA is labeled `TCGA_PAAD` in the saved input.

| Fold | Model | Within-cohort C | TCGA C | TCGA pairs | CPTAC C | CPTAC pairs | Pooled C |
|---|---|---|---|---|---|---|---|
| 1 | LB | 0.666 | 0.655 | 226 | 0.679 | 187 | 0.664 |
| 1 | YFB | 0.671 | 0.664 | 226 | 0.679 | 187 | 0.671 |
| 2 | LB | 0.737 | 0.764 | 216 | 0.704 | 179 | 0.741 |
| 2 | YFB | 0.739 | 0.759 | 216 | 0.715 | 179 | 0.744 |
| 3 | LB | 0.613 | 0.676 | 259 | 0.529 | 191 | 0.621 |
| 3 | YFB | 0.600 | 0.660 | 259 | 0.518 | 191 | 0.612 |
| 4 | LB | 0.667 | 0.638 | 224 | 0.701 | 197 | 0.677 |
| 4 | YFB | 0.679 | 0.643 | 224 | 0.721 | 197 | 0.681 |
| 5 | LB | 0.575 | 0.451 | 193 | 0.736 | 148 | 0.587 |
| 5 | YFB | 0.584 | 0.451 | 193 | 0.757 | 148 | 0.590 |

Five-fold means and paired differences, with SE = SD of the five differences / √5:

| Metric | LB mean | YFB mean | YFB − LB | Paired SE |
|---|---|---|---|---|
| C_within | 0.6516 | 0.6546 | 0.0029 | 0.0044 |
| C_pooled | 0.6580 | 0.6595 | 0.0015 | 0.0029 |
| C_TCGA_PAAD | 0.6367 | 0.6354 | -0.0014 | 0.0042 |
| C_CPTAC | 0.6698 | 0.6780 | 0.0083 | 0.0060 |

Every mean uses all five defined folds. The primary paired difference is 0.0029
with SE 0.0044. These concordances refer to held-out LB fold-in risks and YFB
projection risks, not fitted Lβ for the held-out patients.

All five fitted programs have nonzero loadings in every training fold and
in each full-data fit. Survival activity requires a positive posterior SD and
|E[β_k]| ≥ 1.96 SD(β_k). Counts are:

| Fold | LB nonzero loadings | LB survival-active | YFB nonzero loadings | YFB survival-active |
|---|---|---|---|---|
| 1 | 5 | 1 | 5 | 2 |
| 2 | 5 | 2 | 5 | 2 |
| 3 | 5 | 2 | 5 | 2 |
| 4 | 5 | 3 | 5 | 3 |
| 5 | 5 | 2 | 5 | 4 |
| full | 5 | 3 | 5 | 3 |

The mean number of survival-active programs across folds is 2.0 for LB and
2.6 for YFB. These counts describe the fitted posteriors under the fixed rule,
while held-out concordance evaluates their predictive scores.

## Full-data score comparison

Both full-data fits converged after 52 sweeps. The fitted linear predictors
correlate 0.9977605: LB uses fitted Lβ and YFB uses its normalized YF projection
multiplied by β. This fitted-L comparison is distinct from the fold-in risks
used for held-out LB evaluation.

The correlation between fitted L_k and its own model's projection (YF)_k is:

| Program | LB | YFB |
|---|---|---|
| 1 | 0.9992 | 0.9992 |
| 2 | 0.9954 | 0.9955 |
| 3 | 0.9961 | 0.9961 |
| 4 | 0.9848 | 0.9849 |
| 5 | 0.9946 | 0.9947 |

These high within-fit correlations describe this dataset. They do not make the
two score definitions mathematically equivalent.

### Gene content and cross-model agreement

The full [top-20 table](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/top_genes.csv)
contains each model's rank and signed posterior mean loading. Both models have
five nonzero programs and the same 20-gene set in each corresponding program.
The lists below use LB's magnitude ordering. YFB's full ranked lists are in the CSV.

| Program | Top 20 genes in both models |
|---|---|
| 1 | MYO5A, ZEB2, CCDC80, TCF4, ATP8B2, BNC2, DLC1, FBN1, FERMT2, GFPT2, WIPF1, GAS7, IL1R1, NIN, FSTL1, CD93, SERPINF1, MEF2C, C1S, CYP1B1 |
| 2 | METTL7A, ATP2A3, ACACB, PDCD4, LOXL2, ACSS1, MMP14, AKAP7, LDHA, KIAA1324, P4HA1, LGALS2, ADH1B, C7, TMED6, SLC16A3, CLU, PGK1, GJB2, AQP1 |
| 3 | KIAA1211, SYTL2, PLS1, CORIN, TOX3, EDIL3, CTSE, SLC44A1, RASEF, CDC42BPA, DKK3, CYP2C18, SH3D19, HNF4A, CAPN5, CEACAM6, TSPAN8, RNF128, BCL2L15, COL10A1 |
| 4 | SERPINA5, CTRB2, CPB1, REG1B, CTRB1, CPA1, CELA3A, CEL, CPA2, PNLIPRP1, SLC30A2, PLA2G1B, CTRC, CRP, REG3G, CELA3B, REG1A, REG3A, PKHD1, FGA |
| 5 | UHMK1, ARHGAP5, RAD21, PRKDC, TAOK1, DDX21, RCN3, CHPF, MXRA8, LOXL1, ATP13A3, IQGAP1, EFEMP2, NPDC1, IL6ST, JUND, VPS13C, DDX3X, ZDHHC20, SPON2 |

Every same-index top-20 overlap is 20/20, with Jaccard index 1. Every off-diagonal
cross-model top-20 overlap is zero. All 25 comparisons are saved in
[the overlap matrix table](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/cross_model_top20_overlap.csv).
The [all-pairs score and loading correlations](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/cross_model_program_correlations.csv)
retain the full comparison without imposing a matching assignment.

Overlap with the published DeSurv program lists is identical in LB and YFB:

| Program in both models | DeSurv D1, classical tumor | DeSurv D2, stromal/immune | DeSurv D3, basal-like tumor |
|---|---|---|---|
| 1 | 1 | 10 | 0 |
| 2 | 6 | 0 | 2 |
| 3 | 3 | 3 | 2 |
| 4 | 2 | 1 | 1 |
| 5 | 2 | 2 | 0 |

These compare each fitted program's 20 genes with a published 270-gene list,
with both sets restricted to the 2,218-gene D4 background. The [complete DeSurv
overlap table](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/desurv_top20_overlap.csv)
also gives Jaccard indices and unadjusted hypergeometric p-values. These
content summaries do not select a model or determine the simulation gate.
Full-data programs 2, 4 and 5 are survival-active in both models. The [coefficient
table](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/program_coefficients.csv)
gives their posterior means, SDs and plug-in hazard ratios, along with the
remaining programs and all five folds.

## External scoring, descriptive

Each fixed full-data model was scored once in each external series, with no
refitting, tuning or outcome-based sign change. LB uses fold-in and YFB uses
the training-normalized molecular projection. The S-12 preprocessing was kept,
including within-series z-scores and zeros for missing D4 genes.

The intervals below use 2,000 bootstrap replicates, resampling records within
each series with seed 1. Paired differences use the same records in both arms.

| Series | Records | Events | Missing genes | LB C [95% CI] | YFB C [95% CI] | YFB − LB [95% CI] |
|---|---|---|---|---|---|---|
| Dijk | 90 | 81 | 30 | 0.5976 [0.5310, 0.6668] | 0.5989 [0.5303, 0.6669] | 0.0013 [-0.0104, 0.0132] |
| Moffitt_GEO_array | 123 | 83 | 71 | 0.5326 [0.4529, 0.6071] | 0.5375 [0.4576, 0.6080] | 0.0049 [-0.0055, 0.0155] |
| PACA_AU_array | 63 | 38 | 13 | 0.6718 [0.5796, 0.7622] | 0.6645 [0.5726, 0.7558] | -0.0073 [-0.0203, 0.0044] |
| PACA_AU_seq | 52 | 31 | 27 | 0.7088 [0.6011, 0.8119] | 0.7044 [0.5960, 0.8056] | -0.0045 [-0.0171, 0.0104] |
| Puleo_array | 288 | 181 | 82 | 0.6437 [0.5988, 0.6866] | 0.6445 [0.6002, 0.6871] | 0.0008 [-0.0031, 0.0049] |

Pooled intervals use the donor bootstrap from the S-12 audit, with donors
resampled within their assay-membership groups. All records from one donor
are kept together across both PACA-AU assays. The overlap check reproduced
47 shared donors. Per-series counts describe records, which overlap across
the two PACA-AU series. Larger scores retain the Cox hazard direction throughout.

| Pooling rule | Quantity | Estimate [95% donor-bootstrap CI] |
|---|---|---|
| Equal series | LB C | 0.6309 [0.5876, 0.6704] |
| Equal series | YFB C | 0.6299 [0.5874, 0.6701] |
| Equal series | YFB − LB | -0.0010 [-0.0061, 0.0044] |
| Within-series comparable-pair weighted | LB C | 0.6273 [0.5940, 0.6614] |
| Within-series comparable-pair weighted | YFB C | 0.6283 [0.5958, 0.6618] |
| Within-series comparable-pair weighted | YFB − LB | 0.0010 [-0.0023, 0.0045] |

The five series give nearly identical performance for the two parameterizations.
The paired pooled intervals contain zero under both weighting rules. This
external comparison changes no Task E carry-forward decision. Saved S-12
aggregate tables do not contain aligned patient scores, so the requested score
correlation with that historical K = 7 working model is unavailable. No S-12
model was reconstructed to obtain it. The present models use K = 5 and a common
SVD start, so they are not the S-12 fits.

## Interpretation and interfaces

LB specifies the Cox linear predictor as Lβ. L is the posterior latent activity
of a program, conditional on the molecular factor model and the other programs.
For a new patient, fold-in estimates this activity under the stored molecular
priors and precisions. It uses no new survival outcome. Conditional on the other
activities, exp(β_k) is the hazard ratio for a one-unit increase in L_k on the
fitted latent scale. That increase adds F[,k] to the expected standardized
molecular profile, holding the other activities fixed. Its unit depends on
the fitted loading scale.

YFB specifies the predictor through the molecular projection. Prediction uses
the training normalization and any stored projection center, then multiplies the
projected program scores by β. Conditional on the other projections, exp(β_k)
is the hazard ratio for a one-unit increase in this normalized projection score.
It is a weighted sum of the standardized gene levels. Canonicalization gives
each live YFB projection unit SD on the training data. Raw β magnitudes and
hazard ratios per unit therefore use different scales in the two models. It does not estimate the
latent activity by fold-in. In the factor model, YF ≈ L(F′F) plus noise, so shared
genes can mix the latent activities in the projection. The two models therefore
need not give the same program scores or coefficient units.

Point-Laplace loadings can be positive or negative. Gene lists are ranked by
absolute posterior mean loading, with signs retained in the aggregate table.
Lists contain up to 20 nonzero genes. A dead program has no top genes. The
DeSurv comparator uses the published 270-gene program lists restricted to the
D4 background, so those overlap counts compare sets of different sizes.
Program indices come from a common start but can change meaning during fitting.
All cross-model program pairs are reported without imposing a correspondence.

The study helper takes patients-by-genes matrices and aligned outcome, cohort
and fold vectors. It constructs one scaled matrix and one posterior start for
each pair, fits LB and YFB, then uses fold-in and projection respectively.
Its output includes private fits and scores as well as aggregate metrics.
The runner keeps the private fields under `data/`, writes only aggregates to
its output directory and refuses reused output paths. This prevents accidental
overwriting of a registered run. `survival` supplies concordance and `ebnm`
supplies the existing posterior updates. No new package is introduced.

The aggregate summarizer reads saved fits and per-seed tables. It counts
recovery successes, false-active programs and converged fits for each registered
cell, then applies only the relative exclusion rule in the prompt. It retains
defined-value counts beside descriptive means. Failed seeds cannot satisfy the
binary gate. For training fits it reports coefficients, posterior SDs, activity
counts, predictor correlations and gene overlaps. The existing
`top_n_genes_table()` and `compute_geneset_overlap()` helpers supply gene ranking
and background-restricted overlap calculations. These summaries do not refit
or select a model.

The external scorer reads the fixed full-data checkpoint and ordered training
genes. Each series is prepared as in S-12, aligned to those genes and scored
once. Missing genes are assigned zero after within-series standardization.
Per-series intervals use `code/concordance_ci.R` with explicit fixed Cox
orientation. Pooled intervals resample biological donors, retaining all rows
from each sampled donor across the PACA-AU assays. The equal-series mean and
the comparable-pair-weighted stratified C are reported separately because they
weight the five series differently. Neither external score is used for tuning.

## Limits

The w_S = 10 simulation cells use survival-block weighting, which differs from
the historical S-9/S-14 weighting when w_S differs from 1. Training-fold z-scores
are learned within each cohort, so this CV does not exactly reproduce S-11.
All fits in this comparison use K = 5. S-11 and task A results at K = 7 remain
historical references under their original rank labels. External scoring is
descriptive. Its performance does not determine the Task E carry-forward rule.

The registered S-9/S-14 generator already SD-scales training and validation
cohorts separately. The additional study scaling then uses training statistics
for both. This preserves the historical generator and should not be described
as a simulation pipeline in which every preprocessing statistic is learned
only from training rows. Real-data CV scaling is learned within each training
complement.

## Interpretation for the progress book

At K = 5 on the same D4 z-scores, pooled LB and YFB recover the same top-20 gene
sets in the training data and give closely correlated fitted predictors.
Held-out prediction is similar: LB estimates new latent activities by fold-in,
whereas YFB uses a weighted projection of the gene levels. The CV difference
is 0.0029 with SE 0.0044, and the donor-aware external interval includes zero.
Their coefficients refer to different score units: exp(β_k) is a hazard ratio
per unit latent activity in LB and per unit normalized molecular projection
in YFB. The registered simulation rule excludes LB because of false activity
at weight 1. Both models miss the near-threshold recovery criterion. The
iteration-cap failures and missing warning messages limit the numerical
interpretation, so the comparison does not establish recovery in every regime.

## Evidence and rebuild commands

The scientific fitting code was tested and pushed in `859f1401` before either
mode started. External scoring and its aggregates were pushed in `c320f763`,
with explicit frozen point-estimate evaluation in `1f267690`. The tested
summarizer is `76c83a53`. Source hashes in both fitting settings matched the
files on disk throughout fitting. No fitter or its recorded dependencies
changed while the modes ran. The only Part 2-owned code change was the small
initial-state snapshot explicitly authorized by Andrew before the runs.

- [Simulation settings](../../results/benchmark_sim/outputs/lb_yfb_simulation_2026-11-12/settings.txt),
  [registered design](../../results/benchmark_sim/outputs/lb_yfb_simulation_2026-11-12/design.csv)
  and [per-seed metrics](../../results/benchmark_sim/outputs/lb_yfb_simulation_2026-11-12/simulation_metrics.csv).
- [Training settings](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/settings.txt),
  [fold metrics](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/cv_metrics.csv),
  [paired summaries](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/cv_summary.csv)
  and [program activity](../../results/benchmark_sim/outputs/lb_yfb_training_2026-11-12/program_activity.csv).
- [External settings](../../results/benchmark_sim/outputs/lb_yfb_external_2026-11-12/settings.txt),
  [series concordance](../../results/benchmark_sim/outputs/lb_yfb_external_2026-11-12/external_cindex.csv),
  [paired differences](../../results/benchmark_sim/outputs/lb_yfb_external_2026-11-12/external_paired_diff.csv)
  and [donor-pooled intervals](../../results/benchmark_sim/outputs/lb_yfb_external_2026-11-12/external_pooled.csv).

Rebuild aggregate interpretation tables without refitting:

```sh
Rscript results/benchmark_sim/summarize_lb_yfb_study.R --date=2026-11-12
```

The full-data training checkpoint is
`data/lb_yfb_study_lb_yfb_training_2026-11-12/patient_level.rds`. Per-cell
simulation checkpoints are under
`data/lb_yfb_study_lb_yfb_simulation_2026-11-12/`. Input, fits, scores and process
logs remain under ignored `data/`. Only aggregate CSVs and settings are committed.
The private coordinating run log is `data/lb_yfb_part1_2026-11-12_run_log.json`.

## Independent final review

Final independent review checked the code against the registered plan and
recomputed the summaries from the raw aggregates. It also checked the saved
common starts, coefficients, posterior SDs, signed top-gene weights and all
cross-model overlap values against the private fits. Review found no numerical
or reporting discrepancy in those checks. It identified the unclassified
warnings as a material diagnostic limitation, rather than evidence that the
warnings were benign. A fresh final review approved the completed findings with this limitation explicit.

## Text for the shared status files

The K = 5 matched LB/YFB study finished 104 paired simulation cells and all
12 training fits under point-Laplace L and F, per-program Normal β, D4 genes
and pooled cohort-specific z-scores. The registered relative rule excludes
LB because of false activity in all three S-9 cells at weight 1, where YFB
passes. Neither model passes S-14 at either weight. All S-9/S-14 and training
fits converge. Seven matched-generator fits at weight 10 reach the iteration
cap. The simulation process reports at least 50 warnings whose messages were
not persisted, a material numerical-diagnostics limitation. Training CV C is
0.6516 for LB fold-in and 0.6546 for YFB projection, difference 0.0029 with
SE 0.0044. Fitted predictors correlate 0.9978, and the corresponding top-20
gene sets coincide. External donor-pooled pair-weighted C is 0.6273 for LB
and 0.6283 for YFB, difference 0.0010 with 95% interval −0.0023 to 0.0045.
External scoring is descriptive and changes no selection decision. Fits and
scores remain private. Aggregate evidence, code, tests and the full warning
and convergence limitations are in the dated LB/YFB findings note.
