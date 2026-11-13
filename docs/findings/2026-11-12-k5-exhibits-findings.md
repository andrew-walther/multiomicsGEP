# Reported-fit exhibits reproduced at K = 5

Date: 2026-11-12. Author: Andrew Walther.

## Result

Every exhibit that Chapter 4, the November 13 writeup and the November 13
progress-book chapter take from the reported K = 7 fit now has a K = 5
version built from the working model. The K = 5 fit is the one adopted on
2026-11-12: YFB, point-Laplace L and F, per-cohort z-scores, per-program
Normal β, pooled TCGA and CPTAC (273 patients, 139 deaths), w_S = 1 and the
2,218 genes selected by combined mean and variance rank (top 3,000 per
cohort, intersected). It was refitted through the same runner as the K = 7
fit, on the same folds, and is the same fit as the K = 5 fits of the
training-CV redo and of the LB and YFB comparison. Its loadings agree with
the training-CV redo's K = 5 full fit to 7 × 10⁻⁹, its coefficients agree
with the LB and YFB comparison to four decimals, and its external
concordance reproduces that comparison's YFB values exactly.

The main changes in interpretation are:

- Survival-active programs are 2, 4 and 5 (2, 3 and 5 at K = 7).
  Program 3 drops below the threshold (z 3.0 to 1.5). Program 4, with
  acinar genes at its positive pole and B-cell and interferon genes at its
  negative pole, rises above it (z 1.9 to 2.4).
- The estimated differences between the supervised and two-step fits are small. The training CV
  difference is +0.004 (SE 0.004) and the external pooled difference is
  −0.001 (donor interval −0.010 to 0.008). At K = 7 they were +0.011
  (SE 0.005) and +0.007 (−0.003 to 0.017).
- Only Program 2 meets the main-program rule (Programs 2 and 3 at K = 7).
- The training cutpoint moves from −0.046 to 0.618 training SD, so the
  high-risk group is smaller (128 of 616 validation records, against 297).
- The gene-count curve no longer favors 2,218 genes for the supervised fit.
  The curve is flat to about 8,000 genes and lower above, in both arms.
- All five K = 5 fold fits converge in both arms (3 of 5 at K = 7).
- The survival-weight refit selects w_S = 1. Nested CV gives C = 0.658
  with SE 0.027, so the working fit and its external scores remain unchanged.

## How the exhibits were rebuilt

The runners and builders read the reported fit through fixed paths. An
opt-in environment variable, `EXHIBIT_K`, now selects another rank
(`results/benchmark_sim/exhibit_rank.R`). With it unset, every script keeps
its historical rank, paths and file names. All five Chapter 4 builders were
rerun unset before committing. Their CSVs and LaTeX tables are
byte-identical to the committed ones, and every PDF matches after its
creation date is removed. With `EXHIBIT_K=5`:

- aggregate outputs go to `results/benchmark_sim/outputs/k5_working_2026-11-12/`,
  in subfolders named after the K = 7 source;
- patient-level files go to `data/k5_working_2026-11-12/` (ignored);
- figures and tables are written next to the K = 7 ones in
  `paper/prelim/figures/` and `paper/prelim/tables/` with a `_k5` suffix,
  and SVG conversions are in `outputs/k5_working_2026-11-12/figures_svg/`.

No K = 7 file was overwritten. The K = 5 builder runs skip the simulation
exhibits, which do not depend on the working fit, and skip the numeric
checks against the K = 7 findings. Descriptive program labels are not
assigned at K = 5, so the builders name the programs P1 to P5.

Rebuild everything from the repository root, in this order:

```sh
export EXHIBIT_K=5 OMP_NUM_THREADS=1
Rscript results/benchmark_sim/run_yfb_supervised_real.R d4 2      # 12 fits of the 2,218-gene set
Rscript results/benchmark_sim/run_yfb_supervised_real.R screen 2
Rscript results/benchmark_sim/run_yfb_supervised_real.R size 2
Rscript results/benchmark_sim/summarize_yfb_supervised_real.R
Rscript results/benchmark_sim/score_yfb_supervised_external.R
Rscript results/benchmark_sim/save_yfb_external_patient_scores.R
Rscript results/benchmark_sim/run_prompt6_comparators.R
Rscript results/benchmark_sim/run_prompt7_real.R
MC_CORES=1 Rscript results/benchmark_sim/run_prompt7_ebmf.R
Rscript results/benchmark_sim/run_yfb_supervised_real_refit.R fit 2
Rscript results/benchmark_sim/run_yfb_supervised_real_refit.R summary
Rscript paper/prelim/make_figures_supervised_sim.R
Rscript paper/prelim/make_figures_supervised_programs.R
Rscript paper/prelim/make_figures_supervised_external_km.R
Rscript paper/prelim/make_figure_fixed_cutpoint_km.R
Rscript paper/prelim/make_tables_comparators.R
unset EXHIBIT_K
Rscript results/benchmark_sim/k5_program_correspondence.R
for f in paper/prelim/figures/*_k5.pdf; do
  pdftocairo -svg "$f" results/benchmark_sim/outputs/k5_working_2026-11-12/figures_svg/$(basename "$f" .pdf).svg
done
```

The runners need `PDAC_DATA_ROOT` (or the default local path). The
prepared training matrix `data/yfb_supervised_real_2026-11-04.rds` and its
folds do not depend on K and are shared with the K = 7 run. The DeSurv
comparator is read from its saved fit in `data/prompt6_desurv/`. The
elastic-net, PurIST and EBMF rank-15 and rank-40 comparators are refitted
deterministically, and their external concordances equal the K = 7 run's
to machine precision. `run_prompt7_ebmf.R` now honors `MC_CORES` (default
unchanged), so it ran on one core.

## Exhibit ledger

Paths beginning `k5/` are under `results/benchmark_sim/outputs/k5_working_2026-11-12/`.
Other CSV folders are under `results/benchmark_sim/outputs/`. Files named
`tab_*.tex` are in `paper/prelim/tables/` and `sup_*.pdf` in
`paper/prelim/figures/`.

| # | Exhibit | K = 7 source | K = 5 file(s) |
|---|---|---|---|
| 1 | Training CV table (tab:training-cv) and loading precision text | `yfb_supervised_real_2026-11-04/{fits,summary,programs_D4}.csv`, `tab_training_cv.tex` | `k5/s11/{fits,summary,programs_D4}.csv`, `tab_training_cv_k5.tex` |
| 2 | Gene-count curve (fig:size-curve) and screening text | same folder, `{size_agreement,size_preferred,screen_paired}.csv`, `sup_real_size_curve.pdf` | `k5/s11/{summary,size_agreement,size_preferred,screen_paired}.csv`, `sup_real_size_curve_k5.pdf` |
| 3 | Weight refit text | `yfb_supervised_refit_2026-11-06/*` | `k5/s15/{weight_curve,nested_cv,verdict,programs_full_fit,metrics}.csv` |
| 4 | External table and figure (tab:external, fig:external) | `yfb_supervised_external_2026-11-05/*`, `tab_external.tex`, `sup_external_c.pdf` | `k5/s12/*`, `tab_external_k5.tex`, `sup_external_c_k5.pdf` |
| 5 | Fixed-cutpoint risk groups (fig:km-risk) and text | `prompt7_desurv_style/{cutpoints,group_hr,adjusted_hr,km_monthly,at_risk,enrichment}.csv`, `sup_external_km_fixed.pdf` | `k5/prompt7_desurv_style/` (same names), `sup_external_km_fixed_k5.pdf` |
| 6 | Comparator table (tab:comparators) | `prompt6_comparators/*`, `prompt7_desurv_style/rank_*.csv`, `tab_comparators.tex` | `k5/prompt6_comparators/*`, `k5/prompt7_desurv_style/rank_*.csv`, `tab_comparators_k5.tex` |
| 7 | Program classification (tab:program-rule) | `prompt7_desurv_style/{program_rule,fold_matches}.csv`, `tab_program_rule.tex` | `k5/prompt7_desurv_style/{program_rule,fold_matches}.csv`, `tab_program_rule_k5.tex` |
| 8 | Program table, enrichment table, loading heatmap, signature correlations | `prelim_ch4_2026-11-06/{programs_table,programs_enrichment,programs_subtype_concordance}.csv`, `tab_programs.tex`, `tab_enrichment.tex`, `sup_programs_heatmap.pdf` | `k5/prelim_ch4/` (same names), `tab_programs_k5.tex`, `tab_enrichment_k5.tex`, `sup_programs_heatmap_k5.pdf` |
| 9 | External program associations and per-program KM curves | `prelim_ch4_2026-11-06/external_program_direction.csv`, `tab_external_program_direction.tex`, `sup_external_km_program_{2,3,5}.pdf` | `k5/prelim_ch4/external_program_direction.csv`, `tab_external_program_direction_k5.tex`, `sup_external_km_program_{2,4,5}_k5.pdf`, `sup_external_km_risk_k5.pdf` |
| 10 | K = 5 to K = 7 program correspondence (new) | none | `k5/program_correspondence/{matches,score_correlations,k7_unmatched}.csv` |

`sup_programs_correspondence_k5.pdf` and the `programs_correspondence*.csv`
files in `k5/prelim_ch4/` compare the K = 5 loadings with the historical
α_F = 0 benchmark. Chapter 4 no longer reports that comparison, so they are
by-products only. Use item 10 instead.

Unchanged, because they do not depend on the working fit: all simulation
exhibits (fixed-weight, cross-validated-weight and near-threshold
simulations, the pruning comparison, the permuted-outcome null check, the
post hoc false-active contributions, the simulation checks of the program
rule and the mixed-signal simulation), the DeSurv, elastic-net and PurIST
values themselves, the prior, gene-screen, scaling and rank comparisons
(already reported at their own ranks), and the data counts.

The October YFB/EBMF verification figures are also unchanged. Their
builder, `docs/progress_book/figures/verify_october_yfb.R`, compares the
archived October `supervised_real_check` fit with the saved rank-7
`ebmf_cox_external` fit using historical preprocessing. It writes
`docs/progress_book/figures/2026-11-06/{external-concordance,program-correspondence}.png`.
These are historical diagnostics, and none of the three current document
sources references these images. Substituting the working rank-5 fit would
change the comparison being documented. The current rank-5 external
comparison is covered by items 4 and 6, and its program results by items
7 to 10. No rank-5 copy of the archived October comparison was made.

## Numbers at K = 5

### Training cross-validation (item 1)

Same five folds as the K = 7 table. Concordance is on all pairs within
each held-out fold, with z-scores computed once on all 273 patients. This
is a different estimand from the 0.6546 within-cohort concordance of the
LB and YFB comparison, which uses fold-specific scaling and counts only
same-cohort pairs. Report 0.6546 as the working model's CV concordance and
use the values below only in the supervised-versus-two-step table.

| Fit | Mean fold C (SE) | Pooled out-of-fold C | CV partial log-likelihood | Folds converged | Survival-active |
|---|---|---|---|---|---|
| w_S = 1 | 0.660 (0.027) | 0.654 | −773.6 | 5/5 | 3 |
| Two-step | 0.656 (0.029) | 0.651 | −773.6 | 5/5 | 3 |

Paired fold difference +0.004 (SE 0.004), against +0.011 (SE 0.005) at
K = 7. Every fold fit converged (largest sweep count 151 of 300). The
largest fold coefficient was 0.42 per SD of the score (0.46 at K = 7).
Survival supplied 4.9 × 10⁻⁷ to 1.2 × 10⁻⁶ of the loading precision of
survival-active programs (`share_active_median` in `k5/s11/fits.csv`,
median over active programs in each of the six 2,218-gene fits). Over all live
programs the range is 3.1 × 10⁻⁷ to 9.3 × 10⁻⁷. The supervised and
unsupervised loadings agree at |r| = 1.000 for every program
(`k5/prelim_ch4/programs_correspondence.csv`, rows "w_S=1 vs w_S=0"), as
at K = 7.

### Gene-count curve and screening (item 2)

| Genes | w_S = 1 C (SE) | Two-step C (SE) | Folds converged (each arm) | Summed fit minutes, 6 fits per arm (w_S = 1 / two-step) |
|---|---|---|---|---|
| 2,218 | 0.660 (0.027) | 0.656 (0.029) | 5/5 | 1.5 / 1.9 |
| 3,000 | 0.655 (0.031) | 0.660 (0.031) | 5/5 | 0.9 / 0.9 |
| 5,000 | 0.661 (0.029) | 0.664 (0.031) | 5/5 | 0.6 / 0.6 |
| 8,001 | 0.657 (0.032) | 0.659 (0.033) | 5/5 | 0.7 / 0.8 |
| 12,001 | 0.634 (0.034) | 0.638 (0.031) | 5/5 | 0.9 / 0.9 |
| 17,900 | 0.635 (0.035) | 0.633 (0.034) | 5/5 | 7.2 / 7.5 |

Fit minutes sum the elapsed times of six fits (full fit and five folds)
of one arm. They are neither processor CPU time nor the elapsed time of
the parallel run. All 72
size fits converged, and the whole S-11 rerun (the 2,218-gene set, sizes and screening,
144 fits) accumulated 0.62 summed fit hours on two workers. At K = 7 the supervised fit
fell from 0.663 to 0.643–0.648 above 5,000 genes while the two-step fit
stayed level, so the supervised advantage reversed. At K = 5 the two arms
track each other: both are level to 8,001 genes and 0.018 to 0.026 lower
at 12,001 and 17,900. The largest difference between sizes (0.027 for
w_S = 1, 0.031 for two-step) is about one fold SE (0.027 to 0.035). This
compares unpaired means. No paired size-versus-size difference was
computed. The registered
preferred-size rule (smallest set within one SE of the reference whose
active programs match, |r| ≥ 0.7) picks 8,001 genes for w_S = 1 and 12,001
for two-step. The 2,218-gene set fails the program-match part against the
all-genes fit: the matched active-program correlation is 0.67. The rule
was registered as a description, not a selection step, and the gene
screen remains the 2,218-gene set selected by combined mean and variance
rank, supported by the separate training-only screen comparison.
Screening genes for survival inside each fold gave
gains of +0.008 (SE 0.014), +0.014 (0.016) and +0.023 (0.016) over an
outcome-free set of 250, 500 and 1,000 genes for w_S = 1, each within
about 1.5 paired SEs, and no screened set beat the 2,218-gene fit (best
0.653). Screened fits converged in 7, 12 and 10 of 12 at 250, 500 and
1,000 genes, the rest at the 300-sweep cap.

### Weight refit (item 3)

All 216 checkpoints are present, with one unique fit for each combination
of six weights, six outer levels and six inner levels. No checkpoint
contains a fit error. The grid accumulated 4,383 seconds of fit time.
Of the 216 fits, 189 converged and 27 reached the 300-sweep cap. All six
full-data fits and all 30 folds used for the overall weight curve
converged. The 27 capped fits belong to inner CV within outer training
sets. No fit recorded a relative ELBO decrease.

| w_S | CV partial log-likelihood | Mean fold C (SE) | Folds converged | Median survival precision share among active programs |
|---|---|---|---|---|
| 1 | −773.551 | 0.660 (0.027) | 5/5 | 8.44 × 10⁻⁷ |
| 2 | −773.562 | 0.660 (0.027) | 5/5 | 1.70 × 10⁻⁶ |
| 5 | −773.598 | 0.659 (0.027) | 5/5 | 4.30 × 10⁻⁶ |
| 10 | −773.664 | 0.657 (0.027) | 5/5 | 9.53 × 10⁻⁶ |
| 20 | −773.925 | 0.656 (0.025) | 5/5 | 1.87 × 10⁻⁵ |
| 50 | −774.186 | 0.654 (0.026) | 5/5 | 5.80 × 10⁻⁵ |

The selected weight is 1, which maximizes the summed CV partial
log-likelihood. The range across weights is 0.635. The five nested-CV
outer training sets select weights 1, 1, 1, 10 and 50, respectively.
Their outer-fold concordances are 0.677, 0.742, 0.608, 0.674 and 0.589,
giving a mean of 0.658 and fold SE of 0.027. All five selected outer fits
converged. These concordances use the same preprocessing and all-pairs
fold scoring as the supervised-versus-two-step table above, rather than
the fold-specific scaling and within-cohort scoring of the working-model
CV estimate 0.6546. Scaling was prepared before the nested splits, so
the nested result evaluates weight selection conditional on that matrix.

The selection rule permits inner fits that reach the iteration cap if
their CV partial log-likelihood is finite and no fit error occurred.
The capped fits therefore remain in the nested selection calculations,
as specified by the original procedure. They limit interpretation of the
nested estimate and are not counted as convergence.

The full-data fit at weight 1 converged in 52 sweeps. Its loadings,
activities, coefficients and coefficient second moments are identical
to the saved K = 5 supervised full fit. The maximum difference between
the five weight-1 fold concordances and their S-11 reference values is
4.44 × 10⁻¹⁶. Programs 2, 4 and 5 remain survival-active. The selected
weight is unchanged, so no new external scoring is required and the
working configuration remains w_S = 1.

In `programs_full_fit.csv`, `risk_var_share` means
Var(β_k Z_k) / Var(Σ_j β_j Z_j). These ratios omit covariance terms and
sum to 1.140 across programs. They are not additive shares. The signed
covariance shares reported for program classification in item 7 use a
different definition and sum to one. `std_beta` is the absolute
coefficient per training SD, so use `beta` for its sign.

The historical rank-7 refit also selected weight 1. Its nested C was
0.656 with SE 0.027 and its CV partial log-likelihood range was 0.26.
Replace those numerical descriptions with the rank-5 values above.
The interpretation, that this weight search leaves the unit-weight fit
in place, is unchanged.

Files: `k5/s15/{metrics,weight_curve,reproduce_s11,nested_cv,verdict,programs_full_fit}.csv`.
Rebuild aggregates from the saved checkpoints with
`EXHIBIT_K=5 OMP_NUM_THREADS=1 Rscript results/benchmark_sim/run_yfb_supervised_real_refit.R summary`.
The fitter and predictor snapshot is recorded at commit `0555d299` in
`k5/s15/provenance.txt`. Checkpoints are git-ignored and are not included
in the commit.

### External scoring (item 4)

| Series | n (events) | w_S = 1 C (95% CI) | Two-step C | Paired difference (95% CI) | Slope w_S = 1 (SE) | Slope two-step (SE) |
|---|---|---|---|---|---|---|
| Dijk | 90 (81) | 0.599 (0.530, 0.667) | 0.599 | 0.000 (−0.018, 0.020) | 0.74 (0.38) | 0.68 (0.33) |
| Moffitt | 123 (83) | 0.538 (0.458, 0.608) | 0.538 | 0.000 (−0.024, 0.023) | 0.35 (0.36) | 0.33 (0.34) |
| PACA-AU array | 63 (38) | 0.664 (0.573, 0.756) | 0.666 | −0.001 (−0.021, 0.016) | 0.83 (0.35) | 0.74 (0.30) |
| PACA-AU RNA-seq | 52 (31) | 0.704 (0.596, 0.806) | 0.699 | +0.006 (−0.014, 0.026) | 1.85 (0.48) | 1.60 (0.41) |
| Puleo | 288 (181) | 0.645 (0.600, 0.687) | 0.651 | −0.007 (−0.017, 0.003) | 1.07 (0.21) | 0.97 (0.19) |

Equal-series mean C is 0.630 for w_S = 1 and 0.631 for two-step. The
pooled difference is −0.001 (record interval −0.009 to 0.007, donor
interval −0.010 to 0.008, 568 donors, 47 shared PACA-AU donors), and
−0.002 (−0.011 to 0.007) without PACA-AU RNA-seq. The stratified HR per
pooled SD is 1.39 (1.26 to 1.53) for w_S = 1 and 1.40 (1.27 to 1.54) for
two-step, with donor-clustered intervals 1.24 to 1.56 and 1.24 to 1.57.
The supervised slopes are closer to one than the two-step slopes in three
of five series (Dijk, Moffitt, PACA-AU array) and farther in PACA-AU
RNA-seq and Puleo. Every
w_S = 1 concordance equals the YFB value in
`outputs/lb_yfb_external_2026-11-12/external_cindex.csv`.

### Fixed-cutpoint risk groups (item 5)

The cutpoint maximizes the log-rank statistic on the out-of-fold training
scores of the K = 5 fold fits (each standardized by its fold's training
mean and SD), with at least 20% per group. It is 0.618 training SD for
w_S = 1 (59 of 273 training patients high) and 0.485 for two-step. Applied
unchanged to the external scores of the fixed K = 5 fit:

| Series | Low / high | High vs low HR (95% CI) |
|---|---|---|
| Dijk | 75 / 15 | 1.63 (0.88 to 3.02) |
| Moffitt | 101 / 22 | 0.99 (0.57 to 1.74) |
| PACA-AU array | 50 / 13 | 1.60 (0.77 to 3.33) |
| PACA-AU RNA-seq | 39 / 13 | 3.13 (1.46 to 6.69) |
| Puleo | 223 / 65 | 1.95 (1.39 to 2.72) |
| Pooled (stratified) | 488 / 128 | 1.69 (1.34 to 2.14), donor-clustered 1.31 to 2.18 |

Two-step: 469 / 147, pooled 1.78 (1.43 to 2.23), donor-clustered 1.40 to
2.27. The pooled HR per training SD of the supervised score is 1.54
unadjusted, 1.53 after adjusting for PurIST subtype and 1.44 after
adjusting for DeCAF subtype. In every series the supervised high-risk
group holds a larger share of PurIST basal-like tumors (0.31 to 0.73
against 0.02 to 0.10) and of DeCAF permCAF tumors (0.77 to 1.00 against
0.39 to 0.72). Moffitt is the only series with an HR below one (0.99). It
was also the weakest series at K = 7 (1.06).

### Comparators (item 6)

| Model | Mean C | Slope range | HR per SD (donor-clustered 95% CI) | Supervised minus model (95% CI) |
|---|---|---|---|---|
| Supervised (w_S = 1) | 0.630 | 0.35 to 1.85 | 1.39 (1.24, 1.56) | |
| Two-step (rank-5 fit with w_S = 0 + Cox) | 0.631 | 0.33 to 1.60 | 1.40 (1.24, 1.57) | −0.001 (−0.010, 0.008) |
| DeSurv | 0.614 | 0.45 to 0.67 | 1.43 (1.25, 1.64) | +0.016 (−0.009, 0.043) |
| Elastic-net Cox (2,218 genes) | 0.619 | 0.49 to 1.29 | 1.34 (1.20, 1.49) | +0.011 (−0.008, 0.031) |
| PurIST | 0.599 | 0.20 to 0.98 | 1.25 (1.13, 1.38) | +0.031 (−0.004, 0.066) |
| EBMF rank 15 + Cox | 0.575 | 0.10 to 0.49 | 1.16 (1.04, 1.29) | +0.055 (0.027, 0.084) |
| EBMF rank 40 + Cox | 0.581 | 0.10 to 0.31 | 1.16 (1.05, 1.29) | +0.049 (0.013, 0.087) |

The two-step score is the same factorization at rank 5 fitted with
w_S = 0, followed by Cox regression on its program scores (`twostep` in
`run_prompt6_comparators.R`), not flashier's EBMF. The table builder and
Chapter 4 previously labeled it "rank-K EBMF + Cox" and now say "rank-K
fit with w_S = 0 + Cox" at both ranks. The comparator rows equal the K = 7 table, apart from the two-step row,
which is now rank 5, and the differences, which are now taken against the
K = 5 supervised score. Every interval includes zero except those against
the rank-15 and rank-40 EBMF fits. The rank-5 flashier EBMF fit with Cox
regression, refitted with the rank sweep, has mean C 0.631 and a
difference of −0.001 (−0.010, 0.008), the same as the two-step row. Its
rank-40 fit stopped at flashier's iteration limit, as at K = 7.

### Program classification (item 7)

| Program | Active | Coefficient ratio | Signed share of risk variance | Active folds | Main |
|---|---|---|---|---|---|
| 2 | Yes | 1.00 | 0.591 | 5 | Yes |
| 3 | No | 0.27 | 0.034 | 3 | No |
| 4 | Yes | 0.56 | 0.141 | 3 | No |
| 5 | Yes | 0.62 | 0.234 | 2 | No |

Program 1 (ratio 2 × 10⁻⁸, coefficient zero) is omitted by the table's
rule. Survival supplies at most 2.5 × 10⁻⁶ of any program's loading
precision (2.1 × 10⁻⁶ at K = 7). Only Program 2 is a main program. Programs
4 and 5 meet the coefficient and variance thresholds but are active in
three and two folds. At K = 7, Program 3 also qualified.

### Programs, enrichment and signatures (item 8)

| # | E[β] (SD) | z | Active | Positive pole genes | Negative pole genes |
|---|---|---|---|---|---|
| 2 | −0.37 (0.08) | −4.5 | yes | METTL7A, ATP2A3, ACACB, PDCD4 | LOXL2, MMP14, LDHA, P4HA1 |
| 4 | 0.21 (0.09) | 2.4 | yes | SERPINA5, CTRB2, CPB1, REG1B | SLAMF6, BANK1, AKNA, EMB |
| 5 | 0.23 (0.08) | 2.8 | yes | UHMK1, ARHGAP5, RAD21, PRKDC | RCN3, CHPF, MXRA8, LOXL1 |
| 1 | 0.00 (0.00) | 0.0 | no | ZEB2, MYO5A, CCDC80, TCF4 | RAB25, AP1M2, PKP3, LAD1 |
| 3 | 0.10 (0.07) | 1.5 | no | RPS7, RPL6, RPS2, EEF1D | KIAA1211, SYTL2, PLS1, CORIN |

The two sets with the smallest adjusted p-value at each pole, as in the
chapter's enrichment table (all adjusted p < 0.001). Sets with larger
|NES| but equal adjusted p exist in `programs_enrichment.csv`.

- Program 2: positive pole Bailey ADEX (NES 3.56) and DeSurv D1 classical
  tumor (2.67). Negative pole Hallmark epithelial–mesenchymal transition
  (−2.88) and DeSurv D3 basal-like tumor (−2.49). Same sets as K = 7.
- Program 4: positive pole Bailey ADEX (2.75) and Reactome developmental
  lineage of pancreatic acinar cells (2.42). Negative pole Hallmark
  allograft rejection (−2.40) and interferon-gamma response (−2.14).
- Program 5: positive pole DeSurv D2 stromal-immune (2.78) and Hallmark
  mitotic spindle (2.49). Negative pole epithelial–mesenchymal transition
  (−2.26) and DeSurv D1 classical tumor (−1.89). Same sets as K = 7.

Signature correlations (Spearman, pooled training; range over the five
external series):

- Program 2 with Moffitt basal-like −0.68 (−0.74 to −0.69). With PurIST
  basal probability −0.66 to −0.45 externally. (K = 7: −0.71, −0.78 to
  −0.71, PurIST −0.64 to −0.45.)
- Program 4 with Bailey ADEX 0.73 (0.50 to 0.77) and Collisson
  exocrine-like 0.73 (0.47 to 0.75). With Bailey immunogenic −0.29 (−0.63
  to −0.12). With PurIST basal probability 0.05 to 0.50 externally.
- Program 5 with DeSurv stromal-immune 0.18 (−0.08 to 0.14). With PurIST
  −0.20 to 0.00. Its stromal association remains weak, as at K = 7
  (0.21, −0.03 to 0.23).

Proposed labels for the documents session to accept or change: Program 2
"Classical / basal-like" (unchanged), Program 4 "Acinar / immune",
Program 5 "Mitotic / stromal signals" (unchanged). Program 4's adverse
direction combines a higher acinar (exocrine) signal with a lower B-cell
and interferon signal. Acinar content in a bulk tumor sample often reflects
normal pancreas in the sample, so this label describes the genes and does
not imply tumor biology.

### External program associations (item 9)

HR per cohort-specific SD of each program score, separate unadjusted Cox fits:

| Series | Program 2 | Program 4 | Program 5 |
|---|---|---|---|
| Dijk | 0.64 (0.49 to 0.83) | 1.15 (0.92 to 1.44) | 0.75 (0.60 to 0.94) |
| Moffitt | 0.90 (0.73 to 1.12) | 1.02 (0.82 to 1.27) | 0.90 (0.72 to 1.12) |
| PACA-AU array | 0.80 (0.60 to 1.08) | 1.50 (1.10 to 2.05) | 1.12 (0.82 to 1.51) |
| PACA-AU RNA-seq | 0.72 (0.52 to 1.00) | 1.76 (1.21 to 2.57) | 1.43 (0.95 to 2.15) |
| Puleo | 0.72 (0.60 to 0.85) | 1.21 (1.04 to 1.41) | 0.96 (0.84 to 1.10) |
| Training direction | Protective | Adverse | Adverse |

Direction agrees with training in 5/5 series for Program 2, 5/5 for
Program 4 and 2/5 for Program 5 (K = 7: 5/5 for Program 2, 4/5 for Program
3, 2/5 for Program 5). The per-program KM figures are
`sup_external_km_program_{2,4,5}_k5.pdf`. There is no Program 3 figure at
K = 5 because Program 3 is not survival-active.

### Program correspondence with K = 7 (item 10)

Training score correlation r between program scores of the two fits
(273 patients) and the overlap of the 25 largest loadings at each pole of
the matched pair:

| K = 5 | Best K = 7 match | r | Next best abs r | Top-25 overlap (+ / −) | K = 5 z (active) | K = 7 z (active) |
|---|---|---|---|---|---|---|
| 1 | 1 | 0.995 | 0.45 | 16 / 18 | 0.0 (no) | 0.0 (no) |
| 2 | 2 | 0.993 | 0.45 | 14 / 21 | −4.5 (yes) | −4.8 (yes) |
| 3 | 3 | 0.944 | 0.25 | 19 / 14 | 1.5 (no) | 3.0 (yes) |
| 4 | 4 | 0.924 | 0.30 | 19 / 21 | 2.4 (yes) | 1.9 (no) |
| 5 | 5 | 0.988 | 0.07 | 22 / 23 | 2.8 (yes) | 2.8 (yes) |

The match is one to one and keeps the program numbers. K = 7 Programs 6
(IFI30, LILRB4, FCER1G, CD74, myeloid) and 7 have no K = 5 counterpart.
Their best absolute score correlations with any K = 5 program are 0.45
(Program 1) and 0.24. Neither was survival-active at K = 7 (z 0.0 and
−0.9). Program 6's scores correlate 0.45 with K = 5 Program 1 and −0.30
with K = 5 Program 4.

## Proposed status text

For PROJECT_STATUS.qmd and the progress book (owners decide):

> The exhibits that reported the K = 7 fit have been rebuilt at K = 5 from
> the working model (2026-11-12 K = 5 exhibits findings). Programs 2, 4 and
> 5 are survival-active, and only Program 2 meets the main-program rule.
> At K = 5 the estimated supervised-minus-two-step differences are small in
> training CV (+0.004, SE 0.004) and in the external series (pooled
> −0.001, donor interval −0.010 to 0.008). External concordance matches
> the earlier K = 5 scoring of the LB and YFB comparison exactly.
> The training-only survival-weight search also selects w_S = 1. Its
> nested CV concordance is 0.658 with SE 0.027. Of 216 grid fits, 27
> reached the iteration cap, all within the nested inner fits. Every
> selected outer fit and the selected full-data fit converged.

## Limits

- External scoring is descriptive. Nothing here selects or changes the
  model.
- The supervised-versus-two-step comparison is closer at K = 5 than at
  K = 7. Its sign in training CV (+0.004) and externally (−0.001) differs,
  and both are well within their SEs.
- `cairo_pdf` cannot load on this machine (XQuartz is not installed).
  The figure builders first fell back to `pdf()`, which does not embed
  fonts, for both ranks. They now fall back to macOS's Quartz PDF device,
  which embeds font subsets (Arial in place of Helvetica), before `pdf()`.
  All Chapter 4 figures at both ranks were rebuilt this way. Every font is
  embedded, and each figure's extracted text matches the earlier version
  except for how the minus sign is encoded.
- The program builder's correspondence with the historical α_F = 0
  benchmark is not interpretable at K = 5 and is not reported.
- During the fits, up to four R processes ran at once for a few minutes
  (two fit workers, the comparator run and a builder), one more than the
  three-worker limit.

## Verification and document handoff

The main exhibit outputs and rank-aware runners were committed in
`5982ce76`. Their file list and numbers were sent to the documents session
before the pause. The resumed work adds the S-15 aggregates and completes
this note. No new fit was run, and no file under `code/` was changed by
the exhibit work. The requirement to rerun `tests/run_tests.R` after a
`code/` change therefore does not apply.

A fresh independent review checked the S-15 results directly against all
216 checkpoints before the final commit. It verified the complete grid,
unique keys, absence of fit errors, finite test concordances and inner CV
partial log-likelihoods, five live programs, selected weights, convergence
counts, and the nested mean and SE. The full unit-weight fit exactly
reproduces its S-11 reference. The saved fitter and predictor and the
current `code/` and configuration have no differences from the recorded
snapshot commit. Review corrections distinguish convergence from reaching
the iteration cap, variance ratios from covariance shares, and summed
fit time from CPU time. Sampled external, cutpoint, program-rule,
correspondence, training and screening values agree with the aggregate
CSVs.

For the documents session, replace the remaining rank-7 weight-refit
paragraph and table using item 3 above. The selected weight remains 1.
Update the nested C from 0.656 to 0.658, its SE stays 0.027, the outer
choices are 1, 1, 1, 10 and 50, the likelihood range is 0.635, and the
weight-50 median active-program precision share is 5.80 × 10⁻⁵. Report
the 27 capped inner fits. The full weight curve is in `k5/s15/weight_curve.csv`.
The writeup exhibit swap, rendering, Chapter 4 review, full-prelim rebuild
and status updates remain with the documents session. Those files were
not edited or staged in this continuation.
