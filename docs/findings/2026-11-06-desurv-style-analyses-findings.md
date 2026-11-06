# DeSurv-style analyses of the fixed S-11 models

The [registration](../plans/2026-11-06-desurv-style-analyses-plan.md) fixes the
two S-11 models, the cutpoint search, the program rule and the simulation
design before external evaluation.

## What was estimated

For each reference, one cutpoint maximizes the training OOF log-rank
statistic subject to the minimum group size. Each held-out score was
standardized with the corresponding fold-training mean and SD. External
scores use the full-training mean and SD and the same cutpoint in every
cohort. Kaplan–Meier estimates are displayed on a monthly grid, with counts
at risk at the labeled ticks. The adjusted and unadjusted HRs in the
subtype table use the same eligible rows for each subtype. The local
DeCAF calls use the older permissive `permCAF` and restraining `restCAF`
labels. The enrichment table identifies the call actually tested. The
later publication uses promoting `proCAF` terminology. [DeCAF preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC11118336/)
and [published DeCAF study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12923980/).

Program HRs are marginal, per training SD, and use the stored loading
signs. Signed variance shares allocate each component's covariance with
the total score to that component and sum to one. The fold stability
calculation matches loading directions one to one before applying the
activity and correlation criteria.

The fixed S-11 checkpoints omit E[L²]. For the final-state precision audit,
E[L]² gives a lower bound on the genomics contribution to precision,
which gives an upper bound on survival's share. The precision label is
identified as false where that upper bound does not exceed the registered
threshold. This final-state calculation uses the final fit's risk weights
and is distinct from the within-sweep trace.

The saved S-13 and S-14 checkpoints lack per-program fold identities, so
only the coefficient-and-share portion of the rule can be checked there.
S-13 retains full-fit posterior activity. In S-14 the two-step activity can
be reconstructed by exact-time Cox regression, while the supervised
checkpoints retain coefficients but omit their posterior SDs. Missing
activity is explicit in the CSV. The simulation audit cannot validate the
full three-part rule.

The mixed-signal simulation calibrates the two raw molecular shares using
the S-14 generator and assigns equal standardized prognostic coefficients.
It uses training-scaled log-risk variability and reports gene precision
as the fraction of the top 50 absolute fitted weights in the known planted
support. Molecular recovery and held-out discrimination answer different
questions and are reported separately.

## Reproduce and trace the results

From the repository root, with the S-11 artifacts and `PDAC_DATA_ROOT` set
as for Prompt 6:

```sh
Rscript results/benchmark_sim/run_prompt7_real.R
Rscript results/benchmark_sim/run_prompt7_ebmf.R
Rscript results/benchmark_sim/run_prompt7_simulation.R
Rscript results/benchmark_sim/run_prompt7_simulation.R mixed
```

All reported numbers and exhibits are indexed below. Paths are relative to
[`outputs/prompt7_desurv_style/`](../../results/benchmark_sim/outputs/prompt7_desurv_style/).

| Command above | Files and values |
|---|---|
| Real | `cutpoints.csv` has the frozen cutpoints and OOF group counts and statistics. `group_hr.csv` has external group counts and HRs. `at_risk.csv`, `km_monthly.csv`, `km_w1.pdf` and `km_twostep.pdf` give cohort KM exhibits and counts. |
| Real | `program_hr.csv` and `program_forest.pdf` give the per-program HRs, confidence limits and fit denominators. `adjusted_hr.csv` has same-row unadjusted and adjusted HRs. `enrichment.csv` has subtype counts, proportions and Fisher p-values. |
| Real | `program_rule.csv` and `fold_matches.csv` give all rule inputs, PVE and precision labels, and each fold match. |
| EBMF | `rank_ebmf_fits.csv` has requested and achieved ranks and fit statuses. `rank_external.csv`, `rank_paired.csv` and `rank_pooled_hr.csv` give cohort scores, paired donor contrasts and stratified HRs. |
| Simulation | `s13_rule.csv` and `s14_rule.csv` give checkpoint retention, false-active coverage and top-gene precision. |
| Simulation mixed | `mixed_simulation.csv` gives both planted directions per fit, realized raw shares, matching and recovery, false activity, top-gene overlap, C and fit status. |

The source tables include the contributing sample or fit denominators.
The KM figures are aggregate and the local working objects are kept out of Git.

## Results

- The standardized training cutpoints were −0.046 for supervised YFB and
  0.526 for two-step Cox. The pooled high-versus-low HR was 1.52
  (95% CI 1.25 to 1.85) and 1.69 (1.35 to 2.12), respectively. Both
  models assigned all 616 eligible external records to groups. The
  cutpoints come from `cutpoints.csv` and the counts and HRs from
  `group_hr.csv`. The subtype analyses had complete call coverage in
  the eligible rows (`adjusted_hr.csv`, `enrichment.csv`).
- The rule retained supervised programs 2 and 3. Program 5 met the
  coefficient and variance-share thresholds and was active in only
  three of five folds. No active program met either survival-only
  criterion. The largest active-program upper bound on survival
  precision share was 2.07 × 10⁻⁶ (`program_rule.csv`, `fold_matches.csv`).
- The partial rule retained the planted active program in all 80 S-13
  checkpoint fits and rejected 35 of 36 false-active directions. For
  S-14 it passed all 160 loading matches among 240 fits. In the 60 S-14
  two-step fits with reconstructed activity, it retained 36 matched
  active programs and rejected 61 of 128 false-active directions.
  The remaining 180 S-14 fits do not support an exact false-active
  rejection count (`s13_rule.csv`, `s14_rule.csv`).

| Mixed-signal arm | Strong recovered | Near-threshold recovered | Near-threshold mean top-50 precision | Mean validation C |
|---|---:|---:|---:|---:|
| Two-step | 10/10 | 9/10 | 0.826 | 0.731 |
| A1, weight 1 | 10/10 | 9/10 | 0.858 | 0.735 |
| A1, weight 5 | 10/10 | 2/10 | 0.346 | 0.667 |
| A2, weight 1 | 10/10 | 9/10 | 0.860 | 0.734 |

All 40 mixed-signal fits converged and contributed both planted-direction
rows. The table uses all replicates within each arm from
`mixed_simulation.csv`. The poorer recovery at weight 5 accompanies
poorer validation discrimination in this design. The comparison is a
simulation result under one registered mixed-signal design.

## Fit accounting

The risk-group and program runner reused the full and fold-specific S-11
factorizations. It estimated the two-step Cox model on the full training
data and within each fold. All 12 risk-group HR fits, 42 per-program HR
fits and 48 subtype HR fits returned estimates. There were no external
sample drops for these analyses (`cutpoints.csv`, `fold_matches.csv`,
`group_hr.csv`, `program_hr.csv`, `adjusted_hr.csv`).

The checkpoint audit replayed 60 exact-time two-step Cox regressions for
activity and refit no factor models (`s14_rule.csv`). The mixed-signal
denominator is 40 factor-model fits and 10 two-step Cox regressions.
No fit failed, no replicate was dropped, and every C was defined
(`mixed_simulation.csv`).

All three EBMF fits achieved their requested ranks and returned Cox scores
for all 616 records. The equal-cohort mean C was 0.614 at rank 7, 0.575 at
rank 15, and 0.581 at rank 40 (`rank_external.csv`). The supervised-minus-EBMF
contrasts were 0.0067 (95% donor bootstrap interval −0.0032 to 0.0168),
0.0452 (0.0213 to 0.0707) and 0.0398 (0.0080 to 0.0743), respectively
(`rank_paired.csv`). The fit object has no convergence flag. Rank 40
issued a maximum-iteration warning and its results are reported as
returned (`rank_ebmf_fits.csv`). None of the rank fits threw an error and
no scored patients were dropped. These comparisons do not select a rank.
