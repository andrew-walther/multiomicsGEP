# Chapter 4 numerical sources

The paths headed `outputs/` below are relative to
`results/benchmark_sim/`. Figure axes, table numbers and printed precision
follow the builders listed below. The chapter preserves the frozen comparisons
and distinguishes the two external resampling units.

## Rebuild the figures and tables

Run from the repository root:

```sh
Rscript paper/prelim/make_figures_supervised_sim.R
Rscript paper/prelim/make_figures_supervised_programs.R
Rscript paper/prelim/make_figures_supervised_external_km.R
```

The simulation builder reads saved aggregate CSVs. The program and KM builders
also need saved fits and local PDAC inputs. If the ignored external score cache
is absent, set `PDAC_DATA_ROOT` and run
`Rscript results/benchmark_sim/save_yfb_external_patient_scores.R` before the KM
builder. This cache contains patient-level outcomes and projections and belongs
under ignored `data/`. The chapter embeds generated LaTeX tables and references
vector PDFs in `paper/prelim/figures/`. The ledger below identifies the sources.

| Location | Numbers reported | Source |
|---|---|---|
| Introduction | Five-year pancreatic cancer relative survival, 13.7%; October 15 paired difference −0.003, interval −0.029 to 0.024, margin −0.02 | `SEER2026` in `paper/ssbmf-refs.bib`; `outputs/supervised_real_check/delta_c.csv` at commit `65cf90a`; October 15 entry in `DECISIONS.md` |
| Methods | Point-Laplace molecular priors, Normal coefficient prior, unit projection SD, weights 0 and 1, relative ELBO tolerance 10⁻⁵, maximum 300 sweeps, positive SD and absolute z ≥ 1.96 for survival activity, correlation ≥ 0.7 for recovery, 73 derivation checks | `run_yfb_supervised_real.R`, `run_yfb_survival_weight_p2000.R`, `code/fit_ssbmf_yfb.R`, `config/globals.yml`, `derivations/SSBMF-YFB-model-derivations/verify_derivation.R` and its companion derivation |
| Weight and simulation design (T1) | Fixed and CV weight grids; number of genes, patients, programs, datasets and fits; planted shares; null design; oracle concordance ranges; folds; criteria for recovery, false activity, calibration and convergence | `outputs/yfb_survival_weight_2026-11-04/grid/{design,metrics}.csv`; `outputs/yfb_supervision_strength_2026-11-05/{metrics,null_cells,selected_vs_twostep}.csv`; `outputs/yfb_near_threshold_2026-11-05/metrics.csv`; S-9, S-13 and S-14 registrations and findings in `docs/plans/` and `docs/findings/`; numerical table construction in `paper/prelim/make_figures_supervised_sim.R` |
| PDAC data | TCGA 144, CPTAC 129, total 273 and 139 deaths; five series, 616 records, 568 donors, 47 shared PACA-AU donors; top 3,000 per training cohort, D4 2,218, all 17,900, K = 7; five folds and 2,000 bootstrap replicates | `outputs/yfb_supervised_real_2026-11-04/{folds,summary}.csv`; `outputs/yfb_supervised_external_2026-11-05/{external_cindex,external_paired_diff_donor}.csv`; `docs/findings/2026-11-06-external-cohort-uses.md`; `run_yfb_supervised_real.R`; `score_yfb_supervised_external.R`; `config/globals.yml` |
| Rank selection | KL-2 recovery 17/20 and 19/20 unpruned, 2/20 at both noise levels for either pruning arm; retained median ranks 2 and 1; working rank 7 | `outputs/flat_beta_pruning_rerun_2026-11-06/summary.csv`; `docs/findings/2026-11-06-flat-beta-pruning-rerun-findings.md` |
| S-9 (F1/F2) | At 2,000 genes and weight 1, 10/10 recovered and converged in each cell, mean false activity 0.9, 0.1 and 0.0; weights ≥40 stopped at 300 sweeps; calibration medians 1.45, 2.16 and 3.23; C and precision-share ranges; CV selected weight 1 in 55/60, share-targeted in three, weight 10 in two; earlier reference recovered 0 or 1/10 | `outputs/yfb_survival_weight_2026-11-04/grid/{metrics,summary,cal_slope_exact,cv_chosen,cv_chosen_cells,gate_cells}.csv`; `docs/findings/2026-11-04-survival-weight-p2000-findings.md`. F1 joins exact-time slopes by dataset and arm. |
| S-13 (T6/F12) | All panel values, selected-weight frequencies, 10/10 recovery, differences up to 0.007, calibration medians, null false activity 1.1, joint-posterior comparison within 0.004 | `outputs/yfb_supervision_strength_2026-11-05/{selected_vs_twostep,cv_selected,weight_curves,verdict,null_cells}.csv`; `docs/findings/2026-11-05-supervision-strength-findings.md` |
| S-14 (T6/F13) | All panel values; near-threshold shares 0.05%, 0.1%, 0.2%; strong-survival recovery 5/10 versus 0/10 in the smallest cell, paired differences 0.156, 0.041 and 0.008; 29/30 direction; realistic oracle C about 0.74 and differences 0.002 to 0.038; fixed-weight calibration range 0.04 to 0.65 | `outputs/yfb_near_threshold_2026-11-05/metrics.csv`; `docs/findings/2026-11-05-near-threshold-findings.md` |
| Post hoc false-active contribution | S-14 median variance ratio up to 1.3% at weight 1 in the 0.1% and 0.2% cells; up to 0.008 change in mean C across the two arms in those cells; near-threshold ratios 8.7% to 68.8% | `outputs/yfb_supervision_strength_2026-11-05/false_active_contribution.csv`; `outputs/yfb_near_threshold_2026-11-05/false_active_contribution.csv`; post hoc sections of the S-13 and S-14 findings |
| S-11 (T2/F3) | All training CV C, SE, partial log-likelihood, convergence and active counts; paired ΔC 0.011 (SE 0.005); three active programs; two folds at sweep cap, max coefficient 0.46; precision share 10⁻⁸ to 10⁻⁶; agreement 1.000; feature-size, screening and null comparisons | `outputs/yfb_supervised_real_2026-11-04/{fits,summary,programs_D4,size_agreement,size_preferred,screen_paired,null_cells}.csv`; `docs/findings/2026-11-04-supervised-real-data-findings.md` |
| S-15 | Weight 1 selected; criterion −775.168 to −775.429; nested C 0.656 (SE 0.027) versus 0.663 fixed; at weight 50 median share 7.3 × 10⁻⁵ and agreement to reported precision | `outputs/yfb_supervised_refit_2026-11-06/{verdict,weight_curve,programs_full_fit}.csv`; `docs/findings/2026-11-06-supervised-real-refit-findings.md` |
| S-12 (T3/F4/F9) | Series counts, C, intervals and calibration slopes; mean C 0.620; pooled differences 0.007 (record interval −0.002 to 0.017, donor interval −0.003 to 0.017) and 0.006 excluding PACA-AU RNA-seq; pooled HRs 1.37 (1.24 to 1.51), 1.35 (1.22 to 1.49) and historical DeSurv 1.50 (1.31 to 1.72) | `outputs/yfb_supervised_external_2026-11-05/{external_cindex,external_paired_diff,external_paired_diff_donor,external_pooled_hr}.csv`; DeSurv as cited in the S-12 findings. The primary slope CSV shows four of five estimates closer to one, with PACA-AU RNA-seq at 1.57 versus 1.26. That corrects the five-of-five wording in the earlier summary. |
| Historical rank-matched comparison | Frozen series C and mean 0.627, single-log mean 0.625, paired rank-matched ΔC 0.011 (−0.013 to 0.033) | `results/benchmark_snapshot_2026-10-14/single_modality/external_cindex_ci.csv`; `outputs/kl_fix_rerun_2026-10-30/external_results.csv`; `docs/findings/2026-11-06-orientation-audit.md` and August 20 entry in `DECISIONS.md` |
| Programs (T4/F7/F8) | Seven programs, three active, coefficient means −0.405, 0.218, 0.230 and table SDs and z scores; figure scaling and numbers of genes; Pearson correlations and overlap among top 25 genes, including historical matches 0.623, 0.649 and 0.180 | `outputs/yfb_supervised_real_2026-11-04/programs_D4.csv`; `outputs/prelim_ch4_2026-11-06/{programs_table,programs_correspondence,programs_correspondence_topgenes}.csv`; `paper/prelim/make_figures_supervised_programs.R` |
| Biology (T5/B1–B3) | Enrichment estimates, gene-set eligibility counts and adjusted P values; signature scores summarized by correlations | `outputs/prelim_ch4_2026-11-06/{programs_enrichment,programs_subtype_concordance}.csv`; `paper/prelim/make_figures_supervised_programs.R` |
| External KM and direction (F5/F6/B5) | Risk counts on the displayed 0–60 month time grid; 15 Cox HRs and 95% Wald intervals; direction agreement in 5/5, 4/5 and 2/5 series | `paper/prelim/make_figures_supervised_external_km.R` reads ignored fixed S-12 scores and verifies counts and concordance against `external_cindex.csv`. The public summary is `outputs/prelim_ch4_2026-11-06/external_program_direction.csv`; training signs are in `programs_D4.csv`. |

## Exhibit builders

- `make_figures_supervised_sim.R` reads one row per saved fit or cohort
  comparison, checks design sizes and aggregates by simulation cell.
  It joins the exact-time slopes and pairs the threshold arms by dataset
  before computing differences. Separate source columns supply the
  concordance and calibration intervals. The outputs are vector PDFs and
  LaTeX table fragments. This builder needs no fitted object or patient data.
- `make_figures_supervised_programs.R` reads the saved supervised,
  unsupervised and historical fitted loading matrices. It aligns rows by gene
  name, selects the largest loadings at each pole, and compares the fixed
  columns on their shared genes. Zero historical columns give undefined
  correlations. `fgsea` tests signed loading ranks using the modeled genes as
  the universe. Spearman correlations compare fixed scores with descriptive
  signature scores and saved PurIST probabilities. The five unestimable
  Reactome p-values remain missing in the aggregate output. Gene-set and
  signature summaries describe the fixed programs and play no part in
  model selection.
- `save_yfb_external_patient_scores.R` reads the fixed fits and local
  cohorts, applies the registered preprocessing and reconstructs the
  historical two-step Cox coefficients on training data. It verifies
  series counts, events and frozen concordance before saving the local
  score matrix under ignored `data/`. An existing verified cache is
  preserved. No molecular model is refitted.
- `make_figures_supervised_external_km.R` verifies that the frozen
  projections reproduce the fitted risk score. It then forms within-series
  tertiles with training-fixed direction and uses `survival::survfit()` to
  get curves and risk counts. Separate `survival::coxph()` fits estimate
  descriptive HRs per within-series SD of each active program score.
  Displaying the curves through 60 months does not truncate the Cox
  analyses. Public outputs are curves, counts at risk and aggregate HRs.

## Comparators, risk groups and program rule (added 2026-11-06 review)

Tables `tab_comparators.tex` and `tab_program_rule.tex` are written by
`Rscript paper/prelim/make_tables_comparators.R` and pasted into the chapter.
Text values come from the same files under `results/benchmark_sim/outputs/`:

- Mean C, paired differences, pooled HRs and donor-clustered intervals:
  `prompt6_comparators/external.csv`, `paired.csv`, `pooled_hr.csv`; ranks
  15 and 40 from `prompt7_desurv_style/rank_external.csv`, `rank_paired.csv`,
  `rank_pooled_hr.csv`. Ridge selection: `prompt6_comparators/tuning.csv`.
- Cutpoints −0.046 and 0.526: `prompt7_desurv_style/cutpoints.csv`. Group HRs
  1.52 (1.24 to 1.87), 1.69 (1.32 to 2.16) and group sizes 297 and 144:
  `group_hr.csv` (donor-clustered limits). Unadjusted 1.51 and adjusted 1.48
  and 1.40: `adjusted_hr.csv`. Subtype proportions: `enrichment.csv`.
- Program rule and the 2.1e-6 bound: `program_rule.csv`, `fold_matches.csv`.
  Simulation checks 80/80, 35/36 and 61/128: `s13_rule.csv`, `s14_rule.csv`
  (two-step rows). Mixed-signal 9/10 and 2/10: `mixed_simulation.csv`.
- Supervised and two-step pooled donor-clustered intervals 1.22 to 1.53 and
  1.21 to 1.51, and DeSurv 1.43 (1.25 to 1.64): `prompt6_comparators/pooled_hr.csv`.
- 139 training deaths: `yfb_supervised_real_2026-11-04/folds.csv`.

## Build and review counts

Page counts are from `pdfinfo` and the Chapter 4 and Appendix A boundaries
in `pdftotext -layout` for `bios-dissertation/prelim/full-draft/prelim.pdf`.
The baseline is bios-dissertation commit `a09eb9d`. The standalone and
abstract counts use `project3-ssbmf-draft.pdf` and `abstract.pdf`.
The assembly test count comes from `test_prelim_assembly.py`, and the R suite
count comes from `tests/run_tests.R`. Round 3 counts are in its Chapter 4
`SUMMARY.md` and per-edit dispositions in `independent-review.md`.
