# November 6 writeup numerical sources

The inline R expressions and tables in [writeup.qmd](writeup.qmd) read the
files below. Paths beginning with `outputs/` are relative to
`results/benchmark_sim/`. Method constants come from the cited registration,
configuration or derivation. Study identifiers, program numbers, dates and
table and figure indices label the evidence.

| Location | Reported numbers | Evidence files |
|---|---|---|
| Summary and training CV | Counts of training patients and folds, C and SE for both fits, paired difference and its SE, CV partial log-likelihood, convergence and number of active programs | `outputs/yfb_supervised_real_2026-11-04/{folds,summary,fits,programs_D4}.csv`. The paired SE uses within-fold C differences. |
| Model and fitting | Unit training projection SD; weights including 0 and 1; model rank; relative ELBO tolerance and sweep cap | `code/fit_ssbmf_yfb.R`; `config/globals.yml`; `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`; rank is read from KL-2 `summary.csv`. |
| Data and cohort table | Training 273 patients, 139 deaths; per-cohort counts and external missing genes; 616 external records, 568 donors under the distinct-cohort assumption, 47 shared donors and the PACA-AU RNA-seq duplication | `outputs/yfb_supervised_real_2026-11-04/folds.csv`; `outputs/yfb_supervised_external_2026-11-05/external_cindex.csv`; `paca-au-donor-audit.csv`; verified against `external_paired_diff_donor.csv` in the S-12 directory. |
| Preprocessing, dates and prior external uses | Log transform, top 3,000 ranking rule, single- and double-log gene counts and shared count, historical rank and weight choices | `code/preprocess_desurv.R`; `config/globals.yml`; `outputs/single_log_check/gene_sets.csv`; `docs/findings/2026-11-06-external-cohort-uses.md`; dated entries in `DECISIONS.md`. |
| Simulation designs | Gene counts, planted shares, datasets per cell, candidate weights, survival-activity cutoff and absolute correlation cutoff; S-14 oracle C ranges | `outputs/yfb_survival_weight_2026-11-04/grid/design.csv`; `outputs/yfb_near_threshold_2026-11-05/metrics.csv`; S-9, S-13 and S-14 plans and corresponding `run_yfb_*.R` runners. |
| S-12 external scoring | Series C, intervals, slopes and mean C; pooled +0.007, record interval −0.002 to +0.017, donor interval −0.003 to +0.017, sensitivity without PACA-AU RNA-seq; pooled hazard ratios and intervals; slopes closer to one in four of five series | `outputs/yfb_supervised_external_2026-11-05/{external_cindex,external_paired_diff,external_paired_diff_donor,external_pooled_hr}.csv`. Calibration compares absolute distance to one, matching series by cohort name. |
| October 15 registered check | Registered difference −0.003 (−0.029 to 0.024), margin −0.02; KL-fix difference −0.006 (−0.032 to 0.021); underlying mean C and simulation recovery at both noise levels | `outputs/supervised_real_check/delta_c.csv` at `65cf90a` and the current CSV; `outputs/supervision_diagnostics/c3c4_n400_b3.3/results.csv`; October 15 decision in `DECISIONS.md`. |
| Loading precision, sizes, screening and null checks | Ranges of survival precision share, largest fitted coefficient, matched loading correlation, performance and convergence by size, best-match correlations, paired screening differences and null false activity | `outputs/yfb_supervised_real_2026-11-04/{fits,summary,size_agreement,screen_paired,null_cells}.csv`; `outputs/yfb_survival_weight_2026-11-04/grid/metrics.csv` for the simulated median shares; `outputs/prelim_ch4_2026-11-06/programs_correspondence.csv`. |
| S-9 | Recovery, mean false activity, exact-time calibration medians 1.45, 2.16 and 3.23, convergence and validation C, 55/60 unit-weight CV selections | `outputs/yfb_survival_weight_2026-11-04/grid/{design,metrics,summary,cal_slope_exact,cv_chosen,cv_chosen_cells}.csv`. |
| S-13 | Selected weights by cell and arm, recovery, false activity, validation C, median slopes, null false activity and gate verdicts | `outputs/yfb_supervision_strength_2026-11-05/{selected_vs_twostep,verdict,weight_curves,cv_selected,twostep}.csv`; S-13 registration for the slope band and null rule. |
| S-14 | Recovery, false activity, oracle C, paired C and counts favoring supervision by cell, including 5/10 versus 0/10 and +0.156 at the smallest share with strong survival | `outputs/yfb_near_threshold_2026-11-05/metrics.csv`; the pairing joins on strength, planted share and seed. |
| Post hoc nonmatched contributions | Cell medians of the variance ratio, cell means of full-minus-matched C, 1.3% and 2.9% larger-share maxima by arm, 0.008 maximum absolute mean change, 8.7% to 68.8% at the smallest share | `outputs/{yfb_supervision_strength_2026-11-05,yfb_near_threshold_2026-11-05}/false_active_contribution.csv`; definition in `results/benchmark_sim/summarize_false_active_contribution.R`. |
| S-15 | Chosen weight 1, weights 1 through 50, partial log-likelihood range 0.26, nested fold choices, nested C 0.656 (SE 0.027) versus fixed 0.663, loading precision share 7.3 × 10⁻⁵ at weight 50 and all tabulated fit summaries | `outputs/yfb_supervised_refit_2026-11-06/{verdict,weight_curve,nested_cv}.csv`. |
| KL-2 | Rank-seven working fit; all recovery counts, noise levels, false-active means, returned convergence counts and median ranks by arm | `outputs/flat_beta_pruning_rerun_2026-11-06/summary.csv`; registration `docs/plans/2026-11-06-flat-beta-pruning-rerun-plan.md`. |
| Historical benchmark and corrected input | Per-cohort historical and single-log C, difference intervals, risk-score correlations, mean 0.625 versus 0.627, training CV, rank-matched pooled difference and historical LASSO mean | `outputs/kl_fix_rerun_2026-10-30/{external_results,cv_summary}.csv`; `outputs/ebmf_cox_external/{yfb_vs_ebmf_k7_matched_ci,2026-10-30-orientation-rescore}.csv`. |
| Programs, enrichment and subtype analysis | Program counts, posterior means and SDs, active labels, top genes, enrichment scores and adjusted p-values, unestimable p-value count, training and external signature correlations | `outputs/prelim_ch4_2026-11-06/{programs_table,programs_enrichment,programs_subtype_concordance}.csv`. The reported enrichment row for each pole has the smallest estimable adjusted p-value. |
| Historical matches and external program associations | Signed Pearson loading correlations, top-gene overlaps, per-series hazard ratios and intervals and agreement counts 5, 4 and 2 | `outputs/prelim_ch4_2026-11-06/{programs_correspondence,programs_correspondence_topgenes,external_program_direction}.csv`. |
| Figures | Axes, plot values, risk counts and loading scales of the reused Chapter 4 exhibits | [Chapter 4 numerical source ledger](../../prelim/notes/chapter4-number-sources-2026-11-06.md). The writeup selects the original vector PDF for LaTeX and the conversion for HTML. |

The new results chunks use the saved fit-level summaries directly. They
pair arms by dataset before summarizing concordance, and summarize
nonmatched contributions within each design cell before comparing cells.
These operations preserve the distinction between paired differences,
medians of variance ratios and registered false-active counts. The only
new figure helper selects the appropriate vector file for the active
format and delegates inclusion to `knitr::include_graphics()`.

The rendered PDF has 26 pages, as counted by `pdfinfo writeup.pdf`.
The figure check compares the 11 numbered PDF captions and 11 embedded
HTML SVGs against the figure labels in `writeup.qmd`.
