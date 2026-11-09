# External-cohort use and donor resampling

Dijk, Moffitt, PACA-AU array, PACA-AU RNA-seq and Puleo all contributed to
the historical comparisons that informed the single-modality configuration.
**No series in that panel qualifies as untouched.** The matched ICGC
multimodal set overlaps PACA-AU and has already been examined. The available
records do not establish an untouched external cohort for a final claim.

## Decisions informed by external outcomes

The links point to the decision log at commit `b6ef631`. Earlier choices
are included even when the current survival-weighted fit supersedes them.

| Decision | What the external outcomes did | Record |
|---|---|---|
| Historical β prior | External concordance, alongside simulation, favored point-Normal over point-Laplace. The prior used in the current fit was selected in later work. | [April 24 prior decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L5681-L5687) |
| Initialization constraints | External concordance for TCGA-trained LB, together with training ELBO and retained factors, informed the decision to discard nonnegative and Cox-sign-aligned SVD starts. Restart selection itself used training ELBO. | [May 6 initialization](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L4990-L5038) |
| Merged training and cohort indicators | External comparisons accompanied the structural diagnosis and were used to recommend indicator-adjusted LB as a sensitivity arm and pooled YFB without an indicator as the primary arm. The later cohort-specific loading, β and stratified-hazard comparisons did not displace pooled D4. | [May 22 cohort report](2026-05-22-cohort-lmm-benchmark-report.qmd), [May 25 benchmark](2026-05-25-merged-benchmark-report.qmd), [September 4 cohort-aware panel](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L2193-L2265) |
| Per-cohort standardization and M5 preprocessing | The extended benchmark compared preprocessing, model and cohort-indicator variants on all five external series. Highest mean external C supported M5, per-cohort standardization without rank transformation and pooled YFB. The external result also supported retaining the proposed rank floor. | [May 25 decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L4554-L4633) |
| D4 gene screen and D5 rejection | External C favored per-training-cohort combined mean and variance ranks, intersection and standardization. The D4 versus D2 gain cleared the historical adoption threshold, and adding a cohort indicator in D5 hurt external C. | [May 27 decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L4513-L4550) |
| Starting rank and warm-start/BO alternatives | Rank 7 first came from training CV. July external comparisons informed retaining it over smaller cold starts and over near-tied warm-start and BO alternatives. In August, ranks 5 and 6 were screened out on external C before ranking by ELBO. The wider August sweep used external performance to support the plateau and the rank 7 recommendation. | [July 12 training CV](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L4097), [July 13 parsimony](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L3631-L3668), [August 19 screen](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L2918-L2933), [August 27 sweep](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L2423) |
| Retaining the July methodological changes | The branch review included an external net-benefit gate. The drop in C was accounted for by correcting the mismatch between training and external preprocessing. The correction itself followed the training recipe. | [July 12 correction](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L4267), [July 15 gate](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L3461-L3475) |
| Loading supervision weight α_F and frozen-loading duration N_frozen | The August cold-start sweep and joint grid used external C. The grid did not justify departing from α_F = 0 and N_frozen = 0. These defaults predated that grid. | [August 20 grid and cold-start sweep](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L2975-L3089) |
| Historical EBMF + Cox comparator and headline | External results were reviewed for unpenalized and LASSO Cox arms across the historical EBMF configurations. The rank-matched comparison qualified the earlier claimed advantage. The later manuscript headline used the large-rank LASSO comparator. The LASSO penalty itself was selected by training CV. Issue #6 requested leading with the rank-matched comparison. | [August 20 comparator decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L2568-L2705), [issue #6](https://github.com/andrew-walther/multiomicsGEP/issues/6) |
| October supervised fit disposition | Simulation chose the configuration before external scoring. The registered external criterion failed, and the observed estimate and uncertainty informed the decision to keep the configuration for further development while recording the failure. | [October 15 result and decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L1079-L1136) |
| Single-log correction | The duplicate log was a preprocessing error. A fixed external sensitivity check measured the consequence before the decision to adopt the corrected input for new work. The check showed little change in external concordance and used no competing transformation grid. | [October 26 single-log check](2026-10-26-single-log-check.md), [October 30 decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L288-L309) |

## Evaluation and reporting uses

- The October 2 orientation correction fixes the historical D4 score direction
  from its convention. External outcomes were used to verify the published
  estimates. Earlier scripts did choose score direction separately on external
  outcomes. See the [October 2 decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L1546-L1578)
  and the [orientation audit](2026-11-06-orientation-audit.md).
- The October 16 accuracy review documented overlap and separated the primary
  matched ICGC PDAC set from its nested broader-histology set. The manuscript's
  reporting qualifications followed inspection of the existing analyses.
  See the [October 16 decision](https://github.com/andrew-walther/multiomicsGEP/blob/b6ef631/DECISIONS.md#L336-L359)
  and [preprocessing audit](2026-10-26-preprocessing-audit.md).
- External concordance, calibration, subtype and program concordance, subgroup
  sensitivity and retrospective benchmarking were also reported. Where a
  comparison was used to choose or retain a configuration, it is listed above.
  The fixed November S-12 comparison inherited D4 preprocessing and rank from
  the earlier studies. Its survival weight, feature-set arm and survival-screen
  arm were fixed using S-11 training CV. The [S-12 scoring plan](../plans/2026-11-05-supervised-external-scoring-plan.md)
  documents this separation. The [S-15 refit](2026-11-06-supervised-real-refit-findings.md)
  used training data only and retained the existing fit.

## S-12 donor sensitivity

The runner joins the kept PACA-AU sample IDs from `load_pdac_raw()` to
`icgc_sample_id` and `icgc_donor_id` in the local original metadata. Array and
RNA-seq use the same donor namespace. The remaining cohorts use cohort-prefixed
sample IDs as donor proxies. The resampling frame conditions on each donor's
observed set of series. Each draw carries every record for the donor and both
model scores. Concordance pairs remain within series. The pooled estimand is
the equal-weight mean of the series-specific differences, so the point
estimate is the same as in the record bootstrap.

The scored panel contained 616 records representing 568 donors under the
separate-cohort assumption. The two PACA-AU series shared 47 donors, and
RNA-seq contained 52 records from 51 donors. The overlap and donor total are
in the donor CSV. The RNA-seq counts are also documented in the October
[accuracy dispositions](../../paper/writeup/2026-10-16-writeup/_process/Accuracy_Review_Dispositions_10_16_26.md).
The record count is the sum of `n` for either arm in
`external_cindex.csv` in the same output directory.

| Resampling unit | Pooled ΔC | 95% percentile interval | Aggregate source |
|---|---:|---:|---|
| Record | +0.0071 | −0.0022 to +0.0166 | [Registered record bootstrap](../../results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/external_paired_diff.csv) |
| Donor | +0.0071 | −0.0029 to +0.0173 | [Donor sensitivity](../../results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/external_paired_diff_donor.csv) |

Both intervals include zero. Donor resampling addresses dependence among
records from the same donor. The inferential target still conditions on the
observed series and their membership patterns. The score and outcome records
were processed locally and the output contains aggregate quantities only.
