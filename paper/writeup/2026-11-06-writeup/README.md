# paper/writeup/2026-11-06-writeup/

Version of the project writeup dated November 6, 2026. **Rendered for review.**
The previous version is [`../2026-10-16-writeup/`](../2026-10-16-writeup/README.md), which is not
edited.

This version starts from a new outline, not from the October 16 source. Its centerpiece is the
survival-weighted YFB fit (`fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)`, 2,218 D4 genes, single
log). It leads with the comparison against the same factorization at rank 7 followed by Cox,
reports the pre-registered October 15 null result, and lists every use of the external cohorts
(GitHub issue #6 items 2, 4, 5 and 6). The section plan and numerical sources are below. The [initial plan](_process/plan.md) records the starting scope.

## Reproducibility

The [source](writeup.qmd) evaluates the saved aggregate CSVs and stops if a
required input is missing. The registered October 15 comparison is read
from the CSV at commit `65cf90a`. The aggregate PACA-AU donor audit
transcribes the October 16 [accuracy dispositions](../2026-10-16-writeup/_process/Accuracy_Review_Dispositions_10_16_26.md).
The [number and figure ledger](NUMBER-SOURCES.md) identifies the input files
for every table, figure and result in the text.

From the repository root, with Quarto, R, knitr and a working LuaLaTeX
installation, run:

```bash
quarto render paper/writeup/2026-11-06-writeup/writeup.qmd --to html
quarto render paper/writeup/2026-11-06-writeup/writeup.qmd --to pdf
```

The PDF reuses the Chapter 4 vector figures directly. The HTML embeds SVG
conversions of the same files. To regenerate the SVGs after a Chapter 4
figure changes, run the command below from the repository root.

```bash
for svg in paper/writeup/2026-11-06-writeup/figures/sup_*.svg; do
  stem=$(basename "$svg" .svg)
  pdftocairo -svg "paper/prelim/figures/$stem.pdf" "$svg"
done
```

The writeup needs the committed aggregate evidence and figures. Rendering
does not need patient records, fitted objects or a new analysis. The
[HTML](writeup.html) and [PDF](writeup.pdf) are ready for review. The October
16 version remains the circulated record.

## Section plan and numerical sources

This table preserves the scope of the original section plan. The numerical ledger above reflects the rendered version.

| § | Section | Status | Content | Evidence for every number |
|---|---|---|---|---|
| — | Summary | change; draft complete | Aim, model, lead comparison (supervised vs rank-7 two-step), the October 15 null result, what supervision does and does not do on real data, simulation near the threshold | as in §3.1, §3.2, §3.4, §3.5 |
| 1.1 | Notation | carry, edited | Drop multimodal symbols except where proposed work needs them; keep L/F/β, projection scores | none |
| 1.2 | Model | new | Y = LFᵀ + E; η = (Y E[F])β; survival weight w_S multiplying the Cox term in the F update, β update and objective; canonicalization (sd(YF_k) = 1); point-Laplace L and F, Normal β | `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.pdf`; DECISIONS 2026-11-04 |
| 1.3 | Fitting and prediction | new | Cox working approximation re-expanded at E[η]; stopping rule (relative ELBO 10⁻⁵, 300 sweeps); per-cohort z-scores for new cohorts; frozen orientation; no post-fit sign flip | `config/globals.yml`; DECISIONS 2026-11-06 (sign flip removed) |
| 1.4 | The two-step comparator | new | w_S = 0 loadings, then an unpenalized Cox fit on the training program scores; why w_S = 0 alone gives β = 0 | DECISIONS 2026-11-04 (correction) |
| 1.5 | Choosing the survival weight | new | S-11 fixed w_S = 1 by training CV; the survival-block weight (A1) and CV partial log-likelihood rule from S-13; S-15 applies that rule to the real data | S-11, S-13 findings; `docs/plans/2026-11-06-supervised-real-refit-plan.md` |
| 1.6 | Earlier configuration (α_F = 0) | change, condensed | One paragraph: what the October 16 benchmark was and why it is no longer the method | DECISIONS 2026-10-02, 2026-11-04 |
| 2.1 | Data | change | Training TCGA + CPTAC (273, 139 deaths), single log, D4 (2,218 genes); five external series with n and deaths; 616 records, 568 people; 47 shared PACA-AU donors | `results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/external_cindex.csv` (n, events, missing genes); `single_log_check/gene_sets.csv`; `paca-au-donor-audit.csv` (October 16 accuracy dispositions) |
| 2.2 | Uses of the external cohorts | new (issue #6 item 6) | One table listing every choice the external cohorts informed: K = 7, preprocessing variant, α_F, N_frozen, comparator configuration, plus later scorings (October 15 check, October 26 single-log check, S-12) and the rule each followed; statement that no cohort is untouched | `docs/findings/2026-11-06-external-cohort-uses.md` |
| 2.3 | Simulation designs | new | S-9 (50 and 2,000 genes, planted variance 0.25%, 1%, 5%, weight grid), S-13 (weight chosen by CV, arms A1/A2), S-14 (near-threshold shares 0.05%, 0.1%, 0.2%; strong and realistic survival) | `yfb_survival_weight_2026-11-04/grid/design.csv`; S-13 and S-14 plans |
| 2.4 | Metrics and uncertainty | change | Harrell's C with training-fixed direction; paired bootstrap within cohort; cohort-stratified pooled difference; calibration slope with exact event times; record-level and donor-level resampling | S-12 script and findings |
| 3.1 | Supervised fit against the rank-7 two-step fit | new; lead result (issue #6 item 5) | Training CV 0.663 vs 0.652 (paired +0.011, SE 0.005); external mean 0.620 vs (two-step mean); per-cohort table; pooled +0.007 (−0.002, 0.017) and without PACA-AU seq; pooled HR per SD; calibration slopes | `yfb_supervised_real_2026-11-04/summary.csv`, `fits.csv`; `yfb_supervised_external_2026-11-05/external_cindex.csv`, `external_paired_diff.csv`, `external_paired_diff_donor.csv`, `external_pooled_hr.csv` |
| 3.2 | The pre-registered October 15 supervision check | new (issue #6 item 4) | Registered estimand and margin; ΔC = −0.003 (−0.029, 0.024) misses −0.02; rerun after the KL fix −0.006 (−0.032, 0.021); simulation recovery 17/20 and 19/20 | Git `65cf90a` and current `supervised_real_check/delta_c.csv`; DECISIONS 2026-10-15 (lines about 1041–1100) |
| 3.3 | What survival does to the real-data loadings | new | Survival share of loading precision 10⁻⁸ to 10⁻⁶; w_S = 1 and w_S = 0 loadings agree; gain comes from the β update; feature-set size curve; survival screening; null check | `yfb_supervised_real_2026-11-04/summary.csv`, `screen_paired.csv`, `null_cells.csv`, `size_agreement.csv` |
| 3.4 | Simulation: supervision strength and the weight choice | new; Ch4 figure | S-9 recovery, convergence, attenuation (slopes 1.45–3.23); S-13 CV-chosen weights 1–10, recovery equal to two-step, verdict | `grid/gate_cells.csv`, `cv_chosen_cells.csv`, `cal_slope_exact.csv`, `summary.csv`; `yfb_supervision_strength_2026-11-05/selected_vs_twostep.csv`, `verdict.csv`; Chapter 4 figure F1 |
| 3.5 | Simulation: programs near the detection threshold | new; Ch4 figure | S-14 recovery and paired C by cell; contribution of false-active programs | `yfb_near_threshold_2026-11-05/metrics.csv`, `false_active_contribution.csv`; S-13 `false_active_contribution.csv` |
| 3.6 | Real-data refit at the CV-chosen weight (S-15) | new | Weight curve, w*, nested-CV C, decision under the registered rule; external scoring only if adopted | `yfb_supervised_refit_2026-11-06/` CSVs; S-15 findings file |
| 3.7 | Rank and ELBO pruning | new | KL-2 recovery, returned fit convergence, false activity and retained rank | `flat_beta_pruning_rerun_2026-11-06/summary.csv` |
| 3.8 | The α_F = 0 benchmark under the single log | change (issue #6 item 2) | Single-log refit external C per cohort and mean 0.625 vs frozen 0.627; risk-score correlations; training CV. The October 40-factor LASSO and 20-factor comparators were run only under the double log and are listed as historical, not rerun | `kl_fix_rerun_2026-10-30/external_results.csv`, `cv_summary.csv` |
| 4 | Programs of the supervised fit | change; Ch4 | Seven programs; three survival-active (β, SD, top genes per pole); correspondence with w_S = 0 and with the α_F = 0 programs; enrichment and subtype concordance from Chapter 4 (B1–B4) | `programs_D4.csv`; Chapter 4 T4, T5, F7, F8 |
| 4.x | KM curves by risk tertile | Ch4, if produced | Within-series risk and program tertiles from the frozen scores | Chapter 4 F5 and F6 |
| 5.1 | Interpretation | new | What the supervised fit does (β shrinkage, calibration, small C gain) and does not do (move real loadings); where supervision helps in simulation | as above |
| 5.2 | Limitations | change | K = 7 informed by external C; shared donors and the series-level estimand; no untouched cohort; fold fits at the sweep cap; attenuation in simulation; one generator | as above |
| 5.3 | Continuing work | change | Data-driven supervision strength; M-modality model (proposed, no October numbers); rank selection; comparators from issue #6 item 7 (PurIST, DeSurv, elastic-net Cox, clinical); donor-aware inference and one held-back cohort | DECISIONS 2026-11-04, 2026-11-02; issue #6 |
| A | Model derivation | change | Short summary of the algorithm with a pointer to the full derivation PDF and `verify_derivation.R`; the October 16 Appendix A (α_F = 0 algorithm) is not reproduced | derivation PDF |
| B | Changes since October 16 | new | Double log removed; KL normalizer fix; orientation fix in two comparator scripts; post-fit sign flip removed; exact event times in simulation slopes; multimodal moved to proposed work | DECISIONS 2026-10-30 to 2026-11-06; progress book 11/6 issue table |
| C | Multimodal results of October 16 | carry by reference | One paragraph pointing to the circulated version's multimodal sections; no numbers repeated | `../../2026-10-16-writeup/` |

## Manuscript-planning extension

The expanded source reports the fixed prognostic comparisons, individual
program evidence, and the training-cutpoint risk groups. Appendix D gives
the per-series calibration and risk counts. The two new plots are generated
by `ggplot2` from saved aggregate CSVs at render time. The rebuild commands
above also apply to the expanded source.
