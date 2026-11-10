# README cleanup and content retention

Date: 2026-11-10. Author: Andrew Walther.

The entry page now introduces the model, shows how to fit and predict, summarizes
selected findings and links to the research record. The detailed study instructions
remain available in the guides below. This is a documentation change, with no
change to fitter code, defaults, saved results or scientific decisions.

## Examples and writing approach

Both examples requested by the week plan were read in full through GitHub:

- [DeSurv package README](https://github.com/rashidlab/DeSurv/blob/main/README.md),
  file blob `d3819d4c83d9df80040558a2d1ffc68048309e2d`.
- [DeSurv paper README](https://github.com/rashidlab/DeSurv-paper/blob/public/README.md),
  file blob `e94ab31c9c36e9275dc631852ba42d9a56c6c736`.

The package README provides a short introduction, installation and an executable
example. The paper README distinguishes rendering saved results from recomputing
analyses and gives detailed provenance. The revised entry page uses the first
structure and links to the local results guide for the second. The human-writer
skill was applied while drafting and reviewing the prose.

## Where the original material lives

The original entry page is available at `af6365fe:README.md`. Each section or
separable topic is accounted for here. A retained scientific claim is not a new
endorsement of a historical model.

| Original section or topic | Disposition and retained location |
|---|---|
| Overview, factorization and Cox models, L/F/β/τ definitions | Condensed in the entry page. Full formulas and variational updates remain in the [derivations guide](../../derivations/README.md) and its linked model derivations. |
| LB versus YFB prediction | Kept briefly in the entry page. Fitter/predictor interfaces remain in the [code guide](../../code/README.md). |
| CAVI, EBNM priors, Cox Taylor approximation and second moments | Brief inference description retained. The full explanation remains in the [YFB derivation](../../derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex) and [LB derivation](../../derivations/SSBMF-LB-model-derivations/MF_Derivations_UpdateAlgo_REVISED.tex). |
| Historical predictive benchmark, 2,064 genes, two prognostic programs, C = 0.627, α_F = 0, 616 records and repeated donors | Removed as an entry-page section. A historical summary is retained in the [results guide](../../results/README.md#historical-benchmark-and-workflow-figures). Methods and estimates already appear in [PROJECT_STATUS.qmd](../../PROJECT_STATUS.qmd) and the October writeup. |
| Historical external KM figure and within-cohort tertile caption | Original figure remains at [km_risk_external.png](../../paper/writeup/2026-10-16-writeup/figures/km_risk_external.png) and in the [October writeup](../../paper/writeup/2026-10-16-writeup/writeup.qmd). The results guide links to that version rather than displaying the historical figure as a current result. |
| Multimodal prototype, matched TCGA/ICGC data, prediction/simulation/rank utilities and development status | Condensed to one sentence in the entry page. Details remain in the [code guide](../../code/README.md#multimodal-yfb-expression--methylation-in-progress), [derivations guide](../../derivations/README.md) and current [week plan](../plans/2026-11-09-week-plan.md). |
| Project history, choices and advisor progress | Entry-page links retained to PROJECT_STATUS.qmd, DECISIONS.md, ROADMAP.md and the published progress notebook. |
| October writeup workflow diagrams and rebuild command | Confirmed in the [October writeup README](../../paper/writeup/2026-10-16-writeup/README.md#rebuilding-figures). Added a results-guide pointer. The circulated version was not edited. |
| Accuracy-review dispositions and remaining inference/provenance questions | Confirmed in the [accuracy dispositions](../../paper/writeup/2026-10-16-writeup/_process/Accuracy_Review_Dispositions_10_16_26.md) and October writeup README. Retained a results-guide pointer. |
| Prelim Chapter 4 exhibit commands, required fits/local data, external cache and output locations | Moved the three commands and cache instructions to the [Chapter 4 source ledger](../../paper/prelim/notes/chapter4-number-sources-2026-11-06.md#rebuild-the-figures-and-tables). One entry-page link remains. |
| Installation, core R requirements, root-relative execution | Kept in the entry page. |
| Optional analysis packages and Bioconductor installation | Moved the package table to the [code guide](../../code/README.md#dependencies). Clarified that Rcpp accelerates the canonical loading sweep and that a base R fallback exists. |
| Minimal fit-and-predict example | Retained and updated to explicitly select canonical YFB with point-Laplace molecular priors and Normal β. Standardization now uses training means/SDs for the held-out matrix. The convergence field is `fit$history$converged`. |
| Benchmark commands and output/config paths | Moved to the [results guide](../../results/README.md#historical-benchmark-and-workflow-figures), labeled historical. |
| Full test-suite command | Kept in the entry page. No production code or test file changed. |
| Expression cohort table, platforms, sample sizes, censoring, data path and matched-multiomics loader | Moved intact to the [results guide](../../results/README.md#expression-cohort-reference), with a note distinguishing analysis records from unique donors. |
| Repository structure | Condensed to a linked table in the entry page. Detailed derivation folder descriptions already live in [derivations/README.md](../../derivations/README.md). |
| Experimental normalized-projection fitter, inputs/outputs and example call | Moved interface details to the [code guide](../../code/README.md#experimental-normalized-projection-interface). Scientific failure and non-adoption remain explicit there and in the [execution report](../reviews/2026-10-30-codex-execution-report.md). |
| Experimental explicit variances VBeta/VL/VF and avoiding cancellation | Moved to the code guide, including the recommendation to compute SD as `sqrt(fit$VBeta)`. |
| Seventh shared-loading simulation, earlier numerical attempts and simulation_units.csv | Retained in the [results guide](../../results/README.md#historical-supervision-diagnostics) and existing [output README](../../results/benchmark_sim/outputs/yfb_supervision_2026-10-30/README.md). |
| Registered experimental run/summarize commands and saved-fit replay order | Moved all five commands to the results guide. Registered revision, checkpoint reuse and patient-cache restrictions are retained. |
| Stratified YFB CV, fold-specific centers, within-stratum comparisons and unseen-stratum example | Moved to the [code guide](../../code/README.md#stratified-yfb-cross-validation). The [handoff](../reviews/2026-11-06-stratified-yfb-cv-handoff.md) still records the independent review as open. No review of commit 68e3646 was written. |
| Derivation and writeup links | Consolidated into the linked repository table. The derivations guide links all four model/comparison folders. The writeup index distinguishes the circulated October version from the November version. |
| License and author | Retained in the entry page. |
| Targeted noise/guard controls, rollback and stalled-update failure | Moved to the code guide. No scientific default was adopted. |
| Targeted screen and coefficient-reference commands, checkpoint restrictions and limits | Moved commands to the results guide. Existing [targeted findings](2026-10-30-yfb-targeted-diagnostics-findings.md) and [output README](../../results/benchmark_sim/outputs/yfb_targeted_2026-10-30/README.md) preserve the full inventory and fixed-Gaussian, rather than full-Cox, interpretation. |
| Tracker's model states, diagrams, Completed Results checklist and October methods record | Retained in the [results guide](../../results/README.md#progress-book-exhibits), pointing to the existing [current-state chapter](../progress_book/current-state.qmd) and progress book. Diagram and Quarto rebuild commands retained. |
| November 6 saved-fit and planted-factor figure builders | Moved all three commands and their local-data/saved-table requirements to the results guide. |
| Historical risk-direction rescoring and donor-aware external interval commands | Moved to the [results guide](../../results/README.md#external-orientation-and-donor-intervals). [Orientation](2026-11-06-orientation-audit.md) and [external-cohort](2026-11-06-external-cohort-uses.md) audits already preserve the conventions and evidence. |
| Flat β likelihood behavior and regression command | Moved to the [results guide](../../results/README.md#flat-coefficient-likelihood-regression), linked to the existing rerun registration. |
| Fixed comparator commands, pinned DeSurv and local checkpoint restrictions | Confirmed in [comparator findings](2026-11-06-real-data-comparators-findings.md#reproduce-and-trace-the-results), which already contains both commands and output definitions. Added a results-guide pointer. |
| Risk-group/program commands and their aggregate outputs | Confirmed in [program-analysis findings](2026-11-06-desurv-style-analyses-findings.md), which already contains all four commands. Added a results-guide pointer. |
| November writeup's calibration, stability, fixed-cutpoint and simulation synthesis | Confirmed in its existing [README](../../paper/writeup/2026-11-06-writeup/README.md) and numerical ledger. Added a results-guide pointer. |

The three AGENTS.md pointers that referred to commands in the old entry page
now point to their retained guides. Other repository instructions and circulated
writeup versions were left unchanged.

## Fit-and-predict example

**Purpose.** The example demonstrates the current fit/predict interface with
synthetic data and needs no patient files.

**Logic.** It generates a 150-by-200 molecular matrix from three underlying
programs, creates survival times and events, and reserves the last 50 patients.
`scale()` estimates feature means and SDs on the first 100 patients and applies
those same values to the held-out rows. The fitter estimates five program
columns. Prediction multiplies held-out expression by the fitted loadings and
survival coefficients.

**Choices.** Explicit `canonicalize = TRUE`, `w_S = 1`, point-Laplace molecular
priors and a Normal coefficient prior avoid reliance on the historical fitter
defaults. These are illustrative settings, not a replacement for the registered
research design. Training-only scaling avoids using held-out measurements to
prepare the training matrix.

**Inputs and outputs.** The fitter takes a numeric patient-by-feature matrix
and matching numeric survival-time and binary event vectors. It returns
posterior moments and iteration history. The loadings have 200 rows and five
columns, β has length five, and the prediction vector has length 50.

**Limits.** The synthetic times are generated from latent activities. This
example checks the API, not recovery, calibration or scientific superiority of
YFB. Inspect `fit$history$converged` before interpreting a fit. Real analyses
must follow their registered preprocessing, cohort and validation rules.

**R functions.** Base R generates and scales the data. `ebnm` supplies empirical
Bayes normal-means updates, and `survival` supplies survival-analysis machinery.
The fitter can use its compiled Rcpp sweep or the existing base R fallback.
No dependency or production function was added.

## Verification

The entry page has 107 lines and about 580 words, compared with 381 lines
before cleanup. It contains no internal study identifiers. Local Markdown links,
retained script paths and removal destinations were checked. The example was
executed as a documentation smoke check and compared against direct Y_test F β
projection. It converged in 31 iterations and returned 50 finite, nonconstant
held-out scores. The direct-projection identity passed. The full R suite was
not rerun for this documentation-only change. Independent review checked the
plan, full diff, retention ledger, quick-start API, selected findings, local links
and Markdown anchors. It found no actionable issues. No patient data or research analysis was run.
