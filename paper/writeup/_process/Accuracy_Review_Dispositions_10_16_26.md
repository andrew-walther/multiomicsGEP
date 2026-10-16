# Writeup accuracy-review dispositions

This is the current item-level tracker for [the original review](Accuracy_Review_10_16_26.md).
Andrew requested applying its findings on October16,2026, after the cohort
explanations and workflow figures were added. The original review remains a
dated snapshot; this file records what was changed and what still needs a
separate scientific decision. Reporting a limitation is not the same as
completing the analysis that would resolve it.

Scope: manuscript, derivation appendices, bibliography, evidence inventory
and current documentation. Preserve model code/defaults, frozen numerical
outputs and historical failure records. No new model fit, external scoring,
bootstrap or scientific feature/model selection was performed.

## Item-by-item disposition

| Item | Disposition | Change / evidence / remaining boundary |
|---|---|---|
| C1 | Corrected previously; retained | CPTAC is RNA-seq in Methods and Table2. Both training cohorts are RNA-seq with separate upstream processing. |
| C2 | Corrected | Active cohort audit already carries the dated RNA-seq correction. Inventory S03 now agrees; dated historical records remain historical. The new decision/status entry qualifies the earlier blanket preprocessing assurance. |
| C3 | Corrected previously; retained | Table2 reports TCGA144/75, CPTAC129/64 and pooled273/139 separately. |
| C4 | Corrected previously; retained | ICGC/PACA-AU is one validation source. The67-donor set explicitly contains the primary50. Figures/prose preserve primary PDAC as the target population. |
| C5 | Corrected | Main text and Appendix C distinguish144 training patients/75 deaths from50 primary validation donors/30 deaths. Appendix C records exclusions and the supplied-collection boundary;194 is not a training sample size. |
| C6 | Corrected | Summary, performance Methods, external-results prose/table, KM captions and Limitations use cohort/death records for616/414. Unique-person counts and donor assumptions are explicit. |
| C7 | Documented; alternative analysis open | The internal loader applies log2(x+1), then downstream flags repeat it for TCGA/Dijk/PACA RNA-seq. A read-only reconstruction exactly reproduces all273x2064 prepared training values and frozen gene order. This verifies present code/cache identity, not all historical acquisitions or the scientific impact of a single-log alternative. No transformation changed. |
| C8 | Corrected | Data Methods, workflow caption and Appendix C distinguish shared assay availability from TCGA-only ranking and explicitly identify within-cohort imputed methylation. |
| V1 | Reporting corrected; donor-aware inference open | Methods, Results, Table caption and Discussion label the saved record-level CIs as unclustered historical summaries; they do not establish donor-aware superiority. No replacement CI is invented. |
| V2 | Reporting corrected; donor-aware inference open | Program-curve caption and biological Methods state that cohort stratification does not make repeated donors independent. P-values remain historical/descriptive. The repeated RNA donor also affects that series' record-level inference. |
| V3 | Corrected | Verified47-donor cross-series overlap and the52-record/51-donor RNA series replace the “not determined” statement. Appendix C records the conditional568-person/387-death total. |
| V4 | Corrected | Summary now says conditional training CV showed no consistent gain through17900 genes, without implying expanded-feature external validation. |
| V5 | Corrected | Methods/Table9 identify the fold SE as descriptive and retain the fixed-before-fold preprocessing caveat and overlapping training-set dependence. |
| V6 | Reporting corrected; formal increment open | Saved unconverged.685/.676 estimates remain; the difference is explicitly.0095, with no paired interval. Overlapping marginal intervals are not a test. No claim of methylation increment/equivalence. |
| V7 | Corrected | Historical constant-risk C=.500 is retained as no discrimination and distinguished from the later prospective NA convention. No value is silently replaced. |
| V8 | Corrected | Starting-rank development use of external concordance remains visible in Methods, workflow caption and Limitations. No claim of an untouched final external test. |
| M1 | Corrected | Methods/Appendix B.9 distinguish the mathematical beta-prior pushforward from actual carried state. Beta moments and L/F priors rescale; per-factor Normal beta prior variance is not carried as an object. Cached penalties are not certified as current. |
| M2 | Corrected | Discussion acknowledges the single fixed-q stopping-point observation and its failed output-stability check; general effects on rank/generalization remain unknown. |
| M3 | Corrected | Discussion limits all-fit convergence to the162 evaluated-framework fits in the historical sweep. |
| M4 | Corrected | Methods/Limitations explicitly identify three independent datasets per scenario/noise combination refitted at nine ranks;27 fits are not27 independent datasets. |
| M5 | Corrected | Recovery definitions describe independent best-absolute-correlation matching, possible reuse of estimates, and the code's raw coefficient-sign comparison. No matching metric was changed. |
| M6 | Corrected | Limitations explain nonnegative reconstructed means for signed standardized input and the absence of a single-block feature intercept; no cell-population decomposition is claimed. |
| M7 | Corrected | Comparator described as a40-factor greedy cap plus CV-tuned LASSO Cox, with performance attributed to the whole procedures rather than uniquely to shrinkage/supervision. |
| M8 | Corrected | Main text/Appendix B.11 distinguish terminal component levels, deletion changes and coordinate likelihood fractions. The saved joint components are reported without inferring a deletion mechanism from observation counts alone. |
| B1 | Corrected | Single-program headings use associated-expression terminology; joint labels remain provisional. No cell identity or mechanism is established by these labels. |
| B2 | Corrected | First methylation-share claim, caption and plot axis explicitly use normalized program-wise rank-one squared norms, not additive PVE or survival contribution. |
| B3 | Corrected | The14 unadjusted marginal HRs are exploratory joint-projection associations, distinct from conditional coefficients and methylation-specific incremental value. |
| E1 | References added; acquisition/version gaps explicit | A local bibliography cites cohorts, DeSurv, EBMF/ebnm/flashier, Cox/Breslow/C, fgsea/resources, PurIST, imputation and BH. Appendix C supplies public study IDs and saved software versions. It does not invent missing original-fit environments, exact MSigDB releases or assay-specific acquisition manifests. DeSurv is unpublished author material; no journal/DOI status inferred. |
| E2 | Corrected | S03/S08/S10-S13/S39 reconciled; S74-S77 add corrected-observer/fixed-q evidence, donor accounting, transform/cache audit and bibliography/provenance. Historical untraced established-model.685 is distinct from traced joint-framework.685207. |
| E3 | Corrected | Full-training C in feature-count tables rounds to three decimals. Table2 already clarifies totals and separate cohort deaths; its explicit record/death-record wording requested by Andrew is retained. |
| E4 | Corrected previously; refreshed | `writeup.pdf`/`writeup.html` remain canonical. Annotated pre-edit PDF is preserved under ignored data/writeup_review_backups/; dated circulation copies are not overwritten. |
| E5 | Corrected | Detailed25k arithmetic/agreement/balance/filter and fixed-q accounting moved to Appendix C. Main Results retain failures, unconverged states and no-adoption conclusions with cross-references. |

## Read-only preprocessing verification

The current `load_pdac_raw()` calls the supplied internal loader, which logs
TCGA once. `preprocess_merged_cohorts()` then applies the configured log flag.
The audited training call used the historical settings: combined-rank3000
per cohort, intersection, per-cohort standardization, no rank/quantile transform.
The resulting screened gene order is identical to the frozen benchmark list;
its273x2064 matrix is exactly equal to the baseline columns in ignored
`data/round2_feature_cold_primary.rds` (maximum absolute difference0).
This audit read molecular inputs and produced only aggregate verification;
no patient identifiers, matrix values, outcomes or fitted state are exported.
No production runner was executed and no alternative was fitted.

This new observation qualifies, rather than erases, the earlier dated
statement that the CPTAC platform-label correction implied correct overall
preprocessing. Archive-wide acquisition reconciliation remains open.

## Concrete decisions needed before analytical follow-up

1. **Repeated-donor inference.** Choose the target estimand before computing
   replacement intervals or p-values. Keeping all historical records with
   donor-cluster resampling retains cross-cohort pairs and the original
   record-weighted estimand; selecting one prespecified assay/sample per
   donor changes the population weights and point estimates. Neither is a
   wording-only correction. The selected design must specify within-series
   repeats, cohort stratification, model-paired resampling and the pooled
   program test. Per-series reports can remain separate but do not become
   independent studies by being listed in separate rows.
2. **Preprocessing alternative.** A single-log pipeline needs a prospective
   design, recomputed training selection/scaling, fitting, and appropriate
   validation under the existing scientific restrictions. Do not use external
   outcomes to select the transformation or attach frozen estimates to new
   inputs. The present cache audit is not such a comparison.
3. **Multimodal increment and provenance.** A paired interval should follow
   a defined comparison of stable, selected models. The current two capped
   fits cannot certify reliable added value. Retrieve original acquisition,
   gene-set release and older session manifests when available; current
   installed versions are not substitutes.

These are explicit outstanding scientific/provenance requirements, not
authorization to run prepared panels, change defaults or relax criteria.

## Source verification for references

Existing cohort/method entries were reused from `paper/ssbmf-refs.bib`, with
internal audit notes omitted from the rendered bibliography. Additional
resource/method entries were checked against primary sources:

- [fgsea version3 preprint](https://www.biorxiv.org/content/10.1101/060012v3)
- [Hallmark source article](https://doi.org/10.1016/j.cels.2015.12.004)
- [Reactome citation page](https://reactome.org/cite)
- [KEGG source publications](https://www.genome.jp/kegg/kegg1.html)
- [GO citation policy](https://www.geneontology.org/page/use-and-license)
- [BH source article](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [ebnm source article](https://www.jstatsoft.org/article/view/v114i03)
- [DeSurv author repository](https://github.com/rashidlab/DeSurv-paper)

The author repository's journal label is not independent publication
verification. Study accession URLs follow the already-pinned author data
inventory; EGA confirms the parent pancreatic study identity. Some live
BioStudies pages expose only a client-side shell and GEO imposed a browser
check, so complete assay/version provenance is explicitly not claimed.

## Verification record

- The full repository test suite passes: **710 passed, 0 failed** (final log: `/tmp/writeup-accuracy-final-tests.log`).
- HTML and PDF render successfully; the PDF contains 51 pages. Rendered tables, workflow diagrams, revised modality-share figure, mathematical appendix, diagnostic appendix and references were visually checked.
- All HTML fragment links resolve and all 24 cited bibliography entries appear. All 32 review IDs have dispositions; the source inventory has 77 unique IDs.
- External concordance estimates and intervals agree with their saved source table at displayed precision. Diagnostic results remain preserved in Appendix C.
- The exact current historical expression preprocessing reconstruction matches the frozen 273-by-2064 training matrix, including gene order (maximum absolute difference zero). This is a cache audit, not a preprocessing comparison or a new fit.
- All 10 immutable benchmark manifest entries pass SHA-256 verification. No model code, configuration, test or benchmark snapshot was changed. No new model fit or external score was produced.

The sole plotting-code change replaces the modality-share axis label in
`figures/make_writeup_figures.R`. The existing `labs()` call receives a more
precise text label; it changes neither data nor calculations. The figure was
regenerated from the saved seven-program aggregate CSV, and the plotted
shares were checked against that table. Its output remains the existing
PDF/PNG pair. The label describes normalized rank-one squared norms; it must
not be interpreted as additive explained variance or survival importance.
