# Verification of the second writeup revision

Reviewed on 2026-10-14 at `fbf8317`, on `writeup/draft`. The checkout was clean at the start. This review compares `8772b9b..fbf8317`, the previous `verification.md`, the claimed dispositions, the revised sources, and the delivered HTML and PDF. Only this verification file was written; no report, figure, inventory, code, default, data or saved analytical result was changed. No fitting, commit, push or merge was performed.

## Verdict

**The requested scientific corrections are verified.** The added algebra is correct, the algorithms now fit on single pages, and the exhibit abbreviation checklist passes. No new substantive scientific error was found. Two minor cleanup items are recorded below; neither overturns the corrected mathematics or reported results.

Two requirements remain distinct from that verdict. The main text is **24 pages**, exceeding the 15–20-page target. Also, the program-curve producer still needs local expression data: its provenance disclosure is correct, but a literal requirement to regenerate every figure exclusively from saved fits, scores and aggregates remains incomplete. The previously disclosed activity-KL and cached-prior-KL implementation issues remain method work, not fixes made on this writing branch.

## V1–V13

Locations below refer to the current sources. Each corrected scientific statement was also located in the rendered HTML and PDF, rather than accepted from the disposition alone.

| Finding | Verdict | Verification and evidence |
|---|---|---|
| **V1 — Multimodal precision floor** | **Verified** | B.8 (`../appendices/multimodal.qmd:164–176`, PDF p.37) derives the unrestricted optimum for positive residual sum, explains the unbounded zero-residual likelihood, and gives the evaluated protected computation, `n / max(Q, epsilon_machine)`. The actual fitting path at `code/fit_multimodal_yfb.R:505–509` calls the helper at `:639–645`; the removal score uses it too (`:699–701`). The text no longer attributes the unused update helper's error branch to this fitter. |
| **V2 — Zero-column normalized target** | **Verified** | Main §1.2 (`../writeup.qmd:148`, PDF p.3) and A.1 (`single_modality.qmd:29`, p.25) define the normalized direction as zero when the loading column is zero. They distinguish this mathematical convention from the numerical ridge in the evaluated algorithm and state that the Gaussian coordinate derivation additionally freezes the normalizer. This resolves the positive-probability spike-prior boundary. Minor consistency item N1 below concerns the first displayed sweep formula, not the definition of this target. |
| **V3 — Risk orientation** | **Verified** | A.10 (`single_modality.qmd:287`, p.31) now states the evaluated fit's convention only: higher eta denotes longer survival, and validation uses risk = minus eta. It no longer suggests the current risk-oriented check left that convention in place. Current `fit_cox_on_yf.R:885–906` uses risk-oriented concordance and a flip below 0.5. The saved fit has zero recoverable coefficient variance; reversing its coefficients gives the reported adverse/protective signs. Main §1.2.1 and Figure 2 agree. |
| **V4 — Exact target versus working second moment** | **Verified** | A.5 (`single_modality.qmd:164`, p.28) explicitly assigns the E(Z squared) coefficient to the random-projection working quadratic and calls the exact-target update nonquadratic and nonconjugate. Table 10/A.11 (`:310`, p.32) makes the same distinction. The evaluated plug-in update uses the squared fixed projection. |
| **V5 — Missing derivation steps** | **Verified** | A.5 now matches all coefficients and normalizing constants in the Gaussian product and integrates out the normalized posterior factor. B.6 explicitly expands the reconstruction residual square; B.8 differentiates in precision. B.5 states the prior-only boundary for B.5–B.7, agreeing with activity/loading prior moments at `multimodal_yfb_updates.R:36–49,181–187` and the coefficient branch at `:520–527`. The single-modality precision floor in A.3 is distinguished from these multimodal branches. Detailed algebra checks follow below. |
| **V6 — Algorithm/table pagination** | **Verified** | PDF page images show Algorithm 1 entirely on p.4 and Algorithm 2 entirely on p.7, including inputs, all steps, outputs and prediction. No boxes remain. All 11 tables have their caption, header and final row on one page; no duplicate header precedes a caption. A.11 and B.13 headings now stay with their tables. See the layout audit below. |
| **V7 — Local caption definitions** | **Verified** | Rechecked all 11 tables and five figures in both formats using their new numbering. The formerly missing cohort/assay names, CI, LASSO, CV, PL, SD, KL, model-arm labels and normalization exceptions are supplied locally. The deleted single-modality program table no longer needs a caption. Full checklist below. |
| **V8 — Revisit gaps** | **Verified** | New R1 blocks surround A.3's single-modality floor (`single_modality.qmd:114–118`) and A.4's point-exponential posterior (`:122–141`). All three sources have balanced, properly nested markers: 58 blocks in the main source, 16 in Appendix A and eight in Appendix B. The moved matched-donor Table 7 and multimodal-program Table 9, including their captions, remain enclosed. |
| **V9 — Normalization exceptions** | **Verified** | Table 1 uses the stored D-star; §1.2 avoids literal unit-length wording; Figure 4 uses the stored training norm. The multimodal scale paragraph, Algorithm 2, B.9 and Table 11/B.13 all retain the constant-projection exception. B.9's threshold is square-root machine epsilon, matching `fit_multimodal_yfb.R:190–200`. Moment/prior transformations match `:199–214`. See N1 for a small remaining first-definition inconsistency. |
| **V10 — Classification/removal boundaries** | **Verified** | §1.5 and B.7 protect the coefficient SD denominator by 10 to the minus 12, matching `fit_multimodal_yfb.R:548–560`; Table 9 says a zero point mass is inactive. Methods, Algorithm 2, B.11 and B.13 restrict removal to K greater than 1. B.11 correctly says a deletion “does not decrease” the implemented score, including equality, matching `:586–595`. |
| **V11 — Figure provenance** | **Verified disclosure; saved-input-only requirement incomplete** | The disposition now acknowledges that `../figures/make_writeup_figures.R:95–109` loads local held-out expression to rebuild frozen projections. The script header, README command and S63 disclose that dependency. The figure producers do not fit factor models or Cox regression models. No patient-level cache was introduced. This is a transparent remaining dependency, not fulfillment of the narrower saved-fits/scores/aggregates-only condition. No report correction is required to hide or remove it. |
| **V12 — Supporting single-block comparator** | **Verified** | Comparators (`../writeup.qmd:360`, PDF pp.9–10) and Table 7 (`:550`, p.18) identify unstandardized log2 expression and the 100-sweep cap. The comparator paragraph also identifies its own frozen projections. These agree with `code/run_multimodal_yfb_simulation.R:139–149,159–161`; the real-data runner does not override that cap. It is no longer presented as interchangeable with the principal standardized, pooled fit. |
| **V13 — Wording, references and rounding** | **Verified** | Summary qualifies lower-noise mean-rank recovery and its low-variance exception. Both monitoring links name A.9 and resolve there. The matched-donor difference is “about 0.01,” consistent with the saved 0.0094675. Discussion says “no statistically established advantage” over the rank-matched comparator. The deleted table removes the definitive null label; §4.1 retains the appropriately qualified immune-signature association statement. |

### Earlier incomplete A/B items

- **A1 and A8 / B14:** the zero-column target, working-quadratic attribution and remaining algebra now resolve the earlier incompleteness.
- **A4 and A10:** stored normalization and multimodal constant-column exceptions are explained; N1 is a minor formula-consistency refinement.
- **A11:** marker gaps are closed.
- **B6 and B7:** local definitions and algorithm pagination now pass.
- **B8:** no fitting, numerical continuity, vector/PNG formats and printed readability pass; the local-expression dependency remains as V11.
- **B13:** organization and pagination pass, but the main-text length target remains unmet. No References section or Appendix C has returned.

## Added algebra, step by step

| Passage | Check |
|---|---|
| **A.1, normalized direction** | For a nonzero column, dividing by its Euclidean norm defines its direction. The separate zero-vector branch handles the spike-prior event without division by zero. Freezing this normalizer is an additional approximation needed for linear predictor coordinates; the text explicitly distinguishes it from the random target. |
| **A.5, Gaussian square** | With prior variance v and observation variance s squared, the displayed V and m satisfy 1/V = 1/s squared + 1/v and m/V = x/s squared. The residual constant is x squared/s squared minus m squared/V = x squared/(s squared + v). Thus the quadratic, linear and constant terms all match. |
| **A.5, constants and convolution** | V times (s squared + v) equals s squared times v. This proves the equality of the two products of Gaussian normalizing constants. The product identity therefore holds as a density identity, not just up to proportionality. Integrating the normalized N(b; m,V) factor gives the asserted N(x;0,s squared+v) marginal. The following constrained variance maximizer, posterior first/second moments and zero-variance point-mass boundary are correct. |
| **B.6, reconstruction loading** | Expanding (r minus l f) squared gives r squared minus 2 r l f plus l squared f squared. Factor k is excluded from r, so the stated variational independence gives the linear coefficient using mean r times mean l, and the quadratic coefficient using E(l squared). Multiplication by minus tau/2 and summation yield the stated reconstruction A and B with the correct sign and factor of two. |
| **B.8, precision** | Differentiating n/2 log tau minus tau Q/2 gives n/(2 tau) minus Q/2. For positive Q the root n/Q is a maximum (negative second derivative). At Q = 0 the likelihood is unbounded. The numerical residual floor belongs to the computation and is correctly separated from the mathematical optimum. |
| **B.5–B.7, prior-only branch** | The normal-means reduction requires A positive. At A = B = 0 there is no likelihood information and returning the prior is correct; A = 0 with materially nonzero B is an error. Activity/loading helpers tolerate a tiny numerical B remainder, which does not contradict the mathematical statement. The normal coefficient branch agrees. |

Reinspection also found no error in the unchanged point-exponential half-line completion, truncated-normal first and second moments, point-Laplace negative-half substitution and signed first moment, projection/predictor moments, scale pushforwards or removal-score/Jensen explanation. Generic identities and model-specific approximations remain separated. This review verifies the disclosed algorithms; it does not establish exact ELBO ascent or resolve the known implementation issues.

## PDF and HTML layout

The PDF has **39 pages**. I agree that the main text ends on **p.24**: the title, contents and Summary occupy p.1; Methods through Discussion occupy pp.2–24. Appendix A begins p.25 and Appendix B p.33. Excluding the first page leaves 23 pages, so counting convention alone does not bring it within target.

All 39 page images were inspected in overview, with larger images for the algorithms, shortened biology ../figures/table and Gaussian derivation. Table locations are:

| Table | PDF page | Table | PDF page |
|---|---:|---|---:|
| 1 | 2 | 7 | 18 |
| 2 | 8 | 8 | 19 |
| 3 | 12 | 9 | 22 |
| 4 | 13 | 10 | 32 |
| 5 | 15 | 11 | 39 |
| 6 | 16 | | |

No table spans pages or repeats a header before its caption. Algorithms 1/2 occupy pp.4/7 in full. No orphaned heading was found, including A.11 and B.13. Equations, captions, symbols and exhibit numbering are readable; no clipped equation was found. There is still noticeable space below the text on p.6 before the unbroken second algorithm, and below the material on p.21 before Table 9. That is acceptable for keeping these exhibits together, although it contributes to the length. The shorter final appendix pages are not a main-text length problem.

Live HTML inspection at its normal 1,265-pixel viewport found 11 tables, five loaded embedded figures, 886 rendered MathJax containers, zero MathJax error elements and no broken in-document anchors. Document width equals scroll width. Both plain algorithm lists retain their outputs. This was not a phone-width test.

### Exhibit abbreviation checklist

Each row passes in both delivered formats. Definitions in adjacent report prose were not counted as caption definitions.

| Exhibit | Local definitions/details checked |
|---|---|
| **Table 1** | expr/meth, SD, n/p/K, L/F dimensions, posterior-mean F, stored D-star, intercept/centered-data symbols; collapsed-projection exception in its row. |
| **Table 2** | TCGA-PAAD, CPTAC, RNA-seq, PDAC, ICGC/PACA-AU, IPMN, CpG, 450K and DeSurv; training deaths and matched-donor overlap caveat. |
| **Table 3** | Single-modality model arm, survival coefficient, recovery/activity definitions and two-step arm; undefined cells identified. |
| **Table 4** | C, CI, LASSO, PDAC, PACA-AU, RNA-seq; concatenated versus pooled difference and bootstrap definitions. |
| **Table 5** | BIC, CV, PL, PVE and C; n/p population and dimension, operational classifications and monitoring-statistic status. |
| **Table 6** | SD, C, recovery/activity/convergence; base versus evaluated framework and uncertainty absence. |
| **Table 7** | Expr./meth., C/CI, TCGA, ICGC, PDAC, CpG; comparator settings, event counts, bootstrap and convergence caveat. |
| **Table 8** | C, SD, TCGA, ICGC, PDAC; training/validation deaths, classification and stopping rule. |
| **Table 9** | SM/U/EO, R/S, HR, SD, TCGA, ICGC, PDAC; match definition, shares, population and unadjusted tests. Gene symbols need no acronym expansion. |
| **Table 10** | KL, variational objective, survival weights and stored norms; working-quadratic versus exact-target row. |
| **Table 11** | sd, KL and mean projection; constant-column and K greater than 1 exceptions in rows. |
| **Figure 1** | C and direction, paired standard error, model arms and scenarios. |
| **Figure 2** | C, PDAC, PACA-AU, RNA-seq; risk orientation, grouping, censoring, time and tests. |
| **Figure 3** | SD, true rank, noise multiplier, means/ranges rather than intervals, sample/feature counts. |
| **Figure 4** | PDAC; stored-norm projection, grouping, censoring, time and stratified unadjusted tests. |
| **Figure 5** | TCGA, CpG, SD; share denominator, survival-active asterisk and non-PVE/non-survival interpretation. |

## Cuts, numbers, figures and scope

The deleted single-modality program table mostly repeated §4.1. Program classifications, signs and conditional-versus-marginal interpretation remain in §§1.2.1, 1.3 and 4.1. Biological-method details remain in §2.5. The program labels, enrichment qualifications, weak/nonconsistent associations and nonindependence caveat remain. The reduced Table 1 keeps the four contrasts while shared priors and precision details remain in its caption and Methods.

Table 3 retains the rank-6/20 rows unchanged. The rank-4/12 shared-recovery values in its note agree with saved S18 outputs. **Not every deleted cell is preserved:** rank-4/12 specific-recovery and activity cells, including rank-4 nothing-shared cells, are omitted. This is the deliberate supporting-detail cut, not a numerical discrepancy or removed scientific qualification. The prose retains under-specification, false activity and outcome-oriented-null warnings. No dangling reference to the deleted table or unresolved cross-reference was found.

Condensed Limitations retain coefficient-uncertainty, approximation/KL, convergence/removal, small-cohort, overlap/rank-selection, fixed-fit-bootstrap, screening, simulation and biological-dependence/multiplicity qualifications. All five Continuing work aims remain. Discussion's corrected platform claim agrees with Table 2: three validation cohorts are microarray, while two are RNA-seq, a platform already represented in training.

The two figure-height changes are the **only calculation-script changes**: program curves 3.1 to 2.9 inches and modality shares 3.0 to 2.7, with the same 9-point plot theme. Inputs, transformations, values and labels are unchanged. The program curves still display stratified p below 0.0001 and p = 0.0041; the shares agree with Table 9. PNG dimensions are 1,950 by 870 and 1,950 by 810, respectively, with 300-dpi metadata. Their PDFs and the full report contain no raster image objects. At the printed 6.5-inch width, text, curves, censoring marks and angled program labels remain legible. S63 accurately records the dependency and height changes; no new numerical result was introduced.

The complete revision changes only `paper/writeup/**`, `README.md` and `PROJECT_STATUS.qmd`. Model code, defaults, data and saved analytical results are unchanged. The saved figure calculations independently reproduced during the previous verification remain applicable because the current script diff changes only height arguments; those calculations were not rerun here.

## Minor residual/new issues and exact fixes

1. **N1 [P3, residual consistency] Mention the ridge in the first sweep-normalizer formula.** `../writeup.qmd:148` and `../appendices/single_modality.qmd:38` write the evaluated sweep norm as the plain Euclidean norm, whereas A.10 correctly gives the protected square-root norm-squared plus 10 to the minus 10. The fitter agrees with A.10 (`fit_cox_on_yf.R:643–646`). **Fix:** define the evaluated sweep norm as `sqrt(||mean F from the preceding sweep|| squared + 10^-10)` at both first introductions, or explicitly say those displayed norms suppress the numerical ridge specified in A.10. Keep the random-target zero-direction convention separate. This does not reopen V2 or invalidate the disclosed stored-normalizer interpretation.

2. **N2 [P3, introduced documentation pointer] Correct the disposition-file pointer.** `PROJECT_STATUS.qmd:141` says “dispositions in the same file” after naming `paper/writeup/_process/verification.md`. Claude's dispositions are in `paper/writeup/_process/review_dispositions.md`, under “Verification round.” **Fix:** replace that parenthesis with the explicit latter path. The report itself is unaffected.

No new P1/P2 scientific or layout problem was identified. Length and the V11 dependency remain explicitly unresolved requirements rather than undisclosed defects.

## Remaining length options

Estimates concern the present layout, overlap with reflow, and are not guaranteed page-count reductions. No alternative report was edited or rendered to measure them.

| Option | Estimated saving | What is lost / recommendation |
|---|---|---|
| **Rank 6: remove the base, single-block and two two-step rows from Table 7** | Approximately **0.25–0.5 page** of content, including their comparator-specific caption text; possibly zero or one whole page after reflow. | Loses the compact matched-population comparison showing why the framework differs from the base model and supporting pipelines. The two retained framework rows still answer whether adding methylation helps. **Prefer keeping these rows for the advisor-facing report:** the savings are small relative to their comparative value. If cut, retain a short qualified supporting-comparator sentence and the distinct input/cap disclosure wherever the single-block number remains. |
| **Rank 8: move more Methods detail to Appendices A/B** | Approximately **1–2 main-text pages**, depending on how much algorithm detail moves and how pages reflow. | Move the detailed working-curvature/square-completion exposition and the full removal-score expression/old-prior bookkeeping to the existing appendix discussions. Shorten algorithm steps by referring to appendix update equations rather than repeating A/B formulas. Keep each model's reconstruction, priors, predictor, survival-feedback boundary, prediction, classification, evaluated settings and a short approximation warning in the main text, with both algorithms still unbroken. **Prefer this option first** if the page target is firm. It reduces immediate algorithmic detail but preserves it in the same document. |

These two cuts alone do not reliably reduce 24 pages to 20. A firm 20-page ceiling would require a more compact main Methods presentation or an additional exhibit moved to the appendix; it should be an explicit editorial decision, not an assumed consequence of these estimates. For advisor review, the current 24-page scientific account is coherent; it should not be described as meeting the original length target.

## Checks not completed and practical limits

- No alternative length cuts were applied or rendered, so savings above are estimates.
- No new fitting, bootstrap distributions or numerical analysis was run. This review used code, saved outputs, the verified earlier calculations and the revision diff. It does not reconstruct the historical fitting environment or resolve donor overlap beyond the already documented subset.
- No phone-width HTML or physical paper print was tested; printed readability was assessed from page images at the delivered dimensions.
- The local-expression dependency is verified. The broader assertion that no saved projection file exists anywhere was not independently established by this review; absence of such a file is not needed to document the producer's actual dependency.
- The known activity-KL normalizer and cached-prior-KL behavior remain disclosed implementation limitations; this writing review neither fixes them nor assesses the effects of refitting after a method change.

This file is ready for Andrew's review and remains uncommitted pending permission.
