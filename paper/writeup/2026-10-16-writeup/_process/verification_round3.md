# Verification of the final editorial round

Reviewed at `6bd53cd` on `writeup/draft`, comparing `d80b724..6bd53cd` and the uncommitted change to `docs/plans/Meeting_Notes_Followup_10_2_26.md`. Initial branch/status checks passed: that plan file was the only modified file, with no untracked files. Read the editorial check, final dispositions, full source diff, current sources and rendered outputs. Only this verification file was written. No report, plan, code, figures, defaults, data or saved analytical results were edited; no fitting, commit, push or merge was performed.

## Recommendation

**Conditional go for the interim write-up, after one sentence is corrected.** E3's replacement still overgeneralizes the low-variance simulation result. Qualify it by the lower-noise setting, as specified below, then refresh both renders. No further scientific or layout blocker was found in this editorial round.

Andrew must also accept the **22-page main text** as an exception to the 15–20-page target, or request further cuts. The report is suitable for interim advisor review; it does not resolve the disclosed methodological issues or establish survival-guided discovery in the evaluated single-modality configuration.

## Required corrections E1–E5

| Item | Verdict | Evidence in source and rendered formats |
|---|---|---|
| E1 — Matched-donor caption | **Verified** | `../writeup.qmd`, Table 7 caption, now calls the row the “single-block model” and identifies unstandardized log2 expression and the 100-sweep cap. The HTML caption and PDF p.16 agree with the row label. |
| E2 — Base comparator | **Verified** | Comparators §2.3 introduces four further arms, including the base multimodal model without intercepts/removal, with nonnegative loadings and moment-based prior updates. Present in HTML and PDF pp.7–8; consistent with Table 7. |
| E3 — Continuing-work evidence boundary | **Incomplete** | §5.3 item 1 removes the unsupported general claim about single-modality supervision, but replaces it with an unqualified statement that the multimodal planted low-variance program “was not recovered.” This appears in source, HTML and PDF p.22. Section 3.4 and Table 6 distinguish zero recovery at lower noise from 52% recovery at higher noise (48% with correct sign). See R3-1 below. |
| E4 — Exhibit citations | **Verified** | Notation cites Table 1; A.11 cites Table 10; B.13 cites Table 11. Located in HTML and PDF pp.2,31,39 respectively. Source label checks find no unresolved cross-reference; HTML fragment links all have targets. |
| E5 — First prose definitions | **Verified** | §2.1 introduces PACA-AU as the Australian pancreatic cancer cohort of the International Cancer Genome Consortium (ICGC) and expands CpG at its first prose use. Definitions are present in HTML and PDF p.6, before subsequent prose uses. |

## Methods cuts and preservation

| Item | Verdict | Evidence |
|---|---|---|
| M1 — Relocated algorithms | **Verified** | Algorithms 1 and 2 open Appendices A and B after their introductory paragraphs. Their complete text is identical to the previous main-text versions except Algorithm 1's necessary replacement of `@eq-cox-working` with `@eq-a-working`. Inputs through prediction remain whole on PDF pp.23 and32. HTML contains both complete lists. The main Methods Fitting summaries preserve initialization/update order, working-quantity timing, precision updates, stopping, orientation and the removal check, and link to the algorithms. |
| M2 — Predictor and working approximation | **Verified** | Methods §1.2 retains the posterior-mean fixed projection within a sweep, stored training norms for prediction, and the explicit Poisson-style/not-Taylor qualification. A.1 retains the random target's zero-column convention and the ridge in the sweep-fixed normalizer; A.10 retains the stored final-sweep normalizer and prediction formula. A.2's working equation is now Equation 9. Source, HTML and PDF pp.3,23–25,30 agree. |
| M3 — Removal score | **Verified** | Methods §1.5 retains plug-in survival versus posterior expectation, KL penalties cached against priors at the latest coordinate updates, candidate refitting of intercept/precisions, deletion only if the score does not decrease, checks even on unconverged fits, at least two factors before deletion, and no exact-evidence/current-ELBO guarantee. B.11 contains the relocated score (Equation 36), deletion difference (37), and the same restrictions. Algorithm 2 and §3.5 resolve their forward references to Equation 36. Present in HTML and PDF pp.5,16,32,38. |
| M4 — Removed Approximations subsection | **Verified** | Limitations §5.2 retains the unequal monitoring weights, approximate activity KL, plug-in survival and earlier-prior cached penalties, adding that neither algorithm is exact coordinate ascent on one ELBO or guaranteed to increase its monitor. B.11 retains Jensen's inequality, the full Hessian/covariance correction, the diagonal variance-penalized alternative and its lack of a lower-bound/ascent guarantee. The changed working surrogate and curvature are explained in Methods/A.2/B.2. No scientific qualification was lost by removing the separate subsection. |
| M5 | **Not applied, as disclosed** | Rank-selection cross-validation numbers remain in §1.3. This was not one of the authorized cuts in this round. |
| O1–O4 | **Not applied, as disclosed** | Biological-characterization methods remain in Evaluation §2.5; the peripheral simulation sentence, repeated standardization/caption wording and optional minor wording remain. These are editorial suggestions, not approval blockers. |

All requested source references resolve: `eq-a-working`, `eq-mm-elbo`, `sec-appendix-a-alg`, `sec-appendix-b-alg`, and `sec-appendix-b-removal`. No reference to the removed `eq-cox-working` or Approximations section remains. The delivered HTML has no duplicate IDs or missing fragment targets. PDF equation references identify the correct displayed equations.

## Continuing work

- **Item 1: aim preserved; evidence claim incomplete.** Recovery in most simulated datasets, stability and external concordance remain the stated aims. Candidate explanations are explicitly labelled plausible rather than demonstrated. The first result sentence needs R3-1's noise qualification.
- **Item 3: verified.** Both models, survival and reconstruction, selection without validation outcomes, and replicated simulation evaluation remain explicit.
- **Item 4: verified.** Initialization, loading-prior and feature-count sensitivity remain, as do comparisons against unsupervised factorization, two-step fitting and DeSurv on common data. “Feature counts” condenses the earlier gene/CpG wording without deleting that scope.

These statements agree between source, HTML and PDF p.22. All five aims remain.

## Problems ranked by severity and exact fixes

1. **R3-1 [P2, residual evidence overstatement]: qualify the E3 replacement.** At `../writeup.qmd:600`, Continuing work item 1, replace its first sentence with: **“In the lower-noise multimodal simulation, where strong reconstruction programs accompanied a planted low-variance prognostic program, the latter was not recovered despite survival entering the loading update.”** Keep the following explanations labelled plausible. This restricts the claim to the reported failure setting without contradicting the higher-noise result. Evidence: §3.4, Table 6, inventory S27 and the saved `pruning_comparison.csv` rows for `variant == pruned`, `scenario == low_variance_prognostic_factor`. The 27 lower-noise rows have mean prognostic recall 0; the 27 higher-noise rows have mean recall 0.5185185 and correct-sign recall 0.4814815. These are checks of existing saved outputs, not new fitted analyses. Re-render HTML/PDF and check this sentence and p.22's pagination. The original editorial recommendation itself omitted this qualification; applying that recommendation literally did not make the resulting generalization correct.

No new algebra error, lost scientific qualification, broken reference or layout defect was found in the material moved this round.

## PDF layout and length

Inspected all **39 page images** in overview, with enlarged images of the two relocated algorithms, A.11, B.11 and Continuing work. Main text ends on **p.22**; Appendix A starts p.23 and Appendix B p.32. Excluding the title/contents/Summary page leaves 21 main-body pages, so counting convention alone does not meet 20 pages.

| Table | Page | Table | Page |
|---|---:|---|---:|
| 1 | 2 | 7 | 16 |
| 2 | 6 | 8 | 17 |
| 3 | 10 | 9 | 20 |
| 4 | 11 | 10 | 31 |
| 5 | 13 | 11 | 39 |
| 6 | 14 | | |

Each table's caption, header and final row stay on one page; no stray header appears before its caption. A.11 begins on its own fresh page, with Table 10 beneath its heading and introductory citation. Algorithms 1/2 are complete and unboxed on pp.23/32. No orphaned heading, clipped equation or caption collision was found. Equations 9 and36 are legible, including the relocated removal-score display.

The former biology gap is resolved: Figures 4 and5 share p.19, Table 9 is on p.20, and its interpretation follows on the same page. The space beneath Table 6 on p.14 precedes the unbroken rank-sensitivity figure; p.22 is a shorter final main-text page. These are explainable layout costs rather than missing material or page-sized blank gaps. Remaining appendix space is acceptable for keeping tables intact. Further length cuts remain an author decision, not a scientific correction required by this review.

## Revisit markers

**Verified balanced and properly nested:** 54 blocks in `../writeup.qmd`, 18 in `single_modality.qmd`, and10 in `multimodal.qmd`. Algorithm 1 remains within R1/R2; Algorithm 2 within R3/R2 (and the appendix-wide R3 block); the removal-score equation remains within R3/R4. Their enclosing blocks include the complete algorithms/equation, not just their headings. The new Methods summaries and shortened Continuing work remain marked.

## Plan update and scope

- **Prompt 5: verified.** It no longer points to an evidence-gap table; it names the implementation-plan Revisit list, Discussion's Continuing work and hidden R1–R9 markers.
- **Prompt 7: verified for manuscript conventions.** It specifies source/path/date-free manuscript prose, self-contained captions, appendix algorithms, unbroken tables, both renders, dispositions and approval before commit. Its saved-output-only figure instruction is a future production requirement: the existing V11 caveat still applies to the current program-curve producer, which needs local expression for frozen projections. That pre-existing dependency has not been removed by this editorial round and should be accounted for when executing the future prompt.
- **Plan-only scope: verified.** The entire uncommitted plan diff consists of those two prompt edits; no other section changed. This review leaves that edit untouched and unstaged.
- **Committed scope: verified.** `d80b724..6bd53cd` changes only `paper/writeup/**` and `PROJECT_STATUS.qmd`. No model code, defaults, data, saved analytical results or figure producers/assets changed. The existing table values and figure assets were retained. `git diff --check` passes.

## Checks not completed

- No fresh browser runtime/MathJax or responsive-width test: opening a loopback preview server was denied by the sandbox (`PermissionError: Operation not permitted`). Saved HTML content, complete algorithm lists, IDs and links were checked directly; PDF page images independently establish the mathematical and pagination checks above. The previous round's browser runtime check is not represented as a new check here.
- No full re-derivation of unchanged appendices or repetition of all earlier numeric audits. The current source diff relocates the previously reviewed algorithms and score without changing their algebra; the changed evidence claim was checked against saved simulation rows.
- No re-render, new simulation, fitting, figure regeneration or physical print test was performed. The delivered HTML/PDF were the review targets.

## Resolution after Andrew's authorization

R3-1 is corrected in `../writeup.qmd` using the exact replacement above. HTML and PDF rendered successfully after allowing Quarto access to its rendering cache. The qualified sentence is present in both outputs; only PDF p.22's extracted text changed. Its page image was inspected: the five continuing-work items remain legible and unbroken, with no clipping. Pagination remains 39 total pages / 22 main-text pages; algorithms and tables retain their reviewed pages. HTML fragment targets resolve and Revisit markers remain balanced/nested. No model code, data, figures or saved results changed. The required correction from this round is resolved; accepting the 22-page interim length remains Andrew's decision. No commit, push or merge was made.
