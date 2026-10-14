# Claude Code drafting handoff: current project writeup

## Assignment and authority

On **writeup/draft**, execute **Tasks 2, 3 and 5** of `paper/writeup/implementation_plan.md`. Produce the initial technical report and both included derivation appendices. Aim for a complete initial draft by **Friday, October 16, 2026**. This supersedes the earlier October 14 target; it does not imply independent review or author approval will be complete by then.

Run `git status -sb`, `git branch --show-current`, and `git log --oneline -1` before writing. The branch must be `writeup/draft`. If it differs, stop and ask rather than switching automatically. Work in this checkout only, one session at a time; no worktree or second clone. Preserve unrelated files/changes. Ask before a draft commit, push, or merge; this prompt does not authorize those actions. Check the branch again immediately before any separately authorized commit and stage only your files.

Read project/global instructions, then:

1. `paper/writeup/planning.qmd` and its rendered HTML.
2. `paper/writeup/implementation_plan.md`, including **Review disposition**, **Global constraints**, **Review focus**, **Revisit after feature/supervised-F**, the scope cuts, and evidence-gap table task.
3. `paper/writeup/source_inventory.csv`; use exact paths/configurations/selectors and distinguish inspected code from verified numerical evidence.
4. The Sequencing and handoffs section of `docs/plans/Meeting_Notes_Followup_10_2_26.md` and the current dated decisions in `DECISIONS.md`.
5. `paper/writeup/plan_review_findings.md` as background, with the checked dispositions in the implementation plan taking precedence over the review's unchecked claims.

Andrew requested this finalized handoff after a review of the plan. Do not reopen the planning cycle or start method work. Surface a concrete contradiction or missing source when it prevents a defensible statement; otherwise apply the specified evidence boundaries and continue drafting.

## Deliverables

- `paper/writeup/writeup.qmd`, rendered to `writeup.html` and `writeup.pdf`, including `appendices/single_modality.qmd` and `appendices/multimodal.qmd` in both outputs.
- An article-style main text of approximately 15–20 pages including exhibits: findings summary; Methods with two algorithm boxes; Evaluation framework; single-modality then multimodal Results; Biology; Interpretation with one Continuing work subsection; two derivation appendices outside the main-text page target.
- Updated source inventory tracing each quantitative claim to an exact source path/configuration/selector and recording what was actually verified. Existing source IDs must remain stable.
- Necessary presentation assets only under `paper/writeup/figures/`; reuse existing figures with corrected captions/labels when sufficient.
- Task 5 documentation updated to **draft awaiting independent review**, with render commands, artifact paths, source-check evidence, remaining gaps, and draft limitations. The report does not become ready to share merely because it renders.

Keep the initial report about evaluated methods and existing evidence through October 12; later planning/verification dates can be stated accurately. Do not anticipate M0–M8 results. Planned manuscript authors are Walther, Liu, Rashid per the October 14 decision. Retain attribution to Young et al./DeSurv using existing `Young2026` if a formal citation is rendered; the existing entry is a submitted, unpublished manuscript. Do not infer publication venue/DOI from old report or manifest wording. No literature-review or bibliography-expansion project is required.

## Biological characterization: explain methods before interpreting programs

Include the adopted table below in Evaluation, or a clearly referenced supporting table if layout requires it, with explanatory prose. Preserve all six column meanings. For each biological label in Biology, explain which analyses support it and how their outputs lead to that interpretation. Do not supply only pathway names or top genes.

| Analysis | Input object | Population | Statistic and settings | Multiplicity | What the output supports |
|---|---|---|---|---|---|
| Ranked enrichment (fgsea) | Nonnegative posterior means E[F], over the 2,064 selected genes | Established single-modality fit; one program at a time | Raw loading ranking; `scoreType="pos"`; set sizes 10–500; seed 1. Hallmark, Reactome, KEGG, GO-BP, custom PDAC (Moffitt/Bailey and DeSurv lists). Explain normalized enrichment and leading-edge genes. | BH within each program × collection; reporting padj < 0.10. No stated control across programs/collections. | Association of high-loading genes with predefined sets, not mechanism or survival-driven discovery. S21/S22/S41/S42. |
| Over-representation sensitivity | Top 50, 100, or 150 genes by E[F] | Same selected feature universe | Hypergeometric over-representation via `enricher`; supplied background is the 2,064 selected genes, not the genome. Explain the cutoff sensitivity. | Report the saved adjusted p-values and the within-call family; do not imply global error control. | A threshold-based cross-check using the same gene universe, not independent replication. S21/S41/S42. |
| Subtype concordance | Training latent activities E[L], not YF projections | TCGA, n = 144 | Spearman correlation with PurIST.prob; Kruskal–Wallis comparison across PurIST calls | Saved test p-values are unadjusted; state that scope. | Within-training correspondence to an established subtype axis. S23/S41/S53. |
| DeSurv gene-list overlap | Top 270 genes by E[F]; DeSurv lists | Both lists restricted to the selected shared-gene background | Jaccard similarity and hypergeometric overlap test | Reported overlap p-values are unadjusted | Gene-list correspondence, not a predictive head-to-head comparison. The custom fgsea collection reuses the same lists. S35/S41. |
| External signature association | Mean of available leading-edge genes after within-cohort gene z-scoring | Five held-out cohorts | Univariable Cox HR per unit of that mean signature score; state gene coverage. C2 chooses the lowest-padj signature; C4 subsequently excludes PDAC_custom for genomics-only controls. | Unadjusted association tests; C4 is a post hoc correction | External descriptive associations of the selected signature, not necessarily the full program; not preregistered confirmation. S24/S40/S41/S51/S52. |
| Program directions and survival curves | Marginal projection association for direction labels; frozen Program 3/7 projection scores and overall risk for curves | Direction labels documented in decisions; five external cohorts for KM | Marginal directions differ from conditional joint β signs. Within-cohort tertiles; pooled program curves with cohort-stratified log-rank. Historical D4 risk = −η. | Descriptive grouping/tests; report actual scope | Direction and held-out separation, without clinical cutoffs. Keep these scores distinct from signature scores and E[L]. S01/S36/S49. |
| Multimodal top-gene interpretation | Signed expression loadings: up to top 8 positive and 6 negative genes | Matched TCGA fit at starting rank 7 | Manual interpretation of the output gene lists, not fgsea/ORA | No enrichment test performed | Provisional expression-supported labels; weaker characterization than the single-modality enrichment analysis. S31/S43. |
| Multimodal modality shares | Posterior-mean L and modality-specific F | Same fitted programs | ‖L_k‖²‖F_mk‖², divided by the sum of these rank-one squared norms within modality m | Descriptive | Relative fitted rank-one signal, not additive PVE or survival contribution. S31/S43. |
| Multimodal correspondence | Expression loading columns on shared genes | Joint fit versus expression-only, established D4, or unsupervised flash reference | Saved absolute Pearson correlations and best matches; identify the reference and common-gene count. Saved tables contain only \|r\|. | Descriptive matching | Pattern similarity with direction unavailable; do not imply sign agreement or a one-to-one match from independent best matches. S32/S43/S44. |
| Multimodal marginal survival association | Z_k = Σ_m Y_m F_mk, standardized within cohort | TCGA and ICGC primary PDAC | Univariable HR per one SD of the joint projection, including methylation | 7 programs × 2 cohorts = 14 unadjusted tests | Descriptive marginal association, distinct from joint posterior survival-active classification. S31/S43. |

Explain that custom fgsea sets and gene-list overlap reuse DeSurv lists and the DeSurv-aligned selected universe, so they are not independent validations. C4 is a post hoc correction to the borrowed-signature problem in C2, excluding the entire PDAC_custom collection. Subtype tests are within training. BH correction is within program × collection, not across all programs/collections. The August 21 report says method choices were availability-driven and warrant assessment; carry this as a scientific qualification, not a request to run new analyses.

Distinguish signature HRs per unit mean-z signature score from multimodal HRs per within-cohort SD of the joint projection. The signature helper's `max(C,1-C)` is outcome-oriented and cannot substitute for frozen-risk predictive validation. For the historical D4 cached fit, `make_km_figures.R` uses **risk = -eta**, fixed across validation cohorts; do not negate every current model automatically.

Use Program 6's updated **immune infiltrate with a vascular component** label from the August 21 report, explaining mixed stromal genes where supported. The historical label alone is not an additional independent finding.

## Multimodal characterization: completed versus proposed

Include this distinction explicitly, and give the multimodal program table a **label basis** column (manual top signed genes and named correspondence reference):

| Completed multimodal characterization | Proposed characterization (M8; not completed evidence) |
|---|---|
| Signed top-expression-gene lists and manual provisional labels | Sign-aware/two-sided enrichment of expression loadings; the existing EBMF pathway runner illustrates an alignment convention, not a completed multimodal analysis |
| Within-modality rank-one fitted-signal shares | CpG annotation and genomic-context enrichment |
| Absolute loading correlations with expression-only, D4, and unsupervised flash references | Checks of cell-composition and technical explanations; direction-aware correspondence if later computed |
| Fourteen unadjusted, within-cohort standardized joint-projection Cox associations | Characterization and validation of the final programs after method work |

The saved multimodal correspondence sources contain only absolute r. Report that limitation; do not compute new signed correlations. Different best matches need not define a one-to-one alignment. Molecular shares normalize the sum of squared rank-one reconstruction norms within modality, not the total reconstruction norm including cross-terms and not survival attribution. The signed EBMF reference enrichment runner is an existing methodological example; it is not a completed multimodal enrichment analysis or an authorization to add one.

## Scientific and numerical boundaries

- Describe the **evaluated** single-modality alpha_F=0 algorithm accurately. Survival informs the coefficients, not discovery of F; the plug-in coefficient update does not propagate full loading uncertainty. Keep L, F, normalized ZF, signature scores, and risk scores distinct.
- Reconcile equations and algorithms with the actual evaluated settings and implementation. Appendix A reuses the existing derivations plus a mathematical-versus-implemented table; retain sufficient derivation of the evaluated beta/L/F/tau updates. Appendix B requires full reconciliation of intercept, signed priors, moment terms, scale transformations, objective and removal/refitting.
- Explain approximations: Cox working quadratic, plugin versus expected survival contribution, and the absence of a general exact-bound/monotone-ascent guarantee. Pruning can run after reaching the iteration cap; do not describe these removal comparisons as all based on converged fits.
- Distinguish starting rank, retained molecular programs, survival-active classification, and actual removal/refitting. Explain criterion disagreement and the external-performance evidence behind the current K=7 recommendation. Do not state the historical unverified reconstruction-versus-K explanation as a demonstrated result.
- Use K=40 independently tuned LASSO EBMF → Cox as the principal real-data comparator, with exact S10/S11 sources. The rank-matched K=7 and older K=20 unpenalized comparisons are separate supporting results. Simulation comparators have their own settings. `desurv_comparison` contains our LB/YFB preprocessing comparisons, not a fitted DeSurv benchmark.
- Identify CPTAC as proteomics; distinguish the single-block pooled analysis from matched expression/methylation. Keep the matched 50-donor ICGC set, overlapping PACA-AU RNA-seq cohort, and 67-donor sensitivity set distinct.
- Report full-sweep multimodal results by actual variant/scenario/noise, not an across-variant average. State replication and incomplete real-data convergence. Similar observed C-indices are not an equivalence test; added structure is not automatically added prediction or methylation mechanism.
- Keep prospective supervision, prior/weight, rank-selection, feature-count and biological-characterization plans in **one Continuing work subsection**. State there that alpha_F=0 is not the intended final method. Methods states the evaluated configuration with a short cross-reference; factual qualifications may remain next to affected results.

## Evidence-gap table and future revision markers

Produce the Task 3 evidence-gap table, with six claims:

1. Value over EBMF/NMF alone.
2. Value over two-step EBMF → Cox.
3. Value over DeSurv.
4. K selection accounting for survival and reconstruction.
5. Multimodal added value.
6. Biological understanding.

For each give existing evidence, exact source ID/path/selector, conclusion justified now, and missing evidence. S45/S46 provide EBMF factor/pathway correspondence; S32 also contains the unsupervised flash reference. These are not predictive head-to-heads. No matched DeSurv or NMF predictive benchmark was verified in the reviewed sources; do not imply universal absence beyond inspected evidence or manufacture a comparison. Unsupervised decomposition alone needs an explicit scoring rule for a survival C-index.

Mark **every method-dependent source passage**, equation/algorithm block, table/figure definition and caption using the R1–R9 mapping in the implementation plan:

```markdown
<!-- REVISIT R1 START -->
Current evaluated loading-update passage.
<!-- REVISIT R1 END -->
```

| ID | Passages/exhibits to revisit | Method dependency |
|---|---|---|
| R1 | Single-modality F update, alpha/alpha_F, beta priors and supervision weights, initialization and stability rules; single-modality algorithm; Appendix A's affected update/objective derivations and mathematical-versus-implemented table. | M1–M3 |
| R2 | Notation and definition/normalization of ZF, coefficient units/uncertainty, scale transformations, training/new-patient prediction, corresponding steps/equations in both algorithms and appendices. | M2; M4 |
| R3 | Multimodal loadings/coefficient updates, priors/weights, normalization, objective, factor-removal/refitting algorithm, and Appendix B. | M4; M5 |
| R4 | Starting/effective rank definitions, retained/survival-active classification, K recommendations/criteria, model/notation table, rank/recovery panels (Figure 4 and supporting K sweeps), and their captions. | M5; M6; M8 |
| R5 | Evaluation design and sample/feature tables; preprocessing/feature counts, prior/initialization/convergence settings, robustness/sensitivity descriptions. | M0; M6; M7 |
| R6 | Simulation recovery, merging, false prognostic programs, low-variance supervision evidence and convergence; Figure 1, simulation tables and supporting panels/captions. | M1–M7 |
| R7 | Comparator definitions and prediction rules, single-/multimodal validation text, paired differences, Table 3/Figure 2, Table 4/Figure 5, risk KM curves and supporting comparator exhibits/captions. | M7–M8; effects of M2–M6 |
| R8 | All Biology interpretations and final program identities; characterization methods whose inputs/sign/selection change; label-basis table; Figure 3, Table 5/Figure 6, program HR/KM, enrichment/subtype/overlap, cohort attribution, modality shares, CpG discussion and captions. | M8; effects of M2–M6 |
| R9 | Findings summary, cross-model conclusions, evidence cutoff/code-result provenance, claim-centered evidence-gap table/caption, limitations and Continuing work. | M0–M8 |

Use the appropriate ID(s), balanced markers, and multiple occurrences as needed. Do not mark the whole report once as a substitute. Include summary/conclusions, both appendices, algorithm steps, all affected exhibit captions, and biological methods when their inputs change. List marker coverage by ID in the handoff. Hidden source markers must not interrupt the rendered scientific narrative.

## Scope cuts for the October 16 target

Reuse existing figures with needed caption/label corrections rather than redrawing for polish. Reduce supporting merge diagnostics (S19), coefficient-attribution simulations (S20), and threshold sensitivity (S15) to one verified sentence each where feasible; retain necessary caveats. Reuse and reconcile Appendix A rather than re-deriving unrelated alternatives; prioritize Appendix B's newer equation-level additions.

Cut new model fits, simulations, statistical/enrichment analyses, new signed correlations, CpG annotation, exhaustive historical variants, literature review, bibliography expansion, journal selection/template work, and cosmetic figure redesign. Do not cut biological methods, the completed/proposed distinction, algorithms, needed derivations, source checks, Revisit markers, or essential scientific qualifications.

If a source/quantity is still unverified when drafting ends, **omit that number or state the missing evidence in the gap table**. Do not insert placeholder numbers, run a new analysis to fill the gap, or silently repeat a narrative claim as if checked. In particular, trace the documented established-model C on matched donors to an existing saved source before treating it as verified; if no corroborating source is found, retain an explicit evidence gap.

Change **no model code, defaults, raw/processed patient data, saved analytical results, historical derivations, prelim chapter or manuscript stub**. Do not run existing analysis scripts merely to regenerate a result. Reading saved aggregates and checking their selectors is permitted. Rendering must not require private patient-level inputs or fit a model. Do not install dependencies.

## Verification, status, and next handoff

Use the installed renderer (the standalone Quarto theme compiler has an architecture mismatch):

```bash
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/writeup.qmd --to html
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/writeup.qmd --to pdf
```

Check source claims/settings, output coverage, citations/cross-references, figure/table/math legibility, and the main-text/appendix split. Inspect both outputs and record actual verification, omissions, warnings and limitations. These drafter checks are necessary and are **not Task 4's independent review**.

Update Task 5 documentation to **draft awaiting independent review**, leave Task 4 unchecked, and do not create `review.md` claiming a favorable review. End with a concise explanation of deliverables, key design choices, inputs/outputs and gotchas so Andrew can understand the document. Ask before committing and never push or merge without separate authorization.

**Stop after delivering the initial draft.** Codex next performs the independent review under Prompt 3 of the sequencing section and writes findings to `review.md` without editing the draft. Claude later verifies/fixes those findings and records dispositions under Prompt 4; Andrew approves before merging. No additional reviewer agent is a substitute for that separate Codex handoff.
