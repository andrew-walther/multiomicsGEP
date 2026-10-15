# Project Writeup Implementation Plan

> **Current state (2026-10-14):** the initial draft exists (`0a7cfc9`, styling `a9f4b47`). Codex performs Task 4 using `review_prompt.md`. The **Author feedback on the initial draft** section below takes precedence over earlier requirements in this plan and in `drafting_prompt.md` where they conflict.
>
> **For agentic workers:** Claude Code executes Tasks 2, 3, and 5 on `writeup/draft`. Codex independently performs Task 4 after the draft exists; Claude addresses confirmed findings afterward. The October 14 request authorizes finalizing and committing these planning materials, not drafting the report in this session. Read the Review disposition before the author starts the separate drafting prompt.

**Goal:** Produce an advisor-facing account of the current single-modality and multimodal methodology, statistical performance, and biological findings, suitable for development into a manuscript.

**Architecture:** A standalone Quarto planning page establishes section briefs, evidence boundaries, and exhibits. After its approval, one main Quarto report includes two derivation appendix sources and renders to HTML/PDF using existing aggregate outputs and figures. A CSV source inventory records claim selectors and analytical caveats.

**Tech stack:** Existing Quarto and TinyTeX; Markdown and mathematical notation; existing aggregate CSV outputs and PNG/PDF figures. Any later tabulation code uses the project's existing R dependencies or base R; no new dependency or model-fitting runner is introduced.

**Spec:** Andrew's October 14 instructions to check Claude's review claims, finalize planning, and prepare the handoff; `planning.qmd`; and the October 13 Sequencing and handoffs in `docs/plans/Meeting_Notes_Followup_10_2_26.md`. The latest draft target is October 16, superseding the earlier October 14 target. The source review is preserved in `plan_review_findings.md`.

## Review disposition (2026-10-14)

The review was checked as claims, not applied verbatim. Only source/code/metadata inspection was performed; no fits, tests of biological hypotheses, enrichment, or new correlations were run. Source IDs below refer to the expanded inventory. Detailed numerical extraction is still the drafter's obligation, not implied by inspecting an implementation.

| Finding | Decision | Evidence and amendment |
|---|---|---|
| 1. Single-modality characterization methods | **Adopt in modified form.** | S41/S42 confirm E[F] ranking, `pos`, seed 1, size 10–500, and top-50/100/150 ORA; S53 confirms E[L] versus PurIST. S49 separates frozen projections/KM from signatures; S41/S51/S52 show mean-z signature scoring and the post hoc C4 exclusion. Add the requested methods table and dependence/multiplicity/availability caveats. Do not call the checks independent or the corrected control prespecified. Additionally, the helper's signature HR is per unit of the mean score, not automatically per SD, and its `max(C,1-C)` is outcome-oriented, not frozen-risk validation. |
| 2. Multimodal characterization | **Adopt in modified form.** | S43/S44 confirm top 8 positive/6 negative gene lists, rank-one modality shares, standardized joint-projection HRs, and saved absolute correlations. No fgsea/ORA calls were found in the inspected multimodal R scripts. Add completed-versus-proposed work and label-basis columns. Reject recomputing signed r now: saved tables S31/S32 contain only absolute correlations; label that limitation. S47's signed EBMF enrichment is an example for later methodology, not a completed multimodal analysis or an automatically approved choice. |
| 3. Evidence-gap table/source traps | **Adopt in modified form.** | S48's D1–D5 configurations and saved table compare LB/YFB preprocessing, not a fitted DeSurv comparator. Add distinct S45/S46 for EBMF factor/pathway correspondence and S32's flash-reference selector. The six-claim table is required in Task 3. State that a matched DeSurv/NMF benchmark was not verified in the reviewed sources, rather than asserting universal absence. The dated rank explanation in DECISIONS remains a hypothesis. Correlation supports correspondence, not identical programs or an alpha_F=0 supervision effect. |
| 4. Revisit mechanism/Continuing work | **Adopt.** | The section was absent; the sequencing M0–M8 table explicitly requires it. Add R1–R9 source markers covering methods, algorithms, both appendices, evaluation, results, exhibits/captions, biology, and conclusions. Keep proposed supervision work in one Continuing work subsection; retain factual caveats at affected results and one current-algorithm statement in Methods. |
| 5. Executable handoff/scope | **Adopt in modified form.** | The original Task 4 conflicts with the sequencing handoff. Reassign independent review to Codex, separate Claude's later repairs, and create `drafting_prompt.md`. Use the newly requested October 16 target. Reuse figures and derivations; reduce supporting S15/S19/S20 to concise text when verified. Retain sufficient Appendix A derivations plus a mathematical-versus-implemented table, and fully reconcile Appendix B. Omit unverified quantities or put the missing claim in the gap table. Do not freeze the review's pending-source count: verification status is per claim and must be updated as extraction proceeds. |

| Minor suggestion | Decision and evidence |
|---|---|
| Program 6 naming | **Adopt in modified form.** Use the August 21 report's immune infiltrate with a vascular component as the primary label (S21); discuss mixed stromal genes where supported rather than treating two labels as interchangeable. |
| Historical D4 sign | **Adopt.** S49 uses risk = −η and verifies it against the cached scores. S01 documents the historical sign issue. State the convention for this saved fit, not as a universal rule for current fitters. |
| Pruning after iteration cap | **Adopt.** S08's pruning branch checks `settings$prune && K > 1L`, without a successful-convergence guard. State this beside the unconverged real-data removal comparisons. |
| Multimodal flash reference | **Adopt.** S44 builds references `d4` and `flash`; S32 stores `reference=flash` rows. Extend S32 selectors and describe this as loading correspondence, not predictive superiority. |
| DeSurv citation | **Adopt in modified form.** Keep Young et al. attribution, using existing `Young2026` in S50 if a formal citation is rendered. The entry is `@unpublished`, submitted manuscript; reject inferring publication/PNAS details from old report or manifest wording. This needs no bibliography expansion or literature search. |

## Author feedback on the initial draft (2026-10-14)

Andrew reviewed the rendered draft. These requirements apply to the review (Task 4) and the subsequent revision, and supersede conflicting earlier requirements.

**Principle.** The writeup is a standalone narrative of the current method, results, and biological findings, written so its text and exhibits can move into a journal manuscript. It is not meeting documentation and not an audit of earlier claims.

| # | Requirement | Supersedes or amends |
|---|---|---|
| 1 | Title: *Survival-Supervised Bayesian Matrix Factorization for Prognostic Molecular Programs in Pancreatic Cancer* (prelim title with "Molecular" for "Gene Expression"; "Multi-Omic" is a later option if the multimodal model proves effective). No subtitle. | Draft title and subtitle |
| 2 | Finalized document shows no drafting traces: no "initial draft awaiting independent review" header or status box, no mention of reviewers, agents, the source inventory, or untraced claims. It reads as written by Andrew. | Draft status markings (removed at finalization, after review) |
| 3 | No repository references in the document: no scripts, file paths, functions, arguments, configuration labels, internal labels (D4, P3, flash k5), or source IDs. Source IDs are dropped entirely, not hidden in comments. Software packages may be named. `source_inventory.csv` remains the repository's verification record, kept current but not referenced by the document. | "Trace each claim to an exact path/selector" in the document; tracing now lives only in the inventory |
| 4 | No dated or process references (meeting notes, plans, progress chapters, dated decisions, "September model", "October chapter", "August 21 report"). Describe items by content. | Narrative references to dated sources |
| 5 | Exhibits are self-contained: every acronym and abbreviation defined in captions or table notes (e.g., Table 4); captions give population, metric, orientation and uncertainty without selectors. Noise levels are labelled by value (e.g., low/high noise), not "September". | Caption conventions |
| 6 | Algorithms use numbered, indented formatting without a box that breaks across pages. | Algorithm boxes |
| 7 | Figures are sharp and well sized (vector PDF for PDF, high-resolution PNG for HTML; fonts proportional to body text). Figure 1 (held-out risk-tertile curves) is to be redrawn. Redrawing from saved scores/aggregates is allowed; no refitting. Text stays monochrome; figures keep colour. | "Reuse existing figures" scope cut, for figures of poor quality |
| 8 | Drop the claim-centered evidence-gap table and the claim-by-claim framing. Supported findings go into Results and Discussion; unmet goals go into Continuing work as plain aims. | Task 3 "Claim-centered evidence-gap table" requirement; R9 mapping of that table |
| 9 | Section 5 becomes a manuscript-style Discussion: meaning and implications, limitations, continuing work (α_F = 0 statement kept there). | Interpretation section structure |
| 10 | Remove Appendix C: the characterization methods become a concise methods subsection in Evaluation (scientific caveats retained); the single-modality simulation table and figure move into Results 3.1 with context, or one is dropped. | Characterization table placement; Appendix C |
| 11 | No References section for now; DeSurv attributed in text (Young et al., submitted). | Formal citation via `Young2026` |
| 12 | Reduce length and redundancy; organize as Methods, Evaluation, Results, Biology, Discussion. | Page target unchanged (15–20 main-text pages) |
| 13 | Appendices A and B show complete derivations from the model and objective to the EBNM pseudo-observation $x = B/A$ and standard error $s = A^{-1/2}$, with every intermediate step and no unexplained jumps, and state where the evaluated algorithm departs from the joint model. | Appendix A "reuse and reconcile existing derivations" scope cut |

Revisit markers remain required and hidden; they must be carried through the restructuring. R9 now covers the findings summary, the Discussion (implications, limitations) and Continuing work.

## Global constraints

- Work on `writeup/draft` in the existing checkout. No worktree or second clone; Claude and Codex take turns. Check the branch before each commit and stop if it differs.
- Target a 15–20-page main text, including exhibits; derivation appendices are additional.
- Produce an article-style technical report, with numbered scientific sections, a findings summary, notation, and two algorithm boxes.
- Include no literature review, bibliography development, journal shortlist, personal updates, or fixed-date research timeline.
- Use evidence available through October 12, 2026. Proposed method/sensitivity work remains explicitly prospective.
- Change no fitting algorithms, APIs, hyperparameter defaults, raw data, existing prelim source, manuscript stub, or historical derivation.
- Rendering may read existing aggregate results and figures, but must not refit models, launch simulations, or require patient-level files.
- State actual evaluated settings, not assumed package defaults; preserve training-fixed risk orientation.
- Distinguish the principal regularized comparator from older and rank-matched comparisons, and molecular retention from survival-active classification.
- Explain CPTAC's proteomics platform and distinguish the single-block analysis from matched expression–methylation integration.
- Do not claim that survival shapes single-modality F under alpha_F=0 or that approximated objectives guarantee exact monotone CAVI.
- Preserve the `AGENTS.md` symlink; edit its target `CLAUDE.md` for the concise writeup pointer.
- Ask before commits, merges, or pushes unless Andrew explicitly authorizes the specific action. This October 14 planning commit is explicitly authorized; a draft commit/merge/push is not. Stage only task-owned files. Do not send/publish the report.
- Biological characterization must explain input, population, statistic/settings, multiplicity, and interpretation, with completed multimodal analyses separated from proposed M8 work.
- Report evaluated alpha_F=0 accurately; keep prospective supervision methodology in one Continuing work subsection. Apply Revisit markers at all method-dependent source passages and exhibits.
- Retain Young et al. attribution through existing `paper/ssbmf-refs.bib` entry `Young2026` if formal citations are used. No bibliography expansion; do not describe the submitted manuscript as a verified PNAS publication. Current planned manuscript authors are Walther, Liu, Rashid (DECISIONS.md, 2026-10-14).

## Review focus

1. **Mixed comparator outputs:** selecting a plausible but wrong CSV row must not change the identity of the main baseline. Check configuration, rank, penalty, cohort, and patient alignment.
2. **Theory versus evaluated algorithm:** equations must distinguish normalized plug-in single-modality prediction from the multimodal moment-based coefficient update, and state where survival actually enters.
3. **Overlapping cohorts and selection evidence:** the matched 50-donor set is not an additional independent cohort; external-performance evidence used for rank recommendation must be disclosed.
4. **Incomplete convergence and objective approximations:** describe iteration-cap fits as such, and the pruning objective as the implemented approximation rather than an exact survival expectation.
5. **Biological overinterpretation:** enrichment and unadjusted associations support program characterization, not causal mechanism, novel survival-guided discovery, or demonstrated methylation benefit.

## Revisit after feature/supervised-F

M0–M8 are defined in the sequencing plan. These are dependencies of the current evaluated account, not results expected from future work. Mark **every affected source passage**, equation/algorithm block, table/figure definition, and caption. The same ID can occur at several locations; mark individual blocks rather than enclosing the whole document once. When a new passage depends on the same method choices, apply the relevant existing ID.

```markdown
<!-- REVISIT R1 START -->
Current evaluated F-update passage, equation, or algorithm step.
<!-- REVISIT R1 END -->
```

Use balanced START/END pairs. A block depending on several rows may use nested balanced IDs. Markers remain in Quarto source and are hidden in HTML/PDF; rendered scientific prose remains a coherent description of the evaluated configuration.

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
| R9 | Findings summary, cross-model conclusions, Discussion (implications, limitations) and Continuing work. (The claim-centered evidence-gap table is dropped per Author feedback item 8.) | M0–M8 |

Drafter verification: list marker occurrences by ID and check that each named passage/exhibit above is covered, including captions and mathematical blocks. Codex checks this independently in Task 4. Following approved method results, update each ID from the new source selections; preserve or identify historical results rather than silently replacing their provenance.

## Draft schedule and scope cuts

Claude aims for a complete initial draft by **Friday, October 16, 2026**. This replaces the earlier October 14 draft target. Independent review and author approval still follow; a deadline does not make an unchecked report ready to share.

Reuse existing figures with necessary caption/label corrections instead of redrawing them for polish. Reduce threshold sensitivity (S15), merge diagnostics (S19), and cohort attribution (S20) to one well-supported sentence each where possible; retain a necessary caveat if more space is required. Appendix A reuses and reconciles existing derivations plus a mathematical-versus-evaluated table; Appendix B receives full equation-level reconciliation because its evaluated additions are newer. Do not omit derivations or characterization methods to save time.

Cut new analyses, fitting, simulations, enrichment, new signed correlations, CpG annotation, literature review, bibliography expansion, journal-template work, and exhaustive histories of discarded variants. If a quantity remains unverified at handoff, omit it or state the missing evidence in the claim-centered table; no placeholder numbers. Verifying existing files is required and is not permission to regenerate analytical results.

## File responsibilities

| File | Responsibility |
|---|---|
| `planning.qmd` / `planning.html` | Self-contained review page with outline, source map, evidence assessment, and exhibit plan. |
| `source_inventory.csv` | One source entry per ID; section, use, repository path, selector, configuration, verification status, caveat. Extend it to trace every quantitative claim in the final report. |
| `implementation_plan.md` | The execution checklist, checkpoints, verification contract, and progress record. |
| `writeup.qmd` / `writeup.html` / `writeup.pdf` | Main report; both rendered outputs include the derivation appendices. |
| `appendices/single_modality.qmd` | Single-modality derivations reconciled with the evaluated algorithm. |
| `appendices/multimodal.qmd` | Multimodal derivations with intercepts, signed priors, uncertainty, scaling, and pruning. |
| `figures/` | Only adapted or newly assembled presentation figures derived from saved results. Reuse existing figures by reference where sufficient. |
| `drafting_prompt.md` | Claude handoff for Tasks 2, 3, and 5; target October 16. |
| `review_prompt.md` | Codex handoff for Task 4, including the author feedback to review against. |
| `review.md` | Codex independent draft findings, followed by Claude dispositions; not a drafter self-review. |

## Task 1: Planning page and inventory

**Files:** Create `planning.qmd`, `source_inventory.csv`, and this plan. Render `planning.html`. Add concise pointers/status entries to README, CLAUDE, PROJECT_STATUS, ROADMAP, and DECISIONS.

**Consumes:** Existing code, decisions, derivations, progress chapters, saved result CSVs, and figure assets identified in the inventory.

**Produces:** A page whose section briefs are the authority for drafting after user approval; an inventory whose `source_id` keys identify exact source selections.

- [x] Inspect the current checkout and the current project instructions; confirm no pre-existing local changes.
- [x] Write section briefs for Summary, Methods, Evaluation, Results, Biology, Interpretation, and both appendices, with approximately 18 pages allocated to the main text.
- [x] Record source paths, selectors, verification status, and caveats. Include the principal comparison, full multimodal sweep, and matched validation rows.
- [x] Select the main tables/figure families; add two collapsible previews and explicitly distinguish measured findings from outstanding evidence.
- [x] Render from the repository root:

  ```bash
  quarto render paper/writeup/planning.qmd --to html
  ```

  Expected: exit 0 and `paper/writeup/planning.html`; image/style resources are embedded, and no fitting code executes.

  Local rendering note: the standalone Quarto 1.3.450 theme compiler is x86_64 and fails with `Bad CPU type in executable`. RStudio's already-installed ARM-native Quarto 1.10.18 successfully renders the page using:

  ```bash
  /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/planning.qmd --to html
  ```

- [x] Check every repository-relative link/asset and every inventory source path. Recompute the planning snapshot's means, row counts, and selected values from the named CSVs.
  Expected: all files exist; 54 unique inventory IDs; mean single-modality C rounds to 0.627; principal pooled difference is 0.0259; main single-modality simulation has 540 rows and 15 seeds; full multimodal sweep has 486 rows with 162 per variant; selected matched C values round to 0.6852 and 0.6757 and both selected fits are unconverged.
- [x] October 12 inspection (before the current revision): inspect the rendered page in a browser: title/navigation, tables, mathematical notation, collapsed/expanded previews, and narrow-screen readability. Verify local section anchors in the HTML.
- [x] Update project pointers: finalized planning review and handoff available; report not drafted.
- [x] Claude reviewed the plan; check all findings/minors and finalize revisions as Andrew requested on October 14. See Review disposition.
- [x] Prepare Revisit mapping, claim-centered evidence-gap task, and Claude drafting prompt.
- [x] Andrew inspects dispositions and starts the separate Claude drafting handoff. This session does not write the report.

## Task 2: Methods and derivations (Claude Code)

**Files:** Create `writeup.qmd`, `appendices/single_modality.qmd`, and `appendices/multimodal.qmd`; extend the source inventory.

**Consumes:** Approved section briefs; inventory entries S01–S08, S29–S30, S33–S34, S37–S38; realized settings from the evaluated runners/fits.

**Produces:** Common notation, two algorithm boxes, and derivations that consistently describe the configurations used in Results.

- [x] Establish cross-format HTML/PDF YAML: article layout, 11-point PDF text, one-inch margins, modest line spacing, section numbering, contents, embedded HTML resources. Include the two appendix sources from the main report.
- [x] Define notation/dimensions and write the evaluated single-modality likelihood, priors, normalized projection, working approximation, update order, convergence, prediction, rank recommendation, and factor classification.
- [x] Explain the single-block multi-cohort setting and separate genomic offsets, baseline-hazard strata, and cohort-specific survival coefficients; identify the proteomics training input.
- [x] Write the multimodal model and algorithm, including shared L, signed modality-specific F, feature intercepts/precisions, projection moments, scale canonicalization, and pruning/refitting.
- [x] Reconcile the existing single-modality beta/L/F/tau derivations, EBNM handoff, normalized projection, and objective. Add a compact table separating the mathematical joint model from the evaluated alpha_F=0 algorithm and plug-in score approximation; retain the derivation detail required to understand each evaluated update.
- [x] Derive multimodal intercept/L/F/beta/tau updates, second moments, feature coupling, canonicalization, and removal criterion; distinguish the plugin pruning objective from exact expected Cox likelihood.
- [x] Verify every algorithm step/equation against its implementation and evaluated configuration. Reuse independent algebra checks where they cover the relevant equation; do not change model code to make prose and code agree.
- [x] Render both formats:

  ```bash
  quarto render paper/writeup/writeup.qmd --to html
  quarto render paper/writeup/writeup.qmd --to pdf
  ```

  Expected: both exit 0; both appendices, equations, and algorithm boxes are readable and correctly referenced.
- [x] Check the affected methods, algorithm steps and appendix equations carry their Revisit markers; propose a checkpoint and ask before committing. Continue the author-started drafting assignment without section-by-section approval requests.

## Task 3: Statistical results and biological findings (Claude Code)

**Files:** Complete `writeup.qmd`; update `source_inventory.csv`; add adapted assets only when necessary in `figures/`.

**Consumes:** Approved exhibits; inventory entries S09–S28, S31–S32, S35–S36, S39–S54; the methods from Task 2.

**Produces:** Every numerical claim traces to a source selection, and every exhibit has a design, population, metric, and appropriate qualification.

- [x] Write the cohort/design table with exact sample/event counts, platform identities, transformations, screening, replication, noise settings, comparator definitions, and inference procedures.
- [x] Assemble single-modality simulation results by scenario/rank, distinguishing recovery, merging, prognostic classification, and held-out C; label the simulation comparator accurately.
- [x] Assemble five-cohort external performance with regularized baseline rows from S10 and principal pooled difference from S11; keep S12/S13 as supporting comparisons. Explain the role of external-performance evidence in the rank recommendation.
- [x] Include the biological-characterization methods table from planning.qmd, with columns analysis, input object, population, statistic/settings, multiplicity, and supported interpretation. Explain E[F] ranking, E[L]/PurIST tests, mean-z signatures, frozen projections, gene universes, BH scope, post hoc C4 controls, and overlapping evidence. Include the availability-based method-choice caveat.
- [x] Characterize all four retained programs by connecting their labels to those methods/results. Program 6 uses the updated immune/vascular label. State historical D4 risk = -eta; keep the signature helper's outcome-oriented C-index separate from predictive validation.
- [x] Assemble the full multimodal simulation summaries by variant/scenario/noise/rank; specify that 486 includes all three variants. Do not replace the selected variant's metrics with an across-variant average.
- [x] Assemble matched-donor comparisons and all starting ranks using explicit S28 selectors; include convergence and primary versus broader-histology sensitivity. Trace the documented established-model matched-donor comparison to its saved-score derivation before using it as a verified quantitative result.
- [x] Describe all seven multimodal programs with a label-basis column (top signed genes and named correspondence reference), explicit absolute r where that is all that was saved, within-modality rank-one shares, and 14 unadjusted joint-projection HR tests. Include the completed-versus-proposed list from planning.qmd in the report. No new signed correlations or enrichment.
- [x] Include S32's unsupervised flash-reference rows as correspondence evidence distinct from the expression-only/D4 references; do not treat them as a prediction benchmark. Note that removal checks can run after the iteration cap.
- [x] Write Interpretation with one Continuing work subsection for prospective supervision/rank/feature/characterization work; Methods describes the current evaluated configuration with a brief cross-reference.
- [x] Produce the six-claim evidence-gap table specified below. Omit unverifiable numbers from Results and record missing evidence instead.
- [x] Verify Revisit markers cover all affected passages, algorithms, equations, tables, figures, and captions.
- [x] Verify selected results against CSVs, denominators, comparator identity, orientation conventions, and inventory; render both formats successfully and check no private patient file is a render dependency.
- [x] Propose a full-draft checkpoint commit after verification.

### Claim-centered evidence-gap table (Task 3 deliverable; superseded)

**Superseded 2026-10-14 by Author feedback item 8.** The table below is kept as the record of what the initial draft contained. Its substance moves into the Results/Discussion narrative where supported and into Continuing work where not.

Place one table in Interpretation with columns **claim**, **existing evidence**, **source ID/path and selector**, **currently justified conclusion**, and **missing evidence**. Supply all six rows below. Mark the table and caption R9, plus any relevant row-level dependencies. It reports the present evidence, not an intention to run missing analyses during drafting.

| Claim | Existing source starting points | Boundary/missing evidence |
|---|---|---|
| Value over EBMF/NMF alone | S18/S19 recovery; S45 factor correspondence and S46 pathway concordance; S32 flash-reference rows. | Correspondence/recovery is not survival superiority. An unsupervised decomposition has no C-index without a specified scoring rule. No NMF comparison was verified in the reviewed comparator sources. |
| Value over two-step EBMF → Cox | S09–S13, with K=40 LASSO as the principal real-data comparison; S18 uses a differently configured simulation comparator. | Restrict conclusions to the exact benchmark; finalized-model/robustness comparisons await M6–M8. |
| Value over DeSurv | S35 gene overlap; S48 describes LB/YFB preprocessing comparisons. | No matched predictive head-to-head was verified. The directory name is not evidence of a DeSurv fit. Keep attribution to Young et al.; do not imply superiority. |
| K selection accounting for survival and reconstruction | S14/S15 rank/threshold outputs; S26–S28 multimodal rank/recovery; S01 dated decisions and S08 removal criterion. | Criteria disagree, external performance influenced K=7, real multimodal factors are not pruned, and the dated reconstruction-versus-K explanation remains unverified. Combined selection/replicated validation belongs to M5. |
| Multimodal added value | S28 same-donor comparison; S31/S32 program and modality correspondence. | Similar C is not equivalence; both selected fits are unconverged. Prediction improvement, robust selection, and methylation-specific biological explanation are not demonstrated. |
| Biological understanding | S21–S24/S35/S36/S40–S44/S49/S51–S54. | Explain methods and dependent evidence, post hoc controls, training concordance, unadjusted tests and unsigned matching. Mechanism, final-program validation and CpG annotation/enrichment remain gaps. |

## Task 4: Independent complete-draft review (Codex)

**Files:** Read the draft, appendices, outputs and sources; write `review.md`. Do not edit the draft or analytical results.

**Prompt:** `review_prompt.md`.

**Consumes:** Claude's completed draft; source inventory; Review focus; the Review disposition; Author feedback on the initial draft; methods table, completed/proposed list, and Revisit mapping.

**Produces:** Independent findings ranked by severity, with location, evidence and suggested fix. State limits of checks that cannot be completed.

- [ ] Check every quantitative claim against its path/selector, comparator, cohort size/overlap, risk orientation, denominator and uncertainty.
- [ ] Check equations/algorithms against the evaluated implementation; inspect HTML/PDF and record missing references/assets or layout defects.
- [ ] Check that biological methods are explained, claims follow the actual evidence, methods availability/dependence/post hoc caveats are included, and no new multimodal analysis is presented as completed.
- [ ] Check every Revisit dependency is marked in source, including captions and appendices.
- [ ] Review the draft against every Author feedback item, proposing concrete text or structure for each affected location (Part B of `review_prompt.md`).
- [ ] Write `review.md`, then ask before committing it. Do not mark a drafter's own checks as this independent review.
- [ ] Later handoff: Claude verifies each finding, applies all Author feedback items and the confirmed findings, records fixed/not-confirmed/deferred-with-reason dispositions, and rerenders. Andrew approves before merge. This later repair stage is not part of the initial drafting assignment.

## Task 5: Initial-draft documentation and delivery (Claude Code)

**Files:** Update README, CLAUDE, PROJECT_STATUS, ROADMAP, and DECISIONS; finalize the inventory and review record.

**Consumes:** Initial-draft render/source checks. Independent review dispositions are added only after they exist.

**Produces:** Accessible draft HTML/PDF and an accurate handoff for independent review.

- [x] Update README/CLAUDE/PROJECT_STATUS/ROADMAP/DECISIONS as relevant to **draft awaiting independent review**, with artifact/render instructions. Leave Task 4 unchecked; do not claim readiness to share or invent a review disposition.
- [x] Verify documentation links and the `AGENTS.md` symlink; inspect the final diff for unrelated changes and accidental private-data inclusion.
- [x] Deliver the report and explain the main organizational/mathematical choices, inputs/outputs, and remaining scientific limitations. Ask which parts warrant deeper review, in keeping with the project's post-implementation review practice.
- [ ] Check the branch is writeup/draft and ask before a draft commit; stage only task-owned files. Ask separately before any push/merge. End by handing the draft to Codex for Task 4, not by marking the whole workflow complete.

## Progress record

- 2026-10-12: user approved native execution with an HTML planning-page review before drafting. Task 1 started in the existing checkout; the working tree was clean.
- Task 1: source inventory contains 40 entries; all source paths and Markdown links exist. Selected numerical checks passed: mean single-modality C=0.62668; principal difference=0.0259 (0.0002, 0.0498); 540 single-modality simulation rows, 15 seeds; 486 multimodal fits, 162 per variant; selected K=7 matched values and unconverged flags reproduced.
- Task 1: Quarto HTML rendering passed using RStudio's existing ARM-native installation. HTML inspection found valid anchors/links and embedded image, script, and stylesheet assets; AGENTS.md still resolves to CLAUDE.md. Browser checks confirmed loaded embedded figures, rendered math, functioning previews, white title text on the dark banner, and no document overflow at 390-pixel and 1280-pixel viewports after the narrow-screen table repair.
- Task 1: rendering adjustment is confined to the command used locally; no dependency installation, global setting change, or report-format change was needed. The working command is documented in the page and README.
- Execution choice: this deliverable is documentation, so verification consists of source/claim checks, rendering, and visual inspection rather than new model tests. No model code is changed. The explicit user review checkpoint overrides the skill's general continuous-execution instruction.
- Tasks 2–5 remain pending; no report or derivation appendix has been drafted.
- Task 1 deliverables verified and presented for the required planning-page review. No commit or push has been made. Suggested checkpoint: `Add advisor-facing writeup planning page and source inventory`.

- 2026-10-14: checked Claude's committed plan review against source files. Finalized the plan and handoff on writeup/draft under explicit author authorization; no report drafting or new analyses. Independent review is reassigned to Codex. The native-drafting and initial no-commit statements above describe October 12 history and are superseded by this entry and the current workflow header.

- 2026-10-14 verification: RStudio Quarto rendered the revised planning HTML successfully. Static checks passed for 54 unique inventory IDs and existing source paths, authored Markdown links, HTML anchors/local links, embedded image/script/stylesheet assets, the 10-row six-column characterization table, the four-row completed/proposed table, identical handoff tables, R1–R9 mappings, and the AGENTS symlink. Whitespace checks preserve the CSV's existing CRLF style. No models, statistical analyses or model tests were run.
- Preview recovery: the old localhost server was unavailable. Automatic approval review rejected a repository-wide restart because it could expose unrelated private files. A narrower, approved temporary loopback server now allows only planning.html and its 34 explicit local links, with no directory listing. A fresh browser tab verified the October 14 page, updated methods table (11 rows including header), navigation, and no document overflow at the current 1265-pixel width. The earlier narrow-screen checks apply to October 12; no new narrow-screen claim is made. The preview script is temporary and is not part of the repository deliverable.

- 2026-10-14 (Claude Code, drafting): Tasks 2, 3 and 5 executed on `writeup/draft`. Status: **initial draft awaiting independent review**; Task 4 not started. Deliverables: `writeup.qmd` → `writeup.html`, `writeup.pdf` (38 pages: 21-page main text, Appendices A–C); `appendices/single_modality.qmd`, `appendices/multimodal.qmd`; `figures/make_pruning_rank_figure.R` → `figures/multimodal_pruning_final_rank.png`; source inventory extended from 54 to 62 IDs (33 existing entries' verification status/caveat updated, IDs unchanged).
- Render commands (RStudio Quarto; `latex-auto-install: false`, so rendering installs nothing): `quarto render paper/writeup/writeup.qmd --to html` and `--to pdf`. Both exit 0 with no warnings. The first PDF attempt failed on math in appendix headings and Quarto attempted a TinyTeX package update before failing; auto-install was then disabled. No other dependency changes.
- Source checks: every number in the report was read from the CSV/RDS named in the inventory. Saved per-cohort C-indices for D4 and the K40, K7 and K20 comparators were recomputed from saved risk scores with frozen orientation and matched (S56). Multimodal simulation summaries were recomputed from raw rows (S26), not the summary file. No model was fitted and no analysis script was rerun.
- Output checks: no unresolved cross-references or broken HTML anchors; five figures embedded in HTML; Revisit markers hidden in both outputs and balanced/nested in source (R1 13 blocks, R2 14, R3 5, R4 10, R5 5, R6 8, R7 8, R8 9, R9 6). PDF pages inspected as images; an overflow of code-formatted paths in the evidence-gap table was fixed with `hyphenat`. HTML was checked structurally; a browser view of the local file was not available.
- Remaining gaps and limitations: main text is 21 pages (target 15–20); the established model's matched-donor C (0.685) and the October per-program evidence decomposition have no saved source; the K40 per-cohort paired intervals have no CSV; the pathway "suppression" reading conflicts with the verified orientation (resolved in text; flagged for review); the multimodal `fixed_working_objective` diagnostic omits cross-factor terms (code not changed). No commit has been made.
- 2026-10-14 follow-ups: (1) Andrew requested monochrome text; HTML CSS and PDF link/callout colours set to black, figures keep colour. A pixel check of the PDF finds colour only on the five figure pages; the HTML was checked for the CSS but not viewed in a browser. (2) The monitoring-objective defect noted above was fixed separately on `fix/multimodal-monitor-objective` (`d1f8b77`, 532 tests passing); it does not affect any reported result.
- 2026-10-14 (Claude Code, revision): applied `review.md` (Part A verified against code and saved outputs; all 14 Part B items). Dispositions for every finding are in `review_dispositions.md`. Restructured report (Summary paragraph, Methods, Evaluation, Results, Biology, Discussion; Appendices A–B only), rewrote both appendices as complete derivations, removed source IDs/repository/process references, the evidence-gap table, Appendix C and the bibliography. Figures redrawn from saved outputs by `figures/make_writeup_figures.R` and `figures/make_pruning_rank_figure.R` (vector PDF / 300-dpi PNG); inventory S39/S62 updated and S63 added. Both formats render without warnings; Revisit markers balanced; no table splits across pages. PDF 42 pages, main text 27 (target 15–20; cuts need author judgement). Implementation defects A3/A5 documented, not fixed. No model code, defaults, data or saved results changed.
- 2026-10-14 (Claude Code, verification round): checked and applied `verification.md` V1–V13 and its incomplete A/B items (dispositions appended to `review_dispositions.md`). Required corrections: protected multimodal precision floor, zero-column normalized target, orientation wording, working-quadratic attribution, remaining derivation steps, local caption definitions, Revisit gaps, unbroken algorithms (PDF minipage) and tables (measured needspace, three pinned figures). Low-loss length changes (Codex ranks 1–5 and 7) applied; ranks 6 and 8 left for the author. PDF 39 pages, main text 24 (was 27). No model code, defaults, data or saved results changed.
- 2026-10-14 (Codex, second verification corrections): Andrew authorized N1/N2 from `verification_round2.md`. The first evaluated sweep-normalizer formulas in Methods §1.2 and Appendix A.1 now include the numerical ridge already disclosed in A.10; the project-status disposition pointer is corrected. Review/status pointers are updated. Further length cuts and final approval remain pending; the local-expression figure dependency remains disclosed. No model code, defaults, figures, data or saved results changed.
- Verification: HTML and PDF rendered successfully; both corrected formulas are present, PDF pages 3/25 visually checked, and only those pages' extracted text changed. Pagination remains 39 total / 24 main-text pages, with algorithms and tables on the same pages. Revisit markers remain balanced/nested; `git diff --check` passes.
- 2026-10-14 (Claude Code, final editorial check): wrote `final_editorial_check.md`; applied required corrections E1–E5 and Methods cuts M1–M4 at Andrew's request (algorithms and removal-score display moved to the appendices, Approximations subsection folded into Limitations). Main text 22 pages (PDF 39). Dispositions appended to `review_dispositions.md`. No model code, defaults, data or saved results changed.
