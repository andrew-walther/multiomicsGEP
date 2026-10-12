# Project Writeup Implementation Plan

> **For agentic workers:** Use `superpowers:executing-plans` for native drafting in the current checkout. Steps use checkbox syntax for tracking. The author-selected review of the HTML planning page is a mandatory checkpoint before Task 2. An independent scientific review follows the complete draft, not each section.

**Goal:** Produce an advisor-facing account of the current single-modality and multimodal methodology, statistical performance, and biological findings, suitable for development into a manuscript.

**Architecture:** A standalone Quarto planning page establishes section briefs, evidence boundaries, and exhibits. After its approval, one main Quarto report includes two derivation appendix sources and renders to HTML/PDF using existing aggregate outputs and figures. A CSV source inventory records claim selectors and analytical caveats.

**Tech stack:** Existing Quarto and TinyTeX; Markdown and mathematical notation; existing aggregate CSV outputs and PNG/PDF figures. Any later tabulation code uses the project's existing R dependencies or base R; no new dependency or model-fitting runner is introduced.

**Spec:** The user-approved plan in this chat; `planning.qmd` provides the detailed reviewable structure. If the planning page is revised on review, record the revision here before drafting.

## Global constraints

- Work directly in `paper/writeup/` in the current project checkout; do not create a worktree.
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
- Propose checkpoint commits and obtain explicit confirmation before pushing. Do not send or publish the report.

## Review focus

1. **Mixed comparator outputs:** selecting a plausible but wrong CSV row must not change the identity of the main baseline. Check configuration, rank, penalty, cohort, and patient alignment.
2. **Theory versus evaluated algorithm:** equations must distinguish normalized plug-in single-modality prediction from the multimodal moment-based coefficient update, and state where survival actually enters.
3. **Overlapping cohorts and selection evidence:** the matched 50-donor set is not an additional independent cohort; external-performance evidence used for rank recommendation must be disclosed.
4. **Incomplete convergence and objective approximations:** describe iteration-cap fits as such, and the pruning objective as the implemented approximation rather than an exact survival expectation.
5. **Biological overinterpretation:** enrichment and unadjusted associations support program characterization, not causal mechanism, novel survival-guided discovery, or demonstrated methylation benefit.

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
| `review.md` | Complete-draft scientific/presentation findings, their resolution, and final verification evidence. Create after the draft exists. |

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
  Expected: all files exist; 40 unique inventory IDs; mean single-modality C rounds to 0.627; principal pooled difference is 0.0259; main single-modality simulation has 540 rows and 15 seeds; full multimodal sweep has 486 rows with 162 per variant; selected matched C values round to 0.6852 and 0.6757 and both selected fits are unconverged.
- [x] Inspect the rendered page in a browser: title/navigation, tables, mathematical notation, collapsed/expanded previews, and narrow-screen readability. Verify local section anchors in the HTML.
- [x] Update project pointers to say **planning page available; full report awaiting approval**. Suggest a checkpoint commit covering the verified planning artifacts.
- [ ] **Stop for planning-page review.** Approval of the implementation plan does not substitute for this separately selected checkpoint.

## Task 2: Methods and derivations

**Files:** Create `writeup.qmd`, `appendices/single_modality.qmd`, and `appendices/multimodal.qmd`; extend the source inventory.

**Consumes:** Approved section briefs; inventory entries S01–S08, S29–S30, S33–S34, S37–S38; realized settings from the evaluated runners/fits.

**Produces:** Common notation, two algorithm boxes, and derivations that consistently describe the configurations used in Results.

- [ ] Establish cross-format HTML/PDF YAML: article layout, 11-point PDF text, one-inch margins, modest line spacing, section numbering, contents, embedded HTML resources. Include the two appendix sources from the main report.
- [ ] Define notation/dimensions and write the evaluated single-modality likelihood, priors, normalized projection, working approximation, update order, convergence, prediction, rank recommendation, and factor classification.
- [ ] Explain the single-block multi-cohort setting and separate genomic offsets, baseline-hazard strata, and cohort-specific survival coefficients; identify the proteomics training input.
- [ ] Write the multimodal model and algorithm, including shared L, signed modality-specific F, feature intercepts/precisions, projection moments, scale canonicalization, and pruning/refitting.
- [ ] Derive single-modality beta/L/F/tau updates, EBNM handoff, and objective calculation; identify the implemented alpha_F=0 restriction and plug-in score approximation.
- [ ] Derive multimodal intercept/L/F/beta/tau updates, second moments, feature coupling, canonicalization, and removal criterion; distinguish the plugin pruning objective from exact expected Cox likelihood.
- [ ] Verify every algorithm step/equation against its implementation and evaluated configuration. Reuse independent algebra checks where they cover the relevant equation; do not change model code to make prose and code agree.
- [ ] Render both formats:

  ```bash
  quarto render paper/writeup/writeup.qmd --to html
  quarto render paper/writeup/writeup.qmd --to pdf
  ```

  Expected: both exit 0; both appendices, equations, and algorithm boxes are readable and correctly referenced.
- [ ] Propose a methods/derivation checkpoint commit after verification; continue native drafting without another section-level approval request.

## Task 3: Statistical results and biological findings

**Files:** Complete `writeup.qmd`; update `source_inventory.csv`; add adapted assets only when necessary in `figures/`.

**Consumes:** Approved exhibits; inventory entries S09–S28, S31–S32, S35–S36, S39–S40; the methods from Task 2.

**Produces:** Every numerical claim traces to a source selection, and every exhibit has a design, population, metric, and appropriate qualification.

- [ ] Write the cohort/design table with exact sample/event counts, platform identities, transformations, screening, replication, noise settings, comparator definitions, and inference procedures.
- [ ] Assemble single-modality simulation results by scenario/rank, distinguishing recovery, merging, prognostic classification, and held-out C; label the simulation comparator accurately.
- [ ] Assemble five-cohort external performance with regularized baseline rows from S10 and principal pooled difference from S11; keep S12/S13 as supporting comparisons. Explain the role of external-performance evidence in the rank recommendation.
- [ ] Characterize all four retained programs with gene weights, enrichment, subtype correspondence, overlap, and external association. Explain the sign conventions and exploratory nature of interpretation.
- [ ] Assemble the full multimodal simulation summaries by variant/scenario/noise/rank; specify that 486 includes all three variants. Do not replace the selected variant's metrics with an across-variant average.
- [ ] Assemble matched-donor comparisons and all starting ranks using explicit S28 selectors; include convergence and primary versus broader-histology sensitivity. Trace the documented established-model matched-donor comparison to its saved-score derivation before using it as a verified quantitative result.
- [ ] Describe all seven multimodal programs, modality shares, correspondence, and unadjusted survival associations. Identify CpG interpretation as incomplete and evaluate the evidence for added methylation information separately from added prediction.
- [ ] Write Interpretation and remaining evidence in narrative prose, distinguishing planned feature-count, filtering, weighting, and supervision studies from completed analyses.
- [ ] Verify selected results against CSVs, denominators, comparator identity, orientation conventions, and inventory; render both formats successfully and check no private patient file is a render dependency.
- [ ] Propose a full-draft checkpoint commit after verification.

## Task 4: Complete-draft review and repair

**Files:** Review the full writeup and appendices; create `review.md`; modify the authored sources/assets where substantive findings require correction.

**Consumes:** The draft, source inventory, current code/configurations, this plan, and the five Review Focus conditions.

**Produces:** A report with resolved substantive accuracy/presentation findings and an auditable verification record.

- [ ] Audit each quantitative claim's recorded path and selector; check source rows, model identity, cohort sizes, uncertainty, denominators, scenario grouping, and rounded values.
- [ ] Render both formats. Search for unresolved cross-references, missing assets/citations, and visible editorial placeholders; inspect PDF pages visually for overflow, cropped figures, poor table breaks, and unreadable math.
- [ ] Dispatch one independent scientific reviewer with the full draft, inventory, plan, and Review Focus. Request findings on unsupported claims, omissions, theory/implementation mismatch, biological interpretation, and coherence; do not request implementation of new analyses.
- [ ] Resolve substantive findings, recording their evidence and disposition in `review.md`. If evidence is unavailable, narrow or remove the affected claim and describe the actual evidence gap.
- [ ] Rerender both formats after repairs, repeat the checks affected by changes, and confirm the 15–20-page main-text target without sacrificing needed qualifications.

## Task 5: Documentation and delivery

**Files:** Update README, CLAUDE, PROJECT_STATUS, ROADMAP, and DECISIONS; finalize the inventory and review record.

**Consumes:** Verified report outputs and review disposition.

**Produces:** Accessible HTML/PDF links and sufficient context for a future collaborator to maintain the report.

- [ ] Replace planning-only status with the actual draft/review state; add render instructions and exact artifact locations. Record only completed work as completed.
- [ ] Verify documentation links and the `AGENTS.md` symlink; inspect the final diff for unrelated changes and accidental private-data inclusion.
- [ ] Deliver the report and explain the main organizational/mathematical choices, inputs/outputs, and remaining scientific limitations. Ask which parts warrant deeper review, in keeping with the project's post-implementation review practice.
- [ ] Suggest the final checkpoint commit. Ask about pushing only if there are unpushed commits, and never push without explicit confirmation.

## Progress record

- 2026-10-12: user approved native execution with an HTML planning-page review before drafting. Task 1 started in the existing checkout; the working tree was clean.
- Task 1: source inventory contains 40 entries; all source paths and Markdown links exist. Selected numerical checks passed: mean single-modality C=0.62668; principal difference=0.0259 (0.0002, 0.0498); 540 single-modality simulation rows, 15 seeds; 486 multimodal fits, 162 per variant; selected K=7 matched values and unconverged flags reproduced.
- Task 1: Quarto HTML rendering passed using RStudio's existing ARM-native installation. HTML inspection found valid anchors/links and embedded image, script, and stylesheet assets; AGENTS.md still resolves to CLAUDE.md. Browser checks confirmed loaded embedded figures, rendered math, functioning previews, white title text on the dark banner, and no document overflow at 390-pixel and 1280-pixel viewports after the narrow-screen table repair.
- Task 1: rendering adjustment is confined to the command used locally; no dependency installation, global setting change, or report-format change was needed. The working command is documented in the page and README.
- Execution choice: this deliverable is documentation, so verification consists of source/claim checks, rendering, and visual inspection rather than new model tests. No model code is changed. The explicit user review checkpoint overrides the skill's general continuous-execution instruction.
- Tasks 2–5 remain pending; no report or derivation appendix has been drafted.
- Task 1 deliverables verified and presented for the required planning-page review. No commit or push has been made. Suggested checkpoint: `Add advisor-facing writeup planning page and source inventory`.
