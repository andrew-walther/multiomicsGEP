# Cross-Modality Integration: Pre-Processing Reorder & Batch Effect Diagnostics

## Context

I have two related planning documents on the SSBMF (Survival-Supervised Bayesian Matrix Factorization) pipeline for joint TCGA + CPTAC training. The newer one — meeting notes from a team discussion on 4/27/26 — should drive the work. The older plan was written before that meeting and was never refined or executed; some of its proposed approaches are no longer the right first move given what came out of the meeting. Treat it as background context only.

I'd like you to read both documents and the current pre-processing implementation, then walk me through an integrated plan section-by-section so I can suggest edits before execution.

I will say **"go ahead"** when I'm satisfied and ready to begin execution. Until then, **do not modify any code or files** — we are in planning/discussion mode only.

## Documents to Read

### Primary — drives the plan
- **Meeting notes (4/27/26):** `/docs/plans/Meeting_Notes_Followup_4_27_26.md`

### Secondary — background context only
- **Earlier plan (pre-meeting, never executed):** `/Users/ajwalther/.claude/plans/what-are-the-key-luminous-treehouse.md`
  - Use to fill in detail where the meeting notes are silent
  - Where the two conflict, the meeting notes win
  - Some items in this older plan are no longer the right first approach — flag them rather than carrying them forward

### Current implementation — for code context
- `code/preprocess_desurv.R`

## Files Likely to Change

- `code/preprocess_desurv.R` — pre-processing step reorder
- `results/benchmark_sim/run_ssbmf_benchmark.R` — to consume revised pre-processing
- **New auxiliary script** for unsupervised factor analysis (PCA + EBMF/flashier) — also in `results/benchmark_sim/`

## What the Meeting Covered

### 1. Pre-processing reorder for multi-platform cohorts
Reordering the pre-processing steps when working with cohorts spanning multiple platforms (TCGA + CPTAC). Exact ordering is specified in the meeting notes document.

### 2. Naïve / unsupervised factor selection as a diagnostic
Consider simpler approaches that **omit the joint NMF/Survival objective**:
- **PCA** on the merged training data
- **EBMF** via the `flashier` R library

The goal is to understand whether the current batch-effect issue (all factor coefficients shrunk to zero under the Beta/ARD prior) is caused by the supervised objective itself, or is inherent in the merged data structure regardless of method.

### 3. Exploratory visualizations of subject loadings by cohort
Generate plots that reveal whether any factor's subject loadings stratify by cohort membership (e.g., consistently high values in TCGA subjects, low in CPTAC, or vice versa). This would directly identify factors that are encoding cohort/platform signal rather than biological/prognostic signal.

## Verification Checks

| Check | Current value | Question to answer |
|---|---|---|
| Common genes selected when merging CPTAC + TCGA (cap = 2000) | 838 | Why so few? Does the reorder recover more? |
| Non-zero (Beta-prior) factors selected | 0 — all shrunk to zero | Does the reorder yield any informative factors? |

## Overarching Goals

1. **Diagnose** the cross-platform batch effect — specifically, identify *which* factors (or whether all factors) are explaining cohort-difference signal rather than biological/survival signal.
2. **Test** whether the proposed pre-processing reorder resolves the issue when training on the joined cross-platform cohort.

## Hard Constraints

- ❌ **Do not implement batch correction via `limma::removeBatchEffect()`** — my advisor explicitly advised against this in the meeting due to concerns about generalizability when applying the trained model to new data sources.
- ❌ **Do not start coding** until I say "go ahead."
- ✅ **Engage me with clarifying questions** until we both have a comprehensive shared understanding before finalizing the plan.

## My Preferred Execution Order

Once we've agreed on the plan, I'd like to proceed in these phases, with **stop-and-verify checkpoints** between each:

### Phase 1 — Diagnostic visualizations (existing fit)
Produce plots of the current factor subject-loadings stratified by cohort. Quantify how many factors (and which ones) appear to encode cohort-difference signal vs. biological signal.
**Checkpoint:** review plots together before moving on.

### Phase 2 — Unsupervised factor benchmark (auxiliary script)
*Only if implementable quickly.* Create a new R script in `results/benchmark_sim/` that runs PCA and EBMF (`flashier`) on the merged cross-platform training data with **no survival supervision**. Compare resulting factor structure to what SSBMF produced. Does the batch-effect issue persist under unsupervised methods, or is it specific to our supervised objective?
**Checkpoint:** review comparison before moving on.

### Phase 3 — Pre-processing reorder + re-fit
Implement the pre-processing reorder per the meeting notes. Re-run SSBMF on the merged cohort. Re-check both verification metrics:
- Did `p` (gene count) recover beyond 838?
- Are any factor coefficients now non-zero?

**Checkpoint:** review metrics together before moving on.

### Phase 4 — Evaluate & decide
- If the reorder improves both training fit and external-cohort validation → update `README.md`, `CLAUDE.md`, and methods notes with the methodological change.
- If not → iterate. Surface candidate next steps (the older plan may be useful here, but treat its suggestions as options not defaults) for me to weigh in on.

## What I Need From You Right Now

1. Read the three documents listed above.
2. Walk me through your understanding of the integrated plan section-by-section.
3. Ask clarifying questions wherever the meeting notes are ambiguous, conflict with the earlier plan, or leave implementation details under-specified.
4. **Do not produce a finalized written plan document or touch any code until I say "go ahead."**
