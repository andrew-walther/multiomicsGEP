# November research documents: completion review

Author: Andrew Walther

Date: 2026-11-12

## Scope and completed work

The paused document update is complete. The source update is committed in
`multiomicsGEP` at `ede07ead`. The installed hook synced Chapter 4 to
`bios-dissertation` at `5876e34`. The tested full-prelim rebuild is committed
there at `bfed766`.

- Chapter 4, `paper/prelim/project3-ssbmf.qmd`, now includes the rank-5
  weight search, the methylation cross-count decision, implementation checks,
  two-modality simulations and the descriptive matched-TCGA comparison.
- The November 13 writeup replaces the remaining rank-7 working-fit exhibits
  with rank-5 versions. Its source ledger and README describe the current inputs.
  Historical simulations and prior/input comparisons keep their original ranks.
- The November 13 progress chapter, current-state page and notebook index
  report the completed work. The evidence checklist is 13 of 24, counting
  completed diagnostic results independently of model adoption.
- Current project summaries and the roadmap reflect completed evaluation,
  rather than an ongoing simulation or unfinished weight search.

No patient data, fitted checkpoints, new model fits or analytical choices are
part of this update. The working model remains rank-5 YFB with unit survival
weight. The descriptive joint fit is not validated by simulation and is not adopted.

## Independent review

A fresh reviewer compared the document sources with the aggregate CSVs. The
checks covered program coefficients, genes, activity and classification,
rank correspondence, external counts and donor intervals, training cutpoints,
risk groups, subtype-adjusted hazard ratios, comparator differences, weight
selection and convergence, both simulation reruns, matched-TCGA results and
methylation cross-count tracking. The final reread found no blocking
scientific or numerical issue.

Corrections made during review:

- Replace a stale rank-7 EBMF label and distinguish the historical pruning
  simulation's seven starting columns from the working fit's five.
- Name precision medians and their upper bounds correctly. The active-program
  range in Chapter 4 is 4.9 × 10⁻⁷ to 1.2 × 10⁻⁶ across its saved fits.
- Report small performance differences with uncertainty, without an equivalence claim.
- Distinguish the current training-only selection rule from earlier development
  that used external concordance for the screen and scaling.
- Retain the distinction between the capped expression reference in the joint
  comparison and the separately restarted expression fits in the methylation study.
- Repair the progress chapter's findings link.

## Rendering and verification

Renders use RStudio Quarto 1.10.18 with TinyTeX. Both writeup formats evaluate
the saved aggregate checks. Rendering uses no patient-level input or model fit.

| Artifact | Pages | Verification |
|---|---:|---|
| Standalone Chapter 4 | 52 | PDF rendered, text has no unresolved references, changed result pages visually inspected |
| Prelim-format Chapter 4 | 53 | Installed sync succeeded, source and 33 citation keys checked, joint table visually inspected |
| November 13 writeup | 45 | HTML and PDF rendered, aggregate checks passed, ten vector figure pairs verified, new result page visually inspected |
| Full prelim | 238 | Four chapters, three appendices, 133 bibliography keys and 45 assets verified by the builder. All six assembly tests passed. Joint result page visually inspected |
| Progress book | HTML | Full book rendered, final chapter rerendered, no unresolved-link warnings. Historical generated-page version churn was excluded from the commit |

PDF text checks found no unresolved cross-references or task placeholders.
The seven rank-5 writeup SVGs exactly match their committed exhibit sources.
`git diff --check` passed. No file under `code/` or modular update script was
changed, so the model test suite was not rerun for these document-only changes.
The chapter and writeup source ledgers identify the aggregate files and commands.

## Explanation of the document calculations

**Purpose.** The writeup reads saved aggregate evidence so its displayed numbers
follow the rank-5 results rather than a manually edited copy of the rank-7 text.

**Logic.** The input chunks select the rank-5 output directories. Existing
summary expressions compute the displayed means, paired differences, intervals
and labels. Existing assertions compare recomputed values with saved summaries.
The figure helper chooses a vector PDF for LaTeX or its SVG conversion for HTML.

**Design.** Historical sensitivity studies keep their own saved ranks. The
rank-5 study paths replace only working-fit evidence. This avoids conflating
new rank-5 results with the earlier simulation designs.

**Inputs and outputs.** Inputs are committed CSVs and aggregate figures.
Outputs are the rendered writeup, chapter and progress pages. No new dependency
was introduced. Existing base R and knitr calculations format and check the evidence.

**Limits.** The weight search's preprocessed all-pairs CV estimate differs from
the working model's within-cohort estimate with fold-specific scaling. Its
27 capped inner fits remain in weight selection. The joint comparison retains
three capped expression folds, is not validated by simulation, and differs in
loading-prior divergence accounting between grouped and pooled arms.

## Remaining research

The 10,000-CpG joint fit has not been checked specifically against the
methylation-specific direction described at 100,000 CpGs. A larger joint fit
remains conditional on that result and a modality-balance check. Cohort-stratified
evaluation, scale-robust dead-program handling, stronger loading supervision,
external multimodal validation and final software release remain future work.
These are not unfinished document replacements.
