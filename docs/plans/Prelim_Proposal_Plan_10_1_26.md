# Prelim Project Proposal Plan — Project 3 (SSBMF / multiomicsGEP) (10/1/26)

This is the plan for writing this project's proposal chapter for the prelim (Chapter 4 of the prelim document). It is separate from the research plan in `docs/plans/Working_Plan_10_1_26.md`. That plan produces the results; this one turns them into a chapter. The chapter source should also serve as the starting draft for a later journal manuscript.

**Deadlines:**
- Standalone abstract upload: 11/2.
- Proposal document: 11/16.
- Oral exam: 12/7.

## Where the chapter lives, and how it reaches `bios-dissertation`

**Recommendation: follow the Project 2 pattern.** In Project 2, a single source file in the project repo is the only file that gets edited. A post-commit hook generates the copy in `bios-dissertation`.
- Source: `SpatialCRT/projects/IncidenceDesign/paper/dissertation_chapter/Dissertation_Chapter.qmd`.
- Sync machinery: `tools/sync_to_prelim.sh`, `post_commit_hook.sh` and `prelim_transform.py` in that directory.
- Generated copy: `bios-dissertation/prelim/project-proposals/project2-incidence/draft/`. Its README says "GENERATED — do not edit here."

For Project 3:
1. Source (agreed 10/1): `multiomicsGEP/paper/prelim/` holds the prelim proposal chapter (e.g. `paper/prelim/project3-ssbmf.qmd` plus `figures/`). `multiomicsGEP/paper/thesis/` will hold the later thesis chapter and manuscript drafts.
2. Target: `bios-dissertation/prelim/project-proposals/project3-ssbmf/`. It currently holds only a placeholder `project3-ssbmf.qmd`, which becomes a pointer stub, plus a generated `draft/`.
3. Adapt the Project 2 sync tools rather than writing new ones. Two different sync approaches across the dissertation would be harder to maintain than one.
4. The thesis chapter and journal manuscript are written from the prelim source in `paper/thesis/`.

**Why not a symlink:** a symlink that points from one repo into another breaks on a fresh clone and when the folder layout changes. It also doesn't carry the prelim-class transformations that Project 2's `prelim_transform.py` applies. A generated copy works with both repos' git history.

Settled 10/1: the user agreed to this structure so the two repos' versions stay the same.

## Content constraints carried over from `bios-dissertation`

1. **α_F and the motivation (`bios-dissertation/ROADMAP.md`, 2026-09-10 and 2026-10-01 entries):**
   - The literature review motivates SSBMF as recovering low-variance prognostic programs that unsupervised factorization misses.
   - With `alpha_F = 0`, survival does not reach F, so the current method does not yet show that.
   - Before 11/16, one of the following must be true:
     - (a) the α_F experiments (Working Plan WS4.2) show that supervision changes F, or
     - (b) the Introduction and Methods make only the claims that the `alpha_F = 0` model supports.
   - The abstract's Chapter 4 sentences are already worded for (b) (`bios-dissertation/prelim/abstract/notes/abstract-fact-sheet.md`).
   - What the current model does: with α_F fixed at 0, YFB keeps factors that explain a lot of expression variance and are also associated with survival. It cannot find low-variance prognostic programs, which we suspect exist. A simulation with a planted low-variance prognostic program (Working Plan WS4.2) tests this directly. The goal is to remove α_F rather than tune it.
   - **This is the one research result the proposal's framing depends on.** In the working plan it is optional for 10/2; for the proposal it has a deadline.
2. **Error-in-variables correction** (comment in `templates/skeleton/prelim/chapters/04-project3-ssbmf.qmd`):
   - Methods explains that β's precision uses E[L²] (or E[Z²]) rather than the squared posterior mean, which corrects the attenuation that literature review §1.3.3 motivates.
   - Methods cross-references §1.3.3 and does not repeat its motivation.
3. **Structure:** Amber's DeSurv chapter has Introduction / Materials and Methods (Model, Hyperparameter selection/CV, Simulation studies, Real-world datasets) / Results / Discussion. The skeleton chapter already uses that structure.
4. **Repo writing rules:** write for a biostatistician reading the project cold, with no internal session labels and no mannered prose.

## Outline

### 1. Introduction
- **Motivation:** gene expression programs in PDAC, and prognostic subtypes (basal-like/classical).
- **Gap in the literature:** two-step approaches (unsupervised factorization, then Cox) and DeSurv (penalized NMF + Cox). Missing so far:
  - a Bayesian joint model with uncertainty in the factors and shrinkage priors that choose K;
  - joint modelling of more than one omics modality alongside survival.
- **Contributions:**
  - (i) a joint empirical-Bayes matrix factorization + Cox model (YFB, η = (YF)β) fit by CAVI;
  - (ii) K selection by over-specifying K and shrinking;
  - (iii) external validation across 5 PDAC cohorts and several platforms;
  - (iv) a multimodal extension (expression + methylation) with modality-specific priors and shared patient loadings.
- **Outline of the remaining sections.**
- Reuse from literature review Part 1: summarize and cross-reference, don't copy. Topics 9–11 of `prelim/literature-review/outline.md` are the relevant background.

### 2. Methods
2.1 **Single-modality model (YFB)**
- Likelihood, priors (EBNM point-exponential on L and F; Normal empirical-Bayes prior on β) and the survival link η = Ẑβ, with Ẑ = Y F̄ D⁻¹. Use the definition of ZF written in Working Plan WS4.1.
- CAVI updates: summarize here and put the full derivations in the supplement.
- The error-in-variables point from constraint 2.
- α and α_F: say what they are, and state the α_F status from constraint 1.
- Optional cohort terms: the `cohort_id` offset, `strata_id`, `beta_cohort_id`.

2.2 **Multimodal model**
- Y_m = L F_mᵀ + E_m, with L shared across modalities and a prior family for each modality's F_m.
- η = Σ_m Y_m F_m β.
- The factor rescaling step, and how K_eff is chosen.
- How this differs from the multi-cohort setting: cohorts add columns (more samples), modalities add rows (more features).

2.3 **K selection and hyperparameters:** the K_init consensus (ELBO/BIC/held-out log-likelihood/CV C) followed by `classify_factors()`, plus the 1-SE rule.

2.4 **Simulation designs:** single-modality (multi-cohort) and multimodal scenarios.

2.5 **Data and preprocessing:**
- The 7 PDAC expression cohorts and DeSurv gene selection.
- Matched TCGA/ICGC expression + methylation: transforms, CpG matching, imputation, screening.

2.6 **Evaluation:** external C-index, ΔC with a pooled bootstrap, factor recovery, survival-active classification.

### 3. Results (narrative headings)
Each heading should state a finding. **These are candidates, and each one can stay only if its result holds.** Rewrite or drop them as results come in.
- 3.1 Joint fitting matches or improves external discrimination relative to two-step factorization and DeSurv across five held-out PDAC cohorts. *(Have: mean external C ≈ 0.627.)*
- 3.2 Over-specifying K and shrinking returns two survival-associated programs (basal-like/adverse and classical/protective). *(Have: K_eff = 2, pathway enrichment.)*
- 3.3 Under-specifying K merges true programs, and over-specifying does not. *(Have: 9/4 simulation.)*
- 3.4 Cohort-specific survival coefficients locate where a program's effect is concentrated. *(Have: the CPTAC-specific finding for Program 3.)*
- 3.5 [Depends on α_F] Supervision recovers low-variance prognostic programs that unsupervised factorization misses. *Keep only if WS4.2 shows it.*
- 3.6 [Multimodal, simulation] Modality-specific priors recover shared factors and survival-active factors. *(Depends on Working Plan WS2.)*
- 3.7 [Multimodal, real data] Combining methylation with expression in TCGA: ICGC external C compared with the expression-only model. *(Depends on WS3. Report whatever the result is, including no improvement.)*

**Application:** fold this into Results as the real-data part. Present simulation results first, then the PDAC cohorts. No separate Application section unless the committee format requires one.

### 4. Discussion
- **What the findings mean:** for PDAC subtyping and for joint versus two-step modelling.
- **Limitations:**
  - α tempering and α_F = 0, unless resolved;
  - point-exponential priors don't prune whole factors (multimodal);
  - small matched multi-omics cohorts (TCGA n ≈ 150, ICGC n ≈ 50–67 with survival);
  - mixed histology in ICGC;
  - methylation transform choices.
- **Future work / continuing work:**
  - State plainly that the multimodal extension is in progress.
  - Remaining items: ARD-type pruning or a prior on K, the full 305k-CpG fit on Longleaf, more modalities, and the journal manuscript.

### Supplement
Full CAVI derivations (the single-modality derivation PDFs and `derivations/multimodal_YFB/`), additional simulation tables, and preprocessing details.

## Figures and tables to plan
1. Schematic of the model: shared L, F for each modality, and the survival link.
2. ΔC bar chart against the baselines across the 5 cohorts (already exists; regenerate in the chapter's style).
3. Results on under- and over-specifying K (already exists).
4. Program heatmaps and pathway enrichment for the two programs.
5. Multimodal simulation: K_eff and survival-active recovery against K_init.
6. Real-data table for TCGA → ICGC: C for each model, with n and events.
7. Optional: PCA/t-SNE of L (Working Plan WS4.6).

## Steps
1. Create `paper/prelim/` (and an empty `paper/thesis/`) with a skeleton `.qmd` taken from `templates/skeleton/prelim/chapters/04-project3-ssbmf.qmd`.
   - Verify: it renders on its own.
2. Write out the outline above as section headings with bullet points under each.
   - Verify: one review pass against `PROJECT_STATUS.qmd` and `DECISIONS.md` that every listed claim has a source.
3. Draft Methods first (2.1–2.2), since they are the least likely to change. Reuse the 10/2 progress-book chapter's model-form section (Working Plan WS5 §3).
4. Draft Results 3.1–3.4 from existing results. Add 3.5–3.7 as the working plan produces them.
5. Draft the Introduction and Discussion last, after the α_F question is settled one way or the other.
6. Adapt the Project 2 sync tools (into `paper/prelim/tools/`) to generate `bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/`, and turn `project3-ssbmf.qmd` into a pointer stub.
   - Verify: the generated draft renders in the prelim class.
7. Run a review pass of the full chapter (same approach as Chapter 2) before 11/16.

## Progress log
- 2026-10-01: Plan written. Nothing drafted yet.
- 2026-10-01 (late):
  - Steps 1–2 done, and a first full draft written: `paper/prelim/project3-ssbmf.qmd`, which renders to 17 pages. `paper/thesis/` is created (empty).
    - Introduction, Methods (single-modality and multimodal), and Discussion are drafted.
    - Results 3.1–3.4 report the established single-modality findings. 3.5–3.6 (multimodal) are marked [PENDING].
  - Bibliography: `paper/ssbmf-refs.bib`, entries copied verbatim from the master bib. The bib check passes (24 keys).
  - Sync tools adapted in `paper/prelim/tools/` from Project 2:
    - Chapter 4, `project3-ssbmf` paths.
    - Quarto cross-references are no longer counted as citations.
    - The declarations file is optional.
    - **Not yet installed or run.** Installing the post-commit hook and the first sync, which commits in bios-dissertation, need the user's go-ahead.
  - Author list: to confirm.
- 2026-10-02: An independent review was done (`docs/reviews/2026-10-02_prelim_proposal_review.md`).
  - Its accuracy must-fixes M1–M8 are applied. The main changes:
    - Fair two-step baselines: +0.026 against the K = 40 LASSO baseline, +0.011 against the K-matched one.
    - The "joint fitting" framing is removed; it conflicted with the α_F constraint.
    - Multimodal results are reported as parity with the established model (0.685 on the same 50 donors). The draft says that no factors were pruned on real data and that methylation added little.
    - The DeSurv head-to-head comparison is moved to future work.
  - Should-fixes S1, S5, S6, S7, S10, S12, S13 and part of S14 are also applied. The draft is now 25 pp.

### Remaining from the review (before 11/16)
1. **Structure.**
   - Specific aims for the remaining work, each with a testable criterion.
   - A timeline to the defense.
   - A simulation design table for both studies.
   - A short abstract or summary paragraph at the top.
   - A methods paragraph on the comparators.
   - Prediction and preprocessing for new cohorts.
   - Initialization and convergence rules for both fitters.
2. **Prose.** Rewrite Methods, Results and Discussion as paragraphs; keep bullets only for parameter lists.
3. **Figures.** Relabel the code-style labels:
   - Heatmap columns: P3 is not Program 3.
   - ΔC figure legend: all_shared, hybrid and nothing_shared.
   - Per-cohort β figure: add the sign convention.
   - Multimodal K figure legend.
   - Remove the embedded titles.
4. **Citations.**
   - Fix the Young2026 entry (@misc → @article with full details).
   - Add flashier, the cohort papers (CPTAC, Puleo, Dijk), KNN imputation (Troyanskaya 2001), bi-cross-validation (Owen & Perry 2009), BIC (Schwarz 1978), the 1-SE rule, the bootstrap (Efron & Tibshirani), Weibull survival simulation (Bender 2005), and MOFA / iCluster for Contribution 4.
   - Split the SEER sentence.
5. **Remaining should-fixes.**
   - S2: report false survival-active factors in §3.5.
   - S3: tempered and conditional β intervals; partly done.
   - S8: methods vs runs (K grid, noise, n, p); update when the full sweep finishes.
   - S9: pooling description; done.
   - S11: convergence statement; done.
6. **Labels.** Replace internal labels ("New framework", "YFB" undefined, "the frozen model") with descriptive names. Define EBMF at first use.
7. **After the full sweep and remaining starting ranks finish:** update §3.5–3.6 and remove the [Preliminary] and [PENDING] markers.
- 2026-10-02 (late): The structural items from the review are done.
  - Added: summary paragraph; Specific aims (Aim 1 complete; Aims 2–3 with proposed success criteria); timeline; simulation design table; methods paragraphs on prediction in new cohorts, comparators, and initialization/convergence.
  - Methods, Results and Discussion rewritten as prose.
  - Figures redrawn with descriptive labels (`paper/prelim/make_figures.R`): heatmap columns named; scenario names; sign convention on the per-cohort β figure; no embedded titles.
  - Tables are single-spaced (longtable fix). The draft is now 24 pp.
  - **Still open:**
    - Add the missing citations, which are marked with HTML comments in the source. Drafted, verified entries are in `paper/prelim/notes/missing-references.bib`. Adding them requires editing the bios-dissertation master bib, which needs the user's approval.
    - Confirm the proposed thresholds in the success criteria.
    - Confirm the defense date.
    - Confirm the author list.
- 2026-10-02 (night): The 11 verified references were added to the bios-dissertation master bib (`ce3f6fe`, not yet pushed). DOIs were re-checked against Crossref.
  - Young2026 (DeSurv) is now `@unpublished` (submitted, not yet published).
  - All citation markers in the proposal are replaced with real citations: 37 keys, and the bib check passes.
  - Prelim author: Andrew J Walther; prelim exam on 2026-12-07.
  - Planned manuscript authors: Andrew J Walther, Amber Young, Yusha Liu, Naim Rashid.
- 2026-10-02 (late night): Aligned the introduction with the literature review (cites §1.3.1–1.3.6 instead of restating them; 27 cited keys). Following the author's feedback:
  - Specific aims and success criteria replaced by a Remaining work subsection in Methods.
  - Lists rewritten as prose.
  - Tables and figures converted to Chapter 3's raw-LaTeX float style, with embedded KM titles removed.
  - Sync to bios-dissertation working: `tools/sync_to_prelim.sh`, with a post-commit hook (`tools/install_hooks.sh`) that regenerates and commits Chapter 4 whenever the chapter, its figures, its tools or `ssbmf-refs.bib` change. It never pushes. The prelim draft is 23 pp.
  - The standalone introduction for the journal manuscript is planned in ROADMAP.md.
