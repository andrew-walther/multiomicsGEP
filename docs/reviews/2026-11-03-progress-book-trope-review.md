# Progress book trope review (human-writer, review only)

Date: November 3, 2026. Scope: `docs/progress_book/current-state.qmd`, `index.qmd` and the 15 files in `chapters/`, about 36,600 words. No `.qmd` file was edited. The patterns come from the `human-writer` skill: `SKILL.md`, `voice-profile.md` and the new `ai-tropes.md`. The AGENTS.md rules on documentation audience and mannered prose also apply.

## Question for Andrew before any edits

The dated chapters are meeting records. Which files should be rewritten?

1. **Current files only**: `current-state.qmd`, `index.qmd`, and the chapters from October 26 onward (10-26, 10-30, 11-06). The meeting records stay as they were presented.
2. **Current files plus a mechanical sweep of the older chapters**: replace dashes and semicolons, restore missing spaces, and fix the factual inconsistencies listed below. Leave wording and framing alone.
3. **Everything**, including trope rewrites in the older chapters.

`_2026-10-15-methods-record.qmd` is a separate case. The October 16 chapter includes it, so it renders as part of that meeting. Most of its findings are about formatting and audience, not style (see below). It may need cleaning even under option 1.

## Summary

There are 477 flagged items. About 250 are the two hard rules: 203 prose semicolons and 49 rows of em or en dashes (some rows cover several dashes, about 63 in all). The rest are judgment calls. Most of them are repeated caveats and answers to objections nobody raised. Rhetorical tropes such as showy verbs, slogan headings and economic metaphors are rare.

**Which files most need work**

- **`_2026-10-15-methods-record.qmd` (194).** Lines 1–306 are clean apart from semicolons. From about line 315 on, the file reads like a pasted agent log:
  - spaces are missing between words and numbers ("exactly0", "seed20261026", "riskSD39086.312543");
  - it uses internal labels (Task2, Task5-M3, Round2, rung, ladder);
  - it contains imperative instructions ("Continue the prescribed ladder", "must not be described as completed");
  - it has shouted capitals.
  
  Most of the 45 "Other" rows are these problems. This file breaks the AGENTS.md documentation-audience rule more than any style rule.
- **`2026-09-04.qmd` (65).** About 45 em dashes, including every "§N —" heading. It has several bullets of 60 to 200 words and stacked "not X, not Y" defenses. It also uses "honest read" and "correct and trustworthy".
- **`2026-08-21.qmd` (54).** 27 semicolons. Three caveats each appear three times: the self-directed factor-biology work, L being absent from the survival term, and restart stability.
- **`2026-10-30.qmd` (33).** This is the most recent large chapter, so it matters most for option 1. It has 18 semicolons and six answers to unraised objections ("Nothing was deleted", "Substituting either result silently would conceal..."). "Not additional gate datasets" appears three times and "no external analysis" four times. The attempt history at lines 222–246 is one 25-line block.
- **`2026-08-27.qmd` (27).** 14 dashes. "Disagreement shown, not smoothed over" appears about four times, and "costs nothing" is a repeated economic metaphor.

`current-state.qmd` (6) and `index.qmd` (6) are in good shape. `current-state.qmd` has no dashes or semicolons in prose. Its main issue is the zero-feedback ($\alpha_F=0$) caveat, which appears at lines 21–22, 50, 142, 161 and 207–209.

**Hits per file**

| File | Hits | Main patterns |
|---|---:|---|
| `current-state.qmd` | 6 | Other 3, Repeated caveat 2, Unraised objection 1 |
| `index.qmd` | 6 | Semicolon 4, Not X but Y 1, Long sentence 1 |
| `chapters/2026-11-06.qmd` | 3 | Voice habit 1, Repeated caveat 1, Announcing importance 1 |
| `chapters/2026-10-30.qmd` | 33 | Semicolon 18, Unraised objection 6, Repeated caveat 4, Long sentence 2, Lesson drawn 1, Other 1, Voice habit 1 |
| `chapters/2026-10-26.qmd` | 10 | Semicolon 4, Tutorial heading 2, Showy verb 1, Mannered prose 1, Recap box 1, Repeated caveat 1 |
| `chapters/2026-10-16.qmd` | 14 | Semicolon 10, Repeated caveat 2, Voice habit 1, Other 1 |
| `chapters/_2026-10-15-methods-record.qmd` | 194 | Semicolon 100, Other 45, Voice habit 21, Repeated caveat 6, Unraised objection 6, Not X but Y 5, Mannered prose 4, Showy verb 2, Tutorial heading 2, Long sentence 2, Recap box 1 |
| `chapters/2026-10-02.qmd` | 4 | Semicolon 2, Dash 1, Slogan heading 1 |
| `chapters/2026-09-18.qmd` | 13 | Dash 4, Semicolon 4, Repeated caveat 2, Other 1, Not X but Y 1, Long sentence 1 |
| `chapters/2026-09-04.qmd` | 65 | Dash 30, Semicolon 9, Voice habit 5, Unraised objection 4, Moralized precision 2, Long sentence 2, Announcing importance 2, Repeated caveat 2, Not X but Y 2, Stacked hedge 2, Other 2, Mannered prose 1, Counted package 1, Slogan heading 1 |
| `chapters/2026-08-27.qmd` | 27 | Dash 14, Semicolon 5, Economic metaphor 2, Lesson drawn 2, Repeated caveat 1, Announcing importance 1, Moralized precision 1, Long sentence 1 |
| `chapters/2026-08-21.qmd` | 54 | Semicolon 27, Voice habit 8, Repeated caveat 6, Unraised objection 3, Long sentence 2, Not X but Y 2, Showy verb 1, Other 1, Tutorial heading 1, Announcing importance 1, Personality verb 1, Filler 1 |
| `chapters/2026-08-03.qmd` | 7 | Semicolon 3, Repeated caveat 2, Not X but Y 1, Recap box 1 |
| `chapters/2026-07-15.qmd` | 17 | Semicolon 6, Mannered prose 3, Unraised objection 3, Filler 1, Long sentence 1, Repeated caveat 1, Other 1, Not X but Y 1 |
| `chapters/2026-05-28.qmd` | 10 | Semicolon 3, Other 2, Announcing importance 1, Not X but Y 1, Filler 1, Unraised objection 1, Mannered prose 1 |
| `chapters/2026-05-07.qmd` | 6 | Semicolon 5, Recap box 1 |
| `chapters/2026-04-29.qmd` | 8 | Semicolon 3, Other 2, Showy verb 1, Tour-guide transition 1, Lesson drawn 1 |
| **Total** | **477** | |

**Hits per pattern**

| Pattern | Hits |
|---|---:|
| Semicolon | 203 |
| Other | 59 |
| Dash | 49 |
| Voice habit | 37 |
| Repeated caveat | 30 |
| Unraised objection | 24 |
| Not X but Y | 14 |
| Long sentence | 12 |
| Mannered prose | 10 |
| Announcing importance | 6 |
| Tutorial heading | 5 |
| Showy verb | 5 |
| Lesson drawn | 4 |
| Recap box | 4 |
| Filler | 3 |
| Moralized precision | 3 |
| Economic metaphor | 2 |
| Stacked hedge | 2 |
| Slogan heading | 2 |
| Personality verb | 1 |
| Counted package | 1 |
| Tour-guide transition | 1 |

"Other" covers formatting problems: missing spaces, internal labels, imperative agent instructions, bare links and stale cross-references. "Voice habit" covers items in `voice-profile.md`: contractions in academic register, "very", "significant" in its non-statistical sense, and session jargon.

## Factual inconsistencies found in passing

These need checking against the source. They are not wording fixes, so no rewrite is proposed.

- `chapters/2026-07-15.qmd` lines 108 and 147 say "five independent methods" but name only four: gene-set enrichment, PurIST correlation, external hazard ratios and DeSurv overlap.
- `chapters/2026-05-28.qmd` describes an "18-configuration" benchmark as 5 strategies × 2 models × with or without a cohort indicator, which is 20, unless some combinations were dropped. The phrase "18/5 configurations" at line 59 is also unclear.
- `chapters/2026-08-21.qmd` lines 388–389 and 532–534 say the two genomics-only programs are not yet characterized. Lines 403–410 then report Program 5 as stroma and Program 6 as immune.
- `chapters/2026-09-04.qmd` lines 176 and 180 refer to the "second technical note in §3", but the variance note is the third (line 160).

## How to read the findings

Line numbers refer to commit `674927e`. The progress book had no uncommitted changes at review time. Rewrites keep every number, claim, notation and citation unchanged. Dashes in numeric ranges, compound names (Kullback–Leibler, Newton–Raphson), math, code chunks and figure-chunk options were not flagged. Standard statistical verbs ("converges", "assumes", "rejects") were not flagged. Each judgment-call row is a suggestion. Keep the original wherever it does real work.

# Findings by file

## `current-state.qmd`

Clear and well documented, with no dashes or semicolons in prose. The main problems are the zero-feedback caveat repeated across sections, a few defensive sentences about the progress count, and one imperative sentence in a descriptive paragraph.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 33 | Other | "...research goals. [DeSurv reproduction materials](https://github.com/rashidlab/DeSurv-paper)." (bare link as a sentence; same at 259) | "...research goals. The DeSurv figures are available in the [DeSurv reproduction materials](https://github.com/rashidlab/DeSurv-paper)." |
| 46-47 | Other | "For LB prediction, estimate a new patient’s activities from expression with the training loadings fixed." | "For LB prediction, a new patient's activities are estimated from expression with the training loadings fixed." |
| 51-52 | Other | "Its external discrimination establishes a predictive reference." | "Its external discrimination serves as the predictive reference." |
| 142 | Repeated caveat | "The frozen historical benchmark used zero loading feedback ($\alpha_F=0$) while still fitting survival coefficients." (also stated at 21-22, 50, 161, 207-209) | Cut here, since lines 21-22 and 50 already state it. |
| 165-167 | Repeated caveat | "The benchmark comparisons are descriptive and cannot isolate the effect of survival feedback." | Keep here, and drop the parallel point from line 52 ("so it is distinct from the two-step comparator"), or merge the two. |
| 201-205 | Unraised objection | "The count tracks this checklist, with each line counted once. ... Publication readiness depends on the supervision gate and the final evidence, rather than the count." | "Checked historical or diagnostic items have a saved result, even when it was negative. Open items describe evidence still needed for the survival-supervised method." |

## `index.qmd`

Short and mostly plain. The main issues are three prose semicolons and one long synopsis sentence.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 6 | Semicolon | "discussion are included directly; anything that needs the full statistical detail is linked out..." | "discussion are included directly. Anything that needs the full statistical detail is linked out..." |
| 11 | Not X but Y | "This notebook is the meeting-facing summary, not a replacement for that log." | "This notebook is the meeting-facing summary of that log." |
| 16 | Semicolon | "Chapters are listed in the sidebar in date order; click a date to jump straight to that meeting." | "Chapters are listed in the sidebar in date order. Click a date to jump to that meeting." |
| 22 | Semicolon | "still uses $\alpha_F=0$; a fully validated survival-supervised replacement remains an open research goal." | "still uses $\alpha_F=0$. A fully validated survival-supervised replacement remains an open research goal." |
| 29-33 | Long sentence | "Supervised Bayesian Matrix Factorization (SSBMF) jointly factorizes ... unlike standard unsupervised factorization (PCA, NMF), which has no reason to recover clinically meaningful structure." | "Supervised Bayesian Matrix Factorization (SSBMF) jointly factorizes a gene expression matrix and a patient survival outcome, so the resulting gene expression programs (GEPs) are chosen to be both genomically coherent *and* prognostically relevant. Standard unsupervised factorization (PCA, NMF) has no reason to recover clinically meaningful structure." |
| 42 | Semicolon | "(TCGA-PAAD + CPTAC training; 5 independent cohorts held out for external validation)" | "(TCGA-PAAD and CPTAC for training, with 5 independent cohorts held out for external validation)" |

## `chapters/2026-11-06.qmd`

A careful working chapter with plain prose and well-placed caveats. The few problems are small.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 12 | Voice habit | "- Open Items & Questions" | "- Open items and questions" |
| 105 | Repeated caveat | "The later larger-feature evaluation failed its combined recovery, false-active and convergence gate." (also implied at 82-83 and restated at 126-127) | Keep here and cut the repeat at 126-127 (next row). |
| 126-127 | Announcing importance | "...so recovery, false activity and convergence all matter when evaluating a stronger feedback setting." | "...so a stronger feedback setting should be evaluated on recovery, false activity and convergence together." |

## `chapters/2026-10-30.qmd`

The prose is dense and very defensive. Many sentences answer objections nobody raised ("Nothing was deleted", "rather than described as guaranteed..."). The "not additional gate datasets" and "no external analysis" caveats come back several times, and prose semicolons are common. The attempt history at 222-246 is one long paragraph that is hard to follow.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 61 | Lesson drawn | "**Balancing total loss does not guarantee meaningful loading supervision.**" | "**The survival share of loading precision is now audited directly.**" |
| 82-83 | Other | "...and tests that missing methodological claim directly." | "...and tests directly whether survival can guide the loadings." |
| 86 | Unraised objection | "The local development inventory is retained in full. Nothing was deleted." | "The local development inventory is retained in full." |
| 89-90 | Unraised objection | "The table uses anonymous inventory labels so the meeting narrative remains about methods rather than version-control names." | "The table uses anonymous inventory labels." |
| 123 | Semicolon | "previously logged some cohorts twice; each cohort is now logged once." | "previously logged some cohorts twice. Each cohort is now logged once." |
| 144-146 | Unraised objection | "A training-concordance sign correction is a reporting convention rather than a further likelihood optimization, and is frozen for held-out patients." | "Any sign correction from training concordance is a reporting convention and is frozen for held-out patients." |
| 147 | Semicolon | "...in the saved single-log fit are positive; their minimum is..." | "...in the saved single-log fit are positive, with a minimum of..." |
| 148-150 | Unraised objection | "An estimated prior can still collapse with weak evidence, which is a different issue from incorrectly storing a posterior second moment." | "An estimated prior can still collapse when the evidence is weak." |
| 157-158 | Voice habit | "A feature/event weight alone does not control that ratio." | "A feature or event weight alone does not control that ratio." |
| 205 | Semicolon | "...contains `r controls$n_validation`; the initial rank is..." | "...contains `r controls$n_validation`. The initial rank is..." |
| 222-246 | Long sentence | "The first attempt could not calibrate ... The original gate is retained." (one 25-line paragraph covering seven attempts) | Split into two paragraphs: the numerical repairs (222-238) and the unit-mismatch design correction (238-246). Consider one sentence per attempt, keyed to the attempt column of the table above. |
| 226-228 | Long sentence | "An intermediate attempt rescales normal-means problems, checks both posterior moments, reduces damping and returns the last coherent fit with a failure reason whenever a numerical calculation fails." | "An intermediate attempt rescales normal-means problems, checks both posterior moments and reduces damping. When a numerical calculation fails, it returns the last coherent fit with a failure reason." |
| 254 | Unraised objection | "Substituting either result silently would conceal a numerical problem." | Cut. |
| 262 | Semicolon | "...diagnose numerical behavior; they are not additional successful gate datasets." | "...diagnose numerical behavior and are not gate datasets." |
| 264-265 | Unraised objection | "Decreases are explicitly reported, rather than described as guaranteed coordinate ascent." | "Objective decreases are reported." |
| 268 | Semicolon | "...the relative signal variance; the corresponding prepared variance shares..." | "...the relative signal variance. The corresponding prepared variance shares..." |
| 276 | Semicolon | "...at its saved final state; inactive columns have undefined shares." | "...at its saved final state. Inactive columns have undefined shares." |
| 299 | Semicolon | "...with valid stored uncertainty; their denominators are shown explicitly." | "...with valid stored uncertainty. Their denominators are shown in the table." |
| 303-304 | Repeated caveat | "A loading-recovery result alone does not establish calibrated survival inference or stable joint fitting. The gate keeps those requirements together." (same as 65-67) | Cut, since the Takeaways already state it. |
| 312 | Semicolon | "...registered loading means; they are not additional gate datasets." | "...registered loading means." |
| 312 | Repeated caveat | "they are not additional gate datasets" (also 262, 388) | Cut here and at 388. State it once at 262. |
| 330 | Semicolon | "...end-of-sweep changes; it does not attribute the first loss..." | "...end-of-sweep changes. It does not attribute the first loss..." |
| 335 | Semicolon | fig.cap: "...divided by its absolute starting value; loading norms use log10..." | "...divided by its absolute starting value. Loading norms use log10..." |
| 353 | Semicolon | "...with intersection remain primary; survival pre-screening is reserved for an ablation." | "...with intersection remain primary. Survival pre-screening is reserved for an ablation." |
| 361 | Semicolon | "The wrapper retains the first duplicate symbol; the earlier median collapse..." | "The wrapper retains the first duplicate symbol. The earlier median collapse..." |
| 369-370 | Repeated caveat | "These decisions do not authorize another external analysis." (also 68, 125, 305) | Cut. Line 68 already says no new external cohort was scored. |
| 380 | Semicolon | "...remains under evaluation; the latent formulation has not been adopted..." | "...remains under evaluation, and the latent formulation has not been adopted..." |
| 386 | Semicolon | "Holding τ fixed isolates noise feedback; guarded variants also use..." | "Holding τ fixed isolates noise feedback. Guarded variants also use..." |
| 387-388 | Repeated caveat | "These are mechanism tests, not additional gate datasets." | "These are mechanism tests." |
| 397 | Semicolon | "...explicitly stalled fits; `r sum(targeted$success)` meet all screening conditions." | "...explicitly stalled fits, and `r sum(targeted$success)` meet all screening conditions." |
| 406 | Semicolon | fig.cap: "Each row uses common axes; both vertical axes are log10." | "Each row uses common axes, and both vertical axes are log10." |
| 429 | Semicolon | "...smaller than the generating noise; its exact posterior is correct..." | "...smaller than the generating noise. Its exact posterior is correct..." |
| 435 | Semicolon | "...still depends on posterior means; no fixed-generative-model ELBO claim is established." | "...still depends on posterior means, so it is not an ELBO for a fixed generative model." |

## `chapters/2026-10-26.qmd`

Readable and concrete, with good numbers behind each claim. The problems are a few showy or mannered phrases in the summary, tutorial-style run-in headings, a closing "What this means" that restates the opening summary, and semicolons in the decision list and table cells.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 5 | Showy verb | "The single-modality model with expression only still holds the benchmark..." | "The expression-only single-modality model is still the benchmark..." |
| 19-21 | Mannered prose | "Survival's weak pull on the loadings is therefore a property of the model, from the relative size of the survival and reconstruction likelihoods, and not a coding error." | "Survival therefore has little influence on the loadings because of the relative size of the survival and reconstruction likelihoods, not because of a coding error." |
| 28-29 | Semicolon | "...is acceptable; the criteria that define a working model; the proposed preprocessing alignments, and..." | Make it a bulleted list of the four decisions, or use commas: "...is acceptable, the criteria that define a working model, the proposed preprocessing alignments, and..." |
| 161 | Tutorial heading | "**What it means.**" | "**Interpretation.**" |
| 225 | Semicolon | table cell: "...or any real-data fit; worst change −9.9 × 10⁻¹³" | "...or any real-data fit, with a worst change of −9.9 × 10⁻¹³" |
| 226 | Semicolon | table cell: "...every frozen real-data fit; one cross-validation fold fit stopped at 1,000 sweeps" | "...every frozen real-data fit. One cross-validation fold fit stopped at 1,000 sweeps" |
| 232 | Semicolon | table cell: "Check passes; quickstart runs; vignettes not yet written" | "Check passes, quickstart runs, vignettes not yet written" |
| 275 | Tutorial heading | "**What this means.**" | "**Summary of the rebuild.**" |
| 275-278 | Recap box | "The rebuild fixes the algorithm: the objective is a true bound, every step improves it, fits converge ... It does not fix supervision, and the joint model is dominated by methylation." (restates 15-21) | Cut the restatement and start at "Both problems come from how the likelihood weighs survival against reconstruction..." |
| 279 | Repeated caveat | "How to set that balance is now the central question." (also 27-28) | Cut, since the decisions paragraph at 27-28 already says it. |

## `chapters/2026-10-16.qmd`

Accurate but stiff, with ten semicolons. The same caveats (no validated supervised replacement, the prior-consistency repair is unsettled, methylation benefit not established) appear in the summary, in "Progress since October 2" and again at the end. The included `_2026-10-15-methods-record.qmd` was not reviewed.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 10–14 | Voice habit | "The project write-up brings together the Bayesian molecular-program model, ... Its central question is whether molecular programs can..." | Boilerplate opener. "The project write-up asks whether molecular programs can distinguish tumor and microenvironmental signals while giving useful survival prediction in pancreatic cancer. It covers the model, simulation and cross-cohort prediction, feature sensitivity, and biological interpretation." |
| 23 | Semicolon | "...a basis for the manuscript narrative; their biological interpretation remains provisional." | "...a basis for the manuscript narrative, but their biological interpretation remains provisional." |
| 35, 46, 56 | Other | "the retained Single benchmark" ... "Single external C=0.627" ... "using the Single benchmark settings" | Internal label (AGENTS.md line 58). "the single-modality benchmark" |
| 44 | Semicolon | "...a numerical failure that is repaired experimentally; convergence remains unresolved." | "...a numerical failure that is repaired experimentally, but convergence remains unresolved." |
| 49 | Semicolon | "...for the expression-only framework on the same donors; both fits are unconverged." | "...on the same donors. Both fits are unconverged." |
| 49–51 | Repeated caveat | "No fully validated survival-supervised replacement or incremental predictive benefit from methylation has been established." | Already stated at lines 25–26 and 35–37. Delete. |
| 71 | Semicolon | "The uncertainty intervals overlap; no new feature count was selected." | "The uncertainty intervals overlap, and no new feature count was selected." |
| 83 | Semicolon | "...not all seven starting columns; constant or unmatched programs have no defined correlation." | "...not all seven starting columns. Constant or unmatched programs have no defined correlation." |
| 92 | Semicolon | "...contains **305,352**; its imputation exactly reproduces the frozen 10k matrix." | "...contains **305,352**, and its imputation exactly reproduces the frozen 10k matrix." |
| 95 | Semicolon | "These are descriptive alternatives; no filter or SD cutoff has been adopted..." | "These are descriptive alternatives. No filter or SD cutoff has been adopted..." |
| 106 | Semicolon | "...without changing the likelihood; it does not fix the fit." | "...without changing the likelihood, but it does not fix the fit." |
| 110 | Semicolon | "...while three are weakly matched; risk correlation is 0.245." | "...while three are weakly matched. Risk correlation is 0.245." |
| 129 | Semicolon | "...prognostic roles; fewer programs alone are not evidence of a better model." | Also borderline Lesson drawn. "Each retained program needs a reproducible molecular and prognostic role. Fewer programs alone are not evidence of a better model." |
| 131–133 | Repeated caveat | "An initial prior-consistency repair improves local prior updates, but its fitted outputs still move beyond the existing stability diagnostic when the stopping rule triggers. No new supervised configuration is adopted." | Close restatement of lines 32–35. Delete, or fold "No new supervised configuration is adopted" into the Next priorities sentence. |

## `chapters/_2026-10-15-methods-record.qmd`

The first half (lines 1-306) reads as a clean methods record, with semicolons as the main issue. From about line 315 on, the file turns into a pasted agent log: dozens of missing spaces between words and numbers, internal labels (Task2-3, Round2, M3, rung, ladder, handoff, prompt), imperative instructions to an agent, shouted capitals, and a near-constant "this is X, not Y" disclaimer after each result. No dashes are used as punctuation; every dash is a numeric range or a compound method name.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 6 | Semicolon | "It supports the mathematical/software audit; it is not the concise advisor-meeting update..." | "It supports the mathematical and software audit. It is not the concise advisor-meeting update..." |
| 6-7 | Unraised objection | "it is not the concise advisor-meeting update or a validated final method." | "The concise meeting update is in the summary above, and no method here is validated or final." |
| 6 | Voice habit | "mathematical/software audit" | "mathematical and software audit" |
| 8 | Semicolon | "Sections below describe their own checkpoints; later failures supersede..." | "Sections below describe their own checkpoints. Later failures supersede..." |
| 17 | Semicolon | "Fit at K = 7; 20 seeds per noise level." | "Fit at K = 7, with 20 seeds per noise level." |
| 23 | Semicolon | "...and the update overshot;" | "...and the update overshot." |
| 25-26 | Repeated caveat | "Later sections document nonconvergence in other evaluations." | Cut it. The callout at lines 8-9 already says later failures supersede early convergence statements. |
| 34 | Semicolon | "log hazard ratio 1.10; programs 1 and 2 carry the variance..." | "log hazard ratio 1.10. Programs 1 and 2 explain most of the variance..." |
| 34 | Showy verb | "programs 1 and 2 carry the variance and are not prognostic" | "programs 1 and 2 explain most of the variance and are not prognostic" |
| 60 | Showy verb | "for programs that already carry a nonzero coefficient" | "for programs that already have a nonzero coefficient" |
| 65 | Mannered prose | "so it is not the objective the updates climb" | "so the updates do not increase it" |
| 71 | Semicolon | "after every iteration; β̃_k is then a log hazard ratio per SD." | "after every iteration, so β̃_k is a log hazard ratio per SD." |
| 113 | Semicolon | ": C3 + C4 arms; the two values per cell are noise 0.25 and noise 1." | ": C3 + C4 arms. The two values per cell are noise 0.25 and noise 1." |
| 119 | Tutorial heading | "**How much can survival alone recover here?**" | "**Survival-only reference.**" |
| 124 | Semicolon | "...is open; C5 waits on that decision." | "...is open, and C5 waits on that decision." |
| 152-154 | Semicolon | "- recovery 85% and 95%;" (and the next two bullets) | Drop the trailing semicolons from the bullets. |
| 174 | Semicolon | "...reported α_F = 0 model; cohort-stratified paired bootstrap." | "...reported α_F = 0 model, from a cohort-stratified paired bootstrap." |
| 211 | Voice habit | "refers to reconstruction/predictor changes" | "refers to reconstruction and predictor changes" |
| 212 | Semicolon | "before a complete rank batch; there is no completed result table." | "before a complete rank batch, so there is no completed result table." |
| 214 | Voice habit | "Starting K remains 7 with the handoff configurations." | "Starting K remains 7 with the configurations selected earlier." |
| 216 | Semicolon | "with zero variance; that boundary was implemented literally..." | "with zero variance. That boundary was implemented literally..." |
| 227 | Semicolon | "fit did not converge; lasso has one undefined C..." | "fit did not converge. Lasso has one undefined C..." |
| 230 | Voice habit | "The prompt stopped without altering the design." | "The run stopped without changing the design." |
| 239 | Semicolon | "passed correctly; the baseline is allowed to stop below it." | "passed correctly, but the baseline is allowed to stop below it." |
| 241 | Semicolon | "`...two_step_rank_audit.csv`; no core fitter changed." | "`...two_step_rank_audit.csv`. No core fitter changed." |
| 243 | Voice habit | "The selected M4 multimodal configuration" | "The selected multimodal configuration" |
| 255 | Semicolon | "includes 50 donors; all-donor scoring includes 67 and is secondary." | "includes 50 donors. All-donor scoring includes 67 and is secondary." |
| 256 | Repeated caveat | "No prior, rank or stopping rule was tuned." | Cut. Line 246 already says the failure "is retained without tuning the prior, warm-up, rank or stopping rule." |
| 258 | Semicolon | "before selecting primary donors; the review correction preceded expression scoring." | "before selecting primary donors. The review correction preceded expression scoring." |
| 258 | Voice habit | "in the primary/all analyses" | "in the primary and all-donor analyses" |
| 260 | Semicolon | "would give 0.5; both indicate no discrimination." | "would give 0.5. Both indicate no discrimination." |
| 268 | Voice habit | "The five-prompt evidence, independent branch commits and limitations are gathered..." | "The evidence, branch commits and limitations are gathered..." |
| 268-274 | Recap box | "Protected pruning was not adopted... DeSurv was unavailable. The high-noise comparator table uses only ten of twenty datasets... its constant predictions provide no discrimination." | Keep only the pointer to the consolidated report. Each of these points was stated in the preceding four sections. |
| 273 | Semicolon | "reached its cap with all β moments zero; its constant predictions..." | "reached its cap with all β moments zero, and its constant predictions..." |
| 277-278 | Other | "respectively 17,17,15 at noise .25 and 19,17,17 at noise 1" | "17, 17 and 15 at noise .25 and 19, 17 and 17 at noise 1, respectively" |
| 302 | Semicolon | "independently checked; the current runner parses..." | "independently checked, and the current runner parses..." |
| 310 | Voice habit | "it was not installed or reimplemented under this handoff." | "it was not installed or reimplemented for this comparison." |
| 312 | Semicolon | "did not authorize a configuration change; K remains 7." | "did not authorize a configuration change, so K remains 7." |
| 312 | Repeated caveat | "The protected-pruning assessment did not authorize a configuration change; K remains 7." | Cut. Lines 213-214 already state this. |
| 320 | Semicolon | "`...round2-attempts-log.md`; no empirical result or new configuration is claimed yet." | "`...round2-attempts-log.md`. No empirical result or new configuration is claimed yet." |
| 334 | Semicolon | "strongly prognostic and protected; a different program is deleted." | "strongly prognostic and protected, and a different program is deleted." |
| 342 | Semicolon | "remains valid; K = 7/no-pruning is the approved fallback..." | "remains valid. K = 7 without pruning is the approved fallback..." |
| 344 | Voice habit | "See attempts4–5 in the Round2 attempt history." | "See attempts 4 and 5 in the attempt history (`docs/plans/2026-10-16-round2-attempts-log.md`)." |
| 355 | Other | "single protected arm:" | "Single-modality protected arm:" (same for "multimodal protected arm:" at line 368) |
| 382 | Semicolon | "already variance-corrected; the approved full-state remedy investigates..." | "already variance-corrected. The approved full-state remedy investigates..." |
| 385 | Other | "(finalK = 5,500 iterations,loading correlation .934102,output flagFALSE)" | "(final K = 5, 500 iterations, loading correlation .934102, output flag FALSE)" |
| 387 | Semicolon | "one intermediate cap; the 10/1 reference has no final cap..." | "one intermediate cap. The 10/1 reference has no final cap..." |
| 390 | Semicolon | "reported limitations; no model meets every success criterion." | "reported limitations. No model meets every success criterion." |
| 393-396 | Other | "Continue Andrew's approved Single full-state remedy and subsequently the prescribed Task2 training-only collapse diagnosis... do not tune1.96, rank 3 or the output tolerance." | "The next steps are the approved full-state remedy for the single-modality model, then the training-only collapse diagnosis. The 1.96 threshold, rank 3 and the output tolerance stay fixed." (Agent instructions do not belong in a methods record.) |
| 396 | Semicolon | "not a failed recovery criterion; do not tune1.96..." | Covered by the rewrite above. |
| 401 | Voice habit | "#### Task1-5 completed: full-state refits improve recovery but are not adopted" | "#### Full-state refits improve recovery but are not adopted" |
| 406 | Semicolon | "misses the recovery threshold (15/20); five final fits hit the 300-iteration ELBO cap;" | "misses the recovery threshold (15/20). Five final fits hit the 300-iteration ELBO cap." |
| 407 | Semicolon | "ELBO cap; high-noise median K is 1 at every starting K." | "ELBO cap, and high-noise median K is 1 at every starting K." |
| 422-423 | Other | "seed20261026 at noise .25 with K_init5/7/10 (final K3), and seed20261017 at noise1 with K_init7/10 (final K2). Intermediate caps total25" | "seed 20261026 at noise .25 with K_init 5, 7 and 10 (final K 3), and seed 20261017 at noise 1 with K_init 7 and 10 (final K 2). Intermediate caps total 25" |
| 425 | Semicolon | "`.../attempt5_posterior/`; `adoption.csv` records `adopted=FALSE`." | "`.../attempt5_posterior/`. `adoption.csv` records `adopted=FALSE`." |
| 429-433 | Mannered prose | "The listed Single ladder has no remaining applicable remedy... Continue Task2's separate multimodal diagnosis and prescribed ladder." | "No remedy remains on the planned list for the single-modality model... The multimodal diagnosis continues separately." |
| 432 | Semicolon | "require a new decision; the approved K7/no-pruning fallback remains available..." | "require a new decision. The approved fallback, K = 7 without pruning, remains available..." |
| 437 | Voice habit | "### Task2-1 completed TCGA diagnosis (2026-10-15)" | "### TCGA training diagnosis of the multimodal collapse (2026-10-15)" |
| 439 | Semicolon | "preceded all ten fits; its full suite passed584/0." | "preceded all ten fits, and its full suite passed 584/0." |
| 439-442 | Other | "passed584/0... exit0... took100–102 seconds... took2384–2556 seconds... K7" | Insert the missing spaces: "passed 584/0", "exit 0", "took 100–102 seconds", "took 2384–2556 seconds", "K = 7". |
| 443 | Semicolon | "`.../diagnosis1/`; real-data RDS files are ignored." | "`.../diagnosis1/`. Real-data RDS files are git-ignored." |
| 463 | Semicolon | "reaches exactly0 at outer99/inner1; all beta means and SDs become exactly0..." | "reaches exactly 0 at outer 99, inner 1. All beta means and SDs become exactly 0..." |
| 463-466 | Other | "exactly0 at outer99/inner1... stay0... outer54/inner3 and54/inner4... scale1" | Insert the missing spaces ("exactly 0", "outer 99/inner 1", "stay 0", "54/inner 4", "scale 1"). |
| 466 | Semicolon | "with zero slab weight; this is not a numerical scale-zero diagnosis." | "with zero slab weight." (see next row) |
| 466-467 | Unraised objection | "this is not a numerical scale-zero diagnosis." | "The prior scale stays at 1, so the collapse comes from the slab weight and not from a numerical zero scale." |
| 473-480 | Other | "unchanged500-iteration... exceed1e-4... are0... Normal2.23e-10... last100" | Insert the missing spaces ("unchanged 500-iteration", "exceed 1e-4", "are 0", "Normal 2.23e-10", "last 100"). |
| 478 | Semicolon | "Selected and no-start eta changes are0; fixed/framework eta changes are..." | "Selected and no-start eta changes are 0. Fixed and framework eta changes are..." |
| 480 | Other | "Every arm's last100 ELBO differences is positive" | "Every arm's last 100 ELBO differences are positive" |
| 486 | Semicolon | "not a tagged-code reproduction; no equality or superiority claim follows." | "not a tagged-code reproduction. No equality or superiority claim follows." |
| 485 | Not X but Y | "It is a current-code diagnostic refit, not a tagged-code reproduction" | "The refit uses the current code, so it does not reproduce the tagged benchmark." |
| 489-494 | Mannered prose | "Next causally applicable rung... Continue the prescribed ladder with a committed design..." | "The next planned remedy constrains the refitted shared Laplace slab probability below 1/K (R2)... The remaining remedies follow, each with a committed design and unchanged simulation and CV criteria." |
| 489-492 | Long sentence | "Next causally applicable rung: constrain ... (R2); R1 fixed-prior diagnostic already fails ... has not yet been performed." | Split after "(R2)". Second sentence: "The R1 fixed-prior diagnostic already fails the full-training credible-program and convergence requirements, but its simulation and CV evaluation has not yet run." |
| 490 | Semicolon | "below at1/K (R2); R1 fixed-prior diagnostic already fails..." | Covered by the split above. Also fix "below at1/K" to "below 1/K". |
| 492-493 | Other | "Those unrun parts must not be described as completed." | Cut. This is an instruction to the agent, and the previous clause already says the evaluation has not run. |
| 501 | Semicolon | "has run; there are no new scientific results." | "has run." |
| 501 | Repeated caveat | "there are no new scientific results." | Cut. The previous clause already says no DeSurv fit or prediction has run. |
| 501-502 | Voice habit | "The Round2 handoff records the software pins..." | "The handoff notes record the software pins..." or name the file. |
| 506 | Voice habit | "### Task4 completed five-comparator panel (2026-10-15)" | "### Five-comparator simulation panel (2026-10-15)" |
| 508 | Semicolon | "dd11923 and d353322; full suite599/0 before each." | "dd11923 and d353322, with the full suite passing 599/0 before each." |
| 510 | Semicolon | "could not open; its design and execution_status.md remain in round2/." | "could not open. Its design and execution_status.md remain in round2/." |
| 511 | Other | "completed with exit0:40 datasets,200 method rows" | "completed with exit 0: 40 datasets, 200 method rows" |
| 514 | Semicolon | "by Andrew's explicit decision; no arm or threshold changed." | "by Andrew's explicit decision. No arm or threshold changed." |
| 529 | Semicolon | "20 at least2; high noise2 zero-factor..." | "20 with at least 2. At high noise, 2 zero-factor..." |
| 529-532 | Other | "low noise0 zero-factor,0 one-factor,20 at least2... denominator is20... to10 datasets... in6" | Insert the missing spaces ("0 zero-factor, 0 one-factor, 20 with at least 2", "is 20", "to 10", "in 6"). |
| 532 | Semicolon | "has finite C in6; its undefined scores are retained." | "has finite C in 6, and its undefined scores are retained." |
| 535 | Semicolon | "is not Bayesian credibility; false-active isNA." | "is not a Bayesian credibility measure, so false-active is NA." |
| 537-542 | Other | "caps at300... seed20261026/noise .25... seed20261017/noise1... There are30" | Insert the missing spaces ("at 300", "seed 20261026", "noise 1", "There are 30"). |
| 538 | Semicolon | "seed20261017/noise1; both recover the planted program." | "seed 20261017 at noise 1. Both recover the planted program." |
| 541 | Not X but Y | "It is a comparator report, not adoption." | Cut. Line 514 and the section heading already make clear nothing is adopted. |
| 544 | Semicolon | "retained per row; no package/default was changed to silence it." | "retained per row. No package or default was changed to silence it." |
| 546 | Voice habit | "This completes the bounded Task4 rerun, including rank/denominator reporting. Task2's remedy ladder remains open and is the priority." | "This completes the comparator rerun, including rank and denominator reporting. The multimodal remedies remain open and are the priority." |
| 552-555 | Other | "fold1... Folds2/3... below1... is0 from outer167" | Insert the missing spaces ("fold 1", "Folds 2 and 3", "below 1", "is 0 from outer 167"). |
| 557 | Semicolon | "each coordinate and can revive; an unchanged-fold trace will establish..." | "each coordinate and can revive. An unchanged-fold trace will establish..." |
| 561 | Semicolon | "history exactly; only observation changes." | "history exactly, with only the recording turned on." |
| 586 | Semicolon | "remains recognizable; this is not removal of the program by pruning" | "remains recognizable, and the program was not removed by pruning" |
| 587-590 | Unraised objection | "This comparison is a description of saved fits, not a cross-validation estimate using a full-training loading matrix, an external selection criterion, or proof of benchmark reproducibility." | "This comparison only describes saved fits." |
| 593 | Voice habit | "#### Task2-2 final fixed-prior evaluation: not adopted" | "#### Fixed-prior evaluation: not adopted" |
| 595-596 | Other | "The100-fit local evaluation exited0:80 simulations and20 TCGA CV fits, with0 fitting errors and0 nonempty warning records. All100" | "The 100-fit local evaluation exited 0: 80 simulations and 20 TCGA CV fits, with 0 fitting errors and 0 nonempty warning records. All 100" |
| 609 | Semicolon | "unchanged K7 beds; K selection has not been adopted." | "unchanged K = 7 beds, and K selection has not been adopted." |
| 620 | Semicolon | "C=NA/risk SD0; no finite five-pair difference exists..." | "C = NA and risk SD 0. No finite five-pair difference exists..." |
| 622-629 | Other | "cap at500... are0,1,1,0,0... values2–5... has0 credible programs and a500 cap" | Insert the missing spaces ("at 500", "are 0, 1, 1, 0, 0", "values 2–5", "has 0 credible programs and a 500 cap"). |
| 623 | Semicolon | "also cap at500; reference credible-program counts..." | "also cap at 500. Reference credible-program counts..." |
| 630 | Other | "Thus R1 fails necessary real-data criteria despite the multimodal simulation passes." | "Thus R1 fails necessary real-data criteria despite passing the multimodal simulation criteria." |
| 634 | Semicolon | "after exact trace reproduction; its implementation is still private..." | "after exact trace reproduction. Its implementation is still private..." |
| 638 | Semicolon | "#### Task2-3 result: exact reproduction; loss at outer167/inner1" | "#### Fold trace reproduces exactly and loses all coefficients at outer 167, inner 1" |
| 640-652 | Other | "exited0 in2142.164 seconds... All13... complete500-iteration... cell84... Of2500... sweep831... and1.128647... Programs3 and4... All1670... from831 through2500... revive7/8/2/1/2... throughout2500" | Insert the missing spaces throughout this paragraph. |
| 649 | Semicolon | "have zero means/SDs; there is no later revival." | "have zero means and SDs, with no later revival." |
| 656-657 | Other | "retains K7... risk SD0/C=NA... unchanged500 cap" | "retains K = 7... risk SD 0 and C = NA... unchanged 500 cap" |
| 659 | Semicolon | "(one credible program); this is not a CV row." | "(one credible program). This is a training value, not a CV estimate." |
| 659 | Repeated caveat | "this is not a CV row." | Cut or keep once. The same disclaimer appears at lines 627 and 742-743. |
| 671 | Semicolon | "profile flag; it will not combine the unresolved shared-prior scale question." | "profile flag. It does not also address the unresolved shared-prior scale question." |
| 674 | Tutorial heading | "#### What the profiled-curvature implementation changes (scientific results pending)" | "#### Profiled Cox intercept implementation (results pending)" |
| 678 | Semicolon | "preserving dimensions/names; zero total curvature returns the original input." | "preserving dimensions and names. Zero total curvature returns the original input." |
| 690 | Semicolon | "are per-factor Normal; shared priors fail explicitly..." | "are per-factor Normal. Shared priors fail explicitly..." |
| 693 | Semicolon | "monitor/rule remains unchanged; the multimodal output rule also stays fixed." | "monitor and rule remain unchanged, as does the multimodal output rule." |
| 696-699 | Other | "records80 new simulation fits,10 new candidate CV fits... and10 unchanged... all5 finite folds" | Insert the missing spaces ("records 80", "fits, 10", "and 10", "all 5"). |
| 701 | Semicolon | "the actual model path; separate checks verify..." | "the actual model path. Separate checks verify..." |
| 706 | Voice habit | "### Task2-4 completed: profiled Cox curvature is not adopted (2026-10-15)" | "### Profiled Cox curvature is not adopted (2026-10-15)" |
| 708 | Other | "exit0:102 records,92 new fits and10 unchanged" | "exit 0: 102 records, 92 new fits and 10 unchanged" |
| 722 | Semicolon | "reaches300 iterations; its oscillation is documented above. All40" | "reaches 300 iterations, and its oscillation is documented above. All 40" |
| 723-724 | Not X but Y | "These are the approved recovery/false-positive criteria, not proof that protected K selection is resolved." | Cut. Nothing in the paragraph claims K selection is resolved. |
| 737 | Semicolon | "uses all five folds; no undefined reference is replaced or dropped." | "uses all five folds, and no undefined reference is replaced or dropped." |
| 739 | Semicolon | "formerly collapsed fold1; their credible-program counts are 0,1,0,0,0." | "formerly collapsed fold 1. Their credible-program counts are 0, 1, 0, 0 and 0." |
| 740 | Other | "reaches500 iterations" | "reaches 500 iterations" |
| 742-743 | Repeated caveat | "Frozen training C values are labelled historical training quantities, not CV estimates or external scores." | Cut. Lines 627-628 already say the frozen rows are non-CV references. |
| 754 | Other | "Eight of92 new fits cap" | "Eight of 92 new fits cap" |
| 756 | Semicolon | "the full multimodal fit; another five caps belong to..." | "the full multimodal fit. Another five caps belong to..." |
| 760 | Semicolon | "`training_full_with_benchmark.csv`; the latter two retain frozen rows." | "`training_full_with_benchmark.csv`. The latter two retain frozen rows." |
| 767 | Mannered prose | "Further Task2 work is paused by Andrew, not resolved or ladder-exhausted." | "Andrew paused further work on the multimodal collapse before it was resolved or the planned remedies were exhausted." |
| 771-779 | Other | "tag object0fea14903a75606e495267f2f35712eac10d9fc6 still peels to a352f17ca9cff523caa80756aa2af085e62551e1... No push, protected-branch edit or merge has occurred." | Move the hash and git details to the commit message or DECISIONS.md. Keep one sentence: "The full R suite passes (631 tests), and the benchmark files and write-up are unchanged against the tagged snapshot." |
| 772 | Semicolon | "prior-consistency fixtures); prepared observer syntax parses." | "prior-consistency fixtures). The prepared observer script parses." |
| 784 | Other | "against the10/2 meeting" | "against the October 2 meeting" |
| 791 | Semicolon | "does not pass non-inferiority; the limited gene-count experiment..." | "does not pass non-inferiority. The limited gene-count experiment..." |
| 792 | Other | "screened10k skewness plots" | "screened 10k skewness plots" |
| 795 | Semicolon | "before implementation; no new scientific fit has run." | "before implementation. No new scientific fit has run." |
| 796 | Other | "full suite631 passed,0 failed" | "full suite 631 passed, 0 failed" (same at 811, 855, 870, 934, 952) |
| 802 | Other | "the evidence since10/2" | "the evidence since October 2" |
| 807 | Semicolon | "a methylation-only program; immune structure has substantial methylation loading." | "a methylation-only program. Immune structure has substantial methylation loading." |
| 809 | Semicolon | "K7 itself is not a failure; the roles, stability and redundancy need evaluation." | "The roles, stability and redundancy of the seven programs need evaluation." |
| 809 | Unraised objection | "K7 itself is not a failure" | Cut. No one in the section claims it is. |
| 811 | Semicolon | "Rendered PDF visually checked; full suite631 passed,0 failed." | "The rendered PDF was visually checked, and the full suite passed 631, failed 0." |
| 814 | Voice habit | "### Task5-QC1 result: full-probe methylation variability completed" | "### Methylation variability across all TCGA probes" |
| 816 | Semicolon | "Preregistered in dc69f6d; serial run completed with exit0..." | "Preregistered in dc69f6d. The serial run completed with exit 0..." |
| 817 | Other | "All305492 available TCGA probes" | "All 305,492 available TCGA probes" |
| 833 | Semicolon | "defined sample SD; no probe is constant. There are29414" | "defined sample SD, so no probe is constant. There are 29,414" |
| 835 | Semicolon | "Of these,2130 occur in the screened reference; the unchanged cached imputation contains none." | "Of these, 2,130 occur in the screened reference, and the unchanged cached imputation contains none." |
| 839 | Not X but Y | "These are data descriptions, not survival results." | Cut. The section opening already says no survival fit occurred (line 819). |
| 840 | Semicolon | "strongly enriched for variability; that does not establish the best feature universe." | "strongly enriched for variability, but that does not establish the best feature universe." |
| 841 | Other | "would discard10 of the currently selected probes" | "would discard 10 of the currently selected probes" |
| 843 | Semicolon | "remains unverified; this is an explicit provenance gap." | "remains unverified, which is a provenance gap." |
| 846-848 | Other | "check on21 probes... transformed SD to5.55e-17" | "check on 21 probes... transformed SD to 5.55e-17" |
| 849 | Semicolon | "Input hashes remain unchanged; zero warnings." | "Input hashes remain unchanged, and there were no warnings." |
| 850 | Semicolon | "remains local/ignored; only aggregate tables, plots and... are committed." | "remains local and git-ignored. Only aggregate tables, plots and... are committed." |
| 860 | Other | "### Fixed-K expression curve: preflight complete, initialization decision pending" | Fine as a status heading, but the body (lines 862-870) needs spaces: "contains 17,900 genes on 273 patients", "reaches 2064, 5001, 10000 and 17900", "rank 4 of 7", "suite 651/651". |
| 867-868 | Other | "so pooled Single four classified versus joint seven does not isolate the methylation contribution" | "so comparing four classified single-modality programs with seven joint programs does not isolate the methylation contribution" |
| 874-879 | Other | "completes24/24... through17900... are2064:.6084,5001:.5773,10000:.6088 and17900:.6068... are4/4/4/3... is.505" | Insert the missing spaces ("completes 24/24", "through 17,900", "are 2064: .6084, 5001: .5773, 10000: .6088 and 17900: .6068", "are 4/4/4/3", "is .505"). |
| 877 | Semicolon | "mean external C=.626699; no external cohort was rescored..." | "mean external C = .626699. No external cohort was rescored..." |
| 877 | Voice habit | "Frozen D4 has historical mean external C" | "The frozen benchmark has historical mean external C" (D4 is an internal label) |
| 879 | Semicolon | "to the fresh baseline is.505; dense nonzero loading means..." | "to the fresh baseline is .505. Dense nonzero loading means..." |
| 885-888 | Other | "eligible305352-CpG... frozen10k... retain182934/221142 probes... full144-patient" | Insert the missing spaces ("eligible 305,352-CpG", "frozen 10k", "retain 182,934 and 221,142 probes", "full 144-patient"). |
| 887 | Semicolon | "retain182934/221142 probes; neither is adopted." | "retain 182,934 and 221,142 probes, and neither is adopted." |
| 899 | Semicolon | "riskSD.117264, no warnings/errors; the25k fit stopped..." | "risk SD .117264 and no warnings or errors. The 25k fit stopped..." |
| 899-900 | Other | "The10k fit took963.827s... the25k fit stopped after27.981s" | "The 10k fit took 963.827 s... the 25k fit stopped after 27.981 s" |
| 901-902 | Other | "Preserve NA for the failed fit." | "The failed fit is recorded as NA." |
| 906 | Semicolon | "risk correlation.954; both fits cap500." | "risk correlation .954. Both fits cap at 500." |
| 908-909 | Other | "norm shares58.9%/68.7% and loading correlation.993" | "norm shares 58.9% and 68.7% and loading correlation .993" |
| 911 | Semicolon | "PVE or survival fractions; newer biological labels are not inferred." | "PVE or survival fractions, and no new biological labels are inferred." |
| 915 | Semicolon | "### CpG numeric failure reproduced; manuscript sensitivity update requested" | "### Reproduced CpG numeric failure and requested manuscript update" |
| 917-921 | Other | "unchanged25k cell... after32.218s... outer6/inner1/program1:17 nonfinite... span5734 causes143/144... and17 risk-set sums" | Insert the missing spaces ("unchanged 25k cell", "after 32.218 s", "outer 6, inner 1, program 1: 17 nonfinite", "span 5734 causes 143/144", "and 17 risk-set sums"). |
| 917 | Voice habit | "Task5-M2 repeats the unchanged25k cell" | "The rerun repeats the unchanged 25k cell" |
| 919 | Semicolon | "Cox working weights/numerators; projection moments are finite..." | "Cox working weights and numerators. Projection moments are finite..." |
| 929 | Semicolon | "cannot be estimated yet; the failed25k row stays undefined." | "cannot be estimated yet, and the failed 25k row stays undefined." |
| 931-933 | Long sentence | "The bounded manuscript revision adds those verified results, benchmark-linked tables and the CpG diagnosis without claiming a supervised replacement or changing the preserved benchmark tag/snapshot." | "The manuscript revision adds those verified results, the benchmark-linked tables and the CpG diagnosis. It does not claim a supervised replacement or change the preserved benchmark snapshot." |
| 934-935 | Other | "Additional numeric remedies are separate prospective work, not inferred from this reporting authorization." | "Further numerical remedies need their own approval." |
| 940-942 | Other | "shares are6.30e-11–3.03e-7... shares8.70e-7–1.82e-4... program4 has... shares2.03e-7/9.00e-10... precision,7.27e-5/4.40e-6" | Insert the missing spaces ("are 6.30e-11–3.03e-7", "shares 8.70e-7–1.82e-4", "program 4", "shares 2.03e-7 and 9.00e-10", "precision, 7.27e-5 and 4.40e-6"). |
| 944 | Semicolon | "capped terminal states, excluding priors; they do not identify trajectory causes..." | "capped terminal states, excluding priors. They do not identify trajectory causes..." |
| 949 | Voice habit | "Task5-M3 prospectively registers captured-point numerical verification..." | "A preregistered run checks the log-space calculation at the captured failure point, then fits one unchanged 25k case." |
| 952 | Semicolon | "Fresh full suite684 passed,0 failed; no scientific M3 run before its commit." | "A fresh full suite passed 684 tests with 0 failures. No scientific run took place before its commit." |
| 955-959 | Other | "### Stable-arithmetic CpG checkpoint and requested pause" followed by four blank lines | Remove the extra blank lines. |
| 960-961 | Other | "The SAME process99846 exits0 at its unchanged500 cap... scientific convergence FAILS." | "The same process (99846) exits 0 at its unchanged 500-iteration cap... but the fit does not converge." |
| 962-965 | Other | "runtime2660.496s (44.3min)... riskSD39086.312543, credible7/retained7... tol1e-4... outer500 is312267.6–520103.6... in178/499" | Insert the missing spaces ("runtime 2660.496 s (44.3 min)", "risk SD 39086.312543, 7 credible of 7 retained", "tolerance 1e-4", "outer 500 is 312267.6–520103.6", "in 178 of 499"). |
| 962 | Semicolon | "No warnings/errors; runtime2660.496s (44.3min)." | "There were no warnings or errors, and the runtime was 2660.496 s (44.3 min)." |
| 965 | Semicolon | "is312267.6–520103.6; logPL worsens in178/499 outer transitions." | "is 312267.6–520103.6. LogPL worsens in 178 of 499 outer transitions." |
| 966 | Semicolon | "need not be monotone in joint fitting; this observation does not establish a cause." | "need not be monotone in joint fitting, so this observation does not establish a cause." |
| 967-968 | Other | "High training C with these extreme, unconverged estimates is NOT a predictive gain or meaningful discovery." | "A training C this high from extreme, unconverged estimates does not indicate a predictive gain." |
| 969 | Unraised objection | "Saved terminal history is an aggregate export, not another experiment." | "The saved terminal history is an aggregate export of this run." |
| 971 | Semicolon | "Benchmark context (same144 TCGA training patients; different molecular inputs):" | "Benchmark context (the same 144 TCGA training patients, with different molecular inputs):" |
| 981 | Voice habit | "A1 executes only after authoritative M3 termination, exits0 and verifies unchanged inputs." | "The agreement check ran only after the log-space run finished. It exited 0 and confirmed the inputs were unchanged." |
| 982 | Semicolon | "Same144 patients/common10000 probes; all7 activity pairs match identical indices." | "It uses the same 144 patients and 10,000 common probes, and all 7 activity pairs match identical indices." |
| 983 | Semicolon | "Both fits are capped;25k arithmetic differs." | "Both fits are capped, and the 25k fit uses different arithmetic." |
| 995 | Semicolon | "Risk r=.244887983. Accurate comparison succeeds; no minimum correlation was specified..." | "Risk r = .244887983. No minimum correlation was specified..." |
| 995-996 | Unraised objection | "Accurate comparison succeeds; no minimum correlation was specified and no stability/adoption criterion is invented." | "No minimum correlation or stability criterion was specified in advance, so none is applied." |
| 998 | Not X but Y | "This is provisional saved-fit agreement, not a clean converged count curve." | "This agreement is provisional because both fits are unconverged." |
| 1007 | Semicolon | "The benchmark remains unchanged; no model adopted." | "The benchmark remains unchanged, and no model is adopted." |

## `chapters/2026-10-02.qmd`

The cleanest chapter: short sentences, concrete numbers, no hedging. One em-dash pair, one slogan-shaped heading and a few semicolons.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 37 | Dash | "the recommended fit — **intercept + signed loadings + empirical-Bayes priors + evidence-based pruning** — returns about 3 programs" | "In simulation (486 fits with a known truth of 3 programs), the recommended fit returns about 3 programs from any starting rank. It combines **an intercept, signed loadings, empirical-Bayes priors and evidence-based pruning**. The September version kept them all:" |
| 62 | Slogan heading | "## Methylation adds structure, not prediction" | "## What methylation adds" |
| 108 | Semicolon | "Their risk scores had stabilized; the overall fit had not." | "Their risk scores had stabilized, but the overall fit had not." |
| 196–198 | Semicolon | "rescale programs as in the multimodal model, so the plain sum is stable;" / "warm-start from the α_F = 0 fit;" / "phase survival in gradually;" | Drop the trailing semicolons from the bullet items. |

## `chapters/2026-09-18.qmd`

Formal and precise. The section headings use em dashes, and the point that rank must come from replicated held-out performance rather than active-factor counts is repeated about five times. The `expression--methylation` compounds (lines 5, 167) render as joining en dashes and are left alone.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 7 | Other | "...screening panel confirms that active-factor counts cannot select rank." | Line 100 says the panel is "this one draw", so "confirms" is stronger than the evidence. Consider "shows that, in one draw, active-factor counts do not select rank." (Flagged because it may change the claim. Author's call.) |
| 12 | Dash | "## §1 — Matched multimodal YFB implementation" | "## §1. Matched multimodal YFB implementation" |
| 41 | Dash | "## §2 — Resolving raw-scale convergence" | "## §2. Resolving raw-scale convergence" |
| 78 | Dash | "## §3 — Preliminary K-init diagnostic" | "## §3. Preliminary K-init diagnostic" |
| 136 | Dash | "## §4 — Follow-up questions from the September 4 meeting" | "## §4. Follow-up questions from the September 4 meeting" |
| 29–31 | Not X but Y | "This is the intended supervised multimodal extension, not an unsupervised factorization followed by Cox regression." | "This makes the model a supervised multimodal extension." |
| 82 | Semicolon | "(seed 20260919; 120 training and 120 validation patients; 50 features per modality)" | "(seed 20260919, 120 training and 120 validation patients, 50 features per modality)" |
| 103–106 | Long sentence | "The K-selection runner therefore uses identical simulated training and validation cohorts for every K within a replicate, checkpoints every fit, and will select the smallest K within one standard error..." | "The K-selection runner therefore uses identical simulated training and validation cohorts for every K within a replicate and checkpoints every fit. It will select the smallest K within one standard error of the best mean held-out C-index across independent replicates." |
| 117 | Semicolon | "...and eight union-retained factors; it recovered both true prognostic programs..." | "...and eight union-retained factors. It recovered both true prognostic programs..." |
| 118–120 | Repeated caveat | "The extra nominal survival factors are another reason the final rank must be selected from repeated held-out performance, not coefficient flags alone." | Already stated at lines 7 and 104–109. Delete, or shorten to "This is another case where coefficient flags overstate the rank." |
| 134 | Semicolon | "...its own recovery study; it is not part of the current method." | "...its own recovery study. It is not part of the current method." |
| 142 | Semicolon | "...when PDAC survival signal is weak; point-normal shrinkage can collapse coefficients to zero in small folds." | "...when PDAC survival signal is weak, because point-normal shrinkage can collapse coefficients to zero in small folds." |
| 149–153 | Repeated caveat | "Its initial result shows that convergence alone is insufficient: marginal PVE reporting is not a rank-selection method ... The next evidence must come from paired replicated held-out C-index..." | Fourth statement of the same point (also lines 7, 104–109, 179–180). Point back to §3: "§3 gives the initial result. Rank will be selected by paired replicated held-out C-index and partial likelihood, with the fixed-working objective kept as an optimization diagnostic." |

## `chapters/2026-09-04.qmd`

This is the densest of the three files. It has about 45 em dashes over 33 lines, several 60- to 90-word bullets, and a lot of defensive "not X, not Y" framing. It also uses moralized wording ("honest read", "correct and trustworthy"), internal session terminology, and one wrong cross-reference.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 5 | Dash | "...how each cohort contributes to the learned factors** — confirmed on real data and with simulated ground truth." | "...how each cohort contributes to the learned factors**. Both results hold on real data and with simulated ground truth." |
| 6 | Dash | "...choice of retention threshold** — the same 2 survival-active factors are kept..." | "...choice of retention threshold**: the same 2 survival-active factors are kept..." |
| 7 | Semicolon | "Now have genuine held-out versions alongside them; found and fixed a real orientation bug in the process." | "Now have genuine held-out versions alongside them. Found and fixed a real orientation bug in the process." |
| 8 | Semicolon | "...merges true factors together**, as expected; over-specification remains safe" | "...merges true factors together**, as expected. Over-specification remains safe" |
| 9 | Dash | "...out of scope for now** — no such data exists yet, and the current preprocessing can't represent it structurally." | "...out of scope for now**. No such data exist yet, and the current preprocessing cannot represent it structurally." |
| 15 | Dash | "## §1 — Clarifying the 8/27 results" | "## §1. Clarifying the 8/27 results" (same fix for l. 66, 87, 164, 184) |
| 30 | Dash | "...cross-validated?** No — `compute_joint_ll_bic()` fits on all 273 training patients..." | "...cross-validated?** No. `compute_joint_ll_bic()` fits on all 273 training patients..." |
| 30 | Unraised objection | "That's not wrong, but it isn't what was implied at the meeting." | "That is what was computed, although the meeting discussion implied a cross-validated version." |
| 32 | Dash | "(holding out patient AND gene blocks jointly — an ordinary held-out check doesn't work here, since...)" and "...a \"cross-validated BIC\" — the two correct for..." | "(holding out patient and gene blocks jointly). An ordinary held-out check does not work here, because..." and "...a \"cross-validated BIC\", because the two correct for..." |
| 33 | Dash | "...to a corrected 1,614,079 — small, doesn't change which K_init BIC prefers" / "...not just worked around downstream — the practical consequence is..." / "...flip sign as a result — magnitudes and every qualitative conclusion..." | "...to a corrected 1,614,079, a small change that leaves BIC's preferred K_init unchanged" / "...not just worked around downstream. As a result, every joint survival coefficient..." / "...flip sign as a result. Magnitudes and every qualitative conclusion..." |
| 33 | Semicolon | "...drawn from them are unchanged; see §3." | "...drawn from them are unchanged (see §3)." |
| 33 | Moralized precision | "...is now on a *correct and trustworthy* sign convention for the first time, rather than one requiring a symmetric workaround to interpret safely." | "...now uses the correct sign convention, so no symmetric workaround is needed to interpret it." |
| 33 | Long sentence | "It has since been fixed directly in the model's fitting code ... rather than one requiring a symmetric workaround to interpret safely." (about 60 words, in a bullet of about 200) | Split as in the Dash and Moralized precision rows, and move the bug explanation into a short paragraph of its own. |
| 35 | Unraised objection | "**Six criteria on K_init now, not four — and they still disagree, which is informative, not a problem:**" | "**Six criteria on K_init now, not four, and they still disagree:**" (also removes the dash) |
| 47 | Dash | "**Generalization to genuinely new cohorts — the actual deployment question**" (table cell) | "**Generalization to new cohorts, the deployment question**" |
| 49 | Dash | "...as also preferring K_init=3 — it doesn't." | "...as also preferring K_init=3, which was incorrect." |
| 49 | Voice habit | "...so the three needn't agree, and here they don't." | "...so the three need not agree, and here they do not." |
| 51 | Semicolon | "Each answers a different question; the one closest to the deployment question (external C-index) is the one to trust." | "Each answers a different question. The one closest to the deployment question (external C-index) is the one to trust." |
| 51 | Dash | "...the smallest K_init on that plateau — a parsimony argument, not a claim that the other criteria agree." | "...the smallest K_init on that plateau. This is a parsimony argument, not a claim that the other criteria agree." |
| 62 | Dash | "They didn't — the original fit (restart 1) was already the best..." | "They did not. The original fit (restart 1) was already the best..." |
| 66 | Dash | "## §2 — Are retained factors sensitive to threshold variation?" | "## §2. Sensitivity of retained factors to the thresholds" |
| 68 | Dash | "Neither has ever been checked for sensitivity — until now." | "Neither had been checked for sensitivity before this analysis." |
| 68 | Mannered prose | "— until now." | The dramatic reveal adds nothing. Use the plain statement in the row above. |
| 75 | Semicolon | "...no threshold at all keeps 4 factors; the old default of 0.05 keeps 0..." | "...no threshold at all keeps 4 factors, and the old default of 0.05 keeps 0..." |
| 81 | Dash | "...around the K_init=11/13 dip (§1) — `threshold_vs_kinit.csv` shows..." / "...a stable retained-factor count — the threshold analysis doesn't explain the dip..." | "...around the K_init=11/13 dip (§1). `threshold_vs_kinit.csv` shows..." / "...a stable retained-factor count. The threshold analysis does not explain the dip..." |
| 81 | Not X but Y | "the threshold analysis doesn't explain the dip, it just shows the dip's factor count is itself threshold-sensitive in a way the rest of the K_init range isn't." | "The threshold analysis does not explain the dip. It shows only that the factor count at the dip is threshold-sensitive in a way it is not elsewhere in the K_init range." |
| 83 | Dash | "One remaining fairness question — do factors need their coefficients rescaled ... fixed threshold — checked directly: rescaling doesn't change..." | "We also checked whether factors need their coefficients rescaled by how much their signal varies before comparison to a single fixed threshold. Rescaling does not change which factors are called prognostic here." |
| 87 | Dash | "## §3 — Cohort membership in the survival term: performance and interpretability" | "## §3. Cohort membership in the survival term: performance and interpretability" |
| 89 | Announcing importance | "**The headline result: letting the survival coefficient vary by cohort...**" | "**Letting the survival coefficient vary by cohort...**" |
| 89 | Dash | "...contributes to the learned factors — not a predictive-accuracy gain.**" | "...contributes to the learned factors, rather than a gain in predictive accuracy.**" |
| 103 | Semicolon | "$\eta_i$ is patient $i$'s risk score; $(YF)_{ik}$ is how strongly...; $\beta_k^{c(i)}$ is program $k$'s survival coefficient..." | "Here $\eta_i$ is patient $i$'s risk score, $(YF)_{ik}$ is how strongly..., and $\beta_k^{c(i)}$ is program $k$'s survival coefficient..." |
| 103 | Dash | "...rather than estimated fully independently — what's shared across cohorts is the *amount* of shrinkage..." | "...rather than estimated fully independently. What is shared across cohorts is the *amount* of shrinkage..." |
| 114 | Dash | "...reproduces the previously-reported result exactly — mean external C=0.6267 — confirming..." | "...reproduces the previously reported result exactly (mean external C=0.6267), confirming..." |
| 127 | Dash | "(pooled 95% CI includes zero, and is narrow) — the point estimate very slightly favors..." / "(factorize first, then fit survival separately) — a bootstrap CI is only reported here..." | "(pooled 95% CI includes zero, and is narrow), and the point estimate very slightly favors..." / "(factorize first, then fit survival separately). A bootstrap CI is only reported here..." |
| 127 | Semicolon | "...not proof that the two are equivalent; the honest read is..." | "...not proof that the two are equivalent. The supported reading is..." |
| 127 | Moralized precision | "the honest read is \"no evidence of a performance cost...\"" | "the supported reading is \"no evidence of a performance cost...\"" |
| 127 | Stacked hedge | "...was not detected as statistically significant ... A CI that includes zero means a difference wasn't detected here, not proof that the two are equivalent; the honest read is ... not \"underperforms.\"" | The non-equivalence caveat is stated three ways in one bullet. Keep one: "The pooled 95% CI includes zero and is narrow, so there is no evidence of a performance cost to letting the coefficient vary by cohort." |
| 128 | Dash | "...worse than the two-step alternative too — it visibly competes with..." / "...indistinguishable from the current model too — a clean null, not a failure to report." | "...worse than the two-step alternative too. It visibly competes with..." / "...indistinguishable from the current model too." |
| 128 | Unraised objection | "a clean null, not a failure to report" | Cut, as in the Dash row above. |
| 129 | Dash | "...left out of the average — not an artifact of one small cohort." / "...on a common (pooled) risk set** — `joint_yfb_beta_c` > ..." | "...left out of the average, so it does not depend on one small cohort." / "...on a common (pooled) risk set**: `joint_yfb_beta_c` > ..." |
| 129 | Long sentence | "It does **not** extend to the two arms that use `strata_id` ... and are not comparable to the other three or to each other's C-index ranking on this metric." (about 90 words) | "It does **not** extend to the two arms that use `strata_id` (`joint_yfb_all_c`, and a `strata_only` sub-arm not shown in the table above). Stratifying the risk set by cohort changes the scale of the partial log-likelihood, not just its value. Those two arms' held-out log-likelihoods (around -354.8 vs. -445 to -451 for the other three) reflect that different risk-set structure rather than a better fit, so they are not comparable on this metric." |
| 137 | Dash | "(magnitude 0.030 vs. 0.0004 in TCGA_PAAD) — the current (shared-coefficient) model blends..." / "...from external validation performance — it only shows up once..." / "...was corrected this session — see the fix-plan work above these results — so these values..." | "(magnitude 0.030 vs. 0.0004 in TCGA_PAAD). The current shared-coefficient model blends..." / "...from external validation performance. It only shows up once..." / "Because of the sign-convention fix in §1, these values are stated by magnitude." |
| 137 | Semicolon | "...so these values are stated by magnitude; which cohort dominates which program..." | "...so these values are stated by magnitude. Which cohort dominates which program..." |
| 137 | Repeated caveat | "(The joint coefficients' sign convention was corrected this session ... are unaffected by that fix.)" | Repeats l. 33. Shorten to "Because of the sign-convention fix in §1, these values are stated by magnitude." |
| 137 | Announcing importance | "This is exactly the kind of \"does the model actually handle multiple cohorts sensibly\" evidence that would matter for comparison against methods like DeSurv or plain unsupervised factorization, none of which represent per-cohort effects at all." | "This per-cohort evidence is relevant when comparing against DeSurv or plain unsupervised factorization, neither of which represents per-cohort effects." |
| 137 | Other | "corrected this session — see the fix-plan work above these results" (also "this session" at l. 33 and l. 97) | Internal session terminology (AGENTS.md documentation-audience rule). Use "§1" or a date, for example "corrected on 2026-09-04 (§1)". |
| 154 | Dash | "The value again wasn't a bigger number — the shared coefficient already sat close to the true cohort's value — it was that only..." | "As on the real data, the gain is not a larger coefficient. The shared coefficient was already close to the true cohort's value, and only the cohort-specific version makes the attribution visible." |
| 154 | Repeated caveat | "The value again wasn't a bigger number..." | Third statement of the "value is visibility, not accuracy" point (l. 5, l. 89). The rewrite above keeps it to one sentence. |
| 156 | Counted package | "### Three supporting technical notes" | "### Supporting technical checks" |
| 158 | Dash | "Asserted in the 8/21 notes, never tested — tested now:" | "Pooling was asserted in the 8/21 notes but not tested until now." |
| 158 | Voice habit | "...changes both the sample size AND the gene-selection step..." | "...changes both the sample size and the gene-selection step..." |
| 160 | Dash | "...any feedback from survival at all — its update formula is identical..." / "...or which iterate is selected) — a real, currently-disabled option flagged for future work" | "...any feedback from survival at all. Its update formula is identical..." / "...or which iterate is selected). That option exists in the code but is currently disabled, and it is flagged for future work (`ROADMAP.md`)." |
| 160 | Voice habit | "gene-program-DISCOVERY step" / "driven by the FULL objective" | "gene-program discovery step" / "driven by the full objective" (use italics if emphasis is needed) |
| 164 | Dash | "## §4 — Other findings" | "## §4. Other findings" |
| 166 | Dash | "(This diagnostic is from a single simulated replicate — exploratory, not aggregated..." | "(This diagnostic comes from a single simulated replicate and is exploratory. It is not aggregated..." |
| 166 | Semicolon | "...the rest of this section's numbers are averaged over; treat the specific counts as illustrative..." | "...used for the other numbers in this section, so the counts illustrate the pattern rather than estimate its frequency.)" |
| 166 | Stacked hedge | "exploratory ... treat the specific counts as illustrative of the pattern rather than a precise estimate of its frequency" | One hedge is enough, as in the Semicolon row above. |
| 172 | Voice habit | "the advantage is strongly POSITIVE below K=6" | "the advantage is strongly positive below K=6" |
| 174 | Dash | "(no performance metric — a comparison of which genes each method's programs are built from...)" | "(no performance metric: a comparison of which genes each method's programs are built from...)" |
| 176 | Other | "which is the direct explanation for the second technical note in §3" (also l. 180) | The variance question is the third note in §3 (l. 160). The second note (l. 159) is about `cohort_id`. Fix the cross-reference in both the caption and the body. |
| 180 | Dash | "(hypergeometric overlap p-values as low as 1e-129) — not buried deep in..." / "...not independent samples — read them as..." | "(hypergeometric overlap p-values as low as 1e-129), not buried deep in..." / "...not independent samples. Read them as..." |
| 180 | Not X but Y | "These p-values are descriptive, not inferential: ... correlated by construction, not independent samples — read them as \"this overlap is large,\" not as a formal hypothesis test." | "These p-values are descriptive. Both sets of gene programs come from the same expression matrix, so their gene-selection events are correlated by construction. Read them as showing a large overlap rather than as a formal hypothesis test." |
| 184 | Slogan heading | "## §5 — Where the value-add goes next" | "## §5. Multi-modality extension" (also removes the dash) |
| 186 | Unraised objection | "A prerequisite for that direction, not a decision against it." | "Representing multiple feature types per patient is a prerequisite for that direction." |
| 9, 30, 32, 33, 81, 83, 103, 127, 154, 158, 186 | Voice habit | Contractions in academic prose: "can't", "That's", "isn't", "doesn't", "didn't", "what's", "wasn't" | Expand them ("cannot", "That is", "is not", "does not", "did not", "what is", "was not"). |

## `chapters/2026-08-27.qmd`

A short, readable meeting record, but it uses em dashes throughout (15 places, including every section heading) and semicolons. It also states the "disagreement is shown, not hidden" point several times and ends one section with a moralized "honest summary".

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 12 | Semicolon | "ELBO and BIC both prefer K_init=3; external C-index prefers K_init=12..." | "ELBO and BIC both prefer K_init=3. External C-index prefers K_init=12..." |
| 13 | Dash | "...the presented recommendation — survival-active count is stable at 2 from K_init=7 up — but the full curves..." | "...the presented recommendation because the survival-active count is stable at 2 from K_init=7 up, but the full curves..." |
| 19 | Dash | "...flat across K_init=6/12/20 in every scenario — a clean, well-powered confirmation (15 seeds) that..." | "...flat across K_init=6/12/20 in every scenario (15 seeds). Over-specifying K and relying on ARD performs as well as knowing the true K." |
| 20 | Economic metaphor | "...over-specifying K and trusting ARD costs nothing relative to knowing the true K." | "...over-specifying K and relying on ARD performs as well as knowing the true K." |
| 23 | Dash | "**A secondary question — does over-specifying K_init increase the false-positive rate on non-signal factors — softened under more seeds...**" | "**Whether over-specifying K_init increases the false-positive rate on non-signal factors is less clear at 15 seeds than at 5.**" |
| 24 | Lesson drawn | "...which is itself the finding worth recording." | Cut. The next sentence gives the finding. |
| 25 | Semicolon | "At 5 seeds it looked like a clean improvement; at 15 seeds the improvement..." | "At 5 seeds it looked like a clean improvement. At 15 seeds the improvement..." |
| 34 | Dash | "## §1 — K_init consensus sweep" | "## §1. K_init consensus sweep" (same fix for l. 72 and l. 113) |
| 41 | Dash | "...from the fit's *returned* `EF`/`EBeta` — matching exactly what external prediction does — rather than read from a stale in-loop value..." | "...from the fit's *returned* `EF`/`EBeta`, which matches external prediction exactly, rather than from a stale in-loop value..." |
| 54 | Semicolon | "ELBO and BIC both prefer the smallest K_init tested (3); BIC's complexity penalty..." | "ELBO and BIC both prefer the smallest K_init tested (3). BIC's complexity penalty..." |
| 56 | Dash | "...the plateau from K_init=7 onward is wide — most points in that range sit within a few thousandths..." | "...the plateau from K_init=7 onward is wide. Most points in that range sit within a few thousandths..." |
| 62 | Dash | "...for the same reason it has been kept before — the survival-active count (2) is stable there and beyond..." | "...for the same reason it has been kept before: the survival-active count (2) is stable there and beyond..." |
| 63 | Repeated caveat | "The full sweep is presented alongside it specifically so the disagreement between criteria is visible to the audience, not smoothed over into a single number." | Same point as l. 11 and l. 14-15, and again at l. 104-105. Shorten to "The full sweep is shown alongside it." |
| 68 | Dash | "...and a sensitivity sweep — `beta_threshold` in particular was reverse-engineered ... not drawn from a citable source — and whether a fully Bayesian..." | "...and a sensitivity sweep, and whether a fully Bayesian K-selection alternative ... is worth scoping as a longer-term replacement for the current post-hoc consensus. `beta_threshold` in particular was reverse-engineered from this model's own coefficient scale, not drawn from a citable source." |
| 72 | Dash | "## §2 — Multi-cohort simulation under ARD" | "## §2. Multi-cohort simulation under ARD" |
| 81 | Dash | "(e.g. hybrid-scenario held-out C-index 0.815 vs. 0.815) — a useful confirmation that the rewrite changed nothing about the underlying fitting..." | "(e.g. hybrid-scenario held-out C-index 0.815 vs. 0.815). The rewrite changed only which K values get compared, not the underlying fitting." |
| 85 | Economic metaphor | "Over-specifying K and trusting ARD to prune costs nothing relative to knowing the true K in advance" | "Over-specifying K and relying on ARD to prune performs as well as knowing the true K in advance." (Also repeats l. 20.) |
| 86 | Dash | "...knowing the true K in advance — the core validation this exercise was designed to produce." | "...knowing the true K in advance. This was the main check the simulation was designed for." |
| 86 | Announcing importance | "the core validation this exercise was designed to produce" | Cut, or "This was the main check the simulation was designed for." |
| 88 | Lesson drawn | "### §2.2 — A false-positive-rate trend, and a lesson about seed count" | "### §2.2. False-positive rate versus K_init at 5 and 15 seeds" (also removes the dash) |
| 98 | Dash | "...gives p=0.07 — not conventionally significant." | "...gives p=0.07, which is not significant at the 0.05 level." |
| 102 | Moralized precision | "The honest summary is two-layered: the primary metrics..." | "In summary, the primary metrics..." |
| 103 | Semicolon | "...a well-powered null result with no caveats needed; the false-positive-rate question gives..." | "...a well-powered null result with no caveats needed. The false-positive-rate question gives..." |
| 107 | Long sentence | "This also lines up with a separate discussion with advisors about the 8/21 finding ... independent of whether the false-positive *rate* itself trends up or down with K_init." (about 70 words) | "This agrees with a separate discussion with advisors about the 8/21 finding that ARD can over-count survival-active factors in simulation. A small, non-zero coefficient on a non-signal factor adds little to the risk score in practice. Over-counting in this sense is therefore a smaller concern than it might first appear, whether or not the false-positive *rate* trends up or down with K_init." |
| 113 | Dash | "## §3 — The 8/27 lab meeting deck" | "## §3. The 8/27 lab meeting deck" |
| 116 | Dash | "...from the 6/18 deck as a template, not a diff against it — the loadings predictor (Lβ) is dropped..." | "...from the 6/18 deck as a template rather than a diff against it. The loadings predictor (Lβ) is dropped..." |
| 122 | Semicolon | "...the baseline the June deck used; the headline pooled bootstrap difference..." | "...the baseline the June deck used. The headline pooled bootstrap difference..." |

## `chapters/2026-08-21.qmd`

The model specification and K-selection sections are clear and well glossed. The main problems are about 25 prose semicolons, the process caveat on the factor-biology work stated three times, and a contradiction about whether the genomics-only factors have been characterized.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 12 | Semicolon | "...every $K \ge 7$**; added capacity goes entirely to dead factors..." | "...every $K \ge 7$**. Added capacity goes entirely to dead factors..." |
| 14 | Semicolon | "(C=0.596 vs. 0.627); $K=6$ loses on both counts" | "(C=0.596 vs. 0.627). $K=6$ loses on both counts" |
| 42 | Semicolon | "Finalize committee (Lin, ...); prelim targeted for end of Fall 2026; 8/27 lab meeting preparation (...)" | "Finalize committee (Lin, ...). Prelim targeted for end of Fall 2026. Prepare for the 8/27 lab meeting (...)." |
| 56 | Semicolon | "(Program 7, adverse; Program 3, protective)" | "(adverse Program 7 and protective Program 3)" |
| 73 | Semicolon | "(2026-08-19/20, this cycle; §3 below)" | "(2026-08-19/20, this cycle, see §3)" |
| 79 | Repeated caveat | "**Factor stability across random seeds** (8/3 open item). Refitting each $K$ ... from 15 different starting points ..." | Same result appears in the Summary (l. 19) and almost verbatim in §3.1 (l. 457). Shorten to "**Factor stability across random seeds** (8/3 open item): closed, see §3.1." |
| 81 | Semicolon | "...at every $K$; random restarts never beat the SVD initialization." | "...at every $K$. Random restarts never beat the SVD initialization." |
| 88 | Semicolon | "...cross-cohort hazard ratios); the other factors in the $K=7$ fit..." | "...cross-cohort hazard ratios). The other factors in the $K=7$ fit..." |
| 89 | Repeated caveat | "This was pursued independently rather than at the group's direction, so the approach needs review even though the existing results may turn out to be sufficient." | Already said at l. 38-40 and again at l. 396-401. Keep one full statement (l. 398-401) and shorten here to "**Aligning on the factor-biology work** (§2.5)." |
| 99 | Voice habit | "Frequent pushes to the model repository & writing in separate manuscript repo." | "Frequent pushes to the model repository, and writing in a separate manuscript repository." |
| 101 | Semicolon | "...absorbs genomic platform offsets only; it never carries survival information, and is off by default" | "...absorbs genomic platform offsets only. It carries no survival information and is off by default" |
| 136 | Voice habit | "...get driven toward zero by scale alone (not because there's no real signal)." | "...are pushed toward zero by scale alone, even when real signal is present." |
| 146 | Semicolon | "(standard Cox); it cancels out of the partial likelihood" | "(standard Cox). It cancels out of the partial likelihood" |
| 150 | Showy verb | "...the LB parameterization ($\eta = L\beta$), where $L$ carries both roles." | "...the LB parameterization ($\eta = L\beta$), where $L$ enters both likelihoods." |
| 171 | Semicolon | "...under cross-validation (5/28 entry); `normal` is the prior actually used..." | "...under cross-validation (5/28 entry). `normal` is the prior used..." |
| 173 | Semicolon | "$\tau_j$ has no prior in this model; it is fit by constrained maximum likelihood" | "$\tau_j$ has no prior in this model. It is fit by constrained maximum likelihood" |
| 189 | Repeated caveat | "The survival term conditions on $F$ and the *observed* $Y$ ... not on $L$. This is exactly what removes $L$ from the survival side of the model under YFB." | Already stated at l. 148-150 and again at l. 268. Cut the second sentence: "The survival term conditions on $F$ and the *observed* $Y$ (through $Z\!F$), not on $L$." |
| 198 | Semicolon | "...estimated once per update from all the data; the posterior is the individual belief..." | "...estimated once per update from all the data. The posterior is the individual belief..." |
| 201 | Voice habit | "...since they differ whenever there's posterior uncertainty." | "...since they differ whenever there is posterior uncertainty." |
| 218 | Semicolon | "...between a posterior and its prior; it penalizes posteriors that drift too far" | "...between a posterior and its prior. It penalizes posteriors that drift too far" |
| 220 | Voice habit | "...not conjugate to a Normal prior, so it's handled by a second-order Taylor expansion" | "...not conjugate to a Normal prior, so it is handled by a second-order Taylor expansion" |
| 235 | Semicolon | "(how much information the data carries about this one coordinate; larger $A$ means a more confident estimate)" | "(how much information the data have about this one coordinate, where larger $A$ means a more confident estimate)" |
| 239 | Voice habit | "...is what's actually handed to the EBNM solver" | "...is what is passed to the EBNM solver" |
| 240 | Semicolon | "...identical for $\beta$, $L$, and $F$ below; only what goes into $A$ and $B$ changes." | "...identical for $\beta$, $L$, and $F$ below. Only what goes into $A$ and $B$ changes." |
| 244 | Semicolon | "...and $\eta = Z\!F\,\tilde\beta$; run one step of the Cox Newton–Raphson..." | "...and $\eta = Z\!F\,\tilde\beta$. Then run one step of the Cox Newton–Raphson..." |
| 264 | Semicolon | "(see the priors note above; `point_normal` is available but collapses to $\beta=0$ under cross-validation)" | "(see the priors note above)" |
| 264 | Repeated caveat | "`point_normal` is available but collapses to $\beta=0$ under cross-validation" | Already stated at l. 170-172. Cut, as in the semicolon fix. |
| 284 | Semicolon | "...unused in the recommended configuration; see `derivations/...` for the general case" | "...unused in the recommended configuration. See `derivations/...` for the general case" |
| 295 | Unraised objection | "**Dual-source $F$ was tested, not left unused.**" | "**Dual-source $F$ was tested.**" |
| 298 | Voice habit | "[Full method & results: DECISIONS.md, 2026-08-20 entries]" | "[Full method and results: DECISIONS.md, 2026-08-20 entries]" |
| 306 | Semicolon | "...is the posterior-variance correction; it prevents $\tau$ from being systematically overestimated" | "...is the posterior-variance correction. It prevents $\tau$ from being systematically overestimated" |
| 319 | Semicolon | "...separately from $\alpha$; retired 2026-07-12 as redundant with $\alpha$" | "...separately from $\alpha$. It was retired on 2026-07-12 as redundant with $\alpha$" |
| 342 | Long sentence | "First proposed on 5/28 via cross-validated external C-index, but that justification is superseded this cycle (§3.1): ARD pruning within a large-$K$ fit, combined with ... not as the primary selection criterion." | "First proposed on 5/28 via cross-validated external C-index. That justification is superseded this cycle (§3.1). $K$ is now chosen by ARD pruning within a large-$K$ fit, combined with model fit (ELBO) across candidate $K$'s. External C-index serves only as a floor that rules out a $K$ that fits training data best but generalizes worse." |
| 347 | Semicolon | "...takes this current, ARD/ELBO-validated $K=7$ fit as given; it does not revisit the choice of $K$ itself." | "...takes this current, ARD/ELBO-validated $K=7$ fit as given and does not revisit the choice of $K$ itself." |
| 353 | Repeated caveat | "§3.1 below is where decision 1 is actually re-derived; this table only carries out decision 2 on the result." | Restates l. 346-347 and l. 348-351. Delete the sentence (this also removes a semicolon). |
| 385 | Not X but Y | "This is a pattern, not a diagnosis." | "Reading the list suggests a pattern, but it is not a formal characterization." |
| 386 | Semicolon | "Program 3's list doesn't jump out the way Program 7's does; it needs the formal methods rather than a read-off." | "Program 3's list shows no obvious pattern the way Program 7's does, so it needs the formal methods rather than a read-off." |
| 388 | Other | "The 2 genomics-only programs: no gene list or biological characterization attempted yet" | Contradicts l. 403-410, which reports both genomics-only factors characterized this cycle (Program 5 stroma, Program 6 immune). l. 532-534 repeats the stale claim. Update both to point to the 2026-08-21 enrichment report. |
| 391 | Tutorial heading | "**Factor biology: work already done, and where it stands.**" | "**Factor biology from the July 15 cycle.**" |
| 398 | Not X but Y | "The caveat is process, not evidence." | "The open concern is about process." |
| 398 | Repeated caveat | "This line of investigation was pursued on our own initiative rather than at the group's direction..." | Third statement of this caveat (l. 38-40, l. 89-91). Keep this one as the full version and shorten the others. |
| 405 | Announcing importance | "Reporting it turns out to matter: Program 5 reads as desmoplastic stroma..." | "Program 5 reads as desmoplastic stroma..." |
| 430 | Personality verb | "...give the model more programs than it needs and it drives the useless ones' survival coefficients to zero on its own..." | "...when the model is fit with more programs than it needs, the prior shrinks the survival coefficients of the unneeded ones to zero..." |
| 440 | Semicolon | "**dead**: neither; fully pruned" | "**dead**: neither (fully pruned)" |
| 462 | Long sentence | "$K=6$'s ELBO (−812,248) is actually *worse* than $K=7$'s, not better, so it loses on both fit and generalization: both $K=5$ and $K=6$ generalize far worse ... but only $K=5$ can claim a real training-fit advantage that generalization then fails to bear out." | "$K=6$'s ELBO (−812,248) is *worse* than $K=7$'s, so $K=6$ loses on both fit and generalization. Both $K=5$ and $K=6$ generalize far worse than $K=7$ (C=0.596–0.597 vs. 0.626–0.628 for every $K \ge 7$). Only $K=5$ has a real training-fit advantage, and it does not carry over to the external cohorts." |
| 473 | Voice habit | "...whatever the model is fitting internally isn't generalizable structure, so it shouldn't win on training-fit score alone." | "...whatever the model is fitting internally is not generalizable structure, so it should not win on training-fit score alone." |
| 475 | Semicolon | "...is kept as the simpler model; $K=9$ (2 + 3) is recorded as..." | "...is kept as the simpler model. $K=9$ (2 + 3) is recorded as..." |
| 477 | Semicolon | "\"Start as large as possible\" is not supported here; $K$ in the 7–10 range performs best." | "\"Start as large as possible\" is not supported here. $K$ in the 7–10 range performs best." |
| 518 | Filler | "...the same pattern seen in the July two-step comparison (§1), expected for these cohort sizes and not unusual." | "...the same pattern seen in the July two-step comparison (§1), as expected for these cohort sizes." |
| 534 | Semicolon | "...used on Programs 3/7 (§2.5); this is an open follow-up." | "...used on Programs 3/7 (§2.5). This is an open follow-up." (See also the l. 388 contradiction.) |
| 554 | Voice habit | "A large $K$ value that the two-step method wasn't \"told\" ($K$=20/$K$=40)..." | "A large $K$ chosen without reference to the joint model's answer ($K$=20/$K$=40)..." |
| 576 | Semicolon | "(pure noise); strength $1$ is this project's standard..." | "(pure noise). Strength $1$ is this project's standard..." |
| 580 | Unraised objection | "This is a property of the simulated data, not anything either model estimates:" | "These coefficients are fixed by the simulation design:" |
| 593 | Unraised objection | "This is a compute limit, not a reversal." | "The wide intervals reflect the small number of replicates." |

## `chapters/2026-08-03.qmd`

Short and focused, but the "not a performance claim" caveat about DeSurv appears twice here after already appearing in 2026-07-15.qmd. The "Relevance" section mostly restates Results. The figure caption at line 30 has a semicolon, but it sits inside a chunk option and is not tabled.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 20 | Semicolon | "Nothing changed in the model or pipeline this cycle; this was a validation and comparison cycle building directly on 7/15." | Also repeats line 14. "No model or pipeline changes. This cycle was validation and comparison, building on 7/15." |
| 36 | Semicolon | "...independent of the survival term; the survival term identifies it rather than creates it." | "...independent of the survival term. The survival term identifies it rather than creating it." |
| 37–39 | Not X but Y | "SBMF is a supervised re-weighting of EBMF's own factor structure: a variation on EBMF, not a separate model class that happens to also fit survival." | "SBMF is therefore a supervised re-weighting of EBMF's own factor structure, a variation on EBMF." |
| 44 | Semicolon | "Different model classes converge on the same tumor biology; that convergence is the point here..." | "Different model classes arrive at the same tumor biology." |
| 44–46 | Repeated caveat | "...not a head-to-head performance claim. A matched, same-protocol comparison against DeSurv is still the pending next step for an actual performance claim." | State the caveat once: "A performance comparison still needs a matched, same-protocol run against DeSurv." |
| 50–53 | Recap box | "Two independent checks now agree that the 2-program structure is real biology ... That's stronger evidence than either check alone, since neither depends on the joint model's own objective to succeed." | Keep only the new point: "Neither check depends on the joint model's objective, so together they are stronger evidence than either alone that the 2-program structure is real biology." |
| 58–59 | Repeated caveat | "Builds directly on today's gene-overlap comparison, but as a performance claim, not just a biology-convergence claim." | Third time the caveat appears. "Matched, same-protocol head-to-head against DeSurv (same data, same scoring), extending today's gene-overlap comparison to performance." |

## `chapters/2026-07-15.qmd`

The most conversational chapter. It has six semicolons, several "actually"s and idioms ("under the hood", "knob", "apples-to-apples"), and it defends results against objections nobody raised. One factual inconsistency: "five independent methods" lists four.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 17 | Semicolon | "...on model choice this cycle; the work below is corrections and stress-tests..." | "...on model choice this cycle. The work below corrects and stress-tests that same model." |
| 21 | Mannered prose | "Two things changed under the hood, neither of which changes which model we recommend:" | Also repeats line 16. "Two implementation changes:" |
| 29–30 | Unraised objection | "...from an earlier 0.636 to **`r rec_mean`**, a correction rather than a regression." | "...from the earlier 0.636 to **`r rec_mean`**. The earlier value came from mismatched preprocessing." |
| 41 | Semicolon | "because it's the one you get from a single fit with no extra steps; $K=4$ or $K=5$ get you statistically the same performance..." | "$K=7$ stays the default because a single fit reaches it with no extra steps. $K=4$ or $K=5$ give statistically the same performance but need a warm-start procedure." |
| 37, 44 | Filler | "How many of the 7 factors actually matter..." ... "Does the joint model actually need the survival coupling..." | Drop "actually" in both (also lines 27, 56 of 2026-05-28.qmd). |
| 44–45 | Mannered prose | "Note this is a different knob from the $\alpha$ weighting..." | "This is separate from the $\alpha$ weighting in the Method updates above." |
| 47–51 | Long sentence | "Here, instead, we control the *data-generating process*: in simulation, the true survival coefficients ... baseline:" (about 70 words) | "Here we vary the *data-generating process* instead. In simulation, the true survival coefficients used to generate each patient's survival time are multiplied by a scaling factor that runs from 0 (no true prognostic effect) up to a realistic magnitude. These are generating coefficients, not estimates from either model. At each setting we compare the joint model against a two-step (factorize, then Cox) baseline:" |
| 57–58 | Repeated caveat | "At a scaling factor of 0 (no true survival coefficients), joint and two-step are indistinguishable (both near chance, as they should be, since there's nothing to detect)." | The no-signal point is already in line 49 and the figure caption. "At a scaling factor of 0, joint and two-step are indistinguishable and both are near chance." |
| 60–61 | Unraised objection | "...reverses this, but that's a known factor-collapse failure mode in that specific design, not evidence against survival-coupling in general." | "...reverses this. That design has a known factor-collapse failure mode." |
| 84 | Semicolon | "...which lowers its mean; the joint model stays stable across all 5..." | "...which lowers its mean. The joint model stays stable across all 5..." |
| 102 | Semicolon | "We'd flagged this qualitatively before; this is the first time it's quantified." | "We had flagged this qualitatively before, and this is the first time it is quantified." |
| 108 | Other | "Confirmed by five independent methods (gene-set enrichment, PurIST subtype correlation, external hazard ratios, DeSurv gene-list overlap)" | Four methods are listed but the text says five (line 147 also says "all five"). Add the missing method from the pathway-enrichment report: "...DeSurv gene-list overlap, [ADD: fifth method])". |
| 115 | Semicolon | "(basal-like / squamous / MET-EGFR-driven); Program 3 is protective" | "(basal-like / squamous / MET-EGFR-driven). Program 3 is protective" |
| 116–118 | Not X but Y | "...which is reassuring rather than novel; it tells us the model recovered something real, not that we've found something new." | Also a semicolon. "This is the subtype-survival axis already expected for PDAC. The model recovered known biology, which is a check on the model rather than a new finding." |
| 126–128 | Unraised objection | "This is not a claim that we already beat DeSurv directly. A matched, same-protocol comparison ... to make that claim, one way or the other." | "We have not yet compared it directly with DeSurv. A matched, same-protocol comparison is the next step." |
| 137 | Semicolon | "collapses under cross-validation; we're keeping the normal prior, which is the right call for now..." | "collapses under cross-validation, so we are keeping the normal prior. It needs a sentence of justification if asked." |
| 141 | Mannered prose | "(to make sure a future head-to-head is apples-to-apples)" | "(so a future head-to-head uses the same simulation protocol)" |

## `chapters/2026-05-28.qmd`

Dense with useful numbers, but the bold lead sentences lean toward slogans, and the "Relevance" paragraph answers an objection nobody raised. The configuration counts also do not add up.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 9–11, 19, 59 | Other | "18-configuration preprocessing benchmark ... 5 preprocessing strategies ... across LB and YFB, with and without a cohort indicator" ... "10 of the 12 non-per-platform configurations" ... "18/5 configurations" | Counts do not reconcile (5 × 2 × 2 = 20, not 18), and "18/5" is unclear. Check the benchmark report and state the actual grid. No rewrite proposed until the numbers are confirmed. |
| 19 | Announcing importance | "**Critical finding: preprocessing is not a free choice.**" | "**Preprocessing determines whether the fit collapses.**" |
| 25 | Semicolon | "...separately processed RNA-seq datasets; see the 2026-10-16 correction entry..." | "...separately processed RNA-seq datasets. See the 2026-10-16 correction entry in `DECISIONS.md`." |
| 27 | Not X but Y | "**Cohort indicator helps only in specific configurations, not universally.**" | "**The cohort indicator helps only in specific configurations.**" |
| 28 | Filler | "(+0.020 for LB with joint QN, meaningfully helpful there)" | "(+0.020 for LB with joint QN)" |
| 36 | Semicolon | "by +0.014 (D4 vs. D3); the two effects are additive." | "by +0.014 (D4 vs. D3). The two effects are additive." |
| 39–40 | Other | "The same 2-factor structure appears in every preprocessing configuration tested, which supports it being a real property..." | Restates lines 36–39. Fold it in: "...one protective-direction program every time, which supports the 2-factor structure being a real property of the data rather than a fitting artifact." |
| 50–51 | Semicolon | "the number as understood this week; a train/test preprocessing bug ... revises this down to 0.627; see that entry." | "D4's 0.636 is the number as understood this week. A train/test preprocessing bug found later, on 7/15, revises this down to 0.627 (see that entry). It is kept here as the historical record of what we knew at the time." |
| 56–60 | Unraised objection | "This is the week the recommended configuration actually converges ... It beats the alternatives for reasons we understand: ... not because it happened to score highest in a table of 18/5 configurations." | "The recommended configuration is now YFB, DeSurv-aligned gene selection, per-platform z-standardization, no cohort indicator, $K=7$. Each choice has a reason: per-platform standardization avoids platform-scale collapse, the gene selection matches DeSurv, and YFB scores new subjects exactly instead of through LB's OLS approximation." (Also drops "converges", which reads as the statistical term.) |
| 64–65 | Mannered prose | "the main interpretability task left before the biology story can be told" | "the main remaining interpretability task" |

## `chapters/2026-05-07.qmd`

Clean and specific. Five semicolons, and the "Relevance" paragraph repeats the YFB-parity point already made in Results.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 9 | Semicolon | "Fixed; every external result from this point on..." | "Fixed. Every external result from this point on..." |
| 18 | Semicolon | "Reliable for LB; not yet reliable for YFB (see Open items)." | "It is reliable for LB but not yet for YFB (see Open items)." |
| 44 | Semicolon | "picks $K=3$ (1-SE rule); `point_normal` picks $K=8$" | "picks $K=3$ (1-SE rule). `point_normal` picks $K=8$" |
| 52 | Semicolon | "The sign fix matters retroactively; any external number reported before 5/5..." | "The sign fix matters retroactively. Any external number reported before 5/5..." |
| 53–55 | Recap box | "The bigger result is that merged YFB ... reaches parity with LB (0.644 both) instead of collapsing. This is the first real evidence that the YFB reformulation ... is viable, not just a cleaner derivation on paper." | Repeats line 42. "Merged YFB with per-platform z-standardization reaches parity with LB (0.644 both), so the YFB reformulation from the previous entry works in practice." |
| 62 | Semicolon | "LB's loading update still looks genomics-dominated; survival's influence on $L$ may still be too weak..." | "LB's loading update still looks genomics-dominated, so survival's influence on $L$ may be too weak even after the sign fix." |

## `chapters/2026-04-29.qmd`

Terse, readable meeting record. The main problems are semicolons, the internal "Cluster A" label, and a closing paragraph that turns the result into a lesson.

| Line | Pattern | Original | Proposed rewrite |
|---|---|---|---|
| 8 | Semicolon | "...two separately processed RNA-seq datasets; see the 2026-10-16 correction entry..." | "...two separately processed RNA-seq datasets. See the 2026-10-16 correction entry in `DECISIONS.md`." |
| 18–19 | Semicolon | "...before $L$ and $F$; a short $\beta$-only burn-in; an $\alpha$ warmup schedule; and `normalize_AB`..." | "...failed to rescue $\beta$). The changes were reordering the inner CAVI loop so $\beta$ updates before $L$ and $F$, a short $\beta$-only burn-in, an $\alpha$ warmup schedule, and `normalize_AB`, which rescales..." |
| 34 | Semicolon | "Puleo 0.550 vs. 0.626; only Moffitt improved, 0.541 vs. 0.529)" | "Puleo 0.550 vs. 0.626, with only Moffitt improving, 0.541 vs. 0.529)" |
| 16, 26, 32, 40 | Other | "Four in-model changes ("Cluster A")" ... "The Cluster A changes fixed..." ... "after Cluster A" | Internal session label (AGENTS.md line 58). "Four in-model changes, adopted together..." then "These four changes fixed the training-side collapse..." and "after these changes, relative to the earlier merged fit" |
| 20–22 | Showy verb | "That last one mattered most: ... the survival term carried roughly ... so survival had essentially no say in the loading update without it." | "The last change mattered most. At iteration 1 on the merged cohort, the survival term had roughly $0$–$2\times10^{-4}$ of the weight of the genomics term, so without it survival had almost no influence on the loading update." |
| 27–30 | Other | "...confirmed there's real signal to find. 5 of 20 factors ... so real survival signal was present in the data." | Same point stated twice. "An unsupervised EBMF diagnostic on the same matrix found 5 of 20 factors survival-associated, the strongest at C-index 0.629 ($p=3\times10^{-6}$), so the data do contain survival signal." |
| 32 | Tour-guide transition | "The catch: external validation got *worse*..." | "External validation got *worse* on 4 of 5 cohorts after these changes..." |
| 42–43 | Lesson drawn | "A fit with low training error that doesn't generalize is not progress on the problem we actually care about." | Delete. The previous sentence already states the finding. |
