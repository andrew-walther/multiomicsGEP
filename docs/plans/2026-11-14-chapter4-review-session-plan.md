# Chapter 4 (SSBMF) prelim review: session plan

Written 2026-11-14 for a review session run from `bios-dissertation`. Author: Andrew Walther.
Decisions in Section 3 were settled with Andrew on 2026-11-14. Revised 2026-11-15 with lessons
from the Chapter 3 review (`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch4/plan-edits-proposed.md`,
P-01 to P-18, all accepted).

The chapter is Chapter 4 of the prelim (Project 3). Its rendered draft is
`bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/project3-ssbmf-draft.pdf`
(54 pages including front matter and references).

## 0. Read first: where edits go

- **The only editable source is** `multiomicsGEP/paper/prelim/project3-ssbmf.qmd`.
  The `draft/` folder in `bios-dissertation` is generated. Its pre-commit hook
  (`prelim/full-draft/tools/pre_commit_guard.sh`) refuses hand edits to it.
- Committing the source in `multiomicsGEP` runs a post-commit hook that re-renders and syncs
  the draft (`paper/prelim/tools/sync_to_prelim.sh` does the same by hand; it never pushes).
- Edits Andrew makes in the rendered `.tex` (for example in Lattice) are carried back with the
  `prelim-edit-sync` skill (`prelim/full-draft/tools/prelim_review.py`).
- Render with RStudio's Quarto 1.10, not the `quarto` on PATH (see `bios-dissertation/AGENTS.md`).
  After each render, compare the page count with the committed PDF.
- **Time.** The full Chapter 3 review, run with the same phases on 2026-11-15, took about six
  hours, including both independent reviews and the closing checks. Chapter 4 is shorter, with no
  journal manuscript to revise alongside it, so the whole review (Phases 1 to 6) should take four to
  six hours. Confirm the written-prelim date in Phase 1 (`bios-dissertation/ROADMAP.md` lists
  2026-11-16). The current chapter (`97eaf20e`) remains the version of record until the revised one
  is committed.

## 1. Goals for the session

The chapter should read as a journal- or thesis-style chapter on the current state of the work:
what the project aims to do, the current method and what it has shown, and the open items that
continuing work will address (tuning the survival weight w_S and the extension to several data
modalities).

1. Keep the Introduction / Methods / Results / Discussion form.
2. **Curate, do not collect.** The chapter is not a record of every experiment. Include a result
   only if it is relevant, accurate and changes what the reader understands about the current
   method. Experiments that are still in progress, or that did not change the method, stay out
   or get at most a short statement such as "Efforts to tune the supervision weight $w_S$ by
   ... are in progress." Nothing is lost: content taken out of the chapter is kept in
   `multiomicsGEP/paper/prelim/notes/trimmed-material-2026-11.md` (a new dated section), for
   the thesis or a manuscript.
3. Make sure every method and result that remains is up to date with the latest `multiomicsGEP`
   commits.
4. **Keep the body focused.** Chapter 2 (Project 1) is 51 pages and this draft is 54. The page
   count is a guide, not a limit (D3). The body carries the main argument, and bulky methods
   and results detail moves to Appendix C. Trim the Results most, and remove white space left
   by figures that sit alone on a page.
5. Revise the prose section by section for flow and clarity for readers outside the project:
   remove jargon, internal labels and ambiguous wording or claims.
6. Run two independent reviews (writing and committee) whose output is a list of targeted
   current-text / proposed-text edits that Andrew accepts or rejects by ID.

## 2. Context for a cold session

### Project state as of 2026-11-14 (`multiomicsGEP` at `14af817d`)

- **Working model.** Survival-weighted YFB fit (η = (YF)β, the risk score projects measured
  expression onto the loadings), point-Laplace priors on L and F, per-program Normal prior on β,
  w_S = 1, rank 5 chosen by five-fold training CV, 2,218 genes, 273 TCGA and CPTAC training
  patients with 139 deaths. Evaluated once in five external PDAC series.
- **At w_S = 1 survival barely moves the loadings on real data.** In one run of five-fold
  training CV the supervised fit and the rank-matched two-step fit have nearly identical
  concordance, and their loadings correlate above 0.99999. Fixed weights of 40 or more do not
  converge with 2,000 genes in simulation. Choosing a w_S that lets survival shape the loadings
  is the main open problem.
- **Two supervision studies finished after the chapter was last edited** (commit `31411deb`):
  setting w_S from a target survival precision share
  (`docs/findings/2026-11-13-share-targeted-weight-findings.md`) and a fixed prior on the
  standardized survival coefficients (`docs/findings/2026-11-14-bounded-risk-findings.md`).
  Neither changed the working model. Per Decision D2 they get at most a brief mention.
- **Multimodal extension.** Proposed and in development. The unified fitter reproduces the
  single-modality fit with one modality. Preliminary joint expression and methylation results
  exist on matched TCGA patients. The joint fit is not validated in simulation or in external
  data.

### Chapter as it stands (word counts per section of the source)

| Section | Words | Notes for review |
|---|---:|---|
| Introduction: Background and Motivation | 690 | |
| Introduction: Aims and Contributions | 590 | Four aims, chapter scope, open problems, outline paragraph |
| Methods: Model, Weight and Objective, VI | 910 | |
| Methods: Choosing the Survival Weight | 370 | The tested CV rule plus proposed alternatives; check it still describes the current w_S scheme |
| Methods: Prediction, Simulation, Data, Tuning, Screen/Rank, Priors | 1,650 | |
| Methods: Proposed Extension to Additional Data Modalities | 890 | Longest Methods subsection, and it is proposed work; candidate for trimming |
| Methods: Comparators, Risk Groups, Program Classification | 410 | |
| Results (eight subsections) | 7,200 | Too long. Programs (1,660) and Training CV (1,430) are the longest |
| Discussion: Interpretation of the Results, Limitations, Future Work | 730 | Future Work must name the open items |
| Appendix (two prior comparisons) | 950 | |

### Models for structure and content

- **Andrew's Chapter 2 (Project 1)**,
  `bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.qmd`.
  The finished chapter to match for style and section names: Introduction (Background and
  Motivation, Aims and Contributions), method sections, Simulation Study, Application, Discussion
  (Interpretation of the Results, Limitations, Future Work), short supplementary sections. The
  same heading names are now used in Chapters 2, 3 and 4 (`97eaf20e`). The label
  `{#sec-future-considerations}` is kept so existing references still work.
- **Chapter 2 of Amber Young's dissertation (DeSurv)**,
  `~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Liu)/papers/DeSurv/Young_unc_0153D_24541.pdf`
  (chapter starts on page 22, Appendix A on page 98). **Use it only to judge what kind of
  content belongs in the chapter. Do not copy or paraphrase its text.** What it shows: about
  20 pages of body and a long appendix; a compact Methods (formulation and notation, model,
  hyperparameter selection and CV, simulation studies, datasets); four or five Results
  subsections, each built around one finding; algorithm details, update derivations and the
  convergence argument in the appendix. This is a guide to proportion, not a rule, since the
  prelim is not the final thesis.

### Source documents

| Need | Where (paths relative to `~/GithubProjects/`) |
|---|---|
| Chapter source | `multiomicsGEP/paper/prelim/project3-ssbmf.qmd` |
| Number ledger for every value in the chapter | `multiomicsGEP/paper/prelim/notes/chapter4-number-sources-2026-11-06.md`; rank-5 exhibits in `multiomicsGEP/docs/findings/2026-11-12-k5-exhibits-findings.md` |
| Newest results | `multiomicsGEP/docs/findings/` (November 12 to 14 files) |
| Current status and decisions | `multiomicsGEP/PROJECT_STATUS.qmd`, `DECISIONS.md`, `ROADMAP.md`, `docs/progress_book/current-state.qmd` |
| Material cut in the November revamp | `multiomicsGEP/paper/prelim/notes/trimmed-material-2026-11.md` (also the place to record anything cut now) |
| Earlier write-ups | `multiomicsGEP/paper/writeup/2026-11-13-writeup/` (closest to the chapter; has `NUMBER-SOURCES.md`), `2026-11-06-writeup/`, `2026-10-16-writeup/` (appendices: `single_modality.qmd` derivations, `multimodal.qmd`, `diagnostics_and_provenance.qmd`) |
| Derivations | `multiomicsGEP/derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.pdf`, `derivations/SSBMF-LB-vs-YFB-comparison/`, `derivations/SSBMF-YFB-multimodal-model-derivations/multimodal_YFB_derivation.pdf` |
| Review formats to copy | `bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch1/proposed-edits.md` and `committee-proposals.md` |
| Earlier Chapter 4 review rounds | `bios-dissertation/prelim/full-draft/notes/human-writer-review/ch4/` and `round2`–`round4` |

### Writing rules that apply (from both repos' AGENTS.md and Andrew's feedback)

- Present only the current model against external benchmarks (two-step EBMF→Cox, larger-rank
  EBMF→Cox, DeSurv, elastic-net Cox, PurIST). No comparison with earlier versions of SSBMF.
- Name models by configuration, never by date. No internal labels (S-16, "Task A", "D4").
  Say "2,218 genes" and how they were chosen.
- No "registered", "gate" or "criteria set before the fits" language. Say narratively how each
  choice was made and what the result was.
- False-active programs are expected and have small coefficients. Do not frame them as
  pass/fail. Decisions rest on recovery, convergence and training-CV concordance.
- Say "one run of five-fold CV", never "one fold split".
- Planned checks belong in Methods, not the Introduction.
- Appendix subsections are brief: motivation, "Table X compares...", the finding.
- No mannered prose. Avoid semicolons and dashes as punctuation. Title Case headings. Bold
  matrices and vectors, italic scalars. Table captions above, figure captions below.
- Apply the `human-writer` skill to all drafting and revision.
- No ", so" as a connective (40 in the source on 2026-11-15). When replacing several, vary the
  replacement ("which means", "and therefore", a new sentence) rather than using one phrase
  every time.
- No "whereas" (1 in the source). No sentence-initial "Because" (4 in the source).
- No paragraphs that open with a bold lead-in. Use a `###` heading instead; in Appendix C keep
  these headings out of the table of contents (see D6).
- Define every abbreviation at first use, in the body and again in Appendix C if it is first used
  there. Common ones such as iid need no definition.
- American English throughout: spelling and usage.
- "Fitted" only as an adjective before a noun ("the fitted model", "the fitted coefficients").
  The passive verb is "fit": "the model is fit by maximizing...", never "is fitted by" or "is
  fitted to". The same holds for the active past ("we fit the model") and a verb right after a
  noun ("a line fit to them"). Review A flags every misuse.

## 3. Decisions (settled with Andrew, 2026-11-14)

- **D1. Proposed work.** Proposed methods stay in Methods, in a last subsection after the
  current method is fully described, and are named again in the Discussion's Future
  Work. One "Continuing Work" subsection holds all of the proposed methodology,
  with working parts "Supervision Weight Tuning" and "Extension of SSBMF to Several Data
  Modalities", plus any other key outstanding methodological item or result (D10).
- **D2. The two newest supervision studies.** Include them only if they matter to the reader.
  At most a brief comment, used to explain how supervision is currently set through w_S and why
  choosing w_S remains open. No new tables, figures or appendix sections for them.
- **D3. Length and focus.** 54 pages is a guide, not a hard limit, and the finished document
  can run longer. The goal is focus: the body carries the main introduction, methods, results
  and discussion, and the appendix holds the detail a reader can consult for more depth
  (derivations, the optimization algorithm, full gene lists, per-program survival figures,
  every tuning figure). Trim the Results most, and reduce white space from figures and tables
  placed alone on a page (float placement, sizing, combining panels).
- **D4. Joint expression and methylation model.** Described as continuing work in Methods, its
  preliminary results reported in Results (including the methylation-specific programs), and
  both named in Future Work.
- **D5. Scope and tense.** This is a project-proposal chapter for a project with substantial
  progress. It reads as a thesis chapter on the current state of the work and gives proper
  attention to what is proposed next: tuning the survival weight w_S, the extension to several
  modalities, and any other key outstanding item (D10).
  Completed work is in the past tense. Proposed work, in Continuing Work and Future
  Work, is in the future tense.
- **D6. Appendix C layout.** One numbered list of sections, C.1, C.2, ..., with no "Part" or
  "Supplementary Methods/Results" headings. Methods sections come first, then results sections,
  as in Appendix A of Amber Young's dissertation. That appendix and the CTJ supplement
  (`SpatialCRT/projects/IncidenceDesign/paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex`)
  are references for ordering and kind of content, not templates.
  If any section gets `###` subheadings, check how `paper/prelim/tools/prelim_transform.py`
  numbers them. In the prelim they become C.n.m and would enter the table of contents at depth 2.
  Keep them out as Chapter 3 does: a raw LaTeX line after the appendix heading,
  `\addtocontents{toc}{\protect\setcounter{tocdepth}{1}}`, and a reset to `{2}` at the end of the
  file. Check the table of contents of the rebuilt prelim.
- **D7. Headings and names.** Every heading must make sense to a reader scanning the table of
  contents who has not worked on the project: clear, specific to what the section shows, and
  brief. No internal labels or undefined abbreviations in headings (not "LB", "YFB", "D4",
  "S-9"). Only abbreviations defined in the Introduction (SSBMF, EBMF, PDAC) may appear. Let the
  actual content and findings set each heading. The two fits being compared are the "joint SSBMF
  fit" and the "two-step EBMF and Cox fit", used consistently.
- **D8. Reuse of the write-ups.** Content in `paper/writeup/2026-10-16-writeup/`,
  `2026-11-06-writeup/` and `2026-11-13-writeup/` can be copied or revised into the chapter or
  the appendix wherever it is relevant and current. Check each number against its source before
  reuse, and drop history (dated model names, superseded configurations). Their structure
  (Methods, Evaluation, Results, a section on the programs, Discussion with Continuing Work) is
  also a guide for ordering.

- **D9. Describing w_S as a power on the likelihood (Yusha's comment, 2026-11-14).** Yusha
  questioned the term "generalized-Bayes" on page 5 of
  `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.pdf`, which calls the
  w_S ≠ 1 objective the bound for "a power (generalized-Bayes) posterior". The chapter (The SSBMF
  Model, currently source line 89) says the objective "is a variational bound for an update that
  raises the Cox partial likelihood to the power w_S" and cites Bissiri et al. (2016) and Holmes
  and Walker (2017). Weight tuning is continuing work, not a main part of the proposal, so in
  the chapter:
  - Describe w_S as a weight on the Cox log partial likelihood in the objective, with w_S = 1
    the standard variational Bayes fit that the reported results use.
  - The term "generalized Bayes" (or a power-likelihood reading) may be used only if it is
    defined for the reader and cited, never dropped in without context. Keep it brief within the
    Continuing Work subsection (D1, D10). The point is that tuning w_S is work in progress: some
    approaches have been tried, and settling the weight and its interpretation is an open item.
  - Where the reading is used, say what it does not cover: only the survival block is weighted
    while the expression term keeps weight 1; the Cox partial likelihood is not a full
    likelihood (Bissiri et al.'s loss-based updating may cover this; verify in the paper before
    saying so); the cited calibration work mostly concerns tempering with powers at most 1,
    while w_S > 1 upweights survival; and the variant that weights both the Cox term and the KL
    term for **β** is not a power posterior on **β**.
  - Check that the citations support the exact claim made. Review B is asked to check this
    (Phase 5).
  - The derivation document's wording is outside this review. Revisiting it, and the formal
    framing of w_S, is thesis work after the prelim.

- **D10. What Continuing Work contains.** Proposed methodology and results that a reader can
  understand as stated: the supervision weight and how it will be chosen, the several-modality
  model with its specified form, and any other key outstanding conceptual item or result needed
  to complete the project (for example, choosing the number of programs reliably when
  supervision is stronger, if the session judges it key). No internal study labels (such as
  "S-15") and no software deliverables in Methods or Results. Each item says what will be done
  and how it will be evaluated. For the several-modality model, state the results it is expected
  to give: mainly gains in biological interpretation, from programs shared across expression
  and methylation and from methylation-specific programs, rather than an expected improvement
  in the risk score's concordance. The software (an R package) may be mentioned only in the
  Discussion's Future Work.

- **D11. Math accuracy (added 2026-11-14).** Every equation, symbol and formula in the chapter and
  Appendix C must be correct for the current code and the fits behind the reported results. This
  covers the text as it stands and every edit.
  - Each symbol is defined at or before its first use and keeps one meaning across the body and
    Appendix C (for example, w_S, K, η, the priors on L, F and β, and the C-index).
  - Math from outside work (the Cox partial likelihood, point-Laplace and other ebnm priors, the
    variational objective, Harrell's C-index, the one-SE rule, generalized Bayes under D9) is
    cited to its source.
  - Math that describes this project's own model, updates or prediction is checked against the
    derivation and the code that produced the results: the full YFB derivation
    (`derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`, with its
    `verify_derivation.R`), `code/fit_ssbmf_yfb.R` and the update scripts it calls
    (`update_L_surv_YFB.R`, `update_F_surv_YFB.R`, `update_beta.R`, `update_tau.R`),
    `code/predict_ssbmf_yfb.R` for new-patient risk scores, `code/select_K.R` for rank selection,
    the LB comparison (`derivations/SSBMF-LB-vs-YFB-comparison/`, `code/fit_ssbmf_lb.R`) for the
    alternative parameterization, and for the joint model the multimodal derivation
    (`derivations/SSBMF-YFB-multimodal-model-derivations/`) and `code/fit_multimodal_yfb.R`.
    Settings (w_S, α_F, K, prior families, centering) are checked against the configuration the
    reported fits actually ran with, not against comments or defaults.
  - Where the derivation, the code and the chapter disagree and it is unclear which is right, ask
    Andrew. Never change code or the derivation to match the chapter.
  - An edit that adds or changes math names its source in the Phase 4 approval or the review item.

### Answers from Andrew, 2026-11-15 (Phase 2 round 1)

- Deadline is tonight, 2026-11-15, 11:59 pm.
- Results headings: descriptive and short, not finding-style.
- Continuing Work gets a third item, "Number of Programs Under Stronger Supervision".
- Aims and outline: no software or R package outside Future Work.
- C.12 is kept and rebuilt: the saved rank-5 fits at w_S = 1 and w_S = 0
  (`k5_working_2026-11-12/s11/fits/size_D4_N2218_w{1,0}_f0.rds`) match one to one, with loading
  correlations 0.999994 to 0.9999997, score correlations at least 0.99998 and 49 or 50 of the 50
  largest loadings shared. The rank-5 vs rank-7 comparison goes to C.9.
- C.14: 50 genes per active program split by sign. Signed (bipolar) loadings are noted as a
  limitation: programs are contrasts, not additive parts. The nonnegative point-exponential prior was
  tested (C.10) and lowered concordance.
- Converged preliminary two-modality simulations may go to Appendix C.
- Math: define every symbol (MC undefined rows), clear up the inconsistent rows, cite the uncited
  ones (one-SE rule: Breiman et al. 1984 and Hastie et al. 2009, as used by DeSurv).
- Drop the "73 checks" count from the chapter (not a relevant detail).
- Ties: refit with Breslow if quick, else report. Check 2026-11-15: on the full training data the
  two-step Cox risk scores under Efron and Breslow correlate 0.9999998 and C differs in the fifth
  decimal, so the exhibits are not rerun; the chapter states Efron for the two-step Cox step.
- Rescaling (MC-24): the derivation is right (`SSBMF-YFB-model-derivation-full.tex` L536–543,
  L903–905). Rescaling changes neither the reconstruction, the risk scores, the objective nor the
  survival share. It fixes the units of β_k and keeps thresholds on a fixed scale. The chapter's
  claim about w_S meaning the same for every program is removed.
- Generalized Bayes (the D9 sentence at L89): drop the term unless a citation clearly supports it,
  and cite a source for the weighting itself.
- Multimodal model: say "in development" (not "paused").
- Abstract: fix X-01 (only Program 2 stable across folds; Programs 4 and 5 survival-active in the
  adverse direction) and X-02 (rank already chosen by training CV). Update the fact sheet. Slide
  outline is after-submission work.
- All outdated currency rows are fixed.

### Answers from Andrew, 2026-11-15 (Phase 2 round 2)

- Ties: state Breslow for the SSBMF fit and the two-step Cox step (risk scores correlate
  0.9999998; the reported two-step exhibits used Efron, disclosed in one appendix clause unless
  Andrew strikes it at sign-off).
- Body adds one sentence that at w_S = 1 the joint and two-step programs are the same (unlike
  DeSurv against NMF), and that this may change once the supervision weight is settled.
- Signed loadings make program interpretation harder: a Limitation.
- Content map confirmed with the changes below.

**Confirmed map (2026-11-15).** The Section 3a tables stand with these changes. Methods and
Results headings as in the tables, results headings descriptive. Continuing Work has three
items: Supervision Weight Tuning, Extension of SSBMF to Several Data Modalities, Number of
Programs Under Stronger Supervision. Appendix C: C.1 Model Specification and Prior
Distributions; C.2 Variational Update Algorithm and Derivations; C.3 Alternative
Parameterization Linking Survival and Patient Scores; C.4 PDAC Datasets; C.5 Cross-Validation
Procedure for Model Selection; C.6 Rules for Risk Groups and Program Classification; C.7 Tuning
the Survival Weight via Cross-Validation in Simulation; C.8 Recovery of Prognostic Programs With
Weak Expression Signal (adds the non-matched variance paragraph and the two-program simulation in
one sentence); C.9 Choice of Gene Set, Input Scaling, Number of Programs and Survival Weight (adds
nested-CV weight selection, `fig:size-curve`, and the rank-5 vs rank-7 correspondence); C.10
Comparison of Alternative Priors; C.11 Comparison of the Two Model Parameterizations; C.12
Agreement of Joint and Two-Step Programs (from the saved rank-5 w_S = 1 and w_S = 0 fits); C.13
Preliminary Simulations of the Joint Expression and Methylation Model (two-direction and
three-program simulations; the initial simulation is held out); C.14 Survival Analysis of
Programs in the External Cohorts; C.15 Gene Lists for Each Program. Page budget: body 28–30 pp,
Appendix C 15–18 pp.

### Answers from Andrew, 2026-11-15 (Phase 4, Introduction)

- Introduction items I-01 to I-10 accepted; abstract option (a) (drop rank selection from the
  proposed-work sentence, no compensating cut needed).
- Nothing from the math or currency checks needed settling before Phase 4; their rows are applied
  section by section.
- Rescaling of loadings is not a methods issue for later work.
- Joint vs two-step: the takeaway for Results and the Discussion is that the current supervision
  does not change the factorization (apart from recovering a very low-variance survival program
  in simulation), unlike DeSurv, whose supervision visibly changed its NMF factors.
- Signed-loadings Limitation: state it as a tradeoff of point-Laplace against point-exponential
  priors; point-Laplace was chosen for now because the expression inputs are z-standardized.
- Generalized-Bayes sentence: keep one defined, cited sentence (Bissiri et al. 2016) with its
  qualifier; cite DeSurv (Young2026) for the weighting precedent; Mairal et al. (2008) not used.

### Answers from Andrew, 2026-11-15 (Phase 6)

- Review A: groups A-G and I accepted (A6 with the alternative wording); H accepted except H11
  (state the stability issue without "numerical warnings"); J accepted with claims kept to the
  findings. "Therefore, ..." replaces "... are therefore ...". Define and cite every acronym.
  Drop "working model"; define the projection model (η = YFβ) and the latent-score model (η = Lβ).
- J3/C3: validation cohorts only score the fixed model; the Limitations sentence is outdated.
  J17/J18: correct the comparisons and define each study. J19: add the mixed-signal study.
- Review B: C1-C30 accepted, with C-01 narrowed to the planned comparison, C-14 without Efron,
  C-09 as the full derivation and algorithm in Appendix C.2, and C-30 defining a factor as the
  column pair (l_k, f_k) and a program as the genes the factor describes.
- The derivation appendix makes no reference to Andrew's repositories.
- Table 4.7 should not sit alone on a page. Future Work should add extensions beyond the current
  aims, as Chapters 2 and 3 do.

## 3a. Content map: chapter body and Appendix C (confirmed 2026-11-15; changes in Section 3)

The current chapter has 12 Methods subsections, 8 Results subsections, 10 figures and
14 tables, plus two appendix sections. The map below is a starting proposal and the headings
are drafts. "Now" gives the current subsection or exhibit. Anything held out of the chapter goes
to the notes file.

### Methods (body): 12 subsections become 8

| Proposed subsection | Content | Now |
|---|---|---|
| The SSBMF Model | Joint EBMF and Cox model. The risk score projects measured expression onto the loadings, η = (YF)β. The survival weight w_S and the objective. One sentence on why supervision acts through the loadings (detail in C.3 and C.11) | Model; Survival Weight and Objective |
| Priors and Variational Inference | Priors on **L**, **F**, **β** and τ, stated. The coordinate-ascent algorithm in one paragraph, the survival share defined. Updates and derivations in C.2 | Variational Inference; part of the appendix prior text |
| Risk Scores for New Patients | How a fitted model scores patients not used in fitting: projection of their measured expression, with the training centering. Used for the simulation validation sets and the external PDAC cohorts | Prediction in New Cohorts |
| Simulation Studies | Design and metrics, `tab:sim-design` | Simulation Studies |
| PDAC Training and Validation Cohorts | Two training cohorts and five external cohorts in brief (table in C.4) | PDAC Data and Evaluation |
| Model Selection by Training Cross-Validation | Five-fold training CV, how rank 5, the gene screen and input scaling were chosen, and how w_S = 1 was chosen. Full procedure in C.5 | Tuning Protocol; Gene Screen, Input Scaling and Rank; current-practice part of Choosing the Survival Weight |
| Comparison Methods and Evaluation | Two-step EBMF and Cox, larger-rank EBMF and Cox, DeSurv, elastic-net Cox, PurIST. Concordance, calibration, risk groups. Program classification in brief (rule in C.6) | Comparators, Risk Groups and Program Classification; part of PDAC Data and Evaluation |
| Continuing Work | *Supervision Weight Tuning*: why w_S = 1 leaves survival little influence on the loadings, why large fixed weights fail, the alternatives proposed, and one sentence that tuning efforts are in progress. *Extension of SSBMF to Several Data Modalities*: the proposed model form, with patient scores shared across modalities and a loading matrix **F**_m for each modality m with its own estimated prior, and what remains open. Other key outstanding items per D10. A sentence that K = 5 currently comes from training cross-validation, and that letting the priors prune programs from a larger starting K may be revisited if it proves reliable. Future tense | Proposed parts of Choosing the Survival Weight; Proposed Extension to Additional Data Modalities (shortened) |

Sensitivity to the Priors moves to the appendix (C.10).

### Results (body): 8 subsections become 5, each stating a finding

| Proposed subsection (draft heading) | Main exhibits | Now |
|---|---|---|
| The Supervised Model Recovers Planted Prognostic Programs in Simulation | `fig:sim-s9` and `tab:strength`. One sentence that fixed weights of 40 or more did not converge | Survival Weight in Simulation (first part) |
| Joint SSBMF and Two-Step EBMF and Cox Reach Similar Concordance in Training Cross-Validation | `tab:training-cv` | Training Cross-Validation on PDAC (main result only) |
| Risk Score Concordance in Five External PDAC Cohorts | Concordance of each model's risk score in the external cohorts: joint SSBMF against two-step EBMF and Cox, DeSurv, elastic-net Cox and PurIST. `tab:external`, `fig:external`, `fig:km-risk`, `tab:comparators` | External Validation; Comparison With Other Prognostic Models |
| Biological Characterization of the Prognostic Programs | What each survival-associated program of the joint fit is: leading genes, pathway enrichment, correspondence with published PDAC subtypes. `tab:programs`, `fig:programs`, `tab:enrichment`. Survival curves and gene lists go to C.13 and C.14 | Programs (main part) |
| Preliminary Results for Methylation and the Joint Expression and Methylation Model | Methylation alone in a few sentences, the methylation-specific programs found, `tab:joint-tcga`. Labeled preliminary | Methylation Alone; Preliminary Joint Expression and Methylation Results |

### Appendix C: methods sections first, then results sections

| Section | Content | Source |
|---|---|---|
| C.1 Model Specification and Prior Distributions | The full model and the priors on **L**, **F**, **β** and τ (the algorithm is in C.2) | current Methods and appendix text |
| C.2 Variational Update Algorithm and Derivations | The coordinate-ascent algorithm: updates for **L**, **F**, **β** and τ, the survival share derivation, how each update was checked against the objective | `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.pdf`; 10-16 write-up `appendices/single_modality.qmd`; 11-06 write-up Appendix A |
| C.3 Alternative Parameterization Linking Survival and Patient Scores | The model only: Cox coefficients on the latent patient scores (η = **L**β, new patients scored by fold-in), set against the projection form used here (η = (**YF**)β). Results are in C.11 | `derivations/SSBMF-LB-vs-YFB-comparison/` |
| C.4 PDAC Datasets for Model Training and Validation | Table of cohorts: platform, patients, deaths, role, handling of the two PACA-AU series that share donors | number ledger; 11-13 write-up |
| C.5 Cross-Validation Procedure for Model Selection | Fold construction, training centering, the one-SE rule for the number of programs (K from 2 to 12), gene screen, input scaling, and why the rank is not chosen by pruning | Tuning Protocol; Gene Screen, Input Scaling and Rank |
| C.6 Rules for Risk Groups and Program Classification | `tab:program-rule` and the training cutpoints for risk groups | Programs; 11-06 write-up Appendix D |
| C.7 Tuning the Survival Weight via Cross-Validation in Simulation | `fig:sim-s13`: concordance, calibration and convergence when w_S is chosen by cross-validated partial likelihood. Shows how the weight is chosen | Survival Weight in Simulation |
| C.8 Recovery of Prognostic Programs With Weak Expression Signal | `fig:sim-s14`: recovery and concordance when the planted program explains little expression variance. Shows where supervision adds value over the two-step fit | Survival Weight in Simulation |
| C.9 Choice of Gene Set, Input Scaling and Number of Programs | `tab:inputs-rank` (includes the rank comparison behind K = 5), `fig:size-curve` | Training CV |
| C.10 Comparison of Alternative Priors | Two parts: priors on the scores and loadings (`tab:lf-prior-sim`, `tab:lf-prior-real`, point-Laplace against point-exponential) and priors on the survival coefficients (`tab:beta-prior`, per-program Normal against shared sparse priors) | current appendix |
| C.11 Comparison of the Two Model Parameterizations | Results for the model in C.3 against the working model: training and simulation concordance, program agreement | Latent-Score and Projection Parameterizations; `docs/findings/2026-11-12-lb-yfb-study-findings.md` |
| C.12 Agreement of Joint and Two-Step Molecular Programs | Loading correlations and program correspondence of the joint SSBMF and two-step EBMF fits | `results/benchmark_sim/outputs/k5_working_2026-11-12/program_correspondence/` |
| C.13 Survival Analysis of Programs in the Real-World Datasets | `fig:km-program-2/4/5`, `tab:external-program-direction` | Programs |
| C.14 Gene Lists for Each Factor | See the gene-list question below | `results/benchmark_sim/outputs/k5_working_2026-11-12/prelim_ch4/` |

Headings are drafts, several in Andrew's wording.

**Factor and program.** They refer to the same object. A *factor* is a column of the
factorization of the data matrix **Y**. A *program* is the set of genes (or CpG sites) that load
on that factor. Define both once in the Methods (The SSBMF Model) and use each in that sense:
"factor" for the matrix objects, "program" for the biological gene set.

`tab:strength` (joint SSBMF against the two-step fit in each simulation cell, Andrew's draft
heading "Comparison of Joint Supervised and Two-Step Models in Simulation") moves to the body, in
the first Results subsection.

**Gene lists (C.14), settled.** For each active program, list the 50 genes with the largest
absolute posterior mean loading, split by the sign of the loading (genes that raise the program
score and those that lower it), in a compact multi-column table with the rule stated in the
caption. Earlier exhibits used 20 to 25 genes. Alternative, if the saved fit stores each
loading's posterior probability of being nonzero: list genes with probability above 0.95, capped
at 100. Check the saved fit first and use the alternative only if it is available. The full
loading matrix stays in the repository. The gene-by-program heatmap of the rank-5 fit
(`fig:programs`, `figures/sup_programs_heatmap_k5.pdf`) stays in the body.

**Judgement on inclusion.** The map is not a cap. A result that is valuable and current should
be included, in the body or the appendix, even if it is not listed here. Hold out results that
are still in progress, superseded, or do not change what the reader understands. Watch for
bloat in both places: the chapter and the appendix should not overwhelm the reader or bury the
key results under every detail that was examined. Each appendix section should be one the body
points to for a reason.


## 4. Work plan

Record as you go. After each phase and each of Andrew's decisions, write it into this file:
answers under Section 3, the confirmed map in Section 3a, and a dated **Status** line under the
phase with commit hashes. Commit the plan in `multiomicsGEP`. A compacted or fresh session resumes
from this file and the review files, not from the conversation. (The Chapter 3 session was
compacted several times and lost nothing because of this.)

Model use. The main session (Opus) plans, decides with Andrew, and checks every subagent result
before relying on it. Delegate bounded work to Sonnet subagents (model "sonnet"), each with a
complete brief (goal, files, what done looks like):
- applying approved edits from a sign-off list or review file, item by item, and reporting the
  diff for the main session to check;
- the Phase 3 currency check, citation lookups on Crossref or PubMed, render and page-count
  checks, numbering and cross-reference checks, and searches across files;
- the writing review (human-writer).
Keep on Opus: content-map and restructuring decisions, new or substantially rewritten prose, the
math check, and the committee-member review, since these depend on statistical reasoning.

### Phase 1. Set-up (about 15 minutes)

1. Confirm the deadline.
2. Confirm the source and the generated draft are in sync (`git status` in both repos; the draft
   PDF dated after the last source commit).
3. Record the baseline: source commit hash, draft page count, full-prelim page count, and the
   pages where a float stands alone or leaves a large gap. Known: in the prelim's List of Tables,
   Table 4.8 ("Programs of the working model on the training cohorts (n = 273)") has a run of
   blank space before its dot leader.

**Status 2026-11-15 (Phase 1 done).** Deadline confirmed with Andrew: the written prelim is due
tonight, 2026-11-15 at 11:59 pm (not 11-16). Source `97eaf20e` and draft (bios `a19adb8`) in
sync, both repos clean. Baseline: draft 54 pages, full prelim 240 pages. Draft pages with few
words (floats alone or gaps): 21, 23, 24, 37, 38 (97 words), 39–41, 50, 54. Phase 3 and 3b
subagents launched.

### Phase 2. Outline and curation review (about 45 minutes, with Andrew)

Go through the content map in Section 3a row by row. For each subsection and each figure or
table, decide: keep in the body, shorten, merge, move to Appendix C, or hold out of the chapter.
The test for each result: is it relevant to the current method, accurate as of the latest
commits, and does it change what the reader understands? Answer the questions at the end of
Section 3a. Write the confirmed map, with a rough page budget per section, into the session notes
before any prose edits. Check that the Aims paragraph, the "rest of the chapter" paragraph and
the Discussion all match it.

Record anything held out in `multiomicsGEP/paper/prelim/notes/trimmed-material-2026-11.md`
(a new dated section) so it can be reused in the thesis or a manuscript.

**Status 2026-11-15 (Phase 2 proposal sent to Andrew, awaiting confirmation).** Findings while
checking map sources: (a) `program_correspondence/` compares the rank-5 and rank-7 fits, not the
joint and two-step fits; joint–two-step agreement exists only as loading correlations above
0.99999 (training-CV redo and S-16 findings), so C.12 is proposed to be dropped and the rank-5
vs rank-7 correspondence folded into C.9. (b) The saved rank-5 fit
(`k5_working_2026-11-12/s11/fits/size_D4_N2218_w1_f0.rds`) stores posterior means only (EF,
EBeta, EBeta2), no inclusion probabilities, and every loading is nonzero, so C.14 uses the
top-50-by-|loading| rule. (c) Real-data w_S selection (nested CV) proposed for C.9 with the other
training-CV tuning results.

### Phase 3. Currency check (subagent, runs alongside Phase 2)

A fresh subagent compares every number and status claim in the chapter with its source in the
ledger and findings files, and lists any claim made out of date by commits after `31411deb`. It
only checks what the chapter says. It does not propose adding new results. Output:
`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch4/currency-check.md`,
one row per claim: location, current text, source, status (current / outdated / unsourced),
proposed text. It also compares the chapter's Project 3 claims with
`bios-dissertation/prelim/abstract/abstract-text.qmd` and `notes/abstract-fact-sheet.md`, Chapter 1
Section 1.3 (`prelim/literature-review/draft/`) and `presentation/prelim-slides/OUTLINE.md`, and
flags mismatches only. The abstract fills pages iii and iv exactly, so any edit there needs an
equal cut, and Andrew decides it.

### Phase 3b. Math check (subagent, runs alongside Phase 2)

A second fresh subagent applies D11 to the chapter and Appendix C. Output
`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch4/math-check.md`, one row
per equation or symbol: location, the math as written, its source (citation, or file and line of
the derivation or code), status (verified / cited / wrong / undefined / inconsistent / uncited),
and proposed text. Wrong, undefined or inconsistent rows are fixed in Phase 4, before the
Phase 5 reviews.

### Phase 4. Content edits (section by section, with a written sign-off list)

Apply the target outline and the currency fixes first, so the reviewers read the revised
chapter. For each section, write the proposed edits to `round5/ch4/section-edits.md` in the
Chapter 1 round-5 format: ID, the currency-check and math-check rows each item answers, location
(source line), current text and proposed text. Restructuring is described as moves (from, to,
what changes) rather than full text. Items that change a claim or a number go last. **Stop** for
Andrew's decisions by ID, recorded in a log table at the top of the file. Apply the accepted items
with a Sonnet subagent, check its diff, commit and push in `multiomicsGEP`, and show the
re-rendered section. For each section check:

- **Introduction:** the aims match what the chapter delivers, and the open problems match the
  Discussion.
- **Methods:** the current method first, in the order of the confirmed Section 3a map, in the
  past tense. Model Selection by Training Cross-Validation states plainly how w_S = 1 was chosen.
  Then the single Continuing Work subsection in the future tense: supervision weight tuning, with
  at most a sentence on the newest studies (D2), the multimodal extension with its model form
  (D4), and any other key outstanding items (D10).
  Apply D9 wherever w_S is described as a power on the likelihood.
- **Results:** each subsection opens with its question and states its finding in the first two
  sentences. Detail that does not change a conclusion moves to the appendix or is cut. The joint
  model's results are labeled preliminary.
- **Discussion:** Interpretation of the Results, Limitations, and Future Work, which names each
  open item (w_S tuning, the multimodal extension, and any other item under D10) and points back
  to the Continuing Work subsection. It may also mention releasing the software as an R package.
  The Aims and Methods should not.
- **Appendix C:** build sections C.1, C.2, ... per D6 and the confirmed map, methods sections
  before results sections. Each section is
  referenced from the body and kept brief: motivation, "Table X compares...", the finding.
- **Layout:** fix floats that sit alone on a page or leave large gaps.
- **Math:** apply the Phase 3b rows for the section, and check any new or moved equation under
  D11 before showing it to Andrew.

Commit and push the source in `multiomicsGEP` after each section or group of sections, so the
hook re-renders the draft and each step can be reverted. After the first section commit, and after
any commit that changes labels or cross-references, run
`PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$PATH" python3 prelim/full-draft/tools/build_prelim.py`
and `python3 prelim/full-draft/tools/test_prelim_assembly.py` from the bios-dissertation root
(Sonnet subagent). In Chapter 3 the chapter rendered cleanly but the full build's own checks
stopped in Phase 6 (a cross-reference followed by a colon; fixed in bios `d61df5b`).

**Status 2026-11-15 (Phase 4, Introduction applied).** Phase 2 map confirmed (Section 3). The
Introduction items I-01 to I-09 were applied to `paper/prelim/project3-ssbmf.qmd`; I-10 and the
Chapter 4 fact-sheet update were applied in bios-dissertation. Sign-off log:
`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch4/section-edits.md`.
Next: Methods sign-off list.

**Status 2026-11-15 (Phase 4, Methods list sent).** Methods sign-off list M-01 to M-20 and the
full proposed text (`round5/ch4/methods-proposed.qmd`, which also drafts Appendix C.1 to C.6)
sent to Andrew; awaiting decisions.

**Status 2026-11-15 (Phase 4, Methods applied).** M-01 to M-20 accepted except M-20c (no Efron
disclosure; Breslow stated for both fits). New style rule from Andrew, applied from here on: write
what was done and settled on; no negative statements or decision narration unless a direct
comparison is justified by results or literature. Applied in `b9d24164`; chapter 57 pp, prelim
243 pp, assembly tests pass. Next: Results sign-off list.

**Status 2026-11-15 (Phase 4, Results list sent).** R-01 to R-16 and `round5/ch4/results-proposed.qmd`
(also places Appendix C.7 to C.15) sent to Andrew. Joint vs two-step check committed
(`70057347`). Post-review steps recorded in bios ROADMAP (Ch1 review, Ch1-3 writing-rule sweep,
abstract/continuity, full-prelim review).

**Status 2026-11-15 (Phase 4, Results applied).** R-01 to R-16 accepted, R-15a revised to drop the
DeSurv contrast. Applied in `b8ae070d`; chapter 56 pp, prelim 242 pp, tests pass. Andrew asked
again that Chapter 3 be swept for negative statements after Chapter 4 (bios ROADMAP step 2).
Next: Discussion and Appendix C sign-off list.

**Status 2026-11-15 (Phase 4, Discussion and Appendix list sent).** D-01 to D-06 and A-01 to A-05
(`round5/ch4/discussion-appendix-proposed.qmd`). Gene lists script committed. Young (2026)
Sec. 2.3.4 checked: supports the DeSurv vs NMF factor-difference sentence (D-06a).

**Status 2026-11-15 (Phase 4 done).** D-01 to D-06 and A-01 to A-05 accepted (gene lists at 25 per
pole; a loading cutoff needs per-loading posterior summaries the saved fit lacks). Applied in
`c26171de`; chapter 58 pp, prelim 244 pp, tests pass, Appendix C TOC matches the map (C.1-C.15).
Name settled (Andrew 2026-11-15): SSBMF = "survival-supervised Bayesian matrix factorization".
Gene lists: Andrew asked about longer lists (about 250 per program, as in Young 2026, whose 270 is
the BO-tuned n_top that truncates W for its linear predictor). SSBMF's risk score uses all 2,218
genes, so a 250 cutoff would be arbitrary here; options sent to Andrew.

**Status 2026-11-15 (Phase 5 started).** Review A (writing, Sonnet) and Review B (committee, Opus)
launched in parallel on `c26171de`.

**Status 2026-11-15 (Phase 5 done, awaiting Andrew).** Appendix C.15 rebuilt at 50 genes per pole
for all five programs (Andrew; `0bb0fc8f`). Review A: `round5/ch4/proposed-edits.md`, 137 items
(groups A-J). Review B: `round5/ch4/committee-proposals.md`, C-01 to C-30 and Q-01 to Q-20; 8
verified references appended to `new-citations.bib`. Spot-checked: C-01 interval coverage
(`docs/findings/2026-11-04-survival-weight-p2000-findings.md` L26-27) and C-08 prior-variance
estimate (`code/update_beta.R:350`). Stop for Andrew's accept/reject decisions.

**Status 2026-11-15 (Phase 6 done).** Review A applied (`7cdcc3e7`), Appendix C.2 replaced by the
full derivation and algorithm (`f3e825b8`), Review B applied (`b178539c`). An independent Opus
review of the applied edits found no math errors; its fixes (stopping-rule minimum, defined cells
and block-weighted fit, two-step attribution, v_k definition, undefined "arm", positive
restatements) are in `f3ac90ae`, with the expanded Future Work and Table 4.6/4.7 placement.
Table 4.7's caption was shortened so it shares a page with text (`469ea77a`). Q&A prep has 21
questions with format tags. Standalone chapter 71 pp, synced draft 68 pp, full prelim 254 pp;
assembly tests pass (bios `c58c250`).

**Status 2026-11-15 (Chapter 4 review complete).** Done-when audit run (Sonnet); its small items
fixed (`c73071a9`, `837595e4`), review files given a decisions log and check resolutions (bios).
Andrew's final notes applied in `efb1ca6f`: Future Work citations verified on PubMed (Cao 2021
protein, phospho and glyco data; TCGA-PAAD reverse-phase protein array; PurIST for treatment
choice), a paragraph on a risk scale shared across cohorts, coverage removed from planned work
(open problem only), and one expectation notation. Every Section 5 item passes except the slide
outline, which Andrew will update with the slide revision. Synced draft 69 pp, standalone 71 pp,
full prelim 255 pp; assembly tests pass.

### Phase 5. Independent reviews (two subagents in parallel)

Run both on the revised chapter. Give each the source, the PDF, `PROJECT_STATUS.qmd` and the
writing rules in Section 2, but not this session's conclusions. Neither applies any edit.

**Review A: writing.** Uses the `human-writer` skill: it reads all of
`~/.claude/skills/human-writer/ (SKILL.md, ai-tropes.md, voice-profile.md)` before starting and applies every rule in them. Finds writing tics, typos, choppy sentences
that should be joined, ambiguous referents ("it", "this", "that" with an unclear object),
undefined or uncited terms and abbreviations, jargon an outside biostatistician would not know,
and claims stronger than the evidence. Output `round5/ch4/proposed-edits.md` in the Chapter 1
round-5 format: groups by issue type (A typos and grammar, B ambiguous referents, C flow, D
undefined terms and jargon, ...), and for each item an ID, location (section, source line), the
current text and the proposed text. Claim-changing items go in a separate final group that needs
Andrew's decision. Group items so a whole group can be accepted at once, and put items that need
judgment (rewording a claim, cutting a sentence) in their own groups. For a repeated pattern, such
as the ", so" joins, list them in one group with a different replacement for each. Every proposed
sentence follows the writing rules in Section 2, and a proposal that adds or changes a number
names the file it comes from.

**Review B: committee member.** A biostatistics committee member reading the chapter cold. Finds
inconsistencies between sections, gaps in the argument, unsupported claims, missing comparisons
or citations, errors in notation or numbers, and questions likely at the oral prelim. It should
also flag content that could be cut without loss, flag any equation that is undefined, uncited
or disagrees with the derivation or code (D11), and check whether every description of w_S
as a power on the likelihood is accurate and supported by its citations (D9). Output `round5/ch4/committee-proposals.md` in
the Chapter 1 format: items ordered by how likely they are to come up, each with the problem,
the current text and a proposed revision. Any new reference is checked on Crossref or PubMed and
staged in `round5/ch4/new-citations.bib`, not added to `prelim/references.bib` until accepted.
Proposals that add length need a matching cut. Every proposed sentence follows the writing
rules in Section 2, and a proposal that adds or changes a number names its source file. The output
has two parts: (i) proposed text changes, ordered by how likely they are to come up, with each
proposed cut offered as options (a) and (b) where there is a real choice; (ii) questions likely at
the oral prelim that need no text change, each with a suggested answer.

### Phase 6. Decide and apply

1. Andrew replies by item ID ("accept all A, reject C3, edit B2"). Record decisions in a log
   table at the top of each proposals file.
2. Before applying an accepted item that adds or changes a number or a claim, check it against
   its source file; if it does not match, report it to Andrew instead of applying it. Apply
   accepted items to `project3-ssbmf.qmd`, commit and push in `multiomicsGEP`, and let the hook
   sync the draft. An accepted reference is added to `paper/ssbmf-refs.bib` and
   `bios-dissertation/prelim/references.bib` with identical fields before the commit that cites
   it (the sync's bib check stops on any mismatch).
3. Re-render the standalone draft and the full prelim. Compare page counts with the baseline.
   Check figure and table numbering, cross-references and the reference list.
4. Update `multiomicsGEP/PROJECT_STATUS.qmd`, `ROADMAP.md`, `DECISIONS.md` (one dated entry for
   this review's decisions) and README/AGENTS where the global rules ask, and
   `bios-dissertation/ROADMAP.md` (Chapter 4 review done, chapter and full-prelim page counts).
5. Write Review B part (ii), with a format tag for each question (main slide, backup slide or
   spoken answer), to `presentation/prelim-slides/prelim-qa-prep.md` for the Part III slides, as
   Chapter 3 did in `SpatialCRT/projects/IncidenceDesign/presentation/prelim-qa-prep.md`.
6. Done-when audit: a Sonnet subagent checks each Section 5 item against the files (grep, page
   counts, logs) and reports the evidence for each. The main session reviews the report before
   telling Andrew the work is done.

## 5. Done when

- The target outline is recorded and the chapter follows it.
- The body holds the main argument, with derivations, full gene lists, per-program survival
  figures and tuning detail in Appendix C. The body is shorter than now and has no large float
  gaps.
- Completed work is in the past tense and proposed work in the future tense.
- Every number and status claim matches its source as of `14af817d` or later (the currency
  check has no "outdated" rows left).
- Future Work names every open item, and each one appears in Methods as proposed work.
- Every equation and symbol is defined and either cited or verified against the derivation and
  code (the math check has no open rows).
- Both review files exist, every item has a decision, and accepted items are applied.
- The chapter is in American English, and "fitted" appears only before a noun (no "is fitted",
  "we fitted" or "fitted to" left).
- No ", so" joins, "whereas", sentence-initial "Because" or bold paragraph lead-ins remain, and
  every abbreviation is defined at first use.
- Appendix C subheadings, if any, are kept out of the table of contents.
- The prelim abstract, Chapter 1 Section 1.3 and the slide outline agree with the chapter, or each
  mismatch has Andrew's decision.
- The oral-prelim questions are in `presentation/prelim-slides/prelim-qa-prep.md`.
- The draft and the full prelim render without errors, the assembly test passes, and page counts
  are recorded.
- The done-when audit reports evidence for every item above.

## 6. Prompt to start the session in `bios-dissertation`

> I'm reviewing Chapter 4 (Project 3, SSBMF) for the written prelim. Read
> `../multiomicsGEP/docs/plans/2026-11-14-chapter4-review-session-plan.md` first, then this
> repo's AGENTS.md and ROADMAP.md. Edits go to `../multiomicsGEP/paper/prelim/project3-ssbmf.qmd`,
> never to the generated draft. The decisions in Section 3 are settled. Start with Phase 1, then
> walk me through the Phase 2 review of the Section 3a content map (body versus Appendix C), and
> launch the Phase 3 currency check and the Phase 3b math check as subagents while we do it.

### Full-plan prompt

To run the whole plan, Phases 1 to 6, with stops only where Andrew's input is needed:

```text
I'm reviewing Chapter 4 (Project 3, SSBMF) for the written prelim. Read
../multiomicsGEP/docs/plans/2026-11-14-chapter4-review-session-plan.md first, then this repo's
AGENTS.md and ROADMAP.md. Edits go to ../multiomicsGEP/paper/prelim/project3-ssbmf.qmd, never to
the generated draft. The decisions in Section 3 are settled. Confirm the written-prelim date with
me in Phase 1. The Chapter 3 review took about six hours with the same phases, so plan for four
to six hours.

Work through the whole plan, Phases 1 to 6, in order. Stop and wait for me only at these points:
(1) after Phase 2, to confirm the content map before any prose edits;
(2) after each section in Phase 4, to approve that section's sign-off list before you apply the
    edits, then commit and push them in multiomicsGEP;
(3) after Phase 5, when the two review files are ready for my accept/reject decisions.
Launch the Phase 3 currency check and the Phase 3b math check as subagents during Phase 2
(decision D11: all math defined, and cited or verified against the derivation and code), and run
the two Phase 5 reviews as parallel subagents.

After each phase and each of my decisions, write it into the plan (answers under Section 3, the
confirmed map in Section 3a, a dated Status line under the phase) and commit the plan, so a
compacted or fresh session can resume from the files. Keep the plan's Section 5 "Done when" list
as the finish line, and tell me which items are still open whenever you stop.
```

The session must be able to edit and commit in `multiomicsGEP`. Add that folder to the session if
asked. If the session runs long, split after Phase 4: commit, then start a fresh session with the
same prompt plus "Phases 1 to 4 are done."
