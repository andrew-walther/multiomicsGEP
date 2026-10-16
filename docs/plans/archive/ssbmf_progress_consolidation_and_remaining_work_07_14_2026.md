# SSBMF — Progress Consolidation & Remaining Work (since the 6/18/2026 lab meeting)

**Date:** 2026-07-14
**Purpose:** Two things in one document. Part A is an orientation summary — everything accomplished
since the 6/18/2026 lab-meeting deck was merged (`905279b`), for anyone (human or a fresh Claude
Code session) who needs to get oriented without re-reading 34 commits and two prior plan documents.
Part B is the actual remaining-work plan, written to be executed end to end by a new session.

**Source plans this consolidates:**
- `docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md` (6 phases; Phases 1-3 complete,
  Phases 4 and 6 not started, Phase 5 = same work as this plan's Step 5 below)
- `docs/plans/ssbmf_k_parsimony_followup_plan_07_13_2026.md` (6 steps; Steps 1-4 complete, Steps 5-6
  not started)

---

## Part A — What's been accomplished (orientation summary)

If you're picking this up cold: the short version is that the recommended model configuration
itself has **not changed** since the 6/18 deck (still YFB, DeSurv-aligned gene selection, K=7, no
cohort indicator) — what changed is that its numbers were corrected (0.636 → 0.627 external
C-index) and its K=7 choice was rigorously stress-tested from multiple independent angles, all of
which converged on the same answer: **the model only uses 2 of its 7 factors for survival
prediction (K_eff=2), and this K_eff=2 result is robust and real, not a tuning artifact.**

### Original post-lab-meeting action plan — Phases 1-3 (complete)

| Phase | What it did | Result |
|---|---|---|
| **1** (`DECISIONS.md` 2026-07-12, three entries) | Normalized the genomics/survival objective imbalance; retired the redundant `lambda` parameter (`alpha` is now the sole mixing weight); fixed a train/test preprocessing mismatch (external cohorts were being rank-transformed while training was per-platform z-standardized — the reverse). An initial version also boosted β's own precision, which was later found unjustified and reverted (`boost_beta=FALSE` is the corrected default). | External mean C: 0.636 → **0.627** (K=7, K_eff=2). The decline is attributable entirely to the preprocessing bug fix, a genuine correction kept regardless of its small effect on this one metric. |
| **2** (`DECISIONS.md` 2026-07-12) | Compared the joint YFB model against two two-step baselines (PCA+Cox, EBMF+Cox) across 10 seeds × 4 DGP scenarios. | Joint model's advantage emerges and grows from moderate signal onward in 3/4 scenarios; a 4th (`sparse_synthetic`) reversed the ordering — root-caused to a genuine, understood CAVI factor-collapse vulnerability (near-tied-amplitude, disjoint-support factors can collapse to a degenerate near-zero fixed point), not a survival-coupling bug. This finding is the seed for everything in the K-parsimony follow-up plan below. |
| **3** (`DECISIONS.md` 2026-07-13) | Built `results/benchmark_sim/run_k_parsimony_curve.R`: refit the recommended config at K ∈ {2,3,4,5,7} with a single fresh-SVD fit each, re-ran external validation for each. | K=7 won on every one of 5 held-out cohorts; no smaller K reached within 1 SE. **But** flagged that K=2/K=4's suspiciously fast convergence (7-9 iterations to near-zero β) looked like the Phase 2 collapse artifact, not a genuine capacity limit — this caveat is exactly what the follow-up plan below resolves. |

### K-parsimony follow-up plan — Steps 1-4 (complete)

| Step | What it did | Result |
|---|---|---|
| **1** (`DECISIONS.md` 2026-07-13, "Step 1") | Re-fit K ∈ {2,3,4,5} three ways: fresh SVD (Phase 3's baseline), warm-start from the converged K=7 fit's top-K PVE-ranked columns (new: `code/warmstart_from_fit.R`), and best-ELBO multistart (new: `fit_cox_on_yf_multistart()` in `code/fit_modular_multistart.R`). | **K=4 and K=5 both reach K=7-level performance (0.6270 vs. K=7's 0.6267) via warm-start.** Multistart rescued nothing at any K — SVD-init was always the best-ELBO restart. K=2/K=3 remained short even with warm-start (a genuine capacity floor). Decision rule: OPTIMIZATION-LIMITED → routes to Step 2. |
| **2** (`DECISIONS.md` 2026-07-13, "Step 2") | Tried to make Step 1's fix permanent/general: added a deflation-style init (`code/deflation_init.R`, sequential rank-1 SVD of the residual) to both `fit_cox_on_yf()` and `fit_supervised_mf_modular()`. | **Negative result, analytically understood:** deflation-init is mathematically equivalent to batch SVD-init whenever singular values are non-degenerate (essentially always true for real/most synthetic data) — bit-identical results to fresh SVD on real data. Confirms the collapse isn't about *which* linear-algebra decomposition seeds CAVI; only a warm-start from an already-CAVI-shaped solution (Step 1's approach) escapes it. |
| **3** (`DECISIONS.md` 2026-07-13, "Step 3") | User steered toward this given Step 2's negative result: built joint (K, α) Bayesian optimization (`code/select_k_alpha_bo.R`, `rBayesianOptimization` package, branch `joint-k-alpha-bayesopt`) to test DeSurv's own joint-tuning approach. Included a validity gate (`pick_trustworthy_bo_winner()`) after a dry run caught BO finding a degenerate `alpha≈0, K_eff=0` "winner" that looked plausible on raw score alone. | Found K=8, α=0.7072 → external C=0.6282 — marginally better, but **larger**, not smaller. A real limitation was found and documented: `rBayesianOptimization`'s integer-typed K bounds pre-round K to a coarse grid, causing the search to converge early (effectively ~9 unique points explored, not 23) — smaller K is not conclusively ruled out by this specific run, but Step 3 did not deliver its goal (find a smaller K). |
| **4** (`DECISIONS.md` 2026-07-13, "Step 4") | Synthesized Steps 1-3. | **K=4, K=5, K=7, K=8 are all statistically indistinguishable (spread of 0.0015) and all converge to K_eff=2.** Decision: keep K=7 as the recommended default (single-step, dependency-free reproduction) — K=4/K=5 are validated, documented, more-parsimonious alternatives that require a two-step warm-start recipe for zero measurable performance gain. `config/globals.yml` unchanged. `CLAUDE.md`/`ROADMAP.md` updated with the full, honest picture. |

**Current test suite: 311/311 passing.** All of the above is merged to `main` (see commit range
`905279b..b626deb`, 34 commits). Every step went through TDD, an independent
`superpowers:code-reviewer` pass, and a `DECISIONS.md`/`ROADMAP.md` update with real numbers — no
step's numbers are provisional or need re-verification, only re-reading if you need the detail.

### What was drafted but not yet executed
- `docs/plans/pathway_enrichment_plan.md` (a 12-step implementation plan, decided/reviewed in a
  2026-06-16 planning session) — exists only on a now badly-stale branch `pathway-enrichment-plan`;
  see Part B, Step 5.
- A scratchpad-only start on the `SSBMF-paper` manuscript repo scaffold (uncommitted, in a temp
  clone that will not survive a new session) — see Part B, Item 4.

### Two items intentionally NOT carried into Part B (already satisfied, no follow-on work)
The original plan's triage table (Section "Triage summary") included two items disposed of as
"✓ confirm & document" rather than assigned to a phase: **α default = 1/2** (already the default in
`config/globals.yml`, no action needed) and **h₀(t) specification / per-platform z-std-then-merge
ordering** (Breslow partial likelihood already cancels h₀; per-platform z-std already happens before
merge in `preprocess_desurv.R`). Neither needs tracking below — flagged here only so it's explicit
they weren't dropped, just already done before this document existed.

---

## Part B — Remaining work

Standing conventions (apply to every item below, matching everything done in Steps 1-4 above): TDD
for new functions; dispatch `superpowers:code-reviewer` before each commit; run
`Rscript tests/run_tests.R` after any change to a modular update/fitting script (currently
311/311 — do not proceed if this regresses); update `DECISIONS.md`/`ROADMAP.md`/`CLAUDE.md` with
real numbers as each item completes; one branch per item, cut from `main`, merge once tests+review
are clean; **never push to remote without asking**. `PDAC_DATA_ROOT` is not set as a persistent env
var — export it inline per command that needs real data:
```
export PDAC_DATA_ROOT="$HOME/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Liu)/PDAC_data"
```
Work through items 1-3 in order (each depends on the last); items 4-5 are independent of 1-3 and of
each other, and can be done in any order or interleaved. Only pause for a genuine blocker (a test
failure you can't resolve, a missing credential/access, or a true ambiguity this document doesn't
resolve) — not to re-confirm something already decided here.

### Item 1 — Step 5: Pathway enrichment (K-parsimony follow-up plan)

**Goal:** biological interpretation of the recommended config's two active factors (Program 7 =
adverse, Program 3 = protective — confirmed still correct, since Step 4 kept K=7 as the default and
didn't change which factors are active).

**Before implementing:**
1. The old `pathway-enrichment-plan` branch predates ~8,000 lines of subsequent work (all of Steps
   1-4 above) and is unusable as a base. Retrieve *only* the plan document from it, then delete and
   recreate the branch off current `main`:
   ```
   git show pathway-enrichment-plan:docs/plans/pathway_enrichment_plan.md > /tmp/pathway_plan.md
   git branch -D pathway-enrichment-plan   # old branch, being replaced
   git checkout -b pathway-enrichment-plan main
   ```
   (A copy already exists at `/tmp/pathway_plan_from_branch.md` from this session, in case that
   file survives into the new one — verify it matches before trusting it, since `/tmp` is not
   guaranteed persistent.)
2. Re-verify Program 7/3 are still the correct target factors against the current recommended-config
   fit — Step 4 didn't change K=7 or which factors are active, so this should be a quick
   confirmation, not new analysis, but confirm rather than assume.
3. New dependencies needed (Bioconductor: `fgsea`, `clusterProfiler`, `msigdbr`, `org.Hs.eg.db`) —
   already pre-approved in the retrieved plan's own decisions section; install without re-asking.

**Implementation:** execute the retrieved plan's steps as written. Core deliverable: `fgsea`
ranked-by-weight enrichment (primary method — the point-exponential F prior makes all gene weights
≥ 0, a natural one-sided ranking statistic) plus top-N over-representation analysis (confirmatory
cross-check) on Programs 3 and 7, against custom PDAC gene sets (Moffitt/Bailey/DeSurv subtype
programs) and MSigDB collections; figures F1-F5, tables T1-T4, plus three controls (active-vs-inactive
program comparison, external-cohort robustness, factor-sign-direction sanity check per the
`feedback-yfb-km-sign-suppression` gotcha — joint β sign can invert vs. marginal KM direction, label
factors empirically, not by assumed sign).

**Known friction point:** DeSurv gene-list extraction (for the concordance-with-published-comparator
check) pulls from a PDF supplementary appendix; if extraction is unreliable, fall back to the linked
DeSurv GitHub repo for the same gene lists.

**Verify:** re-derived Program 7/3 gene lists match what's already saved in
`presentation/.../assets/active_factor_genes.csv`; enrichment results are directionally sensible
(adverse program enriched for basal-like/aggressive PDAC signatures, protective for
epithelial/classical); full test suite still 311/311 if any shared code was touched.

### Item 2 — Step 6: Executive summary / progress report (K-parsimony follow-up plan)

**Goal:** a research-team-facing progress report, ready to present, synthesizing everything since
the 6/18 meeting.

**Do this last among items 1-3** — it should incorporate Item 1's pathway-enrichment result if
complete by the time you get here.

**Implementation:** synthesize Part A of this document (Phases 1-3, Steps 1-5) into a report
following `CLAUDE.md`'s documentation-audience convention: biostatistician collaborators reading
cold, framed in terms of model/prior/training-set/metric/result — **not** internal phase-label or
step-number narration (use "Phase 1"/"Step 3" etc. only as lookup keys if truly needed, never as the
primary structure of the prose). The honest headline: the recommended configuration's numbers were
corrected (0.636→0.627) and its K=7 choice was independently stress-tested three different ways
(warm-start, deflation-init, joint Bayesian optimization) — all three confirm the same underlying
result (2 factors doing the survival-relevant work, K=7 chosen for one-step reproducibility rather
than because smaller K underperforms). Check `ROADMAP.md`'s flagged to-do item and the
`project_executive_summary_todo` memory entry for any additional context. **Present a draft for
review before finalizing** — this is the one deliverable meant for an audience outside this repo.

### Item 3 — Original-plan Phase 4: Stratified Cox baseline hazard by cohort

**Goal:** let baseline survival risk differ by cohort (`+ strata(study)`), since the existing
`cohort_id` mechanism only absorbs *genomic* platform offsets, not differing baseline risk across
studies.

**Open decision point, not yet resolved — confirm with the user before writing code (the original
plan lists this explicitly as unresolved, do not silently assume the answer):** the original plan's
"Key decisions / open questions" section asks whether the minimal stratified-partial-likelihood
route (no parametric baseline hazard — h₀ still cancels within each stratum) is sufficient, or
whether a fully parametric baseline (e.g. Weibull) is wanted instead. The plan's own recommendation
is the minimal stratified route, but this was never explicitly confirmed with the user — ask before
implementing, don't just proceed on the recommendation.

**What "done" looks like (success criteria, per the original plan — define before implementing,
don't relitigate):**
- **Correctness anchor:** new unit tests for the stratified partial likelihood pass, and the fit is
  *unchanged* when all samples share one stratum (must reduce exactly to the current model — this is
  the test that catches a subtly wrong implementation).
- **Measurable-or-honestly-reported effect:** external mean C-index with cohort strata is ≥ without,
  or the comparison is reported as neutral/negative — either is an acceptable, reportable outcome,
  not a failure condition.

**Implementation (assuming the minimal stratified route is confirmed above):** extend the Breslow
risk-set computation to be stratified by cohort (separate risk sets per study) — the relevant
functions are `calc_cox_taylor_yf()` in `code/fit_cox_on_yf.R` (currently around line 75) and its LB
counterpart `calc_cox_taylor()` in `code/fit_modular.R` (currently around line 84); verify these line
numbers before editing, they will have drifted since the original plan was written. Baseline hazard
still cancels *within* each stratum — this is the standard stratified Cox partial likelihood, no
parametric h₀ is introduced. Higher effort than items 1-2: this touches the CAVI survival derivation
directly and needs new unit tests built from scratch (TDD, not adapting an existing test file).

**Verify:** new stratified-likelihood unit tests; single-stratum reduction test (the correctness
anchor above); external C-index with vs. without cohort strata on the recommended config.
Independent code review, then commit. Branch name: `stratified-cohort-baseline-hazard` (not
previously assigned in either source plan — pick your own if you prefer a different name, this one
just follows the existing branch-per-phase naming convention).

### Item 4 — Original-plan Phase 6: Manuscript track (separate repo)

**Goal:** stand up paper infrastructure (`github.com/andrew-walther/SSBMF-paper`, confirmed
accessible via `gh`) and begin the introduction. Explicitly independent of the model code — can be
done in any order relative to items 1-3, or by a separate/parallel session.

**Current state:** the repo is a from-scratch template (single "Initial commit", placeholder
`.gitkeep` files only). A scaffold was started once already in a scratchpad clone that will **not**
survive into a new session — redo it:
```
gh repo clone andrew-walther/SSBMF-paper
```

**What was learned starting this once already, worth knowing before redoing it:** the template has
an internal inconsistency — `README.md` describes a LaTeX structure (`paper/main.tex`,
`paper/sections/`, `paper/references.bib`, `pdflatex` build instructions) but `.claude/CLAUDE.md`
describes a different, non-existent Quarto/`_targets.R` structure that doesn't correspond to any
real files in the template. **Go with the LaTeX structure** (matches the actual directories and the
README's own build instructions — the Quarto mention is unedited generic boilerplate) and fix
`.claude/CLAUDE.md` to match reality, not the other way around.

**Implementation:**
1. Fill in `README.md` and `.claude/CLAUDE.md` placeholder text (currently generic `[Title]`/`[Brief
   description]` template text).
2. Build `paper/main.tex` + stub section files (`introduction.tex`, `methods.tex`, `results.tex`,
   `discussion.tex` — all placeholder, no real prose yet) + `paper/references.bib`.
3. **Do not** draft real methods/results prose yet (depends on Item 1/Item 3 landing, and possibly a
   final decision on K=7-vs-K=4 framing from Item 2's executive summary). **Do not** copy benchmark
   scripts/figure generators into this repo yet either — that's a bigger reorg decision than this
   scaffolding pass warrants; revisit once there's real content to organize.
4. **Before running any literature research:** ask the user for their existing reference list
   (Zotero export or similar) — the original plan explicitly calls for prompting for this rather than
   guessing what's already covered. Only then run deep literature research (`deep-research` skill) to
   fill gaps: supervised/guided matrix factorization, prognostic gene programs, PDAC subtyping
   (PurIST, DeCAF, basal/classical), Bayesian EBMF. DeSurv is the primary comparator (memory
   `reference-desurv-paper`) — position SBMF as its Bayesian counterpart (same design/data/scoring,
   differentiated by EB-prior shrinkage, posterior uncertainty, error-in-variables E[L²],
   heteroscedastic τ), framed as a *pending* head-to-head, not a superiority claim.
   **Use the original plan's own Appendix** ("Phase 6 starting reference library & gap analysis," in
   `docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md`) as the starting point rather than
   re-deriving a gap list from scratch — it already has the 18-source Zotero library grouped by
   pillar (unsupervised MF/integration, single-cell integration, supervised/semi-supervised NMF, EB/
   Bayesian MF, Bayesian survival, PDAC/cancer-genomics) and a prioritized gap list (HIGH: PDAC
   molecular subtypes — Moffitt/Collisson/Bailey/PurIST/DeCAF; HIGH: gene-program discovery via
   factorization — Brunet NMF metagenes, cNMF/Kotliar; MEDIUM-HIGH: supervised dimension reduction
   for survival — Bair & Tibshirani supervised PCA, PLS-Cox; MEDIUM: penalized/deep-learning survival
   comparators — glmnet, lasso-Cox, DeepSurv; LOW-MEDIUM: VI review, C-index/Harrell). Confirm this
   list is still current with the user (18 sources as of 2026-07-08 — may have grown since) rather
   than assuming it's unchanged.
5. Commit locally in the clone. **Do not push** without asking first.

### Item 5 — Loose ends

- `results/benchmark_sim/outputs/desurv_comparison/` was found stale during this session (reflected
  a superseded pre-Phase-1-correction `K_eff=4` figure) and has already been regenerated — current as
  of this document (D4: mean C=0.627, K=7, K_eff=2). No further action needed unless it goes stale
  again.
- **Open, undiagnosed bug** (unrelated to everything above): YFB with `prior_beta="point_normal"`
  returns C=0.5 for all K in cross-validation. Documented in `ROADMAP.md`. Not blocking anything in
  this document; investigate only if/when it becomes relevant to other work.
- **Action needed, missed in earlier drafts of this document:** the original plan's Global
  Verification section calls for `Rscript tests/test_real_data_loading.R` (77/77, requires
  `PDAC_DATA_ROOT`, local-only) alongside the main suite — this real-data test suite has not been
  explicitly re-confirmed since Steps 1-4's work. Run it once as part of closing out this document,
  not just the main `tests/run_tests.R`.
- **Backlog items from the original plan, still just backlog (no action expected unless the user
  asks):** a post-commit Codex review hook (dev-experience only) and a RAG literature database
  (useful once Phase 6's research volume grows — see Item 4). Neither is being actioned here; noted
  so they aren't silently forgotten.
- **Possibly actionable now, not merely a loose end — flag prominently to the user:** the original
  plan's backlog also lists "schedule a follow-up lab meeting, once Phase 1-2 results are in hand."
  Phase 1-2 (and Phases 3, plus all of Steps 1-4) are now done — this may be worth raising with the
  user as newly-ripe, not something to action unilaterally.
- **Final whole-branch review, called for by the original plan but not yet done as a discrete step:**
  the original plan's Global Verification section asks for "a final whole-branch review against the
  net-benefit gate before any merge is proposed," in addition to the per-phase reviews already
  completed for each individual step. Every phase/step above has been through its own independent
  `superpowers:code-reviewer` pass at merge time, but no single review has looked at the *entire*
  arc (Phases 1-3 + Steps 1-4) together against the original net-benefit gate (all tests passing,
  external C-index not worse than the 0.636 baseline). Worth doing once, ideally right before or
  alongside Item 2's executive summary, since that's the natural point where the whole arc is being
  synthesized anyway.

### Item 6 — Comprehensive HTML project document ("visual executive summary") — DO NOT START until Items 1-5 above are complete

**Goal (as the user described it):** a comprehensive, standalone HTML document — illustrative text,
figures, and workflow diagrams — describing the current state of the project, the key decisions
that constitute the current version of the method, and the method's capabilities. Purpose: get
everyone (the user, collaborators, the research team) to a common point of understanding of what's
been accomplished and what the current shortfalls/open items are. The user's own framing: "this
could illustrate a demo of the method as a visual executive summary."

**Relationship to Item 2 (Step 6's executive summary):** related but distinct. Item 2 is a written,
text-first progress report (likely a `.qmd`/PDF, per this project's existing documentation
conventions). This item is a richer, visual, HTML/artifact-style deliverable — closer to an
interactive one-pager or demo than a report. Item 2's content (the synthesized narrative of Phases
1-3 + Steps 1-6) is a natural input/starting point for this item, not a duplicate effort — once
Item 2 exists, drawing on it should make this item faster, not redundant.

**Trigger condition:** only take this on once Items 1-5 above (pathway enrichment, executive
summary, stratified Cox baseline hazard, manuscript-track scaffold, loose ends) are complete —
i.e., once the full 6/18-meeting-notes plan and its K-parsimony follow-up are both fully wrapped
up. Do not start this speculatively or early; if picking up this document mid-way through Items
1-5, finish those first and only then return to this item.

**What it should likely cover (subject to discussion with the user when this is picked up — not a
finalized spec):**
- Model overview: the YFB architecture (η = (YF)β), how it differs from the LB variant, and why YFB
  is the recommended parameterization.
- Key decisions timeline, condensed from Phases 1-3 and Steps 1-6 (objective normalization fix,
  joint-vs-two-step validation, the K=7/K_eff=2 finding and its triple independent verification via
  warm-start/deflation-init/joint-BO, the K=7-vs-K=4 parsimony judgment call).
- Capability demonstration: external validation results across the 5 held-out PDAC cohorts, the two
  active gene-expression programs (adverse/protective) and their pathway-enrichment biological
  story (from Item 1), ideally with real or illustrative figures rather than tables alone.
- Honest shortfalls/open items: the YFB `point_normal` K-CV bug, Phase 4 (stratified baseline
  hazard) status, anything else still open at the time this is picked up.

**Likely format:** this project has an `artifact-design` skill and an `Artifact` publishing tool
available — worth considering for this deliverable specifically (a self-contained HTML page,
themed, with inline figures) rather than a Quarto-rendered report, given the user's own "visual"
and "demo" framing. Confirm with the user which they'd prefer before building.
