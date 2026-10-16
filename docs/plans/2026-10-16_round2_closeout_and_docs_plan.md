# Round 2 close-out and documentation restructure: plan

> **For agentic workers:** Use superpowers:executing-plans (or subagent-driven-development) task by task. Steps use checkboxes. Andrew approves each review gate (marked **Gate**) before the next task starts.

**Goal:** Close `round2/supervised` without merging it, keep what it learned in one short report, port only vetted code to `main`, and restructure the working documents so the current method and results are easy to see.

**Architecture:** Round2 is frozen under a tag. `main` receives a findings report, the approved code and records, and slimmer canonical documents. Old long-form versions of `PROJECT_STATUS` and `ROADMAP` are kept once in `docs/archive/`.

**Tech stack:** git, R test suite (`Rscript tests/run_tests.R`), Quarto (RStudio's copy at `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`, with `QUARTO_PANDOC=/opt/homebrew/bin/pandoc`).

**Spec:** Andrew's decisions of 2026-10-16 (this conversation's summary is reproduced in the "Decisions" section below).

## Decisions already made (2026-10-16)

- Round2 is not merged. Finish the writeup edits in progress there, then freeze it.
- Supervision of F, rank selection by pruning, and the methylation and joint feature curves close as documented, unresolved tracks. They may be revisited on new branches off `main`; the untried remedies are a candidate list, not a mandatory sequence.
- The DeSurv head-to-head stays deferred; the optional robustness re-run is dropped.
- Document roles going forward:

| Document | Role |
|----------|------|
| `README.md` | Public front page: YFB method, headline results, how to run |
| `PROJECT_STATUS.qmd` | Current best state, a few pages: model, configuration, results, limitations, open questions |
| `ROADMAP.md` | Full project record. Part 1: history (milestones, findings, approaches tried and dropped, one line each with the `DECISIONS.md` date). Part 2: open work in priority order |
| `DECISIONS.md` | Append-only dated decision log |
| `CLAUDE.md` | About 60 lines: working rules, guardrails, pointers |
| Progress book | New pinned first chapter "Current method and findings", plus one chapter per meeting |

- The LB model (η = Lβ) appears in frequently used documents as one line: an earlier formulation superseded by YFB, because YFB scores new patients by projection without refitting L.

## Global constraints

- Never commit patient data or fitted `.rds` objects derived from it.
- Do not modify `results/benchmark_snapshot_2026-10-14/` contents or the tag `benchmark-interim-writeup-2026-10-14`.
- Documentation changes may be pushed straight to `main`. Code ports go through a branch with the full test suite passing first.
- Follow the repo's writing rules: biostatistician audience, no internal labels (Task2-7, Round 2, Cluster A/B) in prose, no assistant-tool names in deliverables, no mannered prose.

## Review focus

1. A frozen-tag reference must resolve: after tagging, `git show <tag>:docs/plans/Round2_attempts_log.md` works from a fresh clone.
2. A ported file must not depend on an unported file: source each ported file in a clean R session.
3. Shrinking documents must not lose a live pointer: grep `main` for links into removed sections and for files that no longer exist.
4. The pinned book chapter and `PROJECT_STATUS.qmd` must not drift apart: one source, included in the other.
5. The headline C = 0.627 must use a training-fixed risk direction (see Task 4, step 4).

---

### Task 1: Freeze round2

**Files:** none (git only)

- [ ] **Step 1:** Confirm the writeup session on `round2/supervised` has committed: `git status --short` shows no tracked changes from it, and Andrew confirms the writeup revision is complete.
- [ ] **Step 2:** Tag the tip: `git tag -a archive/round2-supervised-2026-10 -m "Round2 supervised work, frozen at close-out"` and `git push origin archive/round2-supervised-2026-10`.
- [ ] **Step 3:** Verify from a fresh clone: `git clone -q --branch archive/round2-supervised-2026-10 <repo> /tmp/x && ls /tmp/x/docs/plans/Round2_attempts_log.md`.

Keep the branch itself until Task 3 is merged; delete it afterwards only on Andrew's word.

### Task 2: Findings report

**Files:** Create `docs/reports/Round2_findings_2026-10.md` (on `main`).

One to two pages, for a biostatistician reading cold. Sections:

1. **Scope:** what the work block set out to test (survival supervision of F, rank selection, feature-count sensitivity, multimodal value).
2. **Results table:** one row per item, with status and a one-sentence result. Rows: simulation comparators; expression feature-count curve; methylation descriptive study; rank selection by pruning; survival supervision of F; methylation and joint feature curves; multimodal value; DeSurv head-to-head. Take the numbers from `docs/plans/Codex_to_Claude_Round2.md` on the tag. Include the external non-inferiority failure (ΔC = −0.006, 95% CI −0.032 to 0.021, required lower limit > −0.02) and the expression curve (CV C 0.608 / 0.577 / 0.609 / 0.607 at 2,064 / 5,001 / 10,000 / 17,900 genes).
3. **Remedies tried and untried for supervision:** a short list taken from the remedy table in the same file.
4. **Numerical findings:** log-space Cox arithmetic repairs the 25k-CpG underflow, but that fit still does not converge.
5. **Where the full record is:** the tag name and the paths of the attempt log and status files on it.

- [ ] **Step 1:** Draft and show Andrew. **Gate.**
- [ ] **Step 2:** Commit to `main` and push.

### Task 3: Port vetted code and records

**Files:** a keep/drop table (in chat or appended to the Task 2 report), then a branch `cleanup/round2-port` off `main`.

- [ ] **Step 1:** Build the table for every file round2 adds or modifies in `code/`, plus `results/benchmark_snapshot_2026-10-14/`, `docs/progress_book/chapters/2026-10-15.qmd` and `2026-10-16.qmd`, and the revised `paper/writeup/`. Columns: file, purpose, needed by the current best state or the next experiment (yes/no), recommendation. Starting recommendations:
  - **Keep:** `cox_working_logspace.R` (numerical correctness, opt-in); `feature_sensitivity.R` (nested gene-set manifests and program agreement, reused by the feature-selection experiment); the benchmark snapshot folder (README and manifest only); both book chapters; the writeup, if Andrew approves it as replacing the interim.
  - **Drop (stays on the tag):** `loading_balance_audit.R`, `multimodal_moment_audit.R`, `loading_prior_fixed_q.R`, `profile_cox_intercept.R`, `cpg_feature_sensitivity.R`, `multimodal_feature_sensitivity.R`, `multimodal_program_agreement.R`, all round2 diagnostic runners, SLURM wrappers and output CSVs.
  - **Inspect:** modifications to `fit_cox_on_yf.R`, `update_beta.R`, `update_F_surv_YFB.R` and the multimodal files. Keep a modification only if it is a bug fix or a default-path improvement, with a test showing the default output is unchanged or corrected.

  **Gate.**
- [ ] **Step 2:** On `cleanup/round2-port`, check out each kept path from the tag (`git checkout archive/round2-supervised-2026-10 -- <path>`) together with its tests.
- [ ] **Step 3:** Source each kept code file in a clean R session; then run `Rscript tests/run_tests.R`. Expected: 0 failures. The pathway tests need the ignored `results/benchmark_sim/outputs/desurv_comparison/*.rds` locally.
- [ ] **Step 4:** Commit, merge to `main` (fast-forward or a PR, Andrew's choice), push.

### Task 4: Restructure the canonical documents

**Files:**
- Create: `docs/archive/PROJECT_STATUS_through_2026-10-16.qmd` (+ `.pdf`), `docs/archive/ROADMAP_through_2026-10-16.md`, `docs/archive/README.md` (two lines: what these are and the date)
- Modify: `PROJECT_STATUS.qmd`, `ROADMAP.md`, `CLAUDE.md`, `README.md`, `code/README.md`, `code/SupervisedMF_Context.md`, `results/README.md`

- [ ] **Step 1:** Copy the current `PROJECT_STATUS.qmd`, its PDF and `ROADMAP.md` into `docs/archive/` under the names above. If round2's versions are newer (after Task 1), archive those instead, from the tag.
- [ ] **Step 2:** Rewrite `PROJECT_STATUS.qmd` as the current-best-state document with exactly these sections: *Model* (YFB, priors, α_F = 0, K = 7 with factor classification); *Data and preprocessing* (cohorts, per-cohort z-standardization, DeSurv-style gene selection); *Results* (external C per cohort and mean, KM figure, programs and biology in brief); *Known limitations* (F not survival-informed; pruning not adopted; multimodal adds no demonstrated predictive gain); *Open questions* (pointing to `ROADMAP.md` Part 2). Target: 3–5 rendered pages. **Gate.**
- [ ] **Step 3:** Restructure `ROADMAP.md` into Part 1 (history in date order, one line per milestone or dropped approach, each citing its `DECISIONS.md` date) and Part 2 (open work, prioritized; the first item is the feature-selection experiment in `docs/plans/2026-10-16_feature_selection_experiment_plan.md`). The history of the old `PROJECT_STATUS.qmd` moves into Part 1 in condensed form. **Gate.**
- [ ] **Step 4:** Check the headline: recompute the five external C values for the frozen D4 fit with the training-fixed direction used by `paper/writeup/figures/make_writeup_figures.R`, and confirm they match the writeup figure (0.634, 0.549, 0.648, 0.657, 0.645; mean 0.627). If any differ, stop and report before editing any document. `oriented_cindex()` (max(C, 1 − C)) chooses the direction on the evaluation data and must not be used for reported results.
- [ ] **Step 5:** Rewrite `CLAUDE.md` to about 60 lines: pointers to the documents in the table above; the existing writing rules; and these guardrails:
  - one experiment per branch off `main`, started from a short plan in `docs/plans/`;
  - detailed outputs git-ignored, with only compact summary tables committed;
  - each branch ends with a findings entry in `ROADMAP.md` Part 1 and an explicit merge-or-discard decision; a discarded branch is tagged `archive/<name>` and deleted;
  - no run-instruction paragraphs in `CLAUDE.md` or `README.md`;
  - reported C uses a training-fixed risk direction.
- [ ] **Step 6:** Cut LB to one line in `README.md` (model table row, benchmark command), `code/README.md`, `code/SupervisedMF_Context.md` (LB update table) and `results/README.md`. The LB code stays; the YFB fitter reuses its update functions.
- [ ] **Step 7:** Render `PROJECT_STATUS.qmd` to PDF; grep `main` for links to removed sections or missing files; commit and push.

### Task 5: Pinned "current state" chapter in the progress book

**Files:**
- Create: `docs/progress_book/current_state.qmd`
- Modify: `docs/progress_book/_quarto.yml` (insert `current_state.qmd` directly after `index.qmd` in `book.chapters`)

- [ ] **Step 1:** `current_state.qmd` holds a title ("Current method and findings"), one line saying it is kept current rather than dated, and the body of `PROJECT_STATUS.qmd` through a Quarto include, so there is one source. Try `{{< include ../../PROJECT_STATUS_body.qmd >}}` with the body split out of `PROJECT_STATUS.qmd` into `PROJECT_STATUS_body.qmd`, and include the same fragment from `PROJECT_STATUS.qmd`.
- [ ] **Step 2:** Render the book and check that the chapter appears second in the sidebar, figures resolve, and `PROJECT_STATUS.pdf` still renders. If an include from outside the book folder fails, keep the fragment inside `docs/progress_book/` and include it from `PROJECT_STATUS.qmd` instead.
- [ ] **Step 3:** Publish with `docs/progress_book/publish_gh_pages.sh` after Andrew reviews the rendered chapter. **Gate.**
