# Handoff prompts for 2026-11-06 (usage-limit window)

> Branch routing updated 2026-11-11: `fix/yfb-supervision` was merged into
> `main` on November 9. Use `main` for new work. Branch instructions below
> preserve the original execution context.

Three independent tasks. Run them one at a time, in this checkout on branch `fix/yfb-supervision`.
Other sessions may be working in this folder, so commit only the files you changed, never stash,
reset or switch branches, and do not merge into main. Follow AGENTS.md: no tool attribution in
code or deliverables, plain literal prose with no semicolons or dashes as punctuation, never
use external cohorts to choose anything, and never commit patient data. A reviewer will check
each task afterwards, so list every file you changed and every number you report with its
source file.

---

## Prompt 1: 11/6 project writeup, phase 1 (prep and data-driven skeleton)

Prepare the updated full project writeup in `paper/writeup/2026-11-06-writeup/`. A previous
session started it and was stopped mid-turn. Run `git log --oneline -10 -- paper/writeup/2026-11-06-writeup`
and `git status` first, and build on what is there. Do not edit the circulated
`paper/writeup/2026-10-16-writeup/`.

1. Read the 10-16 writeup and its `_process/` notes, PROJECT_STATUS.qmd, the newest DECISIONS.md
   entries, ROADMAP.md and `docs/findings/2026-11-0*`.
2. Write a section-by-section plan in the folder README. For every section, mark what carries
   over, what changes and what is new, and name the evidence file for every number.
3. Set up `writeup.qmd` so every number comes from R chunks that read saved CSVs. No hand-typed
   numbers (GitHub issue #6 item 8). `docs/progress_book/chapters/2026-11-06.qmd` shows the
   pattern.
4. The writeup must address issue #6 (`gh issue view 6`):
   - item 2: use single-log numbers;
   - item 4: state the pre-registered October 15 supervision null result plainly (ΔC = −0.003,
     −0.029 to 0.024; DECISIONS.md around line 1031);
   - item 5: lead with the rank-7 EBMF + Cox comparison, not the 40-factor lasso comparator;
   - item 6: list every use of the external cohorts (K = 7, preprocessing variant, α_F,
     N_frozen, comparator configuration) and disclose the 47 shared PACA-AU donors.
5. Draft only the sections that do not depend on the prelim Chapter 4 figures. Chapter 4 is
   being revised separately. Leave clearly marked placeholders for its figures.
6. Commit and push. Stop there. Do not render the final PDF or run the prose review yet.

---

## Prompt 2: GitHub issue #6 item 5 leftover (risk orientation in other runners)

Issue #6 item 5 flagged `max(c_raw, 1 - c_raw)`, which orients a risk score with the
validation outcome. Commit b85d390 fixed the two named scripts. About ten other runners still
use the pattern (for example `results/benchmark_sim/run_desurv_comparison.R:164`,
`run_joint_bo.R:107`, `run_cohort_lmm_benchmark.R:126`).

1. Find every remaining occurrence: `grep -rn "1 - c" results/ code/ | grep -i max`, plus variants.
2. For each one, decide whether it feeds a number reported in a current document (writeup,
   progress book, prelim, DECISIONS.md, findings).
   - If it does, replace it with the training-fixed orientation `frozen_reverse_cindex()`
     (`code/concordance_ci.R`).
   - If it does not, still replace it, unless the script is in a dated archive or legacy
     folder. Leave those unchanged and list them.
3. Write a short table (script, line, reported anywhere?, action) in
   `docs/findings/2026-11-06-orientation-audit.md` and add its README row.
4. If rescoring is possible from saved outputs without patient data, rescore, report whether
   any number changes, and save the rescore CSV. Do not refit models.
5. Run `Rscript tests/run_tests.R` if you touched anything under `code/`. Commit and push.

---

## Prompt 3: GitHub issue #6 item 6 (external-cohort uses and donor-level resampling)

1. Write `docs/findings/2026-11-06-external-cohort-uses.md`, an exact list of every decision
   that looked at external-cohort results: K = 7, the preprocessing variant, α_F, N_frozen, the
   comparator configuration and anything else. For each one give the source (DECISIONS.md line,
   findings file or commit). Name at least one external cohort that played no part in any
   choice, or say plainly that none qualifies.
2. Add donor-level resampling for the pooled external comparison in
   `results/benchmark_sim/score_yfb_supervised_external.R` (S-12). PACA_AU_array and
   PACA_AU_seq share 47 donors. The current bootstrap resamples records. Resample donors so
   that a shared donor's array and RNA-seq records move together. Find the donor identifiers in
   the loaders (`load_pdac_raw()`). This needs the local patient data (`PDAC_DATA_ROOT`), and
   if it is unavailable, stop and say so. Write only aggregate outputs:
   `external_paired_diff_donor.csv` next to the S-12 outputs. Report the new interval beside
   the record-level one.
3. Do not choose anything from these results. Add a findings README row. Commit and push.

---

## Prompt 4: finish the prelim Chapter 4 revamp

Continue the prelim Chapter 4 (SSBMF) revamp from where it paused. Read, in order:
- `docs/plans/2026-11-05-chapter4-revamp-prompt.md`;
- `docs/plans/2026-11-05-chapter4-revamp-plan.md`, especially section 0 (Andrew's notes on title,
  format, flow and figures) and section 9 (the pause point and what is next);
- `git log --oneline -10 -- paper/prelim docs/plans/2026-11-05-chapter4-revamp-plan.md`.

The untracked `paper/prelim/make_figures_supervised_programs.R` and `make_figures_supervised_sim.R`
were left unverified when the work paused. Check them before using them.

The evidence has changed since the draft:
- S-15 (`docs/findings/2026-11-06-supervised-real-refit-findings.md`) kept w_S = 1 on real data.
  The CV curve is flat, and survival's share of the loading precision is below 10⁻⁴ even at
  w_S = 50.
- KL-2 (`docs/findings/2026-11-06-flat-beta-pruning-rerun-findings.md`) shows pruning removes the
  planted program, so the working fit is K = 7 unpruned.
- The donor-level interval for S-12 is −0.003 to +0.017
  (`docs/findings/2026-11-06-external-cohort-uses.md`).
- The 11/6 progress-book chapter (`docs/progress_book/chapters/2026-11-06.qmd`) summarizes all of
  this.

Andrew's decisions (2026-11-06):
- The plan is pre-approved.
- Per-patient external risk scores may be saved under `data/` (git-ignored), never committed.
- The prelim abstract may be updated, but only the supervision sentence or sentences, framed as
  progress on a stated limitation.
- Issue #6 items 4 and 5 apply. State the pre-registered October 15 null result (ΔC = −0.003,
  −0.029 to 0.024), and lead with the rank-7 EBMF + Cox comparison.

Finish the checklist rows, the figures and tables (every number traced to a file), the prelim
rebuild with `test_prelim_assembly.py` passing (report page counts), the Round 3 review and the
ROADMAP and DECISIONS updates in both repos. Then answer "Is this completely done?" with the
completeness list from the prompt file. A reviewer will do the final human-writer pass, so list
every changed file.

---

## Prompt 5: finish the 11/6 writeup (after Prompt 4)

Finish `paper/writeup/2026-11-06-writeup/` from the plan in its README (phase 1 was done in
`53b8aaa` and `f5b91e0`). Fill the Chapter 4 figure placeholders with the figures Prompt 4
produced, keep every number computed from saved CSVs, and add S-13, S-14, S-15, KL-2 and the
false-active contribution analysis from `docs/findings/2026-11-0*`. Render HTML and PDF, and
update `paper/writeup/README.md`. List every changed file and every number with its source for
the review.

---

## Prompt 6: real-data comparators (GitHub issue #6 item 7)

Register before running. Write `docs/plans/2026-11-06-real-data-comparators-plan.md` (add a README
row and an attempts-log entry in `docs/plans/2026-10-16-round2-attempts-log.md`), then implement
and run.

- **Fixed reference fits.** Do not refit them:
  - the S-11 supervised fit (`w_S = 1`, D4 genes,
    `results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/fits/size_D4_N2218_w1_f0.rds`);
  - its two-step counterpart (the `w0` fit followed by Cox on the program scores).

  Score both as in `results/benchmark_sim/score_yfb_supervised_external.R`.
- **Comparators, each trained on the TCGA + CPTAC training patients only:**
  1. PurIST as a prognostic score. PurIST calls are already used for subtype concordance:
     `grep -rn -i purist code/ results/`.
  2. DeSurv head to head, trained on the same training patients. Use the existing setup in
     `results/benchmark_sim/run_desurv_comparison.R` and `code/preprocess_desurv.R`, with
     DeSurv's own preprocessing and settings as published, chosen inside training only.
  3. Elastic-net Cox on the D4 genes (glmnet, α and λ chosen by training CV).
  4. A clinical-covariate Cox model, using whatever stage, grade, age and resection-margin
     fields exist in both training and external cohorts. Report which fields were available.
- **External evaluation.** Per-cohort C with 95% bootstrap intervals, the cohort-stratified
  paired difference against the supervised fit with donor-level resampling (the S-12 donor
  bootstrap in `code/concordance_ci.R`), pooled HR per SD of each risk score stratified by
  cohort, and calibration slopes. The risk direction is fixed on training data
  (`frozen_reverse_cindex()`).
- **Rules.** Choose nothing from external results. Patient-level data and scores stay out of
  Git, and only aggregate CSVs are committed. Write a findings document and its README row. List
  every number with its source file.

---

## Prompt 7: DeSurv-style analyses, risk groups and program classification

Register before running (`docs/plans/2026-11-06-desurv-style-analyses-plan.md`, README row,
attempts-log entry). Use the fixed S-11 supervised fit and the two-step fit. Choose nothing from
external results.

1. **Risk groups from one training-CV cutpoint.** Standardize the risk score on training data.
   Choose one cutpoint by training cross-validation, using the out-of-fold scores (state the
   criterion in the plan, for example the maximally selected log-rank statistic). Apply that
   fixed cutpoint to every external cohort. Show Kaplan–Meier curves with numbers at risk, and
   report the high vs low HR per cohort and pooled stratified by cohort. These replace the
   tertile groups.
2. **Per-program forest plot.** Give the HR per SD of each active program's score in each
   external cohort and pooled, for both the supervised fit and the two-step fit.
3. **Adjusted HRs.** Give the risk score's HR adjusted for PurIST subtype and for DeCAF
   (`grep -rn -i decaf`). If DeCAF calls are unavailable, say so.
4. **High-risk-group enrichment.** Show the proportions of Basal-like (PurIST) and proCAF
   (DeCAF) tumours in the high vs low groups, with Fisher tests.
5. **Unsupervised model at matched and larger rank.** EBMF at K = 7 and at a larger K (for
   example 15 and 40) followed by Cox, scored the same way.
6. **Program classification rule.** For each active program of the supervised fit:
   - the standardized |β| ratio to the largest program;
   - its share of the risk-score variance;
   - its stability across the five S-11 CV folds (active in at least 4 of 5, matched by loading
     correlation).

   The proposed rule: a main prognostic program has a ratio ≥ 0.25, a share ≥ 5% and is
   active in ≥ 4 of 5 folds. Check the rule on the S-13 and S-14 simulation checkpoints
   (`results/benchmark_sim/summarize_false_active_contribution.R` shows how to load them). Report
   how often it keeps the planted program and drops false actives. Replace the 1% PVE cut for
   "survival-only direction" with the survival share of the loading precision, and report
   both labels.
7. **Simulation items.** Add a mixed-signal scenario (one strong and one near-threshold
   prognostic program plus neutral programs) and top-gene precision (the overlap of the top 50
   genes of the matched program with the planted genes), using the S-14 generator and runner
   pattern (`results/benchmark_sim/run_yfb_near_threshold.R`). Use exact times in every Cox
   fit on simulated data (`coxph.control(timefix = FALSE)`).

Write a findings document with figures, add its README row, and list every number with its
source file. Patient-level outputs stay out of Git.
