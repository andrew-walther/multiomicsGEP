# Codex execution prompt: 10/30 realignment plan

Written 2026-10-30. Paste this into Codex from the repository root.

---

You are continuing work on multiomicsGEP, an R project fitting a supervised Bayesian
nonnegative matrix factorization jointly with a Cox model (η = (YF)β). Read these first, in order:

1. `CLAUDE.md` (project rules; `AGENTS.md` is the same file).
2. `docs/plans/2026-10-30-realignment-plan.md`: the plan you are executing. Its last section
   records which Phase C items are already done; start from the first item not marked done.
3. `docs/plans/2026-10-30-Unified-Architecture-Alignment.md`: the model specification for Phase S.
4. `docs/findings/2026-10-30-current-method-a1.md`: why a single survival weight did not work.
5. `docs/progress_book/chapters/2026-10-26.qmd`: the last chapter, for reference only.

## Order of work

1. Finish Phase C (cleanup) if any item is open. It merges to `main` once tests pass.
2. Record the Phase P preprocessing decisions in `DECISIONS.md` and in §5 of
   `docs/findings/2026-10-26-preprocessing-audit.md`.
3. Phase S on a new branch `fix/yfb-supervision` from `main`:
   - S0, the precision-share audit, including how DeSurv balances its loss (the DeSurv package
     is in the isolated library the tests use);
   - S1, the new opt-in fitter `code/fit_yfb_supervised.R` with its tests;
   - S2, the simulation decision test, with arms and Gate S registered in
     `docs/plans/2026-10-16-round2-attempts-log.md` before you run anything;
   - S3, real-data training CV, only if Gate S passes.
4. Phase D throughout. Write the new chapter `docs/progress_book/chapters/2026-10-30.qmd`
   and add it to `_quarto.yml`. Keep the issue tracker as a CSV
   (`docs/progress_book/data/2026-10-30-issue-tracker.csv`) and render it with an R chunk. All
   numbers in the chapter come from saved CSV/RDS files through R chunks, not typed by hand.
   Also update `DECISIONS.md`, `ROADMAP.md`, `PROJECT_STATUS.qmd` and the `docs/` READMEs.

5. Close the loop for review. When you finish, or when you stop on a stop condition, write
   `docs/reviews/2026-MM-DD-codex-execution-report.md` (the date you finish), add its row to
   `docs/reviews/README.md`, and fill in the "Review handoff" section at the end of the plan:
   - for each plan item: done, partly done, or not done, with the commit and the evidence (test
     output, result file, or chapter section);
   - every deviation from the plan and why;
   - every failure, non-converged fit or skipped record, with numbers;
   - the decisions that Andrew must make next.
   Commit and push. Do not merge `fix/yfb-supervision` into `main`; a review session reads it first.

## Chapter rules (the chapter is the main deliverable for the next advisor meeting)

- Open with a short list of takeaways that can be shared as is: what was fixed, what was learned,
  and what is still blocked and why.
- Name no reviewer or advisor. Describe issues, solutions and progress.
- No branch names, file paths or commit hashes in the chapter text.
- Plain statistical prose with no metaphor; see the CLAUDE.md prose rules.

## Project rules

- Run `Rscript tests/run_tests.R` after any change to fitting or update code. Do not report
  work as done while tests fail. If a failure comes from git-ignored inputs, say so and name the file.
- Never touch `results/benchmark_snapshot_2026-10-14/` or the benchmark tag. Write new outputs
  alongside the old ones.
- Do not score external validation cohorts to choose anything. Keep at least one cohort unused
  until the final claim.
- Commit messages: a short imperative title and a body that says what changed and why. No
  co-author or tool-attribution lines. Do not name assistant tools in code or deliverables.
- Fail loudly: report non-converged fits, NaN predictors, all-zero factors and skipped records.
- Run long jobs in the background with a socket cluster, not `mclapply`. Forking segfaulted on
  macOS in the A1 runs.

## Stop conditions

Stop, write a status note in the chapter and commit it, without working around the problem, if:
- tests fail and the cause is unclear;
- the benchmark rerun moves external C by more than 0.02 without explanation;
- Gate S fails. Before you stop, diagnose: per-sweep traces of survival share, β, ‖L‖ and the
  objective; the sweep where instability starts; which piece fails first. Then write a standalone
  problem statement that Andrew can take to another reviewer: what was tried, the numbers, and the
  exact point where it fails. DeSurv solved this balance problem for NMF + Cox, so compare
  against its mechanism explicitly.
- a step would use external cohorts to inform a model choice.

## After Codex: review prompt for Claude

Andrew pastes this into a Claude session at the repository root once Codex has finished:

> Review the Codex execution of `docs/plans/2026-10-30-realignment-plan.md`. Read the plan's
> "Phase C progress" and "Review handoff" sections and the report under `docs/reviews/` that the
> handoff names. Check the work against the plan, not against the report: run
> `Rscript tests/run_tests.R`; read the diff of `fix/yfb-supervision` against `main`; check that
> `code/fit_yfb_supervised.R` matches §4 of
> `docs/plans/2026-10-30-Unified-Architecture-Alignment.md` and the plan's S1 list (bounded
> point-normal β prior, damping, hybrid start, KL from ebnm output, no deletion refits); confirm
> that the S2 arms and Gate S were registered before the run and that the gate verdict follows from
> the saved results; and check that the 10/30 chapter's numbers come from saved files and that it
> follows the chapter rules. Report findings ranked by severity, say whether Gate S is
> passed, failed or not reached, and recommend merge, revise or stop. Then update the plan's
> "Review handoff" section with the verdict, and wait for Andrew before merging anything.
