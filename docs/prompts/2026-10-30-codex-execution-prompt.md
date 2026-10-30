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

5. Close the loop for review. When you finish, or when you run out of diagnosed ideas, write
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

## When something goes wrong: diagnose and solve, don't stop

When you hit a problem, diagnose it and try to solve it; do not just stop. That covers failing
tests, non-convergence, an unexpected benchmark change, a failed Gate S, or a specification piece
that does not work as written. For each problem:

1. Find the root cause before changing anything. Use per-sweep traces, small reproductions and
   comparisons with DeSurv or the rebuilt package. Do not try fixes at random.
2. Make the fix, add a test that fails without it, and rerun what it affects.
3. Record every attempt, failed ones included, in `docs/plans/2026-10-16-round2-attempts-log.md`
   and in the chapter: what was wrong, what you changed and what happened.

Gate S needs one rule to stay honest. If the registered run fails and a fix follows, register
the revised arms as a new, numbered attempt before rerunning them. Report the original failure
next to the new result, and never edit the original criterion after seeing results.

Keep going until the plan is done or you have no remaining idea backed by a diagnosis. In that
case, write a standalone problem statement that Andrew can take to another reviewer: what was
tried, the numbers, the exact point of failure, and your best hypotheses. Then finish the review
handoff.

Hard limits that hold whatever the problem:
- Never use the external validation cohorts to choose anything.
- Never modify `results/benchmark_snapshot_2026-10-14/` or the benchmark tag.
- Never report work as done while tests fail.

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
