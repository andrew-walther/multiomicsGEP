# Prompt: finish the supervision work on `fix/yfb-supervision` (attempts S-9 and S-10)

Written 2026-10-30, revised 2026-11-03. Andrew approved this work. Numbering: S-8 is the
targeted mechanism screen already in the attempts log. This prompt registers S-9 (the existing
survival-weighted fitter at 2,000 genes) and, only if S-9 fails, S-10 (downweighted genomics).
It can be run in Codex or in Claude Code. Paste everything below the line from the repository
root, on branch `fix/yfb-supervision`.

---

Continue on `fix/yfb-supervision`. Read these first:

- `docs/findings/2026-10-15-m6-m7-summary.md`, section B, and
  `docs/plans/2026-10-16-round2-return-note.md`: the October 15 survival-weighted results.
- `docs/findings/2026-10-30-current-method-a1.md`: the survival-weight grid.
- `docs/findings/2026-10-30-yfb-targeted-diagnostics-findings.md` and
  `docs/findings/2026-10-30-yfb-supervision-problem-statement.md`: what S-7 and S-8 found.
- `docs/plans/2026-10-30-realignment-plan.md`: Gate S and the "Problems during execution" rules.
- `docs/prompts/2026-10-30-codex-execution-prompt.md`: project rules and chapter rules, which
  still apply.
- `docs/plans/2026-11-01-gemini-multiomicsGEP-diagnostic-and-tdd-guide.md`: a Gemini guide (unchecked
  advice). Its Issue 1 proposes the downweighted genomics likelihood used in Part B.
- `DECISIONS.md`, the July entry near "Two candidate directions were tried and empirically rejected"
  (around line 4029). Shrinking the genomics precision in the latent-L model drove L, F and β to
  exactly zero within 2–3 iterations.

## What the October 15 survival-weighted fit already showed

(The YFB fitter was renamed from `fit_cox_on_yf()` to `fit_ssbmf_yfb()` on 2026-11-04. Earlier
attempt-log entries and findings use the old name.)

The existing fitter `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` puts a survival term into the
gene-by-gene F update, weighted by w_S. On October 15 it showed three things:

1. Survival did act on the F update, although its share of the F precision was small compared with
   the genomics share.
2. In simulation it recovered a low-variance prognostic program that the unsupervised methods
   missed: 17 of 20 datasets at noise 0.25, against 0 of 20 for EBMF → Cox, EBMF → lasso-Cox and
   PCA → Cox, with mean false-active 1.1.
3. On the five external cohorts its C was slightly lower than the two-step method's, within the
   margin of error: ΔC = −0.003 (95% interval −0.029 to 0.024).

All of that simulation evidence is at p = 50 genes (n = 400). The A1 grid found that on the real
training data, at about 2,000 genes, survival's share of the F precision at w_S = 1 is about
10⁻⁶. The S-7 grid tested the new fitter at p = 2,000, but it did not include the survival-weighted
fit as an arm, although the realignment plan listed it as a baseline. So the most direct open
question is whether the fit that worked at 50 genes still recovers the program at 2,000.

## Order of work

1. **Part A (attempt S-9): the survival-weighted fit on the S-7 grid.** Do this first. It uses
   existing code and answers the open question above.
2. **Part B (attempt S-10): downweighted genomics.** Run it only if S-9 fails Gate S, or if S-9
   passes but leaves a problem that the downweighted objective is derived to fix, such as β
   calibration. Say which reason applies.
3. **Part C: the loose ends on this branch.** Do these whatever the outcome of A and B.

## Part A (attempt S-9): the survival-weighted fit at 2,000 genes

1. **Check for prior runs.** First check whether any registered attempt already ran
   `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` on a p = 2,000 simulation. If one did, report it
   and say whether it settles the question.
2. **Prepare the data for this fitter.** It expects signed, per-cohort z-scored input, while the
   S-7 grid prepares nonnegative input for the new fitter. Generate the S-7 datasets (same seeds,
   p ∈ {50, 2,000}, planted variance 0.25%, 1% and 5%, 10 datasets per cell, the same paired
   unsupervised reference) and give each fitter the preprocessing it was written for. Record how
   the planted program's variance share changes under each preprocessing.
3. **Include a replication cell.** Add the October 15 calibrated bed (p = 50, n = 400, β₃ = 3.3)
   as a replication check. It should reproduce about 17 of 20.
4. **Arms: a range of weights, reported as a curve.** The genomics likelihood has n·p terms and
   the Cox likelihood has about one term per event, so at fixed w_S survival's share of the F
   precision falls roughly as 1/p. Run w_S ∈ {1, 10, 100, 1,000, 10,000} and two scaled arms:
   - w_S = p/50, which keeps the October 15 balance (w_S = 1 at p = 50) as p grows;
   - w_S chosen so the median survival share of the F precision equals 0.2 (computed from the
     paired unsupervised fit before the supervised run, not tuned afterwards).
   For every arm report recovery, false-active count, survival share, β calibration, convergence
   and validation C, and plot each against w_S. A1 found that weights above 100 broke convergence
   with free per-gene noise, and S-8 found that fixed noise prevents the blow-up but that a very
   large Cox weight makes β overconfident. Record non-convergence at large w_S as a result, and
   say whether it is the τ mechanism, the β mechanism or something new.
5. **How the weight would be chosen in practice.** Gate S is read at the best registered arm, which
   uses the planted truth and is allowed only in simulation. Also report the arm that
   training-data cross-validated C would pick in each cell, and whether it is the same one. On real
   data the weight may only be chosen by training cross-validation. A prior on w_S is not a
   substitute: w_S tempers a likelihood, so its posterior is not well defined without an extra
   calibration rule (for example SafeBayes or the Holmes–Walker power-likelihood calibration).
   Name that as a follow-up in the findings, do not implement it here.
6. **Gate S is unchanged.** At p = 2,000, a recovery gain of at least 3/10 over the unsupervised
   fit in at least two variance cells, mean false-active at most 1, and at least 90% of fits
   converged.
7. **Measure along the way:** survival's share of the F precision per fit, and how it changes from
   p = 50 to p = 2,000. If recovery drops at 2,000 genes, show whether the falling share explains it.
8. **Register before running:** an S-9 entry in `docs/plans/2026-10-16-round2-attempts-log.md`
   and `docs/plans/2026-11-03-survival-weight-p2000-plan.md`.
9. **If S-9 passes:**
   - run the realignment plan's step S3 with the passing w_S: real-data training cross-validation
     on the single-log data, training data only, no external cohort;
   - then run the null-data calibration check from Part B, step 3.2 at that w_S.

## How to work: root causes and tests first

Work from the model derivation and the update algorithm, not from symptoms. For each problem,
name the mechanism that causes it and fix that mechanism. A change that only hides a symptom
does not count as a fix, for example clipping, extra damping without a derived reason, or a
preprocessing step that removes the visible failure. The per-cohort z-standardization is the
cautionary example: it stabilized the predictor but put negative values into a nonnegative
factorization and hid the scale problem.

Use test-driven development. From the derivation, write tests that encode its mathematical
properties before writing the implementation: the objective does not decrease, the KL terms are
exact, rescaling invariances hold, and the limiting cases (λ_G = 1; no survival signal) reduce
to the known fits. Watch each test fail before the implementation exists, then make it pass.

## Part B (attempt S-10): downweighted genomics

Run only under the condition in "Order of work". The steps below are unchanged from the 10/30 version apart from their attempt numbers.

## Part B, step 1. Derivation (write before any code)

Write `derivations/yfb_downweighted_genomics/2026-10-30-downweighted-genomics-derivation.qmd`
covering:

1. **One objective.**
   L(q) = λ_G · E_q[log p(Y | L, F, τ)] + E_q[log PL(η)] − KL(L) − KL(F) − KL(β), with
   λ_G ∈ (0, 1] and η the normalized centered projection that the current fitter uses. State this
   as a generalized-Bayes (power-likelihood) objective, and say what the β posterior means under it.
2. **Every coordinate update for this objective.** λ_G multiplies the genomics precision and signal
   in the L and F updates. The Cox terms enter at weight 1 everywhere, β included. Include the
   derivative of the loading normalizer, ∂Z_k/∂f_jk = Yc_j/r_k − Z_k f_jk/r_k². Either derive the
   exact coupled update, or state a minorize–maximize or trust-region step that provably does not
   decrease L(q), with the Cox Hessian bound from the guard (D_i = half the risk-set membership).
   The guard showed that cutting the F step alone can stall, so the update has to move F and β
   consistently.
3. **The KL question.** With λ_G small, the L/F priors weigh relatively more against the
   downweighted likelihood, and this is the likely mechanism of the July collapse. Decide and
   justify one option, then register it:
   - (a) keep the KLs at weight 1;
   - (b) temper the L/F KLs by λ_G, which is the same as a power posterior on the genomics block
     with empirical-Bayes priors fitted to the tempered pseudo-data;
   - (c) another choice you derive.
4. **How λ_G is set.** Set it from training data only, so the median survival share of the F
   precision at the paired unsupervised reference hits the existing targets 0.05, 0.2 and 0.5. This
   is the reciprocal of the earlier γ targets, so the arms are comparable.
5. **τ.** Fixed at the training-only unsupervised values is the primary arm. Adaptive τ is a
   secondary arm.

Check the algebra numerically (finite differences of L(q) against each update) before using it.

## Part B, step 2. Implementation

- Add the objective as an opt-in mode of `code/fit_yfb_supervised.R`, for example
  `genomics_weight` with `cox_weight = 1`. Default results must stay unchanged; add a parity test
  against a pre-change fixture.
- Tests:
  - λ_G = 1 with the Cox term at weight 1 reproduces the existing γ = 1 fit;
  - each update does not decrease L(q) on a small fixture, or the decrease is reported;
  - the KL is exact for the tempered case;
  - an all-zero L, F or β is returned as an explicit failure, never as convergence.
- Reuse `code/yfb_diagnostic_controls.R` (fixed τ, guard) and `code/yfb_joint_beta.R` (exact β
  reference) rather than writing new versions.

## Part B, step 3. Registration (before running anything)

Add attempt **S-10** to `docs/plans/2026-10-16-round2-attempts-log.md` and write
`docs/plans/2026-10-30-downweighted-genomics-plan.md`. Register three parts:

1. **Pilot.** The three targeted cases (0.25%, 1%, 5% planted variance; p = 2,000; the same
   synthetic data and seeds as the targeted diagnostics), at each λ_G target and with fixed τ.
   Continue to the grid only if every pilot fit has:
   - no all-zero collapse;
   - bounded reconstruction (expected error ratio below 1);
   - no sweep-cap stall in more than one of the three cases per arm.
2. **Calibration check.** On null synthetic data (no planted program) at each λ_G, the mean number
   of falsely active programs must be at most 1. This is the direct test of whether keeping
   the Cox weight at 1 fixes the overconfidence.
3. **Full grid.** The S-7 grid unchanged (p ∈ {50, 2,000}; planted variance 0.25%, 1%, 5%; 10
   datasets per cell; the same seeds and the unsupervised paired reference). Gate S is
   unchanged: at p = 2,000, a recovery gain of at least 3/10 over the unsupervised fit in at least
   two variance cells, mean false-active at most 1, and at least 90% of fits converged.
   Report S-7 and S-9 next to S-10.

## Part B, step 4. Run, diagnose, document

- Run the pilot, then the calibration check, then the grid if the pilot passed. Use PSOCK workers,
  one BLAS thread each, in the background.
- The diagnose-and-solve rules apply. Any fix after a failed registered run becomes a new numbered
  attempt, S-11 and onward, registered before it is rerun. The original criterion is never edited.
- If the pilot shows the July collapse, diagnose it (per-sweep ‖L‖, ‖F‖, β, the prior scales and
  the KL terms). Try the other KL option from Part B step 1.3 as a registered attempt before concluding.
- Write `docs/findings/2026-10-30-downweighted-genomics-findings.md`. Update the 10/30 chapter
  and its issue tracker (chapter rules: no names, no branch names, file paths or commits, and
  numbers read from saved CSVs), the attempts log, DECISIONS.md, and the plan's "Review
  handoff" section. Commit and push the branch; do not merge into main.

Hard limits: no external cohort informs any choice, the frozen benchmark is untouched, and nothing
is reported as done while tests fail.

## Part C: loose ends on this branch

- The independent review. Fill in the realignment plan's "Review handoff" section, then run the
  review prompt at the end of `docs/prompts/2026-10-30-codex-execution-prompt.md`, or leave it
  ready for a separate review session.
- Step S3 of the realignment plan (real-data training CV), which was never run because Gate S
  failed in S-7. Run it for whichever fitter passes, if one does.
- Done on 2026-11-04: `main` was merged into this branch (656fec7). If `main` gains new commits
  before this work finishes, merge it again and resolve any conflicts the same way.
- Update the October 30 chapter, the issue tracker, DECISIONS.md and the plan's review handoff
  with the S-9 and S-10 results. Commit and push the branch, and do not merge into `main`.

## Model routing (Claude Code only; ignore in Codex)

Use Opus 5.5 as the orchestrator and Sonnet 5.5 subagents for work that does not need its
judgment.

Opus does:
- the derivation and its algebra checks;
- the diagnoses of failures and the choice of the next attempt;
- registration text and gate verdicts;
- reviews of subagent output before any claim is written down.

Sonnet subagents (the Agent tool with `model: "sonnet"`) do:
- running registered runners and the test suite in the background, and reporting exit status,
  counts and errors;
- turning saved outputs into the summary CSVs and tables;
- mechanical documentation updates: README rows, attempts-log formatting, chapter tables that
  read existing CSVs;
- searches across the repository and its branches.

Give each subagent the exact files, commands and expected output, and check its result against
the saved files before relying on it. Never delegate a decision about the model, a gate or what
counts as a fix.

