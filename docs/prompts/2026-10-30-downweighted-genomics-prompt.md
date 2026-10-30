# Codex prompt: register and run the downweighted-genomics test (attempt S-8)

Written 2026-10-30. Andrew approved registering this test. Paste everything below the line into
Codex from the repository root, on branch `fix/yfb-supervision`.

---

Continue on `fix/yfb-supervision`. Read these first:

- `docs/findings/2026-10-30-yfb-targeted-diagnostics-findings.md`: the results this test builds on.
- `docs/findings/2026-10-30-yfb-supervision-problem-statement.md`: the five hypotheses.
- `docs/plans/2026-10-30-realignment-plan.md`: the "Problems during execution" rules and Gate S.
- `docs/prompts/2026-10-30-codex-execution-prompt.md`: the project rules and chapter rules, which
  still apply.
- `DECISIONS.md`, the July entry near "Two candidate directions were tried and empirically rejected"
  (around line 4029). Shrinking the genomics precision in the latent-L model drove L, F and β to
  exactly zero within 2–3 iterations, because L and F co-adapt bilinearly.

## Question

The targeted diagnostics showed two things. A large Cox weight γ makes the β posterior falsely
precise: about 6.9 of 7 null coefficients are flagged active even under exact inference. Free
feature precisions τ_j also feed a reconstruction blow-up. The test is this: if the genomics
likelihood is downweighted and the Cox likelihood stays at weight 1, does survival shape F while β
stays calibrated and the fit stays bounded and converges?

## Step 1. Derivation (write before any code)

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

## Step 2. Implementation

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

## Step 3. Registration (before running anything)

Add attempt **S-8** to `docs/plans/2026-10-16-round2-attempts-log.md` and write
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
   Report S-7 next to S-8.

## Step 4. Run, diagnose, document

- Run the pilot, then the calibration check, then the grid if the pilot passed. Use PSOCK workers,
  one BLAS thread each, in the background.
- The diagnose-and-solve rules apply. Any fix after a failed registered run becomes a new numbered
  attempt, S-9 and onward, registered before it is rerun. The original criterion is never edited.
- If the pilot shows the July collapse, diagnose it (per-sweep ‖L‖, ‖F‖, β, the prior scales and
  the KL terms). Try the other KL option from Step 1.3 as a registered attempt before concluding.
- Write `docs/findings/2026-10-30-downweighted-genomics-findings.md`. Update the 10/30 chapter
  and its issue tracker (chapter rules: no names, no branch names, file paths or commits, and
  numbers read from saved CSVs), the attempts log, DECISIONS.md, and the plan's "Review
  handoff" section. Commit and push the branch; do not merge into main.

Hard limits: no external cohort informs any choice, the frozen benchmark is untouched, and nothing
is reported as done while tests fail.
