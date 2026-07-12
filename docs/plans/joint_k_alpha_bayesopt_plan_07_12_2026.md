# Joint (K, α) tuning via Bayesian optimization — future consideration

**Date:** 2026-07-12
**Status:** Plan only, not implemented. Low/medium priority — revisit after Phase 2/3.
**Motivation:** `DECISIONS.md` 2026-07-12 ("Fresh K-CV under corrected code") found a confirmed
methodology difference behind SBMF's K=7 vs. DeSurv's k=3: DeSurv jointly tunes `k, α, λ` via
Bayesian optimization with a fixed elastic-net penalty; we tune only K via grid-search CV with α
fixed at 0.5 and no penalty term at all. This plan sketches what a comparable joint-tuning
procedure for SBMF would look like, so the K=7-vs-3 gap can eventually be attributed to genuine
model capacity rather than a difference in search procedure.

---

## What DeSurv does (as understood from the lab-meeting deck comparison table)

Bayesian optimization jointly over `(k, α, λ)`, where `λ` is DeSurv's elastic-net penalty on the
factors/coefficients (a role SBMF doesn't have an analogue for — SBMF's shrinkage comes from
empirical-Bayes adaptive priors, learned per-CAVI-step, not a single tunable penalty strength).
DeSurv's search presumably optimizes a validation metric (their reported HR/SD, or an analogous
concordance/likelihood measure) over this joint space, landing on `k=3, α=0.334` for PDAC.

## What SBMF would need for a comparable procedure

1. **Search space.** `K` (discrete, e.g. 2–10) × `α` (continuous, [0, 1]). No direct analogue of
   `λ` exists (see above) — the EB priors are refit at every CAVI step regardless of `(K, α)`, so
   there's no third continuous penalty knob to add unless we introduce one deliberately (not
   recommended — this would reintroduce something like the retired `lambda`, which Phase 1b
   concluded was redundant with `alpha`).
2. **Objective.** Cross-validated external (or held-out-fold) C-index, matching what `select_K_cv`
   and `select_alpha_cv` already compute separately. A joint objective would need a single CV
   harness that fits YFB at a candidate `(K, α)` pair across folds and returns mean C-index — this
   is mostly assembling existing pieces (`select_K_cv`'s fold-fitting logic, parameterized to also
   vary α) rather than new modeling work.
3. **Optimizer.** Each objective evaluation is expensive (a full 5-fold CAVI fit), so a
   sample-efficient method matters — this is exactly the regime Bayesian optimization is designed
   for (Gaussian-process surrogate + acquisition function, few dozen evaluations instead of a full
   grid). Candidate R packages: `ParBayesianOptimization`, `mlrMBO`/`mlr3mbo`, or `rBayesianOptimization`.
   Would need to pick one, define bounds (`K` as an integer-valued dimension — most BO packages
   handle this via rounding, not natively), and set an evaluation budget (e.g., 30–50 evaluations,
   comparable in cost to the 45-fit grid search already run for the K-only CV).
4. **Validation of the optimizer itself.** Before trusting a BO-selected `(K, α)`, sanity-check it
   against the existing grid results we already have (K-CV at α=0.5, alpha-CV at K=7) — the BO
   result should be at least as good as either marginal optimum, and ideally should recover
   something close to the known-good region if the search space and objective are correctly wired.

## Scope estimate

- New code: a joint-objective wrapper around the existing fold-fitting logic (`code/select_K.R`,
  `code/select_alpha_cv.R` share most of the CAVI-fitting-per-fold machinery already); a BO driver
  script in `results/benchmark_sim/`; a comparison report.
- Effort: medium — most of the expensive part (fold-fitting) already exists; the new work is the
  BO wrapper, search-space/budget tuning, and interpreting/validating the result.
- Risk: BO packages can be finicky with mixed discrete/continuous spaces and expensive,
  noisy (CV-based) objectives; expect some iteration to get a sensible, reproducible search.

## Why this is not Phase 2/3 work

Phase 2 (joint-vs-2-step value-add) and Phase 3 (K/K_eff analysis on the corrected objective) don't
depend on this — they use the existing K-CV/alpha-CV machinery as-is. This plan is specifically
about making the **K=7-vs-DeSurv's-k=3 comparison methodologically fair**, which is a narrower,
lower-urgency question than either of those. Revisit after Phase 2/3 land, or sooner if the
manuscript (Phase 6) needs a defensible answer to "why does SBMF use more factors than DeSurv."

## Open questions to resolve before implementing

1. Do we want a real elastic-net-style penalty analogue at all, or is comparing DeSurv's
   3-hyperparameter joint search against our 2-hyperparameter (`K, α`) joint search (no penalty)
   an acceptable, honestly-caveated comparison?
2. Which BO package, and what evaluation budget is affordable given CAVI fit cost at PDAC scale
   (n≈273, p≈2064)?
3. Should this run once (single BO search, reported as-is) or be wrapped in an outer
   stability check (repeated BO runs / seeds) given CV-based objectives are noisy, as the K-CV
   table itself showed (SE ≈ 0.03 per K)?
