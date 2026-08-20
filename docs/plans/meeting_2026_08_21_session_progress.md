# 8/21 Meeting Prep — Session Progress

Tracks execution of the multi-session plan for the 2026-08-21 advisor meeting. Branch:
`meeting/2026-08-21-prep`. Updated at the end of each phase.

- [x] **Prerequisites** — merged `claude/recent-project-updates-n04uuo` into `main` (PR #3,
      `191c2f8`); reconciled and pushed `meeting/2026-08-21-prep` off updated `main`.
- [x] **Phase 0 — Pre-read** — `docs/plans/Meeting_Notes_Followup_8_3_26.md`,
      `docs/progress_book/`, `code/fit_cox_on_yf.R` (YFB — corrected from the originally-cited
      `code/SupervisedMF_Context.md`, which documents the older LB model), `config/globals.yml`.
- [x] **Phase 1 — Model specification chapter** — `docs/progress_book/chapters/meeting_2026_08_21.qmd`:
      §1 (8/3 action-item follow-up, 2 items flagged still open) + §2 (YFB model spec: generative
      model, per-factor CAVI updates, λ deprecation, gene selection, factor composition table).
      Reviewed by code-reviewer subagent (1 real fix: missing α term in q(β_k); 1 notation cleanup);
      narrative trimmed to a plain specification per feedback. Commits `5642999`, `9135532`, pushed.
- [x] **Phase 2 — Factor classification & K-selection validation** — implemented
      `classify_factors()` (survival-active / genomics-only / dead) in `code/select_K.R` + 2 new
      tests (KCV-T17/T18, `tests/test_select_K_cv.R`; T17 covers all 3 categories, T18 the
      borderline-threshold edge case); fixed stale hardcoded `beta_thresh=0.05` in
      `run_phase1_diagnostics.R` to read from `globals.yml`. Ran Analysis A on real TCGA+CPTAC
      data (D4 preprocessing): a single-init sweep (`run_k_init_sweep.R`, K_init ∈
      {5,6,7,8,9,10,15,20}) followed by a 15-restart best-of-multistart ELBO check
      (`run_k_init_multistart_check.R`, K ∈ {5..10}) after the single-init pass showed K=5/K=6
      with the best training ELBO of the whole grid. Result: K_survival_active=2 from K_init=7 up
      (K=5/6 show 3 but fail external validation, C=0.596-0.597 vs. 0.6256-0.6279 for K≥7 — robust
      to multistart, not a single-init artifact). Among K∈{7,8,9,10} (all pass external
      validation), best-of-multistart ELBO is a near-tie between K=7 (2 survival-active + 2
      genomics-only = 4 total factors) and K=9 (2 survival-active + 3 genomics-only = 5 total).
      **Recommendation: K=7**, parsimony tiebreaker; K=9 documented as a near-tied alternative.
      Real-data suite 88/88, full suite 394/394. See DECISIONS.md 2026-08-19.
- [x] **Phase 3 — Analysis B (ARD K-recovery simulation) + Analysis C (signal-ratio sweep re-run)**
      — Analysis B (`run_k_recovery_sim.R`, 45 fits over known ground truth) found ARD
      **over-counts** survival-active factors (91% of fits over-counted, 0% under-counted). A
      follow-up diagnostic (`run_k_recovery_diagnostic.R`) confirmed the "extra" factors are real,
      genuinely non-prognostic study-specific gene programs picking up spurious small β — not
      fragments of the true signal, and not reliably fixed by `cohort_id`. Added an optional
      `rel_thresh` parameter to `classify_factors()` (2 new tests, KCV-T19) that separates real
      from spurious factors cleanly in the simulation — but confirmed it must NOT be applied to
      the real D4 PDAC fit, whose smaller-β factor (Program 3) is independently validated via
      5-cohort external HR<1 (DECISIONS.md 2026-07-15). Analysis C (`run_signal_ratio_sweep.R`,
      new `K_INIT` param) confirmed the original YFB-vs-EBMF finding holds under ARD pruning, and
      confirmed `cohort_id` only partially reduces Analysis B's false-positive pattern. Full suite
      395/395. See DECISIONS.md 2026-08-20. **Open item flagged for discussion, not resolved:**
      ARD's over-attribution of survival signal when no independent validation exists.
- [x] **Phase 3b — Bootstrap CI on the K=5-vs-K=7 external C-index gap** —
      `run_k5_vs_k7_bootstrap_ci.R` (reused the already-fitted K=5/K=7 models, no re-fitting).
      Result: no single external cohort alone is significant (each individually underpowered,
      n=52–288), but pooled across all 5 (n=616) the gap is significant: K=5−K=7 = −0.0312, 95%
      CI [−0.0602, −0.0017]. See DECISIONS.md 2026-08-20 addendum.
- [x] **Phase 4 — Progress book chapter for the 8/21 meeting** —
      `docs/progress_book/chapters/meeting_2026_08_21.qmd` §3 added (how K was chosen from model
      fit + external validity rather than cross-validation; the simulation stress test and why its
      fix can't be mechanically applied to real data; the bootstrap-confirmed K=5-vs-K=7 gap);
      §1's two carried-over 8/3 open items closed out; §2.5's factor table split genomics-only vs.
      fully-pruned. Fact-checked against DECISIONS.md by a review pass (2 findings, both fixed).
      Quarto book renders cleanly. Commit `48c36c7`, pushed.
- [x] **Phase 4b — Two chapter additions requested for §3, plus a 3-round correction to the
      joint-vs-two-step comparison** — added a gene-level summary of K=7's 4 kept factors (table +
      heatmap, `generate_k7_kept_factors_summary.R`) and a YFB-vs-two-step comparison matching the
      real recommended structure exactly. The comparison went through 3 rounds, each catching a
      flaw in the last: (1) matching K to K=7 was itself confounded — it hands the two-step method
      YFB's own answer about model complexity; (2) checked whether EBMF has a natural self-selected
      K (it doesn't — flashier used every factor offered up to K=40, no early stopping) and found a
      second flaw, an unregularized stage-2 Cox that overfits at large K (training C rose while
      external C fell as K grew); (3) fixed both with a large, YFB-uninformed K + LASSO stage 2
      (`run_ebmf_cox_regularized.R`, `run_ebmf_cox_external.R --k`). **Final result: YFB shows a
      real, pooled-significant advantage against the fairest independent baseline tested (+0.026,
      95% CI [0.0002, 0.0498]), and wins numerically in every one of 6 configurations tried** —
      the same fix applied to simulation (`run_k7_signal_sweep.R`'s `EBMF_fair` arm) reproduces the
      originally-intended pattern (equivalence at zero signal, growing advantage as signal
      strengthens), though not independently significant at only 10 seeds. Full numbers: DECISIONS.md
      2026-08-20 (multiple addenda to the same entry). §3.5 rewritten with the final story. Commits
      `51a6dd2`, `91c5cdd`, pushed. Full suite 395/395.
- [ ] **Phase 5+** — not yet scoped in this session; remaining meeting-prep phases (slides/figures,
      any further sections) follow once this phase lands, and are the user's call. Chapter review
      and editorial revision for clarity/concision is the immediate next step, planned for a new
      session (see handoff prompt).
