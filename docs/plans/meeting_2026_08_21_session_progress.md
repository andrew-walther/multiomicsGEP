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
      Real-data suite 88/88, full suite 394/394. See DECISIONS.md 2026-08-19. Analyses B/C
      deferred to a later phase per `docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md`.
- [ ] **Phase 3+** — not yet scoped in this session; remaining steps of the 8/13 plan (Analysis B:
      ARD K-recovery simulation; Analysis C: signal-ratio sweep re-run) and any further meeting-prep
      phases follow once Phase 2 lands.
