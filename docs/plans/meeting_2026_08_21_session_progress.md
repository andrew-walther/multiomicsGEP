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
- [ ] **Phase 2 — Factor classification & K-selection validation** — implement
      `classify_factors()` (survival-active / genomics-only / dead) and the K-init stability sweep
      on real PDAC data, per `docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md`
      (Step 1 + Analysis A). Answers the two items flagged open in §1: how K=7 was chosen, and
      whether the 2 survival-active / 5 genomics-only split is stable under ARD pruning from a
      larger starting K.
- [ ] **Phase 3+** — not yet scoped in this session; remaining steps of the 8/13 plan (Analysis B:
      ARD K-recovery simulation; Analysis C: signal-ratio sweep re-run) and any further meeting-prep
      phases follow once Phase 2 lands.
