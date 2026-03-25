# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

---

## Key Instructions

- **Active implementation:** `code/Supervised_Bayesian_MF_V2.R` (V2). V1 is archived at `code/legacy/Supervised_Bayesian_MF.R` — do not modify.
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.md` instead.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: 105/105 passing.

## Quick Reference

| What | Where |
|------|-------|
| Full project docs & session log | `PROJECT_STATUS.md` |
| Code quick-reference (math ↔ R) | `code/SupervisedMF_Context.md` |
| **Reusable CAVI fitting function** | `code/fit_modular.R` (factor-wise, canonical) |
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| V2 simulation report | `results/full_sim/simulation_report.qmd` |
| Modular simulation report (block-wise, deprecated) | `results/modular_sim_block/modular_sim_report.qmd` |
| **Factor-wise modular simulation report (canonical)** | `results/modular_sim_factor/factor_modular_sim_report.qmd` |
| Test suite | `tests/run_tests.R` |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
