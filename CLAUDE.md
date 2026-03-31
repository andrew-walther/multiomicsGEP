# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

---

## Key Instructions

- **Active implementation:** `code/Supervised_Bayesian_MF_V2.R` (V2). V1 is archived at `code/legacy/Supervised_Bayesian_MF.R` — do not modify.
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.md` instead.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: 105/105 passing.
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 77/77 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git. Stored locally at `~/Library/CloudStorage/OneDrive-.../PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.

## Quick Reference

| What | Where |
|------|-------|
| Full project docs & session log | `PROJECT_STATUS.md` |
| Code quick-reference (math ↔ R) | `code/SupervisedMF_Context.md` |
| **Reusable CAVI fitting function** | `code/fit_modular.R` (factor-wise, canonical) |
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| V2 simulation report | `results/full_sim/simulation_report.qmd` |
| Modular simulation report (block-wise, deprecated) | `results/modular_sim_block/modular_sim_report.qmd` |
| **Synthetic simulation report (canonical)** | `results/modular_sim_factor/synthetic/factor_modular_sim_report.qmd` |
| **PDAC real-data report** | `results/modular_sim_factor/PDAC/factor_modular_sim_report_PDAC.qmd` |
| **Simulation runner (synthetic + real)** | `results/modular_sim_factor/run_factor_modular_simulation.R` |
| Test suite (core modules) | `tests/run_tests.R` (105/105) |
| Real-data test suite | `tests/test_real_data_loading.R` (77/77, local-only) |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
