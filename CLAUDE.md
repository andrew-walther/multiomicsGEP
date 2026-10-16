# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

This is the single source of truth for agent instructions in this repo. `AGENTS.md` is a symlink to
this file (so Codex, Antigravity, and other tools that look for `AGENTS.md` read the same content) —
edit `CLAUDE.md`, never `AGENTS.md` directly. If a `GEMINI.md` or similar is ever added, symlink it
here too rather than maintaining a separate copy.

---

## Key Instructions

- **Deliverable and code attribution:** Do not name assistant tools or include tool-attribution lines in presentation deliverables or code. Tool names remain appropriate in explicit review handoffs. Andrew requested this on 2026-10-15; it supersedes earlier template author lines.

- **Canonical CAVI loop:** `code/fit_modular.R`. `Supervised_Bayesian_MF_V2.R` is reference-only — do not extend. V1 (`code/legacy/`) — do not modify. (`full_sim/` and `modular_sim_block/` are legacy/deprecated.)
- **Formal benchmark pipeline:** `results/benchmark_sim/` — alpha CV, external validation, DeSurv comparison. Exploratory/development fits lived in `results/modular_sim_factor/` (now archived to `results/legacy/modular_sim_factor/`).
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **Global constants:** `config/globals.yml` — all hyperparameters (lambda, alpha grid, K thresholds, DGP params). Never hardcode values defined here.
- **Multimodal YFB proposal:** `derivations/multimodal_YFB/multimodal_YFB_derivation.{qmd,pdf}` contains the concise nine-page linear joint-YFB derivation; `multimodal_YFB_implementation.qmd` in the same directory specifies future R code, and `multimodal_YFB_technical_notes.{qmd,pdf}` preserves the extensive version; `verify_derivation.R` checks the algebra without fitting a model. Read it before multimodal work. An isolated implementation exists on the `codex/multimodal-yfb` branch (`code/*multimodal_yfb*.R`: fitter, updates, prediction, simulation, K-selection runner, with tests); it is work in progress, not the production single-modality model, and its defaults are not approved. Current status and review decisions belong in `PROJECT_STATUS.qmd` and `DECISIONS.md`.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.qmd` instead.
- **Living documents:** Update `DECISIONS.md` when making any architectural choice (algorithm variant, hyperparameter decision, design tradeoff). Update `ROADMAP.md` when completing a milestone or identifying a new priority.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: all passing (532 as of 2026-10-14 on `fix/multimodal-monitor-objective`).
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 88/88 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git (`**/PDAC_data/` and the local `/data/` folder of TCGA/ICGC copies are git-ignored; never commit them). Stored locally at `~/Library/CloudStorage/OneDrive-.../UNC Dissertation (Liu)/PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.
- **Current model status:** Both LB (`code/fit_modular.R`, η = Lβ) and YFB (`code/fit_cox_on_yf.R`, η = (YF)β) are fully implemented. **Recommended configuration:** YFB × per-platform z-std × DeSurv gene selection × no cohort indicator, K=7 default (K selection is a two-stage ARD framework — K_init from an ELBO/BIC/log-likelihood/CV-C consensus, then `classify_factors()` shrinkage gives K_eff from that one over-specified fit), K_eff=2, mean external C≈0.627 across 5 held-out PDAC cohorts. Optional `cohort_id` (genomics offset), `strata_id` (stratified baseline hazard, off by default — performance-neutral), and `beta_cohort_id` (cohort-specific survival coefficients) parameters exist; see the Quick Reference table below. Full history, benchmark numbers, and the reasoning behind each choice live in `PROJECT_STATUS.qmd`; the dated rationale entries are in `DECISIONS.md` — read those before revisiting any of this, don't rely on this summary alone.
- **Documentation audience:** Write ROADMAP.md, DECISIONS.md, PROJECT_STATUS.qmd, and README.md for biostatistician collaborators reading the project cold — not as implementation logs. Avoid internal session terminology (e.g. "Phase A/B/C", "Cluster A/B", "Session N") in prose descriptions; use those labels only in commit messages or as lookup keys. Describe methods and findings in terms a statistical reader would recognize: model variant, prior, training set, metric, result.
- **No mannered prose:** Say what you mean; don't substitute metaphor or flourish for direct statement (e.g. "a dial worth turning" instead of "a parameter worth varying," "this point earns its keep" instead of "this point still matters"). Applies to all writing in this repo — reports, the progress book, plan/decision docs, commit messages, comments. When a literal phrase is available, use it.

- **Protected-pruning evaluation:** opt-in changes and fixed runners are on `m5/k-selection`; adoption was blocked. The controlling fallback decision is recorded at the top of `DECISIONS.md`.

- **Planted comparator evaluation:** runner and partial tables are on `m7/sim-comparators`; see the 2026-10-15 decision for the zero/one-factor EBMF blocker and incomplete high-noise denominator.

- **Multimodal comparator audit:** `results/multimodal_real/run_supervised_multimodal_comparison.R` records actual two-step ranks and the fixed four-arm evaluation on `m7/multimodal-comparators`; see the dated progress-book and decision entries for results.

## Robustness and comparator follow-up

The five-prompt branch/results index is
`docs/reports/Codex_M6_M7_Summary_10_15_26.md`; read its limitations before using
the tables. E did not adopt protection or a new K. `m6/robustness` contains the
training-only robustness runner; current findings belong in `PROJECT_STATUS.qmd`
and `DECISIONS.md`, not this instruction file.

Round 2 work is cumulative on `round2/supervised`. Its controlling specification is
`docs/plans/Codex_Handoff_Round2_10_16_26.md`; status and attempts are in
`docs/plans/Codex_to_Claude_Round2.md` and `docs/plans/Round2_attempts_log.md`.
The benchmark tag and `results/benchmark_snapshot_2026-10-14/` are immutable.

## Quick Reference

Active pointers only — completed one-off analyses and dated reports are indexed in `PROJECT_STATUS.qmd`,
`ROADMAP.md`, and `DECISIONS.md` instead of here.

| What | Where |
|------|-------|
| Code quick-reference (math ↔ R) | `code/SupervisedMF_Context.md` |
| **Reusable CAVI fitting function** | `code/fit_modular.R` (factor-wise, canonical) |
| Hold-out prediction | `code/predict.R` — `predict_supervised_mf()` |
| Train/test splitting | `code/train_test_split.R` — `stratified_split()` |
| Feature selection | `code/feature_selection.R` — `cox_feature_selection()` |
| K selection | `code/select_K.R` — `auto_prune_K()`, `select_K_cv()` |
| Full ELBO computation | `code/compute_elbo.R` — `compute_ebnm_kl()`, `compute_survival_elbo()`, `compute_normal_kl()` |
| Joint log-likelihood / BIC for K_init sweeps | `code/compute_bic.R` — `compute_joint_ll_bic()`; source after `fit_cox_on_yf.R` |
| Cross-validated held-out log-likelihood | `code/compute_cv_loglik.R` — `cv_survival_loglik()`, `bicv_genomics_loglik()`; companion to `compute_bic.R`, reported separately |
| Cohort-specific survival coefficients | `code/update_beta_cohort.R` — `beta_cohort_id` argument on `fit_cox_on_yf()`. Distinct from `cohort_id` (genomics offset) and `strata_id` (baseline hazard) |
| Cohort F update | `code/update_F_cohort.R` — `update_F_cohort_all()` (Normal conjugate) |
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| **LB benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **YFB benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **K_init consensus sweep** | `results/benchmark_sim/run_k_init_sweep.R` — ELBO/BIC/log-likelihood/held-out-C consensus; `--reuse-cache`/`--reuse-cv-cache` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Pathway enrichment functions | `code/pathway_enrichment.R` — `load_d4_weights()`, `run_fgsea_program()`/`run_ora_program()`, `build_pdac_genesets()`, subtype/cohort/DeSurv-overlap concordance functions |
| **Progress notebook (meeting-facing, one chapter per advisor meeting)** | `docs/progress_book/` — Quarto book, `quarto render` to build; add a new `chapters/YYYY-MM-DD.qmd` per meeting |
| **L-update debugging guide** | `docs/update_L_fix.md` — read before any L-update work |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| **Current working plan and goal hierarchy** | `docs/plans/Advisor_Priorities_and_Execution_10_15_26.md` — restart here; `Codex_to_Claude_Round2.md` in the same folder is the current status summary. `Meeting_Notes_Followup_10_2_26.md` preserves advisor requirements; older schedules are historical. |
| Matched TCGA/ICGC expression + methylation loader | `code/load_multiomics_data.R` — `build_multiomics_cohorts()`; real-data checks `tests/test_real_multiomics_loading.R` (local) |
| Multimodal parsimony framework (intercept, ebnm priors, ELBO pruning) and compiled sweep | `code/fit_multimodal_yfb.R` controls `intercept`/`prior_update`/`prune`/`tau_model`; `code/multimodal_yfb_sweep.cpp` |
| Multimodal simulation and real-data runners | `results/benchmark_sim/run_multimodal_pruning_comparison.R`; `results/multimodal_real/run_multimodal_real_fit.R` + `summarize_multimodal_real.R` |
| Prelim proposal chapter source | `paper/prelim/project3-ssbmf.qmd` → synced to `bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` by the installed post-commit hook (or `paper/prelim/tools/sync_to_prelim.sh` by hand; never pushes) |
| Prelim proposal chapter plan | `docs/plans/Prelim_Proposal_Plan_10_1_26.md` — outline, sync to `bios-dissertation`, α_F framing constraint |
| Advisor-facing writeup (approved interim preserved in the benchmark tag; bounded feature-sensitivity revision authorized 2026-10-15) | `paper/writeup/writeup.qmd` → `writeup.html`/`writeup.pdf`, appendices in `paper/writeup/appendices/`; planning page, source inventory, implementation plan and prompts in `paper/writeup/_process/`. Read `PROJECT_STATUS.qmd` for the review state. |

The opt-in single pruning observer is
`results/benchmark_sim/run_single_pruning_trace.R`; it records process-local
nullcheck states without editing the fitter. Its fixed diagnosis and the rank
panel's Longleaf wrapper are documented in the Round 2 attempt history.
Andrew approved an additional opt-in Task1 remedy after the stage trace:
`prune_preserve_posterior=TRUE` carries the retained full single-modality state
into pruning refits. `run_posterior_k_panel.R` evaluates the unchanged 160-cell
protected single panel against the saved original reference. The completed
experiment improves recovery but fails the unchanged adoption rule (recovery, convergence and rank); see the attempt log. Defaults and the
benchmark remain unchanged.

Task2's fixed TCGA-only diagnostic runner is
`results/multimodal_real/run_round2_collapse_diagnosis.R`. Its opt-in trace and
diagnostic prior controls preserve the default fit. Read the attempt history
before execution or interpretation; training contrasts do not adopt a remedy,
and external cohorts cannot select one.

Task3's pinned DeSurv package is installed in ignored `data/round2_R_library/`;
its preprocessing/projection adapter, driver and Longleaf wrappers are prepared
and tested. Andrew deferred the head-to-head for review. Read the Round2 handoff for author SHAs, installation details and the
raw-loader log-transform concern before changing any preprocessing.

The ten-arm TCGA diagnosis has completed; consult the Round2 handoff for
benchmark-linked results and the open remedy ladder. Andrew deferred Task3's
real-data DeSurv head-to-head for review. Prepared DeSurv commands above must
not be executed without renewed authorization; no real-data DeSurv fit or
external scoring has run. Task4 omits DeSurv as explicitly approved.

Task4's comparator runner now preserves all method rows and full denominators
for EBMF rank0/rank1 and isolated failures. The Round2 attempt log controls
its five-arm rerun; DeSurv is omitted by Andrew's explicit decision.

Task2 R1's fixed-prior evaluation is `run_round2_fixed_prior.R` in the benchmark
folder. Its committed attempt history governs80 simulations and20 TCGA CV fits;
all criteria and defaults remain unchanged, with no external data for selection.

`results/multimodal_real/run_framework_fold_trace.R` is the fixed Task2
framework CV fold1 trace. Read its preregistration before execution. It changes
observation only and must reproduce the saved posterior/history exactly.

The additional approved Task2 profiled-curvature experiment uses
`code/profile_cox_intercept.R` and opt-in `profile_cox_intercept=TRUE`;
its runner is `results/benchmark_sim/run_round2_profiled_cox.R`.
Read Task2-4 in the attempt log and
`docs/plans/Round2_profiled_curvature_proposal.md` before execution.
The default/off path must stay identical to `feature/supervised-F`.
The experiment uses framework Normal beta priors; shared-prior combinations
are explicitly unsupported. The paired-CV rule and C=NA remain unchanged.

A prepared prior-consistency observer is
`results/benchmark_sim/run_round2_prior_consistency_trace.R`, with aggregate
point-Laplace moment helpers and focused tests. Its bounded Task2-5 observation is complete; see the result in the handoff.
Any additional execution requires a separate prospective diagnostic attempt in
`docs/plans/Round2_attempts_log.md`. It repeats one completed simulation and
requires exact saved-fit identity; its adjusted monitor is diagnostic only,
never a replacement convergence or adoption criterion.


The revised advisor-priority plan is approved:
`docs/plans/Advisor_Priorities_and_Execution_10_15_26.md`. Its first report is
`docs/reports/Advisor_Progress_since_10_2_26.pdf`. The next descriptive study,
`results/multimodal_real/summarize_cpg_variability.R`, must follow Task5-QC1's
committed design before execution. It preserves raw TCGA missingness, uses
only the existing training patients, verifies frozen screened probe IDs and
applies no new preprocessing/default/SD cutoff. The local inventory is ignored;
only aggregate plots/tables/provenance are committed. Read its actual result
status in the Round2 handoff, not this pointer.


Feature-study preparation/agreement utilities are in `code/feature_sensitivity.R`.
They preserve the existing cohort ranking and undefined dead-factor correlations;
no feature-count fit is implied by their presence. The cold-primary versus
rank-deficient continuation decision is recorded in the active priority plan
and Round2 handoff. Do not start a performance panel before that decision and
its prospective committed design. No model defaults were changed.


Task2-5 registers the bounded prior-consistency observation before execution.
Its additional native second-moment telemetry confirms point-prior support
conflicts independently of tiny reconstructed-probability differences. Read
its actual run status in the attempt log. It does not repair an update or
change convergence/adoption criteria, and does not by itself authorize feature fits. Andrew subsequently approved cold-primary
judgment; the feature-panel design must still be committed before fitting.

Fixed-K expression sensitivity is `results/benchmark_sim/run_feature_expansion.R`.
Read Task5-F1 before execution: cold-primary is authorized, alpha_F0/K7/default
benchmark controls remain fixed, and patient folds condition on the existing
screen/scaling. No external selection or benchmark replacement is authorized.
Use only ignored local patient inputs/fits and safe aggregate outputs.

`results/benchmark_sim/summarize_feature_expansion.R` builds the fixed-K
performance/agreement figures and safe loading summaries from all24 terminal
fits. Its input patient/fitted RDS remain ignored. Consult the live status and
Round2 handoff for results, limitations and the next methylation preflight.

CpG expansion preparation is `results/multimodal_real/prepare_feature_universe.R`;
read Task5-P1 before execution. Its guards are in `code/cpg_feature_sensitivity.R`.
It must reconstruct the unchanged whole-universe imputation and exact frozen
10k matrix/rank before any fit. Probe-count expansion reaches the full eligible
universe;10k is a reference only. Describe the supplied mean/proportion filters
separately; no SD cutoff or default change is authorized by this preflight.

The Task5-M1 methylation-only count driver is
`results/multimodal_real/run_feature_expansion_multimodal.R`, with a manual
Longleaf array wrapper. Follow its prospective entry before running any model.
It uses the historical framework, fixed K7/no deletion, the exact P1 matrix,
and all five counts through305352. Filters remain a separate comparison;
constant-risk C stays NA and no external cohort selects a count or model.

Saved-program role/agreement audit: read Task5-Roles1 before running
`results/multimodal_real/run_saved_program_roles.R`. Its tested helpers in
`code/multimodal_program_agreement.R` preserve the raw multimodal projection,
require exact patient order, compare loadings within modality and describe
within-block rank-one norm shares (not additive PVE). Both saved full fits
are unconverged; no external scoring, enrichment, label transfer or adoption.
Task5-M1 first10k methylation fit has finished at the500 cap with nonconstant
risk but zero credible coefficients. Preserve this failure in count summaries.

Andrew approved Task5-M2, a bounded identical-settings25k failure diagnostic.
Read its committed attempt before running
`results/multimodal_real/diagnose_methylation_moments.R`. It observes F-update
inputs via process-local tracing, preserves the original error/results and
keeps detailed patient state ignored. It is not a scientific remedy; follow
its actual terminal status before any further count fit or corrective change.

Task2-B1 is the bounded saved-state loading-balance observation at
`results/multimodal_real/audit_loading_balance.R`. Read its prospective log
before execution. The helper preserves separate reconstruction/survival
conditional terms; ratios are descriptive, not a selected weight or an
executed sequential sweep. Prior penalties and cached objective components
are distinct. Andrew subsequently prioritizes fixing25k Cox underflow; this
audit must not replace that requested work.

Andrew requests resolving the25k CpG numeric hangup. Task5-M3 is an opt-in
same-likelihood arithmetic check, not a new model weight or default:
`control$cox_arithmetic="logspace"` uses stable risk-set sums in working and
monitoring calculations. Read the prospective attempt and mathematical note
before `results/multimodal_real/run_logspace_cpg_case.R`. Default/original
objects must remain identical to the protected base; no imputation, caps,
priors or adoption rules change. Numeric repair does not prove convergence.

Task5-A1 prepares the requested provisional10k/25k CpG loading agreement:
`results/multimodal_real/compare_cpg_count_programs.R`. Read its prospective
entry and wait for M3's same live process to terminate. It starts no model;
missing/errored output means comparison unavailable. Record original versus
log-space arithmetic and both convergence states; no minimum-r rule or
biological-label transfer is authorized by this descriptive comparison.


Planning consolidation: use the active priority plan for the bounded next
work block and goal hierarchy, and the concise Round2 handoff for status.
The former long handoff is preserved verbatim at
`docs/plans/Codex_to_Claude_Round2_history_10_15_26.md`; its old live/pending
statements are historical. Planning changes do not resume the paused goal.
Original handoffs retain scientific criteria and provenance; their obsolete
branch/merge sequences are not restart instructions. Read current status
before repeating any diagnostic or invoking a prepared runner.
