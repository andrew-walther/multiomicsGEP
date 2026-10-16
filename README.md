# multiomicsGEP

**Supervised Bayesian Matrix Factorization for Joint Genomics and Survival Modelling**

---

## Overview

This project implements a **Supervised Bayesian Matrix Factorization** model that jointly decomposes high-dimensional genomics data — such as gene expression or DNA methylation matrices — while simultaneously modelling patient survival outcomes via a Cox proportional hazards model.

The core insight is that standard unsupervised factorization (e.g. PCA, NMF) finds latent structure that explains genomic variance, but has no reason to recover factors that are *clinically meaningful*. By supervising the factorization with survival data, this model discovers **Gene Expression Programs (GEPs)** that are both genomically coherent and prognostically relevant.

**The model:**

```
Y (n×p) = L (n×K) × F' (K×p) + E        ← matrix factorization (genomics)
h(t_i)  = h₀(t_i) exp(ηᵢ)               ← Cox proportional hazards (survival)
```

- **L** — patient loading matrix: coordinates of each patient in the latent factor space
- **F** — factor weight matrix: gene loadings that define each program
- **β** — survival coefficients: prognostic weight of each factor
- **τ** — feature-specific noise precision

Two parameterizations of the linear predictor η are implemented and benchmarked:

| Model | Linear predictor | Key property |
|-------|-----------------|--------------|
| LB    | η = Lβ          | Factor scores learned jointly with survival; `code/fit_modular.R` |
| YFB   | η = (YF)β       | Predictor computed directly from observed expression; eliminates train/test mismatch in projection; `code/fit_cox_on_yf.R` |

Inference is performed via **Coordinate Ascent Variational Inference (CAVI)**, where each variational update is solved as an **Empirical Bayes Normal Means (EBNM)** problem — promoting sparsity in both F and β through point-normal or normal priors.

---

## Multimodal YFB (in progress)

An extension of YFB to matched gene expression and DNA methylation (TCGA/ICGC),
developed on the `codex/multimodal-yfb` branch.
[PDF](derivations/multimodal_YFB/multimodal_YFB_derivation.pdf) and
[editable Quarto source](derivations/multimodal_YFB/multimodal_YFB_derivation.qmd).
The nine-page document derives linear joint YFB with separate loading priors,
the diagonal Cox working approximation, and L/F/β/Tau updates. The linked
[R implementation specification](derivations/multimodal_YFB/multimodal_YFB_implementation.qmd)
maps the derivation to the isolated implementation. The runnable matched-block
fitter, frozen predictor, simulation generator, and K-selection runner are in
`code/fit_multimodal_yfb.R`, `code/predict_multimodal_yfb.R`,
`code/simulate_multimodal_yfb.R`, and `code/run_multimodal_yfb_simulation.R`.
The multimodal module is not the production single-modality PDAC model; its
replicated K-selection study is still in progress.

Options added on 2026-10-01 (defaults keep the original behaviour):

- **Loading priors per modality.** `control$prior_F` sets point-exponential,
  point-Laplace or Normal, with the signed priors fit by `ebnm`.
- **Factor pruning, which together selects K.**
  - a per-feature intercept (`intercept`);
  - marginal-likelihood prior fits (`prior_update = "ebnm"`);
  - an ELBO nullcheck that removes factors (`prune`), with the Cox partial
    log-likelihood as the survival term.
- **Faster sweep.** The per-feature loading sweep is compiled C++
  (`code/multimodal_yfb_sweep.cpp`), with an automatic fallback to R.

Matched TCGA (training) and ICGC (validation) expression + methylation cohorts are
built by `code/load_multiomics_data.R`, which reads a local, git-ignored `data/`
folder. The real-data runner is `results/multimodal_real/run_multimodal_real_fit.R`;
the simulation comparison of fitting variants is
`results/benchmark_sim/run_multimodal_pruning_comparison.R`.
The original extensive derivation is preserved as technical notes in the same folder.
Run the base-R mathematical checks with
`Rscript derivations/multimodal_YFB/verify_derivation.R`.
Rebuild the PDF with
`quarto render derivations/multimodal_YFB/multimodal_YFB_derivation.qmd --to pdf`;
rendering also executes the checks and stops on assertion failure.

---

## Project writeup

The advisor-facing writeup is being organized in `paper/writeup/`. The
[HTML planning page](paper/writeup/_process/planning.html) covers the proposed Methods,
Results, Biology, and derivation appendices, with a
[source inventory](paper/writeup/_process/source_inventory.csv) and
[implementation checklist](paper/writeup/_process/implementation_plan.md).
The planning review has been checked and the plan finalized (October 14).
**Interim writeup, approved and merged to `main` (October 14).**
The report is `paper/writeup/writeup.qmd`, with derivation appendices in
`paper/writeup/appendices/`, rendered to [HTML](paper/writeup/writeup.html) and
[PDF](paper/writeup/writeup.pdf). It is organized as a manuscript: Methods,
Evaluation, Results, Biological findings and Discussion, with complete
derivations in Appendices A (single-modality) and B (multimodal). The review is
[`review.md`](paper/writeup/_process/review.md) and the disposition of every finding is
[`review_dispositions.md`](paper/writeup/_process/review_dispositions.md); Codex's
verification of the revision is
[`verification.md`](paper/writeup/_process/verification.md). The second verification is
[`verification_round2.md`](paper/writeup/_process/verification_round2.md); its two minor
corrections are applied. A final editorial round
([`final_editorial_check.md`](paper/writeup/_process/final_editorial_check.md),
[`verification_round3.md`](paper/writeup/_process/verification_round3.md)) moved the algorithms
into the appendices; the 22-page main text is accepted for this interim version.

Render the report from the repository root (the RStudio Quarto works on this
computer; rendering reads only saved aggregates, figures and the bibliography,
and installs no LaTeX packages):

```bash
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/writeup.qmd --to html
```

```bash
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/writeup.qmd --to pdf
```

The report's figures are drawn from saved fits, scores and aggregates (no
fitting), as vector PDF and 300-dpi PNG at printed size. The survival-curve
figures need the local held-out expression data:

```bash
Rscript paper/writeup/figures/make_pruning_rank_figure.R
```

```bash
PDAC_DATA_ROOT=data/PDAC_data Rscript paper/writeup/figures/make_writeup_figures.R
```

Rebuild the planning page from the repository root:

```bash
quarto render paper/writeup/_process/planning.qmd --to html
```

If the standalone installation reports `Bad CPU type in executable`, use the
already-installed RStudio Quarto (verified locally):

```bash
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/_process/planning.qmd --to html
```

The page uses existing aggregate results and figures; rendering does not fit
models or require patient-level data. Source links are relative to this
repository, while its displayed content and image resources are self-contained.

## Repository Structure

```
multiomicsGEP/
├── README.md  LICENSE (MIT)
├── CLAUDE.md                ← agent instructions (AGENTS.md is a symlink to it)
├── PROJECT_STATUS.qmd/.pdf  ← full project documentation and development log
├── DECISIONS.md             ← dated architectural and analytical decisions
├── ROADMAP.md               ← prioritized next steps and completed items
├── config/globals.yml       ← single source of truth for hyperparameters
├── code/
│   ├── fit_modular.R        ← canonical CAVI loop for LB (η = Lβ)
│   ├── fit_cox_on_yf.R      ← YFB (η = (YF)β); predict_cox_on_yf.R for hold-out scoring
│   ├── update_*.R           ← modular CAVI updates (β, L, F, τ, and cohort/YFB variants)
│   ├── compute_elbo.R  compute_bic.R  compute_cv_loglik.R   ← ELBO, BIC, held-out log-likelihood
│   ├── select_K.R  select_alpha_cv.R  select_k_alpha_bo.R   ← K and alpha selection
│   ├── preprocess_desurv.R  feature_selection.R  train_test_split.R  predict.R
│   ├── pathway_enrichment.R  concordance_ci.R
│   ├── *multimodal_yfb*.R   ← multimodal YFB fitter, updates, prediction, simulation (in progress)
│   ├── Supervised_Bayesian_MF_V2.R   ← monolithic V2 reference (do not extend)
│   ├── SupervisedMF_Context.md       ← math ↔ code quick reference
│   └── legacy/              ← V1 and early scripts (archived)
├── tests/                   ← run_tests.R plus one test file per module
├── results/
│   ├── benchmark_sim/       ← formal benchmark pipeline: LB/YFB runners, K and alpha CV,
│   │                          external validation, DeSurv and EBMF comparisons
│   ├── figures/  tables/    ← per-cohort outputs
│   └── legacy/              ← retired simulation generations
├── derivations/             ← corrected CAVI derivations (MF_UpdateDerivations/), per-update
│                              derivations (qB, qL, qF, qTau), multimodal_YFB/
├── docs/
│   ├── progress_book/       ← Quarto book, one chapter per advisor meeting
│   │                          (https://andrew-walther.github.io/multiomicsGEP/)
│   ├── plans/               ← working plans (current: Working_Plan_10_1_26.md)
│   ├── reports/             ← dated analysis reports
│   └── *.qmd/.pdf/.html     ← per-update walkthroughs, PDAC data audit
├── presentation/            ← lab-meeting decks (latest: walther_lab_meeting_08_27_2026/)
├── longleaf_setup/          ← UNC Longleaf SLURM scripts
├── paper/                   ← manuscript draft (multiomicsGEP_manuscript.qmd)
└── data/                    ← local TCGA/ICGC copies; git-ignored, never committed
```

---

## Quickstart

### Prerequisites

```r
install.packages(c("survival", "ebnm"))
```

### Run the Formal Benchmark (current entry point)

```bash
# LB model benchmark (η = Lβ; alpha mixing CV-selected per training set)
Rscript results/benchmark_sim/run_LB_benchmark.R

# YFB model benchmark (η = (YF)β; Cox-on-YF reformulation)
Rscript results/benchmark_sim/run_YFB_benchmark.R
```

Outputs go to `results/benchmark_sim/outputs/`. Reports are versioned by date in `docs/reports/` — see `desurv_alignment_report_05_27_26.pdf` for the most recent external validation results.

### Run the Real PDAC Analysis

For real PDAC data, use `fit_supervised_mf_modular()` in `code/fit_modular.R` directly, or use the archived exploratory runner `results/legacy/modular_sim_factor/run_factor_modular_simulation.R`. The active benchmark runners (`run_LB_benchmark.R`, `run_YFB_benchmark.R`) are the canonical entry points for formal evaluation.

**Available datasets:**

| Dataset | Platform | n | Censoring |
|---------|----------|---|-----------|
| TCGA_PAAD | RNA-seq | 144 | 48% |
| CPTAC | Proteomics | 129 | 50% |
| Dijk | RNA-seq | 90 | 10% |
| Moffitt_GEO_array | Microarray | 123 | 33% |
| PACA_AU_array | Microarray | 63 | 40% |
| PACA_AU_seq | RNA-seq | 52 | 40% |
| Puleo_array | Microarray | 288 | 37% |

**Note:** PDAC data files are stored locally (not in git). The default path is
`~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Liu)/PDAC_data`.
Override with `PDAC_DATA_ROOT` (e.g. `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC` on Longleaf).

### Apply to Real Data

The recommended entry point for real data is `fit_supervised_mf_modular()` in `code/fit_modular.R`:

```r
source("code/fit_modular.R")   # also sources update_L/F/beta/tau.R automatically

res <- fit_supervised_mf_modular(
  Y      = your_matrix,    # numeric matrix: n patients × p genes (pre-normalised, column-centred)
  time   = your_time,      # numeric vector: survival/censoring time
  status = your_status,    # integer vector: 1 = event, 0 = censored
  K      = 5,              # number of latent factors (select via cross-validation)
  max_iter = 300,
  tol      = 1e-3,
  verbose  = TRUE
)

# Access results
res$EL     # n×K posterior mean patient loadings
res$EF     # p×K posterior mean factor weights (GEP signatures)
res$EBeta  # K posterior mean survival coefficients
res$EBeta2 # K posterior second moments (uncertainty)
res$Tau    # p noise precision per feature
res$history$rmse        # RMSE per iteration
res$history$elbo_proxy  # genomics ELBO per iteration
```

---

## What to Expect

On the simulated benchmark:

| Metric | Expected |
|--------|----------|
| Reconstruction RMSE | Converges near **1.0** (true noise SD) |
| ELBO proxy | **Non-decreasing** across iterations |
| β sign recovery | **(+, −, +, −, 0)** — matches ground truth |
| C-index (modular) | ~**0.86** on held-out tertiles |
| Factor 5 (β=0) | Correctly shrunk toward zero |

---

## Mathematical Background

The model is derived and documented in two companion documents (both in `derivations/MF_UpdateDerivations/`):

| Document | Description |
|----------|-------------|
| **`MF_Derivations_UpdateAlgo_REVISED.pdf`** | Full corrected CAVI derivations with all 8 errata (R1–R8) from the original working document resolved |
| **`MF_V2_Companion.pdf`** | 17-page companion that walks through every equation and maps it to the exact R variable in V2.R |

The key mathematical concepts are:

- **Mean-field CAVI:** The variational posterior factorises as q(L,F,β) = q(L)q(F)q(β), enabling coordinate-wise updates.
- **EBNM sub-problems:** Each coordinate update has the form: given pseudo-observations x with noise s, estimate a sparse prior g and return posterior moments. Solved via the `ebnm` R package.
- **Cox Taylor expansion:** The non-conjugate Cox likelihood is linearised via a 2nd-order Taylor expansion into a weighted Gaussian form, enabling EBNM updates for L and β.
- **Error-in-variables correction:** The β precision uses E[l²] (not l̄²), preventing survival coefficients from overfitting to uncertain loadings.
- **Gauss-Seidel updates:** Within each CAVI iteration, updates to earlier factors are immediately used when computing later ones — equivalent to block coordinate descent with sequential incorporation of new information.

---

## Version History

| Version | File | Status | Notes |
|---------|------|--------|-------|
| V1 | `code/legacy/Supervised_Bayesian_MF.R` | Archived | Original implementation; 6 known algorithmic issues |
| V2 | `code/Supervised_Bayesian_MF_V2.R` | Reference | Monolithic; all V1 issues corrected (A1–A6); kept for comparison |
| Modular | `code/fit_modular.R` + `update_*.R` | ✅ **Current** | Factor-wise Gauss-Seidel CAVI; tested (`Rscript tests/run_tests.R`); recommended for all new work |

**V2 improvements over V1:**

| ID | Fix |
|----|-----|
| A1 | True Gauss-Seidel CAVI (z_no_k from current EL/EBeta inside k-loop) |
| A2 | Orthogonalisation behind `orthogonalize=FALSE` flag (default off) |
| A3 | Numerical floors `pmax(..., 1e-10)` on all EBNM precision inputs |
| A4 | Dual convergence: both ΔL and Δβ must fall below tolerance |
| A5 | ELBO proxy (genomics log-likelihood) tracked per iteration |
| A6 | `refresh_taylor` flag for optional per-factor Taylor recomputation |

---

## Project Status

The single-modality model is implemented, tested, and benchmarked. The recommended
configuration uses the YFB linear predictor (η = (YF)β) with DeSurv-aligned gene
selection (genes ranked jointly by mean expression and variance within each platform,
top 3,000 per cohort, about 2,064 after intersection). Trained on 273 PDAC patients
(TCGA + CPTAC), it identifies two active prognostic programs, one associated with worse
survival and one with better survival, with mean external concordance 0.627 across five
independent PDAC cohorts (RNA-seq, microarray, proteomics).

Test suite (2026-10-01, `codex/multimodal-yfb`): 509 passing.

**Current work** (prioritized in `docs/plans/Working_Plan_10_1_26.md`; see `ROADMAP.md`):
the multimodal YFB extension, matched TCGA/ICGC preprocessing, and whether survival
supervision can be made to inform the factor matrix F (the `alpha_F` question).
Full history: [`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd); decision rationale: `DECISIONS.md`.

---

## License

MIT (see `LICENSE`). PDAC and TCGA/ICGC data are not included and keep their own terms.

---

## Author

Andrew Walther — May 2026

## Chapter 4 proposal

The canonical chapter is `paper/prelim/project3-ssbmf.qmd`. Its installed
post-commit hook renders and commits the generated Chapter 4 copy in
bios-dissertation; rerun that repository's `prelim/tools/build_prelim.py`
to update the full prelim. Continuing work is described as methodological
aims and expected results without a fixed-date timeline. A substantive
revision will follow the advisor-meeting plan and verified additional results.

## Protected-pruning evaluation

`results/benchmark_sim/run_k_selection_panel.R` contains the fixed simulation design;
`run_k_selection_real_1se.R` is conditional on a positive adoption record in
`outputs/k_selection/adoption.csv`. Run from the repository root with `Rscript`.
The assessment is blocked on the difference between multimodal output convergence and the
required ELBO stopping criterion; protection and a new starting rank were not adopted.
See the 2026-10-15 decision and progress-book entry before running these scripts.

## Planted comparator experiment

Run `Rscript results/benchmark_sim/run_supervision_comparators.R 6` from the root for the
fixed 40-dataset socket-cluster design (flashier, glmnet, ebnm, survival, yaml required).
The first run stopped at high-noise EBMF fits retaining zero or one factor, which cannot
support the requested lasso arm. No workaround was applied. Review `partial_summary.csv`,
`dataset_failures.csv` and `results.csv` in `outputs/supervision_comparators/` under the
benchmark directory; the partial high-noise comparison is conditional on successful fits.
Undefined C values and nonconvergence remain explicit.

## Multimodal comparator audit

Run `Rscript results/multimodal_real/run_supervised_multimodal_comparison.R` from the
repository root to audit cached retained ranks and fit the fixed K = 7 comparison. It uses
`data/processed_multiomics_tcga_icgc.rds`, writes each completed fit and comparison table to
`results/multimodal_real/outputs/supervised/`, and stops on a fit error. This is a potentially
long serial real-data run. Joint fits use raw log2 expression and beta-scale methylation;
the selected expression-only arm uses platform-wise z-standardization on the same genes.
ICGC is standardized once across all donors, then primary-PDAC rows are selected by ID.
Resume reuses this experiment's own checkpoints; remove them before changing inputs or code.

The cached ceilings 6, 7 and 8 all retain the same six-factor EBMF fit, explaining the
identical two-step C-index. The existing real-data reporting scripts now expose the
retained two-step rank instead of reporting it as unavailable. No fixed-rank constraint
has been added to the greedy EBMF procedure.

## Supervised robustness and comparator follow-up (2026-10-15)

The consolidated [five-prompt report](docs/reports/Codex_M6_M7_Summary_10_15_26.md)
records separate branches, test counts, tables and limitations. Protected pruning
was not adopted; K_init=7 and the handoff configurations remain controlling.
On `m6/robustness`, run `Rscript results/benchmark_sim/run_supervised_robustness.R 6`
from the repository root with the existing local PDAC data available. Outputs go to
`results/benchmark_sim/outputs/supervised_robustness/`: all simulated starts,
ELBO-selected starts, convergence, preserved old/new training CV, and feature-count
sensitivity. The unchanged original CV runner also refreshes its own CSVs. Existing
backup CSVs are protected from overwrite; archive a prior output folder before
rerunning. These experiments do not score external cohorts. Run
`Rscript tests/run_tests.R` for the full suite (550 passing on this branch).

The recorded D experiment saved all tables but exited 1 after a trailing parse
error; see the report and `execution_status.md` before treating it as a clean run.

## Round 2 integration

`round2/supervised` integrates the five Round 1 branches and passes 559 tests.
See `docs/plans/Codex_to_Claude_Round2.md` for current progress and
`docs/plans/Round2_attempts_log.md` for designs committed before experiments.
The approved benchmark snapshot and tag remain unchanged.

The Round 2 K panel runs with
`Rscript results/benchmark_sim/run_k_selection_panel.R 6`; it refuses to overwrite
an existing attempt and writes per-cell checkpoints beneath
`results/benchmark_sim/outputs/k_selection/round2/attempt3/`. Its fixed design and
adoption criteria are recorded in the attempt history before fitting.

To observe the unchanged single-modality pruning transitions, run
`Rscript results/benchmark_sim/run_single_pruning_trace.R` after its preregistration
commit. It records the first prescribed seed at both noise levels in a fresh
`outputs/k_selection/round2/attempt4_trace/` folder and refuses overwrite.
The local 800-cell rank panel can also be dispatched one cell per Longleaf array
with `results/benchmark_sim/round2_k_selection.slurm`; no jobs are submitted by R.

The approved full-posterior refit experiment uses
`Rscript results/benchmark_sim/run_posterior_k_panel.R 6`. It first requires the
complete original 800-cell table, then fits the 160 protected single-modality
cells in a fresh `attempt5_posterior/` folder and compares them with the saved
reference and mean-only rows. It changes only the pruning-refit initialization.

For this opt-in experiment, `fit_cox_on_yf(..., canonicalize=TRUE,
prune_preserve_posterior=TRUE)` returns `posterior_state` and carries retained
states into pruning refits. The state contains n-by-K score moments, p-by-K
loading moments, K beta moments, p feature precisions, persistent loading priors
and the shared beta prior. `yf_subset_posterior_state()` keeps the same program
indices in every field. `posterior_init=state` resumes from these quantities;
default pruning refits still initialize from L/F means. Preserving second moments matters
because projection uncertainty enters the beta precision, and the retained beta
sets the first Cox working predictor. This path currently supports canonical
fits without cohort offsets or frozen-F iterations. It changes neither the
model nor its adoption criteria and is not an approved benchmark replacement.

The completed 160-cell full-state panel is not adopted: recovery improved but
recovery, convergence and rank requirements did not all pass. Consult the Round2
attempt history before another experiment; existing output folders refuse overwrite.


For the fixed TCGA-only collapse diagnosis, run
`Rscript results/multimodal_real/run_round2_collapse_diagnosis.R all` only after
its preregistration commit. The ten arms write separate checkpoints to
`results/multimodal_real/outputs/round2/diagnosis1/`; the existing TCGA-only cache
must be present at ignored `data/round2_tcga_training.rds`. `control$trace=TRUE`
observes multimodal beta priors/moments, pseudo-observations and convergence
quantities; `beta_prior_fix=TRUE` and `beta_prior="shared_normal"` are diagnostic
opt-ins. `trace_supervision=TRUE` observes Single's actual prior refits as well.
None changes defaults or the convergence criteria. See the attempt log for the
fixed contrasts and the Longleaf wrapper for individual cells.


DeSurv software preparation uses the isolated ignored library
`data/round2_R_library/`. Add it with
`.libPaths(c(normalizePath("data/round2_R_library"), .libPaths()))` before loading
DeSurv. The pinned package/paper SHAs and installation notes are in the Round2
handoff. `results/benchmark_sim/desurv_author_adapter.R` wraps the authors'
rank preprocessing and frozen top-gene X'W beta projection; its focused tests
run with `Rscript -e 'source("tests/test_desurv_author_adapter.R");report_results()'`.
Andrew deferred the head-to-head scientific comparison for review.

The prepared head-to-head entry is `results/benchmark_sim/run_desurv_h2h.R`:
`all 2` runs the reported-parameter protocol locally; `prepare` creates the
training cache for `round2_desurv_seeds.slurm`, and `finalize` requires all100
seed statuses before consensus/external reporting. `bo 30` is the full author
training-only search (`round2_desurv_full_bo.slurm`), using the unchanged pinned
paper `R/bo_helpers.R` copied to ignored `data/round2_desurv_paper_bo_helpers.R`.
The runner checks the package/loader/helper versions, preserves the benchmark,
and refuses existing outputs. Execute only after its preregistration commit
and the preceding Task2 stopping point.

The ten-arm TCGA diagnosis has completed; consult the Round2 handoff for
benchmark-linked results and the open remedy ladder. Andrew deferred Task3's
real-data DeSurv head-to-head for review. Prepared DeSurv commands above must
not be executed without renewed authorization; no real-data DeSurv fit or
external scoring has run. Task4 omits DeSurv as explicitly approved.

The Round2 comparator rerun is
`Rscript results/benchmark_sim/run_supervision_comparators.R 6`. It uses the
unchanged five methods with explicit EBMF rank0/rank1 rules and per-arm errors,
writing a fresh `results/benchmark_sim/outputs/supervision_comparators/round2/`
folder. Read the preregistered design before executing; old outputs are never
overwritten. DeSurv is omitted as Andrew approved.

Task2 R1's required simulation/CV evaluation is
`Rscript results/benchmark_sim/run_round2_fixed_prior.R 6`, writing fresh
`results/benchmark_sim/outputs/round2/fixed_prior/`. Read its preregistration
before running. It evaluates the fixed shared prior only, with80 simulations
and20 training-CV fits, preserving all failures and frozen benchmark rows.
`prepare`/`finalize` support `round2_fixed_prior.slurm`; do not submit jobs here.
This evaluates a remedy; the full-TCGA diagnostic failures remain necessary
adoption failures even if another table looks favorable.

The Round2 framework-fold diagnostic is
`Rscript results/multimodal_real/run_framework_fold_trace.R`. It repeats the
unchanged failed TCGA CV fold with observational tracing and checks exact
posterior/history identity against the saved fit; it does not change scoring.

The experimental opt-in `profile_cox_intercept=TRUE` profiles the arbitrary
constant from the frozen Cox working quadratic. It weighted-centers only the
survival design; reconstruction and raw prediction use the original matrices.
It supports per-factor Normal beta in the multimodal fitter, and canonical
per-program Normal beta without cohort/strata effects or pruning in Single.
Defaults are unchanged. Run its preregistered training-only evaluation with
`Rscript results/benchmark_sim/run_round2_profiled_cox.R 6`; `prepare` and
`finalize` support the92-task Longleaf wrapper. This is an unadopted experiment;
read the Round2 attempt log and profiled-curvature proposal first.

A prepared prior-consistency observer is
`results/benchmark_sim/run_round2_prior_consistency_trace.R`, with aggregate
point-Laplace moment helpers and focused tests. It has **not been run**.
Before executing, preregister and commit a separate diagnostic attempt in
`docs/plans/Round2_attempts_log.md`. It repeats one completed simulation and
requires exact saved-fit identity; its adjusted monitor is diagnostic only,
never a replacement convergence or adoption criterion.


The approved full-CpG descriptive study is
`Rscript results/multimodal_real/summarize_cpg_variability.R run`, after its
Task5-QC1 preregistration commit. It describes all raw TCGA probes on the
established training patients, preserving missingness and reporting per-probe
means/SDs, predeclared filter overlaps, beta/arcsine SD histograms and a
mean-versus-SD density plot. It does not impute, fit a model or choose a filter.
Fresh safe aggregates/plots go under
`results/multimodal_real/outputs/feature_sensitivity/cpg_qc/`; the probe inventory
RDS stays local/ignored. Input and benchmark identities/hashes are verified.


`code/feature_sensitivity.R` provides preparation and agreement utilities for
the planned fixed-K expression study. Source it, then use
`feature_count_manifest(cohort_matrices,genes,targets,baseline_top_n)` on
training log-scale matrices, or
`feature_program_agreement(old_fit,old_genes,new_fit,new_genes,old_Y,new_Y)`
on aligned Single-block YFB fits. The latter preserves dead-factor NA/unmatched
cases and uses shared-gene loading correlations and normalized projections.
These utilities do not fit models. See the approved plan and input preflight:
the frozen D4 activity rank4/7 makes its prescribed warm inverse singular;
Andrew approved cold starts for the primary curve; its prospective design must
still be committed before fitting.


Task2-5's bounded prior-consistency audit uses
`Rscript results/benchmark_sim/run_round2_prior_consistency_trace.R` only after
its prospective entry is committed. It repeats one saved Single simulation,
requires exact entire-fit identity, and records native returned moments beside
reconstructed prior-density calculations. The observer returns the original
update unchanged. An F-prior-only monitor adjustment is diagnostic, not a new
ELBO/convergence rule or a model remedy. Fresh output only; see the attempt log.


The advisor progress report is available as both
`docs/reports/Advisor_Progress_since_10_2_26.html` (self-contained) and `.pdf`.
Render HTML with RStudio's Quarto on PATH and
`QUARTO_PANDOC=/opt/homebrew/bin/pandoc quarto render docs/reports/Advisor_Progress_since_10_2_26.qmd --to html`.
The completed Task2-5 observation reproduced its saved fit exactly; aggregate
results are in the Round2 handoff. It did not establish a model repair.

The fixed-K expression study is `results/benchmark_sim/run_feature_expansion.R`:
run `Rscript results/benchmark_sim/run_feature_expansion.R prepare` once, then
replace `prepare` with `all` (serial24 fits) or a cell number1-24. It uses only
training cohorts and the frozen benchmark settings, with four nested feature
counts and five fixed conditional-CV folds plus full training. Its committed
prospective attempt controls execution. Local inputs/fits are ignored; safe
CSV tables preserve failures and include a labelled frozen training reference.
`feature_survival_metrics()` in `code/feature_sensitivity.R` keeps constant-risk
concordance undefined and uses the existing Cox partial-likelihood calculation.

Terminal feature-study CSVs preserve empty success diagnostics as text. R
allocation peak MB uses the gc table's final column across R versions. The
first fit's original invalid memory value is retained separately and omitted
from the corrected scientific row; fitting and scores are unchanged.

After all24 cells finish, run
`Rscript results/benchmark_sim/summarize_feature_expansion.R` for aggregate
performance, the fixed-K curve, and baseline/preceding-count program agreement.
It requires all cells and preserves undefined quantities; no model is selected.

The completed expression feature study has24/24 converged fits through17,900
shared genes; conditional CV means .608/.577/.609/.607. Safe tables and two
curves are in the cold_primary output folder. Dense nonzero loading means do
not establish gene-level selection; no feature count/default is adopted. The
advisor report is updated in both HTML and PDF.
