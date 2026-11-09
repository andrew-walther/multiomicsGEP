# Chapter 4 revamp: supervised YFB model as the centerpiece

Written 2026-11-05 from a bios-dissertation planning session (Round 3 prose review of
Chapters 1–3). Status: **chapter assembled for review (2026-11-06)**. The evidence table, checklist and
order of work were updated on 2026-11-06; section 7 records what changed and why. Follow-up prompt:
[`2026-11-05-chapter4-revamp-prompt.md`](2026-11-05-chapter4-revamp-prompt.md).

Chapter source: `paper/prelim/project3-ssbmf.qmd`. The post-commit hook syncs it into
`bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` on any branch. Edits go on
`fix/yfb-supervision`; mirroring to `main` is optional, since that branch will be merged later
(Andrew, 2026-11-05).

**Deadline.** Revamp done by about 2026-11-07, so Andrew has about a week to review and edit before
the written prelim is due on 2026-11-16. Any checklist item that can't be produced by then is
written as planned evaluation, not left as a gap. The oral exam slides (2026-12-07) may show newer
results than the chapter.

**When to start.** After Andrew finishes the current `fix/yfb-supervision` work (S-11 to S-13
line). ROADMAP.md lists this plan as the next item.

## 0. Andrew's notes on the current chapter (2026-11-05)

These are requirements for the revamp, from Andrew's read of the synced prelim PDF.

- **Title.** Change it to "Survival-Supervised Bayesian Matrix Factorization for Prognostic
  Molecular Program Discovery in Pancreatic Cancer" (YAML `title:` in `project3-ssbmf.qmd`; the
  prelim prints it in capitals). It currently reads "... for Prognostic Gene Expression Programs in
  Pancreatic Cancer". Also check the slide outline and the abstract fact sheet for the old title.
- **Format.** Match Chapters 2 and 3: Introduction / Methods / Results (simulation, then the real
  PDAC data application) / Discussion.
  - Blend completed work and remaining proposed work across the whole chapter. Each section says
    what has been done and what the full project will add.
  - The chapter must read as connected prose with transitions between paragraphs and sections,
    not as an outline or a set of separate paragraphs. Each paragraph should follow from the one
    before it.
  - Open with `# Introduction`, not the bold, indented "**Summary.**" paragraph at line 42.
    Chapters 2 and 3 have no summary section; its content moves into the Introduction.
- **No bold paragraph leads.** Bold text at the start of an indented paragraph looks out of
  place.
  - **Default fix:** plain text. The paragraph's first sentence names its topic, as in
    "Survival in the single-modality loadings" (§4.2.8), which becomes an ordinary paragraph about
    that topic.
  - **When to use a heading:** only where the content needs one for navigation. Use a numbered
    `###` heading, or an unnumbered one (`{.unnumbered}`) if a number would be noise.
  - **Current places to fix:**
    - §4.2.1 Single-modality model (line 74, "Initialization and convergence");
    - §4.2.2 Multimodal model (line 88, "The parsimony framework");
    - §4.2.7 Data and preprocessing (lines 129 and 131, "Single-modality application" and
      "Multimodal application");
    - §4.2.8 Remaining work (lines 141, 143 and 145).
  - §4.2.8's content moves into the Methods subsections for proposed work (section 2 below).
- **Figures.** Every figure must be clear and legible at printed size. Rebuild any figure that is
  fuzzy. Problems in the current PDF:
  - Fig 4.1 (`fig:km-risk`, `single_km_risk_external.png`, width 0.75): too small and fuzzy.
  - Fig 4.2 (`fig:programs`, `single_programs_heatmap.png`, width 0.42): far too small; gene
    names (y-axis) and program names (x-axis) can't be read.
  - Fig 4.3 (`fig:km-programs`, `single_km_programs_external.png`, width 0.8): a bit too small
    and hard to read.

  Rules for every new or rebuilt figure:
  - Write vector PDF from R (`ggsave(..., device = cairo_pdf)` or `pdf()`), or PNG at ≥ 300 dpi
    only when vector isn't practical.
  - Size the figure in R at its printed size, about 6.5 in wide for `\linewidth`, so that text is
    ≥ 8 pt after placement. Don't shrink a large canvas with `width=0.4\linewidth`.
  - Heatmaps with gene labels use the full text width and enough height for one readable label
    per row. Show fewer genes per program rather than shrinking the font.
  - Multi-panel KM plots: at most two or three panels per row at full width, with readable
    risk-table text.
  - Check each figure by viewing the rendered prelim PDF page, not only the image file.
  - **Final figure review:** before the chapter goes to Andrew, review every figure in the
    rendered prelim for legibility, size, sharpness, axis and legend labels, and caption match.
    Most figures will be new (supervised fit), but the review covers all of them, kept or new.
    Record the result in the checklist.

## 1. Why the chapter changes

The current chapter (synced 2026-10-30, about 4,860 words) presents the frozen 2026-10-14 fit with
α_F = 0, in which survival never enters the loading update. DECISIONS 2026-10-02 says α_F = 0
can't be the final method. As of DECISIONS 2026-11-04, the preferred single-modality model is the
survival-weighted YFB fit:

- **Fitter:** `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1, prior_LF = "point_laplace",
  prior_beta = "normal")` (`code/fit_ssbmf_yfb.R`; runner
  `results/benchmark_sim/run_yfb_supervised_real.R`, lines 59–61). The code default for
  `prior_LF` is point-exponential, so the chapter must name **point-Laplace** for L and F and a
  **Normal** prior for β. The current chapter and the slide outline say otherwise.
- **Predictor:** η = (Y E[F])β. Under canonicalization each program is rescaled so that
  sd(Y F_k) = 1.
- **Objective:** E_q log p(Y | L, F, τ) + w_S × the Cox partial log-likelihood re-expanded at
  E[η], with Var(η) from q(β) and q(F), − KL. w_S multiplies the Cox term in the F update, the β
  update and the ELBO. Under canonicalization w_S = 0 also gives β = 0, so the unsupervised
  comparator is the w_S = 0 loadings followed by a separate Cox fit ("two-step", `w0cox`).
- **Rank and data:** K = 7 fixed, no pruning. D4 gene set (2,218 genes, single log, per-cohort
  z-scores), no survival screen. Training is TCGA-PAAD + CPTAC (273 patients, 139 events). w_S is
  chosen by 5-fold training CV C.

**Evidence available now** (every number from the findings files named):

| Result | Value | Source |
|---|---|---|
| S-9 simulation, 2,000 genes | planted program recovered in 10/10 datasets in every variance cell at w_S = 1; CV picks w_S = 1 in 55/60; every w_S ≥ 40 fails to converge | `docs/findings/2026-11-04-survival-weight-p2000-findings.md` |
| S-9 calibration (corrected 11-05) | slope 1.45–3.23 at 2,000 genes (attenuated, not overconfident); w_S = 10 closer to calibrated at larger variances; w_S ≥ 100 overfits | same, and DECISIONS 2026-11-04 (correction) |
| S-11 training CV, D4 | C 0.663 (SE 0.024) vs two-step 0.652 (SE 0.026); paired +0.011 (SE 0.005); 3 survival-active programs; 2 of 5 folds stop at the 300-sweep cap | `docs/findings/2026-11-04-supervised-real-data-findings.md` |
| S-11 caveat | survival supplies 10⁻⁸–10⁻⁶ of the loading precision; w_S = 1 and w_S = 0 loadings agree at \|r\| = 1.000 | same |
| S-11 size and screen | D4 preferred; survival screening adds no clear gain; null check passes (false-active 0.1–0.5 per cell) | same |
| S-12 external, 5 cohorts | C 0.584 / 0.537 / 0.656 / 0.681 / 0.644 (mean 0.620); paired gain over two-step +0.007 (95% CI −0.002 to 0.017); calibration slopes 0.36–1.57 (Puleo 0.99), closer to 1 than two-step in every cohort | `docs/findings/2026-11-05-supervised-external-scoring-findings.md` |
| S-12 pooled HR (added 11-05) | 1.37 (1.24–1.51) per SD for w_S = 1, 1.35 (1.22–1.49) for two-step; DeSurv reports 1.50 (1.31–1.72), not paired | same, `external_pooled_hr.csv` |
| S-13 supervision strength (11-05) | CV partial log-likelihood picks w_S 1–10, rising with signal (1–2 weak and 0.25% cells, 5 at 1%, 10 at 5%); never 20 or 50, where fits stop converging. At the chosen weight recovery equals the matched two-step fit (10/10 everywhere); C equal or +0.003 to +0.007. Slopes 1.24–1.45 vs two-step 0.96–1.16. Neither arm (A1 block weight, A2 joint q(β)) passes every registered criterion; A2 adds nothing over A1 | `docs/findings/2026-11-05-supervision-strength-findings.md`, DECISIONS 2026-11-05 |
| S-14 near the detection threshold (11-05) | planted share 0.05%, strong survival: w_S = 1 recovers 5/10 vs 0/10 for two-step, validation C +0.156; 29 of 30 strong datasets favor supervision. Primary recovery criterion fails (false-active 2.5, shared by the two-step fit at 2.8). Realistic survival (oracle C ≈ 0.74): gains +0.002 to +0.038. Fixed w_S = 5 overfits (slopes 0.04–0.65) | `docs/findings/2026-11-05-near-threshold-findings.md` |
| False-active contribution (11-06, post hoc) | where the program is recovered, other active programs carry ≤ 1.3% of risk-score variance and change C by ≤ 0.008; at 0.05% share they carry 9–69% (pieces of the planted direction; 69% two-step vs 8.7% supervised, strong survival) | both S-13 and S-14 findings, added sections |
| October 15 supervised real-data check (issue #6 item 4) | pre-registered ΔC vs the α_F = 0 model −0.003 (−0.029 to 0.024); margin −0.02 not met for lack of precision | DECISIONS 2026-10-15 (lines ~1040–1080) |
| Fair historical comparator (issue #6 item 5) | rank-7 EBMF then Cox, not the 40-factor LASSO pipeline; in the new chapter the matched comparator is the w_S = 0 two-step fit at K = 7 | issue #6; S-11/S-12 |
| Derivation (11-06) | unified YFB derivation with w_S, Cox surrogate, q(β), q(L), q(F), balance of precisions, M-modality extension; 73 numerical checks | `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.{tex,pdf}`, `verify_derivation.R` |
| LB vs YFB at matched weights (11-06, exploratory) | YFB w_S = 1 equals its two-step fit (test C 0.783 vs 0.782, amplitude 0.5); LB at α = 0.5 gains recovery over its two-step (0.828 vs 0.664) but stays below YFB in C | `docs/findings/2026-11-06-lb-yfb-matched-comparison-findings.md` |
| S-15 real-data refit with CV-chosen w_S (registered 11-06, running) | adoption rule: w* > 1, sound full fit, nested-CV C ≥ 0.639; if adopted, one external scoring with the S-12 procedure | `docs/plans/2026-11-06-supervised-real-refit-plan.md`; results pending in `yfb_supervised_refit_2026-11-06/` |

**What the chapter must not claim.**
- That survival shapes the real-data programs. It doesn't at w_S = 1. The gain comes from the joint
  β update with its shrinkage prior.
- That the external gain is established. The pooled interval includes zero.
- Any comparison with the 2026-10-14 benchmark (mean 0.627) as paired. That benchmark used an older
  pipeline.

Making survival move the loadings (a data-driven supervision strength) is the main open method
question. It belongs in the chapter as proposed methodology.

**Coordination.** The S-11 session (configuration review, external scoring) is still running.
Before producing anything, the follow-up reads the latest `DECISIONS.md`, `docs/findings/` and
`git log`. It doesn't rerun, retune or rescore what that session owns, and it doesn't touch its
output folders (`yfb_supervised_real_2026-11-04/`, `yfb_supervised_external_2026-11-05/`). New
outputs go to new dated folders.

## 2. Target outline

The outline matches Chapters 2 and 3: Introduction, Methods, Results (simulation, then the PDAC
application), Discussion. Proposed work sits inside the section it belongs to, not only in Future
Directions. Headings are descriptive, as in Chapter 3; the current Results headings are claims.

There are no bold run-in paragraph openers. The current chapter has 8, at lines 42, 74, 88, 129,
131, 141, 143 and 145. Use a `###` subheading where the content needs one, otherwise a plain topic
sentence. The bold "Summary." paragraph at line 42 goes: Chapters 2 and 3 open directly with the
Introduction, and its content moves into the Introduction's last paragraph.

Folder rule (`paper/prelim/AGENTS.md`): continuing work is described as additional methodology
and expected evidence, with no fixed-date timelines. Proposed methods and expected results are never
presented as findings.

1. **Introduction**
   - The problem: prognostic expression programs in PDAC, and the two-step pipeline (unsupervised
     factorization, then Cox on plug-in scores).
   - The gap: survival doesn't inform discovery, and score uncertainty is ignored.
   - What SSBMF does.
   - Aims and outline, stated as completed aims and proposed aims. *Proposed in this section:* the
     full project's aims (data-driven supervision strength, M-modality model, rank selection).
2. **Methods**
   1. Data and notation; the YFB model (Y = LFᵀ + E; η = (YF)β; Cox partial likelihood).
   2. Priors and empirical Bayes: point-Laplace on L and F, Normal on β, ebnm fits.
   3. Survival weight w_S and canonicalization: definition, and why canonicalization makes w_S
      comparable across programs. Also that w_S = 0 gives β = 0 and how the two-step comparator
      is built.
   4. Variational inference and the ELBO: sweeps, Cox re-expansion at E[η] with Var(η), stopping
      rule (tolerance 10⁻⁵, at most 300 sweeps).
   5. Choosing w_S by training CV; prediction in new cohorts (per-cohort z-scores, frozen
      orientation).
   6. Simulation design (S-9): planted low-variance prognostic program, 50 and 2,000 genes,
      variance cells, w_S grid, recovery, false activity, calibration, CV choice. Table from
      `yfb_survival_weight_2026-11-04/grid/design.csv`.
   7. PDAC data and preprocessing: training and external cohorts, single log, D4 gene set,
      evaluation metrics (C, paired bootstrap, calibration slope).
   8. *Proposed:* **Data-driven supervision strength.** Survival barely moves the real loadings at
      w_S = 1. Candidates come from DECISIONS 2026-11-04: a calibrated likelihood power
      (SafeBayes or Holmes–Walker tempering), a downweighted genomics term (S-10), or normalized
      loss terms as in DeSurv. This subsection also states the acceptance criteria.
   9. *Proposed:* **M-modality model with modality-specific priors on F.** Y_m = L F_mᵀ + E_m for
      m = 1..M, with a shared L, a program- and modality-specific point-Laplace prior g_{F_mk}
      fitted by empirical Bayes, and η from the stacked projections. M = 1 is the
      single-modality model. Expression and methylation are the planned case. Source:
      `derivations/SSBMF-YFB-multimodal-model-derivations/multimodal_YFB_derivation.qmd` (+ technical notes). It's
      written as derived but not yet implemented in this form.
   10. *Proposed:* **Rank selection.** K is fixed at 7 now. Evidence-based pruning (flashier-style)
       under the supervised objective, and the replicated K-selection criterion.
   11. *Proposed:* **Comparators still to run:** DeSurv head-to-head (`results/benchmark_sim/outputs/desurv_h2h/BLOCKED.md`;
       deferred 10-15), and the LB (η = Lβ) vs YFB comparison planned in DECISIONS 2026-11-04.
3. **Results**
   1. Simulation: recovery, false activity, calibration and convergence across w_S, and the CV
      choice of w_S (S-9). *Planned evaluation:* the same grid under the data-driven supervision
      variants.
   2. Training CV and feature-set size (S-11): CV C against the two-step fit, the size curve,
      screening.
   3. External validation (S-12): five-cohort C, paired difference, calibration slopes; KM by
      risk tertile.
   4. Programs and biology: the 7 programs, the 3 survival-active ones, top genes, enrichment,
      subtype and stroma concordance, and correspondence with the α_F = 0 programs.
      *Planned evaluation:* whether stronger supervision recovers programs the w_S = 1 fit misses.
   5. *Planned evaluation:* M-modality and rank selection. Expected evidence and decision rules
      only, with no October multimodal numbers.
4. **Discussion**
   - What the supervised fit does: better-calibrated β, a small gain in discrimination. What it
     doesn't do: move the loadings.
   - Limitations: K fixed; sweep cap; external gain interval includes zero; tertile risk groups
     are within-cohort.
   - Continuing work, tied back to Methods 2.8–2.11, plus transferable risk classifications
     (ROADMAP 2026-10-16) and the R package.

**Material leaving the chapter** goes to `paper/prelim/notes/trimmed-material-2026-11.md`, dated,
with the reason:
- the α_F = 0 framing and its results (tab:external 0.627 vs 0.601, fig:km-risk, fig:km-programs,
  fig:programs, fig:delta-c, fig:percohort), unless a supervised counterpart replaces them;
- the October multimodal parsimony results (fig:mm-pruning, tab:mm-external, the "Multimodal
  model" methods text). Multimodal is proposed work (Andrew, 2026-11-05).

The 2026-10-14 benchmark may still be cited once as the historical reference. It is labelled as
such, per DECISIONS 2026-11-04 "Presentation".

## 3. Results checklist (cross-reference with Andrew's in-progress work)

The completion column gives the saved output or the planned evaluation. Source and remaining columns record the scope at registration. Paths are relative to the multiomicsGEP root.

### Tables

| # | Item | Replaces | Status | Data source | Script | Remaining | Andrew started? | Done |
|---|---|---|---|---|---|---|---|---|
| T1 | S-9 simulation design | tab:sim-design | partial | `results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid/design.csv` | typed from design.csv | write LaTeX table | | [x] `paper/prelim/tables/tab_sim_design.tex` (embedded in the chapter) |
| T2 | Training CV: C, SE, paired difference, folds converged, survival-active count | (new) | partial | `results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/summary.csv`, `fits.csv` | `results/benchmark_sim/summarize_yfb_supervised_real.R` | LaTeX table | | [x] `paper/prelim/tables/tab_training_cv.tex` (embedded in the chapter) |
| T3 | External C in 5 cohorts, two-step C, paired difference, calibration slope | tab:external | partial | `results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/external_cindex.csv`, `external_paired_diff.csv` | `results/benchmark_sim/score_yfb_supervised_external.R` (S-11 session; don't rerun) | LaTeX table; historical column labelled not paired | | [x] `paper/prelim/tables/tab_external.tex` (embedded in the chapter) |
| T4 | Programs: β, SD, z, survival-active, top genes at each pole, short label | (in-text list) | partial | `results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/programs_D4.csv` | — | table + biological labels (needs B1–B3) | | [x] `paper/prelim/tables/tab_programs.tex` (embedded in the chapter) |
| T6 | S-13/S-14 supervision-strength summary: CV-chosen w_S, recovery and C against the matched two-step fit, near-threshold cells | (new, 11-06) | partial | `yfb_supervision_strength_2026-11-05/selected_vs_twostep.csv`, `yfb_near_threshold_2026-11-05/metrics.csv` | — | LaTeX table | | [x] `paper/prelim/tables/tab_supervision_strength.tex` (embedded in the chapter) |
| T5 | Enrichment summary for the survival-active programs | (new) | missing | all-patient fit `fits/size_D4_N2218_w1_f0.rds` (`state$EF`, `genes`) | adapt `code/pathway_enrichment.R` (written for the old D4 fit) | run, then table | | [x] `paper/prelim/tables/tab_enrichment.tex` (embedded in the chapter) |

### Figures

| # | Item | Replaces | Status | Data source | Script | Remaining | Andrew started? | Done |
|---|---|---|---|---|---|---|---|---|
| F1 | Simulation: recovery, false activity, survival share, calibration, convergence and C against w_S | fig:delta-c | exists (progress book) | `yfb_survival_weight_2026-11-04/grid/metrics.csv`, `cal_slope_exact.csv` | `docs/progress_book/figures/2026-11-06/survival-weight-curves.png` (written 11-05 04:00, after the slope correction; its plotting script wasn't found by name, so locate it) | confirm it uses `cal_slope_exact.csv`; restyle for the chapter | | [x] `paper/prelim/figures/sup_sim_weight_s9.pdf` |
| F2 | CV choice of w_S across datasets | (new) | partial | `grid/cv_chosen.csv`, `cv_chosen_cells.csv` | — | small figure or fold into T1/text | | [x] weight frequencies reported in T1 and simulation results |
| F3 | Real-data CV C against feature-set size | (new) | exists | `yfb_supervised_real_2026-11-04/size_elbow.png` (= progress book `real-data-size-curve.png`) | `summarize_yfb_supervised_real.R` | restyle for the chapter | | [x] `paper/prelim/figures/sup_real_size_curve.pdf` |
| F4 | External C by cohort with intervals, w_S = 1 vs two-step | tab:external visual | partial | `external_cindex.csv` | progress book `external-concordance.png` is the **October** fit | new plot from S-12 CSVs | | [x] `paper/prelim/figures/sup_external_c.pdf` |
| F5 | KM by risk tertile in external cohorts | fig:km-risk | missing | needs per-patient external risk scores; S-12 saved only aggregates | extend `score_yfb_supervised_external.R` to save risks (same fixed model, no refit) | approved 11-06: per-patient scores saved under `data/` (git-ignored) only; aggregates in Git | | [x] `paper/prelim/figures/sup_external_km_risk.pdf` |
| F6 | KM by program-score tertile for the survival-active programs | fig:km-programs | missing | same as F5 | same | approved 11-06 (as F5) | | [x] `paper/prelim/figures/sup_external_km_program_2.pdf`, `paper/prelim/figures/sup_external_km_program_3.pdf`, `paper/prelim/figures/sup_external_km_program_5.pdf` |
| F7 | Program heatmap (loadings of top genes, 7 programs) | fig:programs | partial | `fits/size_D4_N2218_w1_f0.rds` | adapt `paper/prelim/make_figures.R` (uses `load_d4_weights()` for the old fit) | new figure | | [x] `paper/prelim/figures/sup_programs_heatmap.pdf` |
| F8 | Correspondence of w_S = 1 programs with the α_F = 0 benchmark programs (and with w_S = 0) | (new) | partial | f0 fits for w1 and w0; `results/benchmark_snapshot_2026-10-14/single_modality/` | progress book `program-correspondence.png` is the October fit | new figure or table; links the new programs to the old Programs 3/7 | | [x] `paper/prelim/figures/sup_programs_correspondence.pdf` |
| F9 | Calibration: external slopes, w_S = 1 vs two-step | (new) | partial | `external_cindex.csv` (`cal_slope`, `cal_se`) | — | small dot plot or a T3 column | | [x] `paper/prelim/figures/sup_external_c.pdf` (calibration panel) |
| F10 | Per-cohort β (TCGA vs CPTAC) for the survival-active programs | fig:percohort | missing | needs per-training-cohort fits | — | moved to planned evaluation (11-06): needs per-cohort refits, outside the chapter's scope | | [x] planned evaluation |
| F12 | S-13 weight curves: CV choice of w_S, C and calibration against w_S, with the two-step reference | (new, 11-06) | partial | `yfb_supervision_strength_2026-11-05/weight_curves.csv`, `cv_selected.csv` | progress-book `s13-weight-curves.png` from the parent session's untracked `plot_supervision_strength.R` (do not edit; restyle in the prelim script) | vector PDF at print size | | [x] `paper/prelim/figures/sup_sim_weight_s13.pdf` |
| F13 | S-14 near-threshold: recovery and paired validation C, w_S = 1 vs two-step, by planted share | (new, 11-06) | partial | `yfb_near_threshold_2026-11-05/metrics.csv` | progress-book `s14-paired-c.png` (same script) | vector PDF at print size | | [x] `paper/prelim/figures/sup_sim_near_threshold.pdf` |
| F11 | Model diagram (Y, L, F, β, w_S) | (new) | exists | `docs/progress_book/figures/two-step-vs-joint.svg`, `supervision-states.svg` | `build_diagrams.py` | optional; convert to PDF/PNG if used | | [x] optional diagram omitted. Methods gives the predictor and objective |

### Biology

| # | Item | Status | Source / method | Remaining | Andrew started? | Done |
|---|---|---|---|---|---|---|
| B1 | Top genes and pole labels for the 3 survival-active programs (hypoxia/glycolysis vs normal pancreas; ribosomal/proliferation vs classical epithelial; third adverse program) | partial | `programs_D4.csv`; the findings file's gene lists | confirm labels against B2–B3 | | [x] `paper/prelim/tables/tab_programs.tex` (embedded in the chapter) |
| B2 | Pathway enrichment (MSigDB Hallmark, Reactome/KEGG via fgsea or clusterProfiler) | missing | adapt `code/pathway_enrichment.R` to `state$EF` of the f0 fit | run and summarize | | [x] `results/benchmark_sim/outputs/prelim_ch4_2026-11-06/programs_enrichment.csv` |
| B3 | Subtype and stroma concordance (Moffitt basal-like/classical, Moffitt normal/activated stroma, Collisson, Bailey, PurIST) | missing | correlate program scores with published signature scores in training and external cohorts | find or write the signature-scoring code; no repo code found by grep for "PurIST" | | [x] `results/benchmark_sim/outputs/prelim_ch4_2026-11-06/programs_subtype_concordance.csv` |
| B4 | Correspondence with the α_F = 0 Programs 3 (protective) and 7 (adverse) | partial | F8 | report \|r\| and shared top genes | | [x] `results/benchmark_sim/outputs/prelim_ch4_2026-11-06/programs_correspondence_topgenes.csv` |
| B5 | Direction of the survival-active programs out of sample | missing | per-cohort Cox on program scores (from F5/F6 data) | needs per-patient external scores | | [x] `results/benchmark_sim/outputs/prelim_ch4_2026-11-06/external_program_direction.csv` |

## 4. Other documents the revamp affects

- **Chapter 1 (bios-dissertation `prelim/literature-review/draft/literature-review-draft.qmd`).**
  These passages describe Chapter 4 in α_F = 0 terms and must be rewritten with the revamp. Round
  3 (2026-11-05) fixed only their prose; exact lines are in
  `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/ch1/` ("Passages for the
  Chapter 4 revamp"):
  - preamble paragraph 3, "The third project, presented in Chapter 4 ...";
  - §1.3.4 "Joint Survival-Supervised Factorization and Its Gaps", closing paragraph;
  - §1.3.6 rank-selection paragraph citing Young2026;
  - §1.3 Summary, closing paragraph ("Chapter 4 presents this framework: a single-modality model
    with external validation ... and a multimodal extension ... with first results on real data").
    The multimodal part is now proposed work.
  - Two optional structure changes deferred from Round 3:
    - fold the end of §1.3.5 (the rank paragraph near line 809) into §1.3.6, which repeats it;
    - add one opening sentence on outcome types to §1.3.8, whose method paragraphs share one
      template.
  - The file's header comment ("Section 1.3 does not describe SSBMF's current implementation
    state") still applies. The literature-review framing about recovering low-variance programs
    must match what the supervised fit shows: survival doesn't move the real loadings yet.
- **Prelim abstract** (`bios-dissertation/prelim/abstract/abstract-text.qmd`; claims sourced in
  `notes/abstract-fact-sheet.md`). The standalone abstract was uploaded 2026-11-02 and fills pages
  iii–iv exactly. Its Chapter 4 sentences (mean 0.627, "slightly higher" than the two-step,
  multimodal first results) no longer describe the chapter. **Andrew decides** whether the
  prelim's copy changes. Any change needs an equal cut and a fact-sheet update.
- **Slides** (`presentation/prelim-slides/OUTLINE.md`, Draft 1, 2026-11-03):
  - the backup slides B1–B9 on the supervised fit move into the main deck;
  - the β prior is Normal, not point-normal;
  - the slides may go beyond the chapter.
- **bios-dissertation full prelim:** `prelim/full-draft/assembly_manifest.json` lists the Chapter 4
  figures. Update it for the new figures, then rebuild (`build_prelim.py`,
  `test_prelim_assembly.py`).
- **ROADMAP / DECISIONS** in both repos: record the revamp, the trimmed material, and the
  supervised model as the chapter's centerpiece.
- **Round 3 human-writer review of Chapter 4** runs last, on the revamped chapter. Use the round 3
  brief in `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/WRITER-BRIEF.md`.

## 5. Decisions

**Settled by Andrew (2026-11-05):**
- **S-11 session:** it is still running. The revamp works around it (section 1).
- **Branch:** `fix/yfb-supervision`; mirroring to `main` is optional.
- **Multimodal:** proposed work, described as planned methodology with modality-specific priors on
  F. The October multimodal results leave the chapter.
- **Scope of the chapter:** it shows the state at writing time; the slides may go further.

**Settled by Andrew (2026-11-06, relayed from the parent session in his words):**
- Plan approval is given in advance. Stop only if the plan departs from section 0 or the evidence
  changes the chapter's main claim.
- Per-patient external risk scores may be saved under `data/` (git-ignored), never committed.
- The prelim abstract may change, limited to the supervision sentence or sentences.
- The finished chapter gets a human-writer review by one or two independent reviewers who did not
  write it.

**For Andrew during the follow-up (answered 2026-11-06 as above; F10 moved to planned evaluation):**
1. F5/F6/B5 need per-patient external risk scores from the fixed model. Is saving them (no refit,
   no retuning) consistent with the "score once" rule, and should the S-11 session produce them?
2. Should the prelim abstract copy change (section 4)?
3. Should F10 (per-cohort β) stay in scope, or move to planned evaluation?

**Checked 2026-11-05:**
- **Prior wording:** point-Laplace on L and F, Normal on β (S-11 runner lines 59–61).
- **Score uncertainty:** the old chapter's aim "score uncertainty in the single-modality
  coefficients" may be partly met already, since the YFB objective uses Var(η) from q(F). The
  follow-up confirms this in the β update before writing it either way.

## 6. Order of work (target: done by about 2026-11-07)

1. Read the latest S-11 session state (section 1). Settle decisions 1–3 with Andrew.
2. Produce T1–T5, F1–F9 and B1–B5 into `paper/prelim/figures/` and a new
   `paper/prelim/make_figures_supervised.R`. Use new dated output folders for anything computed.
   Tick the checklist as each item lands.
3. Move outgoing text to `paper/prelim/notes/trimmed-material-2026-11.md`.
4. Rewrite the chapter to section 2, using the human-writer skill while drafting. Every number
   traces to a file in section 1 or the checklist.
5. Commit the chapter on `fix/yfb-supervision`. Let the hook sync, then rebuild the full prelim
   in bios-dissertation and compare page counts.
6. Rewrite the Chapter 1 passages (section 4) in bios-dissertation.
7. Run the Round 3 review on Chapter 4 (writer, independent reviewer, apply, report).
8. Update ROADMAP/DECISIONS in both repos.

## 7. Update 2026-11-06: what changed since this plan was written

The `fix/yfb-supervision` work after 2026-11-05 changes the evidence, not the chapter's main
claim. The claim stays: the supervised YFB fit predicts as well as or slightly better than the
matched two-step fit, with better-calibrated coefficients, but at w_S = 1 survival does not move
the real-data loadings. What changed:

- **S-13** shows that the survival weight can be chosen from the data. Five-fold CV partial
  log-likelihood under the survival-block weight picks 1 to 10, rising with signal, and never a
  weight that breaks convergence. At the chosen weight supervision equals the matched two-step fit
  on the S-9 cells. "Choosing w_S by training CV" is now partly completed work, not only proposed
  methodology (outline 2.5 and 2.8). New items T6 and F12.
- **S-14** gives the first simulation where supervision changes discovery: near the detection
  threshold the supervised fit recovers programs the two-step fit misses and raises validation C.
  This is the chapter's evidence that supervised discovery is possible, with the qualification
  that the gain is small under realistic survival strength. New item F13.
- **False-active contribution** (post hoc): false-active programs barely affect the risk score
  where the program is recovered. The chapter reports contribution beside the count.
- **The Gate S reference** was the S-7 unsupervised fit on uncentered input. The matched
  two-step fit recovers the S-9 programs too, so the chapter compares against the matched
  two-step fit throughout, and describes S-9 as showing that supervision does not break
  recovery, not that it causes it.
- **S-15** (running) may replace w_S = 1 on the real data. If it adopts w* > 1, the chapter's
  real-data model is that fit and the external results come from its one scoring. Otherwise the
  S-11/S-12 fit stands. Chapter text and figures are written so either outcome drops in.
- **Issue #6 items 4 and 5** are now chapter requirements: the October 15 pre-registered null
  result is stated plainly, and the fair comparator is the rank-matched two-step fit.
- **Derivations** were reorganized; the unified YFB derivation
  (`derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.pdf`) is the Methods
  reference, including its M-modality section.
- **LB vs YFB** (exploratory) and **projection-only centering** (opt-in, nonnegative input) are
  development results. The chapter mentions the LB comparison only as proposed work.
- **Score uncertainty.** The YFB objective uses Var(η) from q(β) and q(F). The chapter states
  that the coefficient update accounts for loading uncertainty, after checking the β update in
  `code/fit_ssbmf_yfb.R` and the derivation.

## 8. Numbered plan as executed (2026-11-06)

1. Update this plan (sections 1, 3, 5, 7). Done 2026-11-06.
2. In parallel subagents, produce into `paper/prelim/figures/` with
   `paper/prelim/make_figures_supervised.R` (vector PDF, sized at print width):
   a. simulation and CV figures and tables from saved CSVs: T1, F1, F2, F12, F13, T6, T2, F3;
   b. external tables and figures from saved CSVs: T3, F4, F9;
   c. per-patient external scoring of the fixed fit with the S-12 procedure (scores in `data/`
      only), then F5, F6 and B5;
   d. programs and biology from the all-patient fit: T4, F7, F8/B4, B1, B2/T5, B3.
3. Move outgoing text to `paper/prelim/notes/trimmed-material-2026-11.md`.
4. Draft the chapter to the section 2 outline, with section 7's changes, using the human-writer
   skill. Every number traces to a file.
5. Rewrite the Chapter 1 passages and the abstract's supervision sentence in bios-dissertation.
6. When S-15 lands, put its result in (and, if adopted, its external scoring).
7. Commit on `fix/yfb-supervision` (hook syncs), update `assembly_manifest.json`, rebuild the
   prelim, run `test_prelim_assembly.py`, report page counts, review every figure on the
   rendered page.
8. Round 3 human-writer review: writer subagent with the round 3 brief, one or two independent
   reviewers, apply, report.
9. ROADMAP and DECISIONS in both repos.

## 9. Paused 2026-11-06 about 04:35 (usage limit), where to resume

Paused at Andrew's request until the usage reset. State:

- **Done:** step 1 (this plan update, commit `cc025b8`); step 3 (outgoing text in
  `paper/prelim/notes/trimmed-material-2026-11.md`); Chapter 1 passages rewritten in
  bios-dissertation (`prelim/literature-review/draft/literature-review-draft.qmd`, cut text in its
  `trimmed-material.md`), including Round 3 structural notes 1 and 2. The β update was checked:
  the canonical path uses E[z²] = z̄² + V^Z (`code/fit_ssbmf_yfb.R` line ~1106), so loading
  uncertainty does enter the coefficients.
- **Chapter draft:** full first draft of every section in
  `paper/prelim/notes/chapter4-revamp-draft-2026-11-06.qmd` (body only; the YAML header with the
  new title still has to be put on). Placeholders to fill: `%%TABLE-T1%%`, `%%TABLE-T6%%`,
  `%%TABLE-T2%%`, `%%TABLE-T3%%`, `%%S15%%`, `%%KM-RISK%%`, `%%PROGRAMS%%`. Not yet run through the
  human-writer pass or a number check.
- **Stopped mid-way (subagents killed, output unverified, not committed):**
  `paper/prelim/make_figures_supervised_sim.R`, `paper/prelim/make_figures_supervised_programs.R`,
  `results/benchmark_sim/save_yfb_external_patient_scores.R`, and
  `results/benchmark_sim/outputs/prelim_ch4_2026-11-06/`. Check, finish or redo step 2 from these.
- **Next:** step 2 (figures/tables, checklist ticks), fill the draft, abstract supervision sentence
  (only the final Chapter 4 sentence of `abstract-text.qmd`; other Chapter 4 abstract sentences
  are out of date too and need Andrew's ruling: "non-negative factorization", "independently
  tuned two-step pipeline", multimodal "selects the correct number of programs"), S-15 result,
  then steps 7–9.


## 10. Evidence carried into the revised chapter

The chapter uses the S-11 fit and the S-12 evaluation. S-15 retained the
unit survival weight and KL-2 supported the unpruned rank-seven fit. The
external interval is shown with and without donor resampling. The primary
slope table gives four of five series closer to calibration, with PACA-AU
RNA-seq as the exception. The nonmatched-program variance analysis is
post hoc and includes programs below the survival-activity threshold.
See the [numerical source ledger](../../paper/prelim/notes/chapter4-number-sources-2026-11-06.md).

F10 remains planned evaluation. F11 was optional and is omitted. Earlier
chapter material remains in
[`trimmed-material-2026-11.md`](../../paper/prelim/notes/trimmed-material-2026-11.md).

The assembled prelim passed all six tests in `test_prelim_assembly.py`.
The full PDF went from 212 to 227 pages. The Chapter 4 body went from
18 to 31 pages and its standalone preview from 23 to 36 pages. The
standalone abstract remains two pages. All 11 vector figures were reviewed
on the assembled pages for legibility and caption agreement. The Round 3
record is in `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/ch4/`.
