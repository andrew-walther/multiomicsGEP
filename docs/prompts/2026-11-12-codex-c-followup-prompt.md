# Codex prompt: task C follow-up (expression reruns and methylation factor characterization)

Date: 2026-11-12. Author: Andrew Walther.

Paste into one Codex session. Part 1 runs first. Part 2 uses Part 1's result for the
point-exponential expression fit. Andrew approved both parts on 2026-11-12, including a
1,000-sweep cap for Part 1 and installation of the Bioconductor 450K annotation package for
Part 2.

---

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or change
branches. Other sessions share this checkout and edit DECISIONS.md, PROJECT_STATUS.qmd, the
attempts log and the week plan. Before each edit to those files, `git pull` and re-read
the file. Make small targeted edits, commit and push straight away. If a push conflicts,
rebase and redo the edit. Never force-push, stash, reset or switch branches. Stage only
your own paths and never use `git add -A`. Write detailed commit messages with no tool
attribution.

Read these first:

- AGENTS.md
- the DECISIONS.md entry "Task C: methylation alone adds no prognostic program at K = 5"
- `docs/findings/2026-11-12-methylation-single-modality-findings.md`, especially
  "Rerun at K = 5" and "Would more iterations help?"
- `results/multimodal_real/run_methylation_single_modality.R`
- `results/multimodal_real/run_methylation_iteration_extension.R`
- the aggregate CSVs in
  `results/multimodal_real/outputs/methylation_single_modality_2026-11-12/K05/`

Everything here is at **K = 5**. Label K in every table, file name and sentence.
Patient-level matrices, fits and scores stay under the ignored
`data/methylation_single_modality_2026-11-12/K05/`. Commit only aggregate CSVs: they may
contain CpG IDs, gene symbols and correlations, but never patient IDs or patient-level
values. No ICGC scoring, and no adoption of any rank, prior, setup or CpG count.

## Part 1: rerun the 8 capped K = 5 expression fits at 1,000 sweeps

At K = 5, 8 of the 12 expression fits reached the 300-sweep cap:

- PL folds 1, 3 and 4;
- the PE full fit (fold 0);
- PE folds 1, 3, 4 and 5.

As a result, most paired methylation-minus-expression differences are computed against
capped expression states. The PE expression-correlation step of the methylation-only rule is
also unassessed.

1. Write `results/multimodal_real/run_expression_iteration_extension_k5.R`, following
   `run_methylation_iteration_extension.R`. Use a separate private directory
   (`.../K05/iteration1000/`) and copy the K = 5 frozen source and design into it, checking
   their hashes. Change **only** `cavi$max_iter` to 1000. `c_one()` resets the seed, so each
   fit restarts from the original initialization rather than a warm start. Select K = 5
   through the runner's `c_set_rank(5)` / `--K=5` mechanism. Do not edit
   `config/globals.yml` or any file under `code/`.
2. Before using a refit, check that its first 300 ELBO values match the saved 300-sweep fit
   (relative error ≤ 10⁻¹⁰). A mismatch is a hard failure.
3. Run the 8 fits with `OMP_NUM_THREADS=1` and at most 3 concurrent fits. Each 300-sweep fit
   took about 1 minute.
4. For each fit, report:
   - status, iterations, and the final relative ELBO change against the tolerance;
   - the correlation of its risk scores with the 300-sweep iterate;
   - matched program correlations with the 300-sweep iterate;
   - held-out concordance under `frozen_reverse_cindex()`.

   A fit that still caps at 1,000 is reported as not converged. Do not raise the cap again
   or tune anything.
5. Recompute the expression reference from the refits and the four already-converged fits:
   mean CV C for PL and PE, and every paired methylation-minus-expression difference and SE
   (SD of the five fold differences / √5).
   - Show the 300-sweep and 1,000-sweep versions side by side. Do not overwrite the
     300-sweep tables.
   - If the PE expression full fit converges, complete the PE expression-correlation step.
   - Neither converged PE methylation full fit has a survival-active program, so no PE
     program can meet the methylation-only rule either way. State this.
6. Write `expression_iteration1000.csv` and `cv_summary_iteration1000.csv` under `K05/`, and
   add focused synthetic tests to `tests/` following `tests/test_methylation_k5_rerun.R`.
   The tests cover the ELBO-prefix check, the cap override being the only config change,
   and the paired-difference arithmetic. They should fail if the function returns a
   constant.

## Part 2: what the methylation factors are, and how they relate to expression

Survival activity is not required here. The questions are:

- What does each methylation factor capture?
- Are any methylation factors not represented in the expression factors?
- Are specific CpGs associated with specific genes across the two data types?

This part is descriptive and adds to task C. It is not a registered study, and nothing is
adopted from it.

**Fits.** Use the converged K = 5 full fits (fold 0) on the 144 matched TCGA patients:

- methylation PL at 10k, 25k, 50k, 100k and all 182,934 CpGs;
- methylation PE at 100k and all CpGs;
- expression PL;
- expression PE, only if Part 1 got it to converge.

Use the saved posterior means. Do not refit any full fit.

1. **Overlap and separation.** For each methylation fit against the expression fit in the
   same setup:
   - Compute the full 5 × 5 matrix of absolute patient-score correlations, using the
     canonical YF scores already in the checkpoints.
   - Compute the canonical correlations between the two 5-dimensional score spaces. They
     measure how much of the methylation score space expression spans.
   - Call a methylation factor **methylation-specific** when its largest |r| with any
     expression factor is below 0.5 **and** it is stable:
     - it is matched in at least 3 of 5 converged CV folds with |r| ≥ 0.5 to its full-fit
       factor (use the existing one-to-one matching in `fold_program_matching.csv`);
     - it recurs across at least two CpG counts (largest |r| ≥ 0.7 between full-fit factor
       scores at different counts).
   - Report every factor's values, not only the ones that pass. Before writing results,
     state the thresholds as fixed in the findings note.
2. **Annotation of each methylation factor.**
   - Install `IlluminaHumanMethylation450kanno.ilmn12.hg19` from Bioconductor (Andrew
     approved this dependency). Record its version in the package-versions CSV.
   - For each factor's top 100 CpGs by |E[F]|, record:
     - CpG-island relation (island, shore, shelf, open sea);
     - promoter region (TSS1500, TSS200, 5′UTR, 1st exon) against gene body;
     - the UCSC RefGene symbol;
     - sign and size of the loading.
   - Compare island and promoter fractions with all eligible CpGs.
3. **Biology of each factor.** On the 144 patients:
   - Correlate each methylation factor's scores with every D4 gene's expression, then run
     fgsea on that ranked list. Use `run_fgsea_program()` and the gene sets from
     `code/pathway_enrichment.R` (Hallmark plus the PDAC sets, including Moffitt and Bailey
     where available).
   - State that the gene universe is the 2,218 D4 genes.
   - Report the same enrichment for the expression factors so the two can be compared.
4. **CpG–gene associations.**
   - For every eligible CpG annotated to a D4 gene's promoter region (TSS1500, TSS200,
     5′UTR, 1st exon), compute the Spearman correlation of β with that gene's expression
     on the 144 patients.
   - Summarize the distribution: fraction negative, and the numbers below −0.3 and below
     −0.5.
   - Compare pairs whose CpG is in a methylation factor's top 100 against all promoter
     pairs.
   - For each matched methylation–expression factor pair (highest |r|), count how many of
     the methylation factor's top-100 CpGs are annotated to genes among the expression
     factor's top 100 genes by |E[F]| (gene loadings, `EF`). Compare that count with a permutation null (1,000
     random draws of 100 eligible CpGs, seeded).
5. **Tumor purity.** Methylation factors often reflect purity or immune infiltration.
   - Look in `PDAC_DATA_ROOT` for TCGA PAAD purity estimates (ABSOLUTE or the TCGA PAAD
     marker paper's purity calls). If they exist, correlate each methylation and
     expression factor with purity and report it next to step 3.
   - If they are missing, say so and stop that step. Do not download anything without
     asking Andrew.
6. **Outputs.** Under `K05/factor_characterization/`:
   - `score_correlations.csv`, `canonical_correlations.csv`, `factor_stability.csv`;
   - `top_cpgs_annotated.csv`, `annotation_enrichment.csv`, `fgsea_by_factor.csv`;
   - `promoter_cpg_gene_correlations.csv` (aggregate per pair);
   - `cross_modality_overlap.csv`, `purity_correlations.csv` if available.

   Write the analysis as a new mode or companion script that reads saved checkpoints only.
   Add synthetic tests: a planted shared factor must be called non-specific and a planted
   methylation-only factor specific; the annotation join must handle CpGs that map to
   several genes; the permutation null must be seeded.

## Reporting

- Append a section "K = 5 expression reruns and factor characterization" to the C findings
  note. Cover:
  - Part 1, with the 300- and 1,000-sweep numbers side by side;
  - the per-factor table: closest expression factor and |r|, stability, specific or not,
    leading annotation and enrichment, and purity |r|;
  - the CpG–gene results.

  Write in plain prose for a biostatistician. Name factors as "methylation factor k at N
  CpGs". Do not interpret a factor as biology unless its enrichment or annotation supports
  it.
- Run `Rscript tests/run_tests.R` and your focused tests. All must pass. Report the counts.
- Get an independent review from a fresh reviewer who sees the code, this prompt and the
  outputs, but not your conclusions. It should recompute the paired statistics, the
  correlation matrices and the specific/not-specific calls. Apply real corrections.
- Then update the shared status files (pull first, small edits):
  - DECISIONS.md: a short addendum to the task C entry, descriptive, nothing adopted;
  - an attempts-log entry;
  - PROJECT_STATUS.qmd;
  - the week-plan C line, with commits.
- Finish with a "Text for the documents" section in the findings note. It should give the
  exact numbers, K = 5, the caveats (n = 144, pre-fold screening and imputation, D4 gene
  universe, any fits that still capped) and the CSV paths, so the documents session can
  fill the Chapter 4, writeup and progress-book methylation placeholders.
- Also say which CpG count you would carry into task D's joint fit (step 4), and why. This
  is a recommendation, not a decision.

Stop and report if a fit errors, an ELBO prefix fails to reproduce, the annotation package
cannot be installed, or anything here contradicts the saved design. Do not work around
problems silently.
