# Round 2: current status and review handoff

Updated2026-10-15. This is the current summary, replacing the former layered
checkpoint preamble. The complete earlier handoff is preserved verbatim in
[the chronological history](Codex_to_Claude_Round2_history_10_15_26.md).
Historical statements there describe their own checkpoint, not current status.

## Start here

1. Read `CLAUDE.md` and this summary.
2. Read [the active priorities and goal hierarchy](Advisor_Priorities_and_Execution_10_15_26.md).
3. Read [the original Round2 specification](Codex_Handoff_Round2_10_16_26.md)
   for unchanged criteria, remedy ladders and per-task rules.
4. Use [the attempt log](Round2_attempts_log.md) for exact designs/results and
   the chronological history for detailed implementation explanations.

**State:** Goal paused at Andrew's request; planning updates do not resume it.
Branch `round2/supervised`; no push, merge or protected-branch edit authorized.
No fitting job remains running. Scientific/report checkpoint `0c0c717`;
restart-plan checkpoint `3e2757d`. No new configuration is proposed or adopted.

**Recommendation:** Keep the north star as the project objective, but use a
smaller active goal when resuming: produce a review-ready diagnosis and one
prospective next experiment, then regroup before another scientific fitting
panel. The active plan defines its deliverables and subsequent milestones.
Do not mark the existing scientific goal complete just to replace it with a
smaller goal. Its older state description incorrectly says external performance
is unchanged; the failed Single noninferiority result below is authoritative.

## Benchmark and current write-up

- Immutable tag: `benchmark-interim-writeup-2026-10-14`, commit
  `a352f17ca9cff523caa80756aa2af085e62551e1`.
- Immutable snapshot: `results/benchmark_snapshot_2026-10-14/`; all10 SHA256
  entries passed before/after the latest scientific/report checkpoint.
- Single D4: alpha_F=0, K7, mean external C=.626699 over5 cohorts.
- Joint framework: ICGC primary C=.685, matching expression-only on the
  same donors; not evidence of incremental prediction.
- `paper/writeup/writeup.qmd` and HTML/PDF now contain Andrew-authorized
  bounded sensitivity/diagnostic revisions. The approved interim is preserved
  in the tag; current manuscript edits do not replace the benchmark model.
- Meeting-facing account: `docs/reports/Advisor_Progress_since_10_2_26.qmd`
  and HTML/PDF (13 PDF pages); current larger write-up42 PDF pages. Both render.

## Current task status

| Task | Status and key result | Main checkpoints; tests at checkpoint |
|---|---|---|
| 0 Integration | Complete; all five Round1 branches integrated into this branch. |2131c30, eb6b290, f2c2a85, bf3e050, b9a3952;559/559|
| 1 K selection | Not adopted. Protected/full-state refits still fail original recovery, convergence and rank requirements. Retain approved K7/no-pruning fallback. |06a18d5,95b8cc4;573/573|
| 2 Supervision/collapse | Unresolved. Exact traces identify prior absorption and nonconvergence; fixed prior fails. Profiling removes constant-risk predictions but all five multimodal CV fits cap500, four lack credible programs. |5a56f66 and subsequent profile result in history; b502795 prior audit;1692431 balance audit; latest suite684/684|
| 3 DeSurv head-to-head | Installed/pinned adapter and runner prepared; Andrew explicitly deferred real-data comparison. No renewed execution authorization. |f37d646;593/593|
| 4 Simulation comparators | Complete:40 datasets,200 rows,0 errors; all denominators and rank0/rank1 cases retained. Supervised planted recovery17/20 and19/20. DeSurv omitted by Andrew. |dd11923,d353322,b941d53;608/608 at results|
| 5 Feature sensitivity | Expression and descriptive CpG work complete; methylation/joint curves partial. Stable-arithmetic25k case still caps; correlations provisional. See tables below. |64f4920,aff0d03,47d4a4a,73c4f0b,15f61b6,0c0c717;684/684 latest|
| 6 Robustness cleanup | Optional clean-process rerun unrun. Earlier saved tables retained with documented nonzero runner exit after script mutation. |Historical limitation retained; lower priority|
| 7 Write-up housekeeping | Complete:14 process files moved, live pointers fixed, write-up and planning renders verified. |639351b;564/564|

Historical test counts are not current suite totals. Latest full suite has
684 passing tests and0 failures. Software tests do not establish the scientific
success criteria. See the history/attempt log for full commits, seeds and arms.

## What is solved, and what remains scientifically open

**Implemented or diagnosed:** positive-SD protection excludes dead point masses;
full posterior state can be retained in pruning refits; comparator failures
stay in denominators; a failed CV fold reproduces exactly; stable risk-set
arithmetic repairs the25k underflow crash. Expression feature expansion and
full-probe descriptive variability studies are delivered.

**Still open:** a defensible, convergent survival-supervised replacement;
valid rank selection; CpG/joint feature-count sensitivity; demonstrated
multimodal survival or additional biological value; final-program characterization.
Experimental code does let survival update F. The evaluated Single benchmark
still uses alpha_F=0, so coefficient fitting does not make it the final
survival-guided program-discovery method.

The profiled panel recovers Single18/20 and20/20 and multimodal16/20 and19/20
at the two fixed noise levels. Convergence requirements still fail. Single
training CV .624776 versus .623985 passes its paired one-SE gate, but full
training caps300 and a simulation caps. Multimodal CV mean=.563678; all five
candidate risks finite/nonconstant, all cap500. The reference has an undefined
fold C, so the prescribed five-finite-pair comparison cannot pass. C remains
NA; Andrew rejected substituting.5. No Round2 external score has run.

Earlier selected supervised Single external C=.620725 versus .626699:
ΔC=-.005974,95% interval[-.031708,.020662]. This FAILS the original lower-limit
>-.02 requirement. It is not external equivalence or preserved performance.
No further external cohort is used to select a remedy or feature count.

## Key feature-count results

### Expression: completed fixed-K7, alpha_F0 sensitivity

Same273 TCGA+CPTAC training patients; cold-primary authorized after frozen
activity rank4/7 made the prescribed warm inverse singular. Screening/scaling
precede folds: conditional CV, not a nested raw-data pipeline.

| Shared genes | Training C | Mean conditional CV C | CV SE | Full retained |
|---:|---:|---:|---:|---:|
| Frozen2064 |.656307|Not recomputed|—|4|
| Fresh2064 |.656259|.608404|.046196|4|
|5001|.615915|.577289|.026584|4|
|10000|.616977|.608789|.039651|4|
|17900 all shared|.619149|.606791|.018562|3|

All24 fits converge; no warnings/errors. Median common-gene loading r with
baseline=.297/.333/.566; risk r=.655/.638/.505 as count rises. Similar C does
not establish stable biology or gene-level feature selection. No count adopted.
Safe tables/figures: `results/benchmark_sim/outputs/feature_expansion/cold_primary/`.

### CpGs: inventory complete; fitting instability remains

Observed-value mean/SD inventory covers305492 raw TCGA probes/144 patients.
The eligible model universe has305352 probes; unchanged whole-universe
imputation exactly reproduces frozen10k values, order and ranking.
Within that eligible universe, mean[.1,.9] retains182934 probes; the95%-extreme
proportion rule retains221142. Neither filter or any SD cutoff is adopted.
Plots: `results/multimodal_real/outputs/feature_sensitivity/cpg_qc/`;
input checks: sibling `universe_preflight/`. Patient inputs/fits stay ignored.

| Full TCGA fit | Genes | CpGs | Training C | Credible | Convergence |
|---|---:|---:|---:|---:|---|
| Frozen joint |3000|10000|.637768|1|No,500 cap|
| Frozen expression-only |3000|0|.643559|1|No,500 cap|
| Original methylation-only |0|10000|.612671|0|No,500 cap|
| Original methylation-only |0|25000|NA|NA|Numeric execution error|
| Same25k, stable arithmetic |0|25000|.983854|7|No,500 cap|

These are contextual training rows with different inputs, not matched CV or
incremental methylation benefit. The original25k error occurs at outer6:
global exponential shifting underflows143/144 weights and17 risk-set sums
are zero; projection moments are valid. Imputation is not the identified cause.
Opt-in log-space arithmetic preserves the same likelihood/working approximation
and fixes nonfinite calculations. The same case then completes in44.3 minutes,
but riskSD39086 and logPL=-40177, with final relative risk/reconstruction
changes=.00589/.0622 against.0001 tolerance. Its high training C is not a
validated gain. Numeric completion does not resolve optimization instability.

Provisional10k/25k comparison on common10000 probes:

| Program | Activity r | Loading r |
|---:|---:|---:|
|1|.982|.945|
|2|.878|.922|
|3|.968|.925|
|4|.513|.118|
|5|.764|.038|
|6|.235|-.087|
|7|.941|.923|

Risk r=.245. Both fits cap and arithmetic differs; not a clean converged
count effect. No labels transferred or new correlation threshold chosen.
Outputs: `results/multimodal_real/outputs/feature_sensitivity/logspace_case/`
and `cpg_count_agreement/`, including frozen contextual benchmark rows.
Counts50k/100k/all, methylation CV, filter and joint count comparisons remain
unrun. Old source fingerprints must not be bypassed after arithmetic changes.

## Other evidence relevant to the next decision

- Saved joint-role audit matches all seven frozen/newer profiled activities
  (r=.832–.999); risk r=.954. Program4 is the sole credible program in both;
  program7 remains methylation-only. Both cap500. Modality norm shares are
  not additive PVE or survival contributions; no new biology labels inferred.
- Saved-state loading-balance audit finds median survival precision fractions
  6.30e-11–3.03e-7 and absolute-numerator fractions8.70e-7–1.82e-4 among
  non-point-mass programs. These are isolated conditional likelihood terms,
  excluding priors, not proof of trajectory causes or a selected mixing weight.
- Prior-consistency tracing reproduces the saved fit and shows support/scale
  concerns. Partial diagnostic monitor adjustments do not cure the late
  oscillation. Do not treat a stale cached objective as a validated common ELBO.

## Next handoff and stop rule

The next deliverable is a bounded mathematical/evidence review with a specific
prospective next experiment and the applicable ladder/approval identified.
Regroup before any fitting panel. Keep original thresholds, model-specific
convergence rules, benchmark/default identity, one-change attempts and training-
only selection. Retain DeSurv deferral and no push/merge authorization.
Write verified findings into the meeting report and current write-up; final
configuration/biology claims retain review and Andrew's approval gates.

## Provenance and history

- [Active plan and goal hierarchy](Advisor_Priorities_and_Execution_10_15_26.md)
- [Original Round2 rules/criteria](Codex_Handoff_Round2_10_16_26.md)
- [Exact prospective attempts and results](Round2_attempts_log.md)
- [Unabridged former handoff](Codex_to_Claude_Round2_history_10_15_26.md)
- [Original advisor notes](Meeting_Notes_Followup_10_2_26.md)

The history is preserved verbatim at consolidation. Do not append new current
status there; update this summary and the attempt log instead. New code still
requires tests, README/CLAUDE updates and the project's teach-back conventions.


Consolidation verification: original handoff preserved byte-for-byte; new local
Markdown links resolve. Fresh full suite684 passed,0 failed. Documentation
only; no fit, benchmark, manuscript content or goal object changed.
