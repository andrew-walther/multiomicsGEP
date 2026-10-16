# Round 2: current status and review handoff

Updated2026-10-16. This is the current summary, replacing the former layered
checkpoint preamble. The complete earlier handoff is preserved verbatim in
[the chronological history](Codex_to_Claude_Round2_history_10_15_26.md).
Historical statements there describe their own checkpoint, not current status.

**Latest execution checkpoint (2026-10-16):** Andrew separately approved Task2-7.
Afterc94c182 it completes in42.867s with exact entire-fit identity, complete
late50 density coverage and no warnings. Latest full suite694/0. Finite late
posterior steps show0/350 losses for each of F/L/beta, but F-prior refits
lose the reference objective50/350 times and refresh18/50. Joint L/beta
comparisons have additional unresolved losses; see the corrected output's
VALIDITY.md and attempt log. This is not certification of a common implemented
objective or a cause of M3. Recommendation: derive/review fixed-q loading-prior
EM with coherent KL before damping/weights; no remedy is implemented. Stop and
regroup before that new decision. The invalid Task2-6 remains preserved:
its observer omitted reconstruction cross-products; the corrected independent
multi-program test detects that defect. The October16 chapter/report includes
both the invalid first attempt and the corrected result with its limits.

## Start here

1. Read `CLAUDE.md` and this summary.
2. Read [the active priorities and goal hierarchy](Advisor_Priorities_and_Execution_10_15_26.md).
3. Read [the original Round2 specification](Codex_Handoff_Round2_10_16_26.md)
   for unchanged criteria, remedy ladders and per-task rules.
4. Use [the attempt log](Round2_attempts_log.md) for exact designs/results and
   the chronological history for detailed implementation explanations.

**State:** The bounded equation/evidence review and approved corrected observer
are delivered. Stop and regroup before implementing a remedy. The north star
remains unmet; completing this review does not resume the former broad goal.
Branch `round2/supervised` remains local; only rendered gh-pages publication
is authorized. No merge or protected-branch edit is authorized. No fitting
job remains running. Preceding CpG/report checkpoint `0c0c717`, planning base
`f8c895b`, corrected-observer preregistration `c94c182`. No configuration adopted.

**Recommendation:** Keep the north star as the project objective. The corrected
observer supports reviewing prior-consistent loading-prior maximization before
damping/weights; regroup before implementing that new proposal. The active plan
defines subsequent milestones.
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
694 passing tests and0 failures. Software tests do not establish the scientific
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

The review and corrected observer below supply the evidence for the next
prior-consistency proposal. Regroup before remedy implementation or any
fitting panel. Keep original thresholds, model-specific
convergence rules, benchmark/default identity, one-change attempts and training-
only selection. Retain DeSurv deferral and no push/merge authorization.
Write verified findings into the meeting report and current write-up; final
configuration/biology claims retain review and Andrew's approval gates.

## Bounded supervision review — 2026-10-16

**Review recommendation and subsequent status:** The review recommended a
complete-cycle observer before another remedy. Andrew approved it and then
separately approved the corrected replay after the first observer error.
The corrected result is summarized above and in Task2-7/VALIDITY.md. The
design below is retained to explain what was observed, not to authorize a
repeat. It does not establish a cause for every cap. R2 and R3 remain untried.

### Equations and evidence

For one loading, with the Cox quadratic and other moments fixed, the current
sequential sweep uses

\[
A_j=\tau_j\sum_i E[L_{ik}^2]+E[\beta_k^2]y_j^\top W y_j,\qquad
B_j=\tau_j R_j^\top E[L_k]+E[\beta_k]y_j^\top h_{-k}
 -E[\beta_k^2]y_j^\top W E[Z_{k,-j}].
\]

The last term removes the other features of the same program. It is present
in `multimodal_yfb_update_F_mk()` and the canonical Single caller; the earlier
omission is already corrected, not a new diagnosis. The helper refreshes
projection means/variances after each coordinate. For profiled curvature,
weighted-centered survival covariates implement
\(H_c=W-ww^\top/\sum w\); reconstruction residuals stay on their original
scale. These statements concern the fixed quadratic, not exact Cox ascent.

| Unresolved issue / established observation | Supported cause or limitation | Remaining uncertainty | Minimal falsifiable check; ladder relationship |
|---|---|---|---|
| Shared-prior exact collapse: selected TCGA arm loses slab mass at outer99/inner1, all beta moments at inner2; no-start at54/3 then54/4. | A returned delta0 shared prior forces every subsequent beta posterior to delta0. R1 prevents this absorption but does not validate fitting. | Why marginal refitting selects the boundary, and whether avoiding it yields credible, convergent programs. | R2, fixed lower slab probability1/K, remains untried; first reconcile its shared-prior units. Do not infer R2 succeeds from R1. |
| Per-factor Normal fold1 reaches all-zero beta at167/1 while loading program4 remains recognizable. | Its scalar variance MLE is max(0,x²−1/A); all pseudo-z magnitudes fall below1. This differs from an absorbing fixed shared prior: earlier revivals occurred. | Whether local curvature, changing loadings or prior estimation drives the crossing. | Freeze the same state and independently verify A/B, including uncertainty; profiling has already tested the constant-direction part and did not deliver convergence. |
| Both canonical models can cap with nonzero risk; Single seed20261026 oscillates. | Fixed-W scalar algebra does not ensure ascent after refreshing the local surrogate. W−H is PSD at the expansion point only. Five inner sweeps are a cap, not evidence that the working problem is solved. | Which transition first loses a coherent objective: coordinate sweep, hyperparameter refit, scale change or working refresh. | Proposed complete working-cycle observer below. A lower Cox logPL alone cannot diagnose a joint-update defect. |
| Returned F prior and posterior can have incompatible support; cached KL remains finite. | Posteriors use the incoming prior; the subsequent marginal prior fit uses sequential pseudo-data, then returns a different prior without updating those posteriors/KLs. A product of these changing pseudo-data marginals is not the fixed-q hyperparameter objective. | How much this changes the optimization trajectory. The exact audit has3 early native support conflicts and none late; F-only monitor adjustment retains17 negative late increments. | Account for L, F and beta distributions/priors together in the proposed observer. Fixed-q EM is a possible later repair, not an executed remedy or sufficient explanation of oscillation. |
| Canonicalization preserves fitted products but does not transform a carried shared beta prior. | (L,F,beta)→(cL,F/c,c beta) transforms second moments and molecular priors, but the helper has no g_beta argument. KL to an unchanged shared prior need not be invariant. | Magnitude/trajectory effect; scales near extinction were close to1. Fresh per-factor variance estimation is a different case. | Transform q and its prior together in an independent density/KL identity test; a common shared prior generally becomes factor-specific under unequal c. R2 needs a declared parameterization before implementation. |
| M2 original25k arithmetic fails with finite eta and valid projection variance; M3 completes numerically but caps. | One global exponential shift zeros143/144 weights and17 risk sums. Log-space risk sums remove that demonstrated arithmetic failure. | Cause of the remaining expanding predictors/reconstruction movement; separation, scale drift and surrogate instability remain hypotheses. | Use saved M2 state for a local check with the already validated log-space arithmetic; do not rerun M3 or call high training C a gain. |
| Terminal survival precision fractions are6.30e-11–3.03e-7 in four capped fits. | These are conditional likelihood terms excluding priors, evaluated independently per coordinate. A beta point mass gates both terms; a zero mean with positive second moment need not. | Cumulative/sequential influence, early trajectory, and whether unwanted feature redundancy is present. | Inspect A/B and prior changes in their own units; no fraction threshold or inverse-p weight follows. Survival weighting is the last rung. |

Working-step validity has a concrete mathematical limit: for two subjects,
one event before the second subject is censored, eta0=(-a,0), and a move to
eta=(0,0), write p=1/(1+exp(a)). The diagonal quadratic gives
Q=−log(1+exp(a))+(1−p)a−p a²/2, which approaches0 as a grows; the true new
logPL is−log2. Thus the local quadratic can exceed the likelihood. This is
an algebraic counterexample to a global-minorizer claim, not evidence that
this particular move caused a saved failure. Profiling removes the constant
direction but does not supply a global ascent guarantee.

Scale freedom also cannot explain away M3: a consistent gauge transformation
preserves eta and reconstruction, whereas both still change at the cap.
Measurement scaling is separate: rescaling Y and its reconstruction loading
by d changes raw YF by d². Keep the original preprocessing fixed; this is no
authorization to normalize a modality or change the predictor.

**Evidence anchors:** [corrected Single derivation](../../derivations/cox_on_YF/canonical_ZF/canonical_ZF_exact_F_update.qmd),
[multimodal derivation](../../derivations/multimodal_YFB/multimodal_YFB_derivation.qmd),
[technical notes](../../derivations/multimodal_YFB/multimodal_YFB_technical_notes.qmd)
(Cox expansion, coupled prior estimation, KL and scale sections),
[profile proposal](Round2_profiled_curvature_proposal.md),
[arithmetic note](Cox_Logspace_Arithmetic_10_15_26.md), and attempt-log
Task2-1/2/3/4/5, Task2-B1, Task5-M2/M3/A1.
Implementation references: `code/multimodal_yfb_updates.R` (F A/B, signed
prior refit and returned KL), `code/fit_multimodal_yfb.R` (inner loop,
canonicalization, damping and monitors), `code/fit_cox_on_yf.R` (canonical
Single caller), `code/profile_cox_intercept.R`, `code/cox_working_logspace.R`.
The technical-note formulas are mathematical references, not proof of code
identity: actual controls select signed priors/intercepts; the current
multimodal helper has no historical working-weight floor.

Direct saved-evidence links by issue: [prior support summary](../../results/benchmark_sim/outputs/round2/prior_consistency_seed26/summary.csv),
[exact fit identity](../../results/benchmark_sim/outputs/round2/prior_consistency_seed26/identity.csv),
[loading-prior trace](../../results/benchmark_sim/outputs/round2/prior_consistency_seed26/loading_prior_trace.csv),
[TCGA cap components](../../results/multimodal_real/outputs/round2/diagnosis1/cap_diagnosis.csv),
[loading-balance interpretation](../../results/multimodal_real/outputs/round2/loading_balance/interpretation.txt),
[cached objective components](../../results/multimodal_real/outputs/round2/loading_balance/saved_objective_components.csv),
[M2 cause summary](../../results/multimodal_real/outputs/feature_sensitivity/moment_diagnosis/root_cause_summary.csv),
[M3 terminal case](../../results/multimodal_real/outputs/feature_sensitivity/logspace_case/case.csv),
and [M3 history](../../results/multimodal_real/outputs/feature_sensitivity/logspace_case/terminal_history.csv). These preserve
the3 support conflicts, exact original fit identity and M3 terminal failure.
No new patient-level computation or scientific fit was needed for this review.

### Remedy ledger — do not call the ladder exhausted

| Step | Status | Evidence / prerequisite |
|---|---|---|
| R1 fixed shared point-Laplace beta prior | Tried; failed adoption |80 simulations+20 CV fits. Single convergence19/20 at low noise and CV delta−.015761 versus SE.006968 fail; multimodal full/CV fits cap, no credible coefficients. |
| R2 refitted shared prior with pi≥1/K | Untried; applicable to shared absorption | No registered R2 panel. It does not address historical Normal collapse, and its shared-prior/canonicalization units must be specified. Not ruled out by R1. |
| R3 per-factor Normal beta with supervised start only | Untried exact combination | Historical framework **without** supervised start was diagnosed, not this one-option R3 arm. No simulation/CV adoption evaluation of the exact combination is recorded. |
| Shared Normal warm-up / no supervised start contrasts | Diagnostic only; tried | Ten-arm TCGA study: shared Normal has near-zero risk/no credibility; removing the start still absorbs. Not substitutes for R2/R3. |
| Longer warm-up | Untried | Warmup10 was used; no longer value tested. A tuned constant is not free: the protocol treats it as a weight requiring a prospective selection rule. |
| Staged fit from a converged framework | Currently inapplicable to saved capped starts | Full framework and all five framework CV fits cap500. A cap cannot be relabelled convergence; no longer cap or altered stopping rule is authorized. This does not rule out staged fitting if a genuine converged start later exists. |
| Damping beta/F updates | Untried; evidence differs by model | Single late oscillation is documented; earlier multimodal diagnostic ELBO increments were positive. Existing damping blends means and variances, not posterior densities, and carries undamped KL/prior outputs. A valid objective account is needed; no arbitrary damping value chosen. |
| Joint/collapsed (F_k,beta_k) update | Untried | Corrected scalar F update is not a joint update. Derive in cox_on_YF and check numerically before a prospective fit. |
| Survival weight {1,2,5,10} | Untried; last resort | Training-CV selection only after applicable preceding steps/decisions. Tiny terminal fractions do not justify promoting it. |
| Profiled working intercept | Additional approved remedy; tried, failed adoption | Single recovery18/20,20/20; multimodal16/20,19/20. One Single simulation caps; all five multimodal CV fits and both full candidates cap. Multimodal paired gate undefined. |
| Log-space risk-set arithmetic | Numeric repair demonstrated, opt-in | M2/M3 exact failing point repaired;25k still caps500 with riskSD39086/logPL−40177. Not a supervision remedy or default adoption. |
| Prior-consistency / loading-balance observations | Completed diagnostics, not remedies | Task2-5 and B1. Fixed-q loading-prior maximization and a fully coherent current-state objective remain unimplemented. |

### Complete-cycle diagnostic design, subsequently approved and executed

**Question:** After accounting for all current priors, does the remaining
late Single oscillation originate within the frozen-quadratic updates or
only when the approximation is refreshed? The completed F-only audit cannot
answer this. Another repeat of its three early support conflicts would add
little; the new information must be a complete transition decomposition.

**Exact single change:** Add process-local observation of the full
current-state working objective, including L/F/beta KL, at update boundaries.
Return every original update unchanged. No prior, posterior, curvature,
initialization, damping, scale convention, convergence rule or fit control
changes. This is a proposed diagnostic extension, not a remedy rung or a
new convergence monitor. It needs a separate committed prospective attempt
and Andrew's approval before its observational replay.

**Fixed arms:** (a) saved Task2-4 cell021 as the unmodified reference;
(b) exactly one observer-enabled replay of that same Single simulation,
seed20261026/noise.25, n400/validation120,p50,K7, per-factor Normal beta,
point-Laplace L/F, no supervised start, profile_cox_intercept TRUE,w_S1,
no pruning, original300 cap/ELBO rule. Reuse the existing generator and
controls, not a new test bed. Require exact identity of the entire fit and
history, as in Task2-5. Freeze source/reference hashes before registration.
Do not add a multimodal or25k arm; transfer to their instability is unproven.

**New observation:** Capture the actual posterior distributions, not only
their moments, using the incoming prior and each coordinate's own A/B at
the instant of update. Retain these distributions across updates and transform
them with the existing canonicalization. For every cycle, separately evaluate
(i) posterior-coordinate changes with priors and working terms frozen,
(ii) prior refits at fixed posterior, (iii) canonicalization with corresponding
prior/distribution units, and (iv) the working refresh at unchanged terminal
q. Include L, F and beta, reconstruction, uncertainty and KL; do not reuse
cached normalizers against a different prior. Record actual carried priors
separately from a mathematically transformed reference prior; if a prior is
not retained by the fitter, identify that gap rather than silently supply one
and claim a common implemented objective. Include tau changes separately.
An initialized distribution that cannot be reconstructed from saved moments
must remain unevaluated until that coordinate has actually been updated;
do not infer its entropy from two moments. Preserve the complete record and
mark coverage explicitly, especially the final50 cycles.

Use an anchored quadratic
Q=logPL(eta0)+u'(eta−eta0)−(eta−eta0)'Hc(eta−eta0)/2; its expectation also
subtracts tr(Hc Cov(eta))/2. Its additive constant changes on refresh;
the code's constant-free fixed_working_objective cannot be compared across
different anchors as a common objective. Keep this local working objective,
the refreshed approximation and the actual Cox logPL as distinct columns.
Even with coherent KL, this remains an expected local quadratic, not an exact
evaluation of E_q logPL or a proven global lower bound. No computed diagnostic
replaces the original history or stopping rule.
Keep infinite KL/support conflicts; do not subtract infinities or treat a
finite partial objective as complete. No external patient/outcome input.

**Required tests before execution:** independent dense-quadratic A/B and
expected-objective calculations on correlated-feature fixtures; coordinate
ascent at fixed priors/W (including the within-program residual subtraction);
independent mixed-distribution KL/quadrature for point-Laplace and
point-exponential q, analytic Normal KL, and old/new prior-density identities;
spike/no-spike/zero-information and incompatible-support cases; density/KL
invariance when q and its prior are rescaled together; correct anchored
constant and uncertainty at refresh; observer-on/off full-object identity.
Use existing prior-audit numerical conventions (quadrature1e-8,
moment/unit checks1e-12, analytic/native probability1e-7); these are software
checks, not scientific thresholds. Run the full suite and protected-base
default-identity tests before the prospective code/log commit.

**Decision and stopping rule:** Register/commit first, replay once into a
fresh aggregate-only folder, then stop. Identity failure, incomplete late-cycle
distribution coverage or inconsistent density reconstruction blocks causal
interpretation; preserve the error/partial output with no retry or new case.
If a frozen-prior coordinate reduces the independently evaluated fixed
working objective beyond the tested numerical agreement, derive/fix that
update before trying damping. If loss occurs at prior refits, prioritize a
consistent hyperparameter step; for signed F the fixed-q candidate is
pi_new=S/p, rate_new=S/T, S=sum Pr_q(F≠0), T=sum E_q|F|, with explicit
support/boundary handling. If fixed-working steps pass and loss appears only
on refresh, that supports investigating the working-step strategy (including
the damping rung), not an automatic choice of a damping value. If nonfinite
KL prevents localization, report the support obstruction rather than claim
an ascent diagnosis. These are prospective branches, not reported results.
No outcome alone proves a cure for M3. Return to Andrew before implementing
any remedy; no automatic follow-on panel.

This review is the interim handback required by the six-attempt protocol;
diagnostic entries do not establish ladder exhaustion. The proposed observer
is not a new remedy rung and is not automatically authorized as a continuation
past that limit. Andrew must decide on this bounded extension.

**Original scientific criteria remain binding:** This diagnostic cannot pass
an adoption gate. Any later remedy needs a separate committed design on
the original two models, seeds20261016:20261035, noise.25/1, n400/validation120,
p50 per modality, K7, original controls/caps and references: recovery≥16/20
at each noise, null-active mean≤reference+.5, every fit converged by its
model's rule (Single ELBO; multimodal output changes); median rank3 applies
when rank adoption is evaluated. Training evaluation retains the same144
TCGA patients/cache, seed42 event-stratified five folds, training-only Single
scaling and conditional screening/imputation caveat. Require all five finite
paired C values and mean delta≥−sd(delta)/sqrt5, full-training nonconstant
risk and≥1 credible coefficient with positive SD, and report every cap.
The existing reference fold C=NA prevents that paired gate passing as
written; it remains NA. Changing it needs Andrew and prospective criteria,
not a silent reference substitution. Only a final selected configuration
reaches the permitted external report (lower95% deltaC limit>−.02), with
frozen benchmark rows; external cohorts cannot choose a remedy.

**Decision after execution:** Review the next opt-in fixed-q point-Laplace
F-prior step (pi=S/p, rate=S/T) with KL against the returned prior. It addresses
a demonstrated local non-ascent step, not every remaining failure. A future
implementation/panel needs its own decision, algebra tests and committed
single-change design. Retain R2/R3 for the distinct absorption problem after
resolving shared-prior units. No new reference, relaxed convergence or default
change follows from this diagnostic. The broader method remains unresolved;
the failed external Single result and all capped fits stand.

**Review verification (2026-10-16):** Fresh `Rscript tests/run_tests.R`
completed with exit0:684 passed,0 failed. All24 local Markdown links/anchors
resolve; `git diff --check` passes. Only this handoff and the active plan
changed; no code, configuration, tests, results or manuscript changed.
All10 benchmark SHA256 entries pass before and after review; the immutable
tag still peels to `a352f17ca9cff523caa80756aa2af085e62551e1`.
The process check before the suite found no active R fitting process.
No scientific fit, external scoring or remedy execution occurred. Commit
on `round2/supervised` only, with branch rechecked immediately beforehand.

## Provenance and history

### Observer implementation walkthrough

Purpose: observe the existing Single algorithm without changing its fitted
state. The inputs are the original update arguments/results and the active
simulation call; outputs are aggregate transition records and exact-fit
identity checks. No dependency was added. Base-R densities and matrix products
are used so the calculation can be checked against independent algebra.

| Function/block | Logic, inputs/outputs and design limitation |
|---|---|
| `cycle_log_prior(q,g)` | Sum expected log density under the mixed spike/continuous measure. A support conflict gives−Inf, never a finite clipped penalty. |
| `cycle_posterior(A,B,g,mean,second,log_ml)` | Recover slab/absolute moments and E log q from the actual normal-means normalizer. Native moments are retained; entropy is not guessed from those moments. Normal, point-Laplace and point-exponential families have separate density formulas. |
| `cycle_kl(q,g)` | E log q−E log g in the same units, with infinite KL preserved. This is a declared-prior calculation, not certification that g is carried by the fitter. |
| `cycle_scale(q,g,c)` | Scale means/absolute moments by c and second moments by c². Continuous entropy changes by−Pr(nonzero)log(c); the spike has no Jacobian. Transform the reference prior together; actual carried-prior status is labelled separately. |
| `cycle_prior(fitted_g,family)` | Read the actual ebnm single-slab distribution without adding an empirical threshold. Multiple continuous slabs fail explicitly rather than being silently collapsed. |
| `cycle_likelihood(state,work)` | Evaluate Gaussian reconstruction and the anchored expected profiled Cox quadratic. Reconstruction uses squared mean residual plus posterior variance, including all cross-program products. The initial observer missed those products and is invalidated; the corrected two-program fixture detects that error. |
| `cycle_snapshot(frame)` | Copy the current Y/moment/precision references from the active call; no assignment is made into the fitter frame. The bounded runner uses no cohort augmentation. |
| `cycle_observer()` and its `record` collector | Keep the density/prior ledger across updates. Record incomplete initialization coverage, finite/nonfinite components and native moment agreement. Transformed transient L/beta references do not become implemented priors. |
| Local `wrapper` | Invoke the original function with unchanged arguments, examine its result, and construct a separate proposed snapshot before the caller assigns it. Separate posterior, prior, tau, scale and refresh events; return the original result unchanged. |
| Observer `run(callback)` | Temporarily install wrappers in the existing global fitting environment and restore the original functions with on.exit, including errors. The software fixture checks the entire returned fit/history. |
| Driver | Reuse fixed cell021 controls/generator, refuse existing output, retain errors/warnings, verify exact identity/hashes, and save detailed rows only under ignored data/. Compact per-iteration/event summaries retain nonfinite counts. The original summaries have a validity warning and must not diagnose the fitter. |

The corrected reconstruction identity is
E[(Y−LFᵀ)²]=(Y−EL EFᵀ)²+EL2 EF2ᵀ−(EL²)(EF²)ᵀ, elementwise squares.
The log-density identity is E log q=E log g+E log N(x;theta,s²)−log m(x).
Both matter: a correct entropy calculation cannot rescue an incorrect
likelihood expectation. Eight focused tests now cover the mathematical
pieces and observation identity; the full suite passes692/0. Scientific
convergence and adoption remain independent requirements.

### Repository footprint and meeting publication

Against local origin/main before this work,166443 added lines comprise144461
CSV lines (87%) and7937 R source/test lines. Six traces contribute105006 CSV
lines. This is not160k new lines of model code. No evidence was deleted.
Review archival/storage separately with provenance and a restore procedure;
do not merely discard traces needed to support the diagnoses. The new detailed
20700-row record remains ignored; Git receives compact summaries and the
validity warning.

The October16 chapter adapts the advisor report's tables/figures into the
existing book headings/navigation and adds the mathematical review and
invalid diagnostic status. Native RStudio Quarto renders all13 pages; the
old standalone Quarto's bundled Sass has incompatible architecture. An older
chapter's broken supporting-report link now points to a packaged resource.
All108 landing/adjacent/new-chapter local links/assets resolve; search indexes
the new chapter. The report HTML/PDF is synchronized. Publication targets
only the rendered gh-pages site, preserving its history and leaving main and
feature/supervised-F untouched.

Final publication verification (2026-10-16): the corrected-evidence refresh
is live at `https://andrew-walther.github.io/multiomicsGEP/chapters/2026-10-16.html`.
GitHub Pages reports built for gh-pages commitb5a413a3b562ac8176cb325284ab5e64a1a101aa.
Only rendered content was pushed, preserving publication history; this source
branch was not pushed or merged. All108 local landing/adjacent/chapter
references and the search entry verify; the live chapter and all three figures
were checked. Both report formats render; the fresh result-checkpoint suite
passes694/0. All10 benchmark hashes remain unchanged. Stop and regroup before
any proposed prior remedy or additional fitting panel.

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
