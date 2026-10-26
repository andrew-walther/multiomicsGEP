# Repository size and test review

October16,2026. Baseline: local origin/main `a352f17ca9cff523caa80756aa2af085e62551e1` versus
round2/supervised at f3fea03. This is the accumulated committed branch diff,
not176,342 uncommitted lines. The checkout was clean before this work block.
Binary files are excluded from text line totals; no history is rewritten.

## What the added lines contain

| Category | Added lines | Share | Files |
|---|---:|---:|---:|
| Saved aggregate/diagnostic CSVs |145,198|82.3%|1,687|
| Production R under code/ |1,578|0.9%|16|
| R experiment runners |4,739|2.7%|35|
| R tests |2,508|1.4%|30|
| Rendered book/documents |10,772|6.1%|19|
| Plans, source documents and other text |11,547|6.5%|86|
| **Total** |**176,342**|**100%**|**1,873**|

Rounding prevents percentages summing exactly to100. There are283 deleted
text lines at this baseline. Generated figures/PDFs are additional binary
artifacts, not R code. The full suite's64 active test files contain11,504
source lines, including older tests outside this branch's additions.

**Largest concentration:** six multimodal per-iteration trace CSVs each have
17,501 lines, together105,006 lines (59.5% of all additions). They occupy
28,044,473 bytes uncompressed; an in-memory gzip estimate is4,194,779 bytes
(85.0% smaller). This estimate changed no file. The six are
five diagnosis1 traces (framework, selected, no_start, normal, fixed) and
framework_fold_trace/fold1_trace.csv, all under results/multimodal_real/outputs/round2/.

## Test review: what is and is not redundant

Parsed all64 active files:706 run_test expressions, **zero exactly duplicated
parsed test bodies**. A simple text matcher recognizes695 calls; parsing the R expressions
finds all706. Three tests share the title `custom A_floor is respected`
and two share `all outputs finite`, but they target different modular updates.
The names are repetitive; the coverage is not identical. Exact-body inspection
cannot prove the absence of semantic redundancy, and this is not a full mutation
or coverage study.

Targeted manual review of likely overlaps:

| Tests compared | Finding | Recommendation |
|---|---|---|
| beta/L/F precision-floor and default-divisor checks | Same contract in separate update implementations | Keep; one can regress while the others pass |
| LB and YFB multistart checks | Similar restart metadata, different fitter paths | Keep integration coverage; shared fixtures are optional later |
| Single/multimodal constant-risk feature metrics | Different public helpers/callers and error contracts | Keep until an intentional helper consolidation; do not delete tests first |
| Prior-consistency audit versus fixed-q production tests | Audit helper reconstructs old behavior; production helper changes priors/KL | Keep independent quadrature and support checks for both |
| Fixed-prior versus profiled protocol NA-C checks | Same scientific restriction plus different design/reporting routes | Keep: forbid missing-fold deletion or C=.5 substitution |
| Working-cycle density/algebra versus path guards | Math regression versus protecting invalid/corrected outputs | Keep while that observer remains callable |

No reviewed test is justified for deletion solely by the count706. Many
are small numerical identities, boundaries or regression checks, not706
scientific fits. The existing706-test suite passed in103.75seconds in this block;
timing includes fitting fixtures and is machine-specific. A quick developer subset may be useful later, but required full validation
before commits remains unchanged.

## Prioritized cleanup ledger

| Priority | Candidate | Benefit / risk | Disposition |
|---|---|---|---|
|1|Compress the six bulky diagnostic traces, with byte hashes and documented CSV.gz readers|Large checkout/diff reduction; consumers and historical replay paths must keep working|Recommended follow-up after the methodological checkpoint; do not delete evidence |
|2|Package completed diagnostic evidence separately from active fitting code|Makes the repository easier to navigate; needs a provenance index and consumer audit|Plan an archive/export boundary after review, not a broad deletion now |
|3|Remove generated HTML from selected development-only docs|Reduces generated diffs; may break existing share/read workflows|Decide distribution conventions first; retain the public book and requested HTML/PDF deliverables |
|4|Consolidate duplicated metric/fixture patterns|Small code saving; risks changing numerical/error contracts|Low priority; separate refactor with behavior tests if pursued |
|5|Drop repetitive-looking tests or old runner scripts immediately|Small line saving; can remove independent regressions or reproducibility|Not recommended on present evidence |

Compression reduces checkout size and displayed text additions. It does not
remove objects already in Git history; do not promise a smaller historical
repository without a separately authorized history rewrite. The immutable
benchmark snapshot/tag is excluded from every cleanup option. No patient
inputs, fits, failed attempts or regression evidence were deleted.

## Decision

Keep scientific/regression tests now. For material size reduction, prefer
lossless compression of the six trace files over deleting tests. Before any
conversion, map readers, preserve exact decompressed-byte SHA256 values and
verify existing summaries reproduce. No compression, deletion, broad refactor,
new quick-test default or test-suite policy change was executed in this audit.


Follow-up in this block: four focused tests were added for the bounded case's
sole control change, independent density calculation, known tail-error detection
and correct recovery/sign/false-program reporting. Full suite710/0; these
protect the new numerical check rather than duplicate the audited706 cases.
New case outputs are compact aggregates; detailed fit/update telemetry stay ignored.
