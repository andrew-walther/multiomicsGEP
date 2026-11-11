# Branch review findings

Date: 2026-11-11

Author: Andrew Walther

This is a proposed disposition, pending Andrew’s approval. No branch or tag was changed. The snapshot compares all refs with main at `47690aa2d519d82015748b59b1775664f5e53fc8`. It covers 35 local branches and 17 origin-tracking branches in 38 rows. A live `git ls-remote --heads origin` check confirmed that all 17 remote names and tips match the tracking snapshot. The initial query failed on DNS, but the subsequent live check succeeded. Recheck remote tips immediately before any cleanup.

Each row lists every local and remote-tracking branch ref, combining matching names only when their tips match. Dates are tip committer dates. “Ahead” counts commits reachable from the ref but not main. Zero means its tip is an ancestor of main, not that its scientific conclusions were adopted. References list exact branch-name mentions in DECISIONS.md (D), AGENTS.md (A), and ROADMAP.md (R). Other evidence is summarized in the purpose column. Where only a tip subject is given, no associated plan or decision was established in this review. Delete means propose deleting the named refs only after approval. Keep-reference means first preserve unique history under an approved archive tag. Unfinished means retain the branch while its remaining use is decided.

| Branch refs | Tip | Date | Ahead | Merged into main | Proposed | Purpose / associated decision | D/A/R references |
|---|---|---|---:|---|---|---|---|
| local `cavi-deflation-init` | `3e91f60e` | 2026-07-13 | 0 | Yes | Delete | K-parsimony deflation negative result (ROADMAP, DECISIONS) | D 3874; R 517, 602 |
| `origin/claude/recent-project-updates-n04uuo` | `ee86a063` | 2026-08-13 | 0 | Yes | Delete | Add factor classification and K-selection validation plan (2026-08-13) | None |
| local `cleanup/before-next-step` | `2ac8c574` | 2026-10-30 | 0 | Yes | Delete | Rerun the D4 benchmark after the KL fix, with single and double log | None |
| local `cleanup/sign-variance-audit` | `4feec243` | 2026-10-30 | 0 | Yes | Delete | Verify hazard orientation and posterior coefficient variance | None |
| local `codex/multimodal-yfb`, `origin/codex/multimodal-yfb` | `6b8efe48` | 2026-10-12 | 0 | Yes | Delete | Multimodal prototype; week plan task D replaces extension of branch | A 46 |
| local `codex/multimodal-yfb-docs` | `f63effec` | 2026-09-17 | 0 | Yes | Delete | Archive September progress notebook review | None |
| local `cohort-cols-L`, `origin/cohort-cols-L` | `b8f0d57b` | 2026-05-25 | 0 | Yes | Delete | Add caffeinate sleep prevention to long-running job steps in merged benchmark plan | None |
| local `desurv-gene-selection`, `origin/desurv-gene-selection` | `f049e623` | 2026-05-30 | 0 | Yes | Delete | Ignore Quarto render cache directories (*_cache/) | None |
| local `feature/supervised-F` | `2bebbd44` | 2026-10-15 | 0 | Yes | Delete | Keep the benchmark snapshot copies out of git; the manifest verifies them; local and tracking tips differ | D 348, 550, 955, 991, 1055, 1525 |
| `origin/feature/supervised-F` | `275a30a1` | 2026-10-15 | 0 | Yes | Delete | Set the order of survival-supervision approaches and the C3-C5 thresholds; local and tracking tips differ | D 348, 550, 955, 991, 1055, 1525 |
| local `fix-L-update-beta-cycle` | `e74a6224` | 2026-04-29 | 0 | Yes | Delete | L-update and beta cycle (DECISIONS) | D 5682 |
| local `fix/current-method`, `origin/fix/current-method` | `bde88ff8` | 2026-10-30 | 0 | Yes | Delete | A1 larger survival-weight assessment | None |
| local `fix/yfb-supervision`, `origin/fix/yfb-supervision` | `d635bab9` | 2026-11-09 | 0 | Yes | Delete | Merged 2026-11-09; week plan section 0c | None |
| `origin/gh-pages` | `476f7315` | 2026-11-09 | 27 | No | Unfinished | Published progress-book deployment history; retain service branch | D 532, 550 |
| local `joint-k-alpha-bayesopt` | `c5069511` | 2026-07-13 | 0 | Yes | Delete | Joint K/alpha search (ROADMAP, DECISIONS) | D 3780, 4087; R 537, 920 |
| local `lab-meeting-06-18`, `origin/lab-meeting-06-18` | `d8870bbb` | 2026-06-16 | 0 | Yes | Delete | Add DeSurv comparison backup slide (accurate to Young et al. PNAS 2026) | None |
| local `m5/k-selection` | `ac572f6e` | 2026-10-15 | 0 | Yes | Delete | Protected-pruning assessment; adoption blocked (AGENTS, DECISIONS) | D 953, 1051, 1102; A 62 |
| local `m6/robustness` | `2729b577` | 2026-10-15 | 0 | Yes | Delete | Training robustness and feature sensitivity (AGENTS, findings index) | A 95; R 1206 |
| local `m7/desurv-h2h` | `54b3df88` | 2026-10-15 | 0 | Yes | Delete | DeSurv prerequisite missing (findings index) | None |
| local `m7/multimodal-comparators` | `33a716a1` | 2026-10-15 | 0 | Yes | Delete | Multimodal comparator audit (AGENTS, findings index) | A 66 |
| local `m7/sim-comparators` | `5ebf47f9` | 2026-10-15 | 0 | Yes | Delete | Planted comparator evaluation, incomplete denominator (AGENTS, 2026-10-15 decision) | A 64 |
| local `main`, `origin/main` | `47690aa2` | 2026-11-11 | 0 | Yes | Unfinished | Current working branch | Generic main mentions not counted as branch-specific references |
| local `meeting-synthesis-4-30` | `ac5b2158` | 2026-05-05 | 0 | Yes | Delete | Rewrite ROADMAP and README for biostatistician readability | None |
| local `meeting/2026-08-21-prep`, `origin/meeting/2026-08-21-prep` | `6fc44174` | 2026-08-21 | 0 | Yes | Delete | Meeting preparation and pathway enrichment (ROADMAP) | R 368, 385 |
| local `meeting/2026-08-27-prep` | `aa3d0042` | 2026-08-27 | 0 | Yes | Delete | Clean up stale K=20 fallback string after independent deck review | None |
| local `merged-benchmark`, `origin/merged-benchmark` | `88ed06d2` | 2026-05-27 | 0 | Yes | Delete | Add DeSurv alignment implementation plan and ROADMAP entry | D 4719, 4795; R 698 |
| local `multi-cohort-sim`, `origin/multi-cohort-sim` | `fa1b46ce` | 2026-06-14 | 0 | Yes | Delete | Multi-cohort simulation (ROADMAP) | R 1112 |
| local `objective-normalization` | `6719e2c6` | 2026-07-12 | 0 | Yes | Delete | Document Phase 1 completion: updated D4 metrics, K_eff open item, session log | D 4129, 4285, 4335 |
| local `phase1-platform-norm` | `3f10a1b4` | 2026-05-06 | 2 | No | Keep-reference | Early canonical YFB platform-normalization result | None |
| local `phase2-init-constraints` | `9460bf8c` | 2026-05-06 | 1 | No | Keep-reference | Initialization constraints; tip records discard verdict | None |
| local `phase3-followup-warmstart` | `099e3e5b` | 2026-07-13 | 0 | Yes | Delete | K-parsimony warmstart result (ROADMAP) | D 3944; R 503 |
| local `phase3-multistart` | `07e3be1b` | 2026-05-06 | 0 | Yes | Delete | Phase 3 (multi-initialization): implement wrapper, sweep, verdict=Discard | None |
| local `phase4-k-cv` | `c9b74659` | 2026-05-06 | 0 | Yes | Delete | Phase 4 K-CV results: normal→K=3, point_normal→K=8 (spike artifact); recommend K=5 | None |
| local `pre-merge-synthesis-05_06_26`, `origin/pre-merge-synthesis-05_06_26` | `e5875425` | 2026-05-20 | 0 | Yes | Delete | Fix YFB K-CV sign orientation; delete stale benchmark artifacts | None |
| local `prelim/part3-slides`, `origin/prelim/part3-slides` | `b80991d3` | 2026-11-03 | 9 | No | Keep-reference | Draft already incorporated on main by `bc22bf24`; two outline files have newer derivation paths on main. Preserve old tip only as historical provenance, then delete refs after approval | Main commit `bc22bf24`; path update `ce16573c` |
| local `rebuild/unified-ssbmf`, `origin/rebuild/unified-ssbmf` | `41f8b9f1` | 2026-10-26 | 17 | No | Keep-reference | Unmerged rebuild and balance panels; AGENTS and realignment plan identify reference; week plan task D supersedes extension | D 384, 407; A 17 |
| local `round2/supervised`, `origin/round2/supervised` | `2fbb78f7` | 2026-10-24 | 0 | Yes | Delete | Round2 diagnostics; archive/round2-supervised-2026-10 tag | D 401; R 1184 |
| local `stratified-cohort-baseline-hazard` | `92d49593` | 2026-07-15 | 0 | Yes | Delete | Address independent review of stratified Cox partial likelihood | None |

## Unique commits

These lists enumerate the commits counted above. The two `feature/supervised-F` tips differ, but both are on main. `origin/gh-pages` is generated deployment history and should stay available for the published book. `prelim/part3-slides` has nine non-ancestor commits, but its slide content is already present on main. Commit ancestry alone does not establish missing content. The unified rebuild is an unmerged reference, not an approved production model.


### `phase1-platform-norm`

- `3f10a1b4 2026-05-06 Canonical YFB merged result: K=3, no_rank, per_platform_standardize`
- `20a123ea 2026-05-06 Resolve YFB merged beta->0 collapse: remove rank transform, K=3 selected`

### `phase2-init-constraints`

- `9460bf8c 2026-05-06 Phase 2 (init constraints): implement, test, sweep, verdict=Discard`

### `prelim/part3-slides`

- `b80991d3 2026-11-03 Align the purpose and test rows on the aims slide`
- `fe23e8e5 2026-11-03 Use the spare space on the model-status and aims slides`
- `f1da648f 2026-11-03 Lengthen the Part III script to about 20 minutes`
- `85844442 2026-11-03 Apply the trope review to the progress book`
- `0a004785 2026-11-03 Discourage prose semicolons and dashes in AGENTS.md`
- `cc3abc4a 2026-11-03 Mirror the Part III prelim deck into bios-dissertation on commit`
- `29274ff1 2026-11-03 Build the Part III prelim deck from the approved outline`
- `863d9be8 2026-11-03 Add trope review of the progress book`
- `674927ec 2026-11-03 Draft the Part III prelim slide outline for review`

### `rebuild/unified-ssbmf`

- `41f8b9f1 2026-10-26 Add the per-arm row files of the balance panels`
- `fb0d83d6 2026-10-26 Regenerate ssbmf Rd files after the derivation section renumbering; ignore run logs`
- `2d881d38 2026-10-26 Report the first results of the ssbmf rebuild`
- `99ea88cd 2026-10-26 Expression + methylation real-data results (rebuild Task 17)`
- `c16aac0c 2026-10-26 Record the Part 4 balance results on the low-variance cells`
- `e3096c31 2026-10-26 Add the ssbmf README with the quickstart, and the balance summary script`
- `24c34f80 2026-10-26 Expression-only real-data results; ebnm failures become rejected steps`
- `0611e62e 2026-10-26 Let the panel runner run the Part 4 balance arms`
- `cbae5952 2026-10-26 Record the w = 1 simulation panel results`
- `bacf0c29 2026-10-26 Render the ssbmf derivation and document its numerical safeguards`
- `865fd853 2026-10-26 Add the balance rules of rebuild Part 4 and the real-data runners`
- `7fa210b2 2026-10-26 Register the w = 1 simulation panel design and runner (rebuild Tasks 9-10)`
- `37c63b7e 2026-10-26 Handle programs that switch off; report the survival share of F precision`
- `1ae7f826 2026-10-26 Test the ssbmf updates, fitter, methods, simulator and CV (rebuild Tasks 6-9)`
- `3d48df5a 2026-10-26 Implement the ssbmf core and verify update stationarity (rebuild Tasks 3-5)`
- `90628e30 2026-10-26 Write the ssbmf model specification and check the Cox bound (rebuild Tasks 1-2)`
- `376f337e 2026-10-26 Scaffold the ssbmf package (rebuild Task 0)`

### `origin/gh-pages`

- `476f7315 2026-11-09 deploy: d635bab9ecf4081524037500e320c8bbf03a1313`
- `784e1488 2026-11-09 deploy: bee2ed7a6f875dc347eb4a829e9040e12ccc2e4b`
- `4449c6a0 2026-11-09 deploy: 39dbdb890a3f7e0ddaa03bdf8bdca3c3f6252abe`
- `b1215fc2 2026-11-09 deploy: 052a580bf0a4afd7b9cfa7bb0e10c26f00d9822e`
- `7746f450 2026-11-06 deploy: 8bbaeece7dea970f0f1787c78ec928b09538cf61`
- `75d49d65 2026-11-06 deploy: 7e43e2765ea8a867e8bd555c9caf5fb48bbee24d`
- `9b96b386 2026-11-06 deploy: f564bc944932338ea692c2821f0cad44b3a913de`
- `60b72c91 2026-11-06 deploy: d2e80977e3e36bf50051f37e0ae835f0ebc45da7`
- `2984605c 2026-11-06 deploy: 1facb21e9ec28f3dd7b050a7f30a7bd8bde3f2c0`
- `7243a0e6 2026-11-04 deploy: bea5a1d26bae7ad4b3d5aef6bfe2d6404897795d`
- `0648b7ea 2026-11-03 deploy: 3d700231d906584ab78eaa5e123ca423258f8d52`
- `0974fc5d 2026-11-03 Publish planted-factor recovery and comparator figures`
- `466fbc52 2026-11-03 Publish November 6 joint YFB verification and simulation figure`
- `ed2fa0a8 2026-11-03 Publish checklist re-verification caveat`
- `463474a3 2026-11-03 Publish read-only open decisions section`
- `d5463d24 2026-11-03 Publish matrix comparison of sequential and joint fitting`
- `54759d90 2026-11-03 Publish reviewed comparison and consistent sidebars from 4dee2ce`
- `b463dce2 2026-11-03 Publish reviewed current project state from bf4b6bb`
- `64b17b40 2026-11-02 deploy: 0d0dbcbf1b03f57868f2d7a92e2e0864e397edca`
- `eb8afe73 2026-11-02 deploy: 01c2ef10645ef0c179b4345c231b25d8128fdc59`
- `aa4dc79f 2026-10-30 deploy: f0828e538178cccb9dcba600cfc0322c5b749f04`
- `3d18bd59 2026-10-30 deploy: a8498df259ff2335b3f19d5ebb53a6dfa6bbb6e5`
- `0a6e5b2f 2026-10-30 deploy: 509823a2970baa8e61bba2e9021a8410eed74d1a`
- `2e6d2677 2026-10-29 deploy: 95e65c663700c3c27a3535c908c449f4f3dc2f01`
- `da8ed9c1 2026-10-25 deploy: bf1526217ae022f42f12aa8fb77dfbc4d902bf13`
- `4e0f13c9 2026-10-25 deploy: 223c87bd7520b74a35fe92158767461455ff232a`
- `b8e07784 2026-10-24 Publish progress book to GitHub Pages (2026-10-25T01:35:49Z)`

## Instruction and link follow-up

The initial scan found branch instructions in five historical prompts/plans:

- `docs/prompts/2026-10-30-downweighted-genomics-prompt.md`
- `docs/prompts/2026-10-30-codex-execution-prompt.md`
- `docs/prompts/2026-11-06-codex-handoff-prompts.md`
- `docs/plans/2026-11-05-chapter4-revamp-prompt.md`
- `docs/plans/2026-11-05-chapter4-revamp-plan.md`

Each now has a prominent routing note directing new work to main, while retaining
the original branch instructions as historical execution context. The active
ROADMAP URL in `docs/progress_book/current-state.qmd` now points to main.

The AGENTS pointers for the merged multimodal prototype, protected-pruning and
comparator evaluations, and training-only robustness runner now identify the
code as merged into main. Scientific limitations and non-adoption remain explicit.
The pointer to `rebuild/unified-ssbmf` remains valid as an unmerged reference.
Dated decisions and meeting records retain branch names as provenance. Any archive
tag or branch deletion remains a separate approval decision.

## Pointer updates and review

The active AGENTS pointers now identify the merged code on main, preserving the
scientific limitations. The current-state chapter's ROADMAP link now points to
main. The five prompts/plans listed above have a branch-routing note directing
new work to main while preserving their original execution context. No historical
commit reference or scientific decision was rewritten. The chapter source was
updated without rendering or publishing the book.

Independent review verified that all 52 local/tracking refs appear exactly once
in the 38 rows and match their pinned tips and ahead counts. Its request to
qualify old branch instructions in historical prompts/plans was addressed with
the routing notes. The live remote heads were checked separately against the
17 tracking refs. No branch cleanup, archive tag or branch deletion was performed.
The table is awaiting Andrew's approval.

## Slide-draft comparison after Andrew's clarification

Andrew asked to ensure the draft is merged onto main. Main commit `bc22bf24`
(November 3, “Add the Part III prelim deck and its bios-dissertation mirror”)
already incorporated the draft. A tree comparison shows identical deck source,
rendered deck, figures and slide tools. Only `OUTLINE.md` and `outline.html`
differ under `presentation/prelim-slides/`: main uses the newer derivation path
from `ce16573c`. The Chapter 4 sync tools on main also contain later rendering
and author-header corrections. No missing slide content was established.

The proposed disposition is therefore keep-reference until an archive tag is
approved, then delete the old branch refs. This would preserve the side history
without reverting newer paths or tools. No merge, tag or deletion was performed.
Andrew's goal of having the draft on main is already met. The full branch table
still awaits approval.
