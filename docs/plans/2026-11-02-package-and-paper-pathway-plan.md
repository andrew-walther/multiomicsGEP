# From development to package and paper

Written 2026-11-02. This plan says when work moves from this repository into the `SSBMF`
package and the `SSBMF-paper` manuscript, and what moves, so that neither is put off until
"everything is done".

## The three repositories and the role of each

| Repository | Role | What belongs there |
|---|---|---|
| `multiomicsGEP` (this one) | Development workspace and lab notebook | Experiments, diagnostics, registered attempts, failed ideas, progress book. It may stay messy and is never published as is. |
| `SSBMF` | The R package | Only code that the working method needs, with tests, help pages and a vignette. One main fitting function, `ssbmf_fit()`, plus prediction, simulation and factor-classification helpers. |
| `SSBMF-paper` | The manuscript | Scripts that call the package to produce every simulation and real-data result, and a Quarto manuscript whose numbers, tables and figures are computed from those saved results. |

The package never depends on the paper or on this repository. The paper depends on a pinned
package version. Patient data is in none of the three: it stays in the git-ignored local `data/`
folder (and later in a shared lab location), and the paper repository has the scripts and the
instructions for obtaining TCGA, CPTAC and the validation cohorts.

## How code moves from development into the package

Develop and test ideas here. When a piece is shown to work, re-implement it cleanly in `SSBMF`
with tests, and check that it reproduces the development result on the same data. From then on,
results that use that piece come from the package, not from this repository's code.

## Stages and their triggers

Each stage starts when its trigger is met, not when the previous stage is "finished".

**Stage 1. Working single-modality supervision (now, in this repository).**
Diagnose and fix the YFB model's supervision failure. The next test is attempt S-9. Fall back to
η = Lβ with fold-in prediction only if YFB makes no progress.
*Trigger to move on:* a fitter passes Gate S in simulation and holds up in training
cross-validation on the PDAC data (convergence, survival-active programs that are not all false
positives, no worse cross-validated C than the unsupervised fit).

**Stage 2. Start the package (single modality).**
In `SSBMF`, write `ssbmf_fit()` for one modality from the written specification. Do not copy files
wholesale. Bring over only the update equations, the Cox working quantities, the priors and KL
terms, prediction, factor classification and a simulator. Add:
- unit tests, including parity tests showing that `ssbmf_fit()` reproduces the development fitter's
  results on the same simulated data;
- help pages for exported functions and one vignette on synthetic data;
- `R CMD check` with no errors or warnings, and a GitHub Actions workflow running it.

Tag the first usable version (v0.1.0).
*Trigger to move on:* v0.1.0 is tagged and passes check.

**Stage 3. Start the paper's simulation results.**
In `SSBMF-paper`, install the tagged package (record versions with `renv`), and write the
simulation scripts against `ssbmf_fit()`. Port the registered simulation designs from this
repository as scripts, not their old outputs. The manuscript's methods section can be written now
from the specification.
*Trigger to move on:* the simulation results regenerate from a clean checkout.

**Stage 4. Real-data results through the package.**
Training cross-validation, the single external validation, comparators (unsupervised EBMF + Cox,
DeSurv, PurIST, elastic-net Cox) and biological characterization, all calling the package. Hold
out at least one validation cohort from every model choice until the final claim.

**Stage 5. The M-modality model.**
Develop the multimodal extension here as before, with M = 1 as the single-modality case. When it
works, add it to `ssbmf_fit()` (for example, a list of data blocks), release a new package version,
and rerun the paper's analyses on that version.

**Stage 6. Release and submission.**
- Freeze package and paper versions, and archive both on Zenodo for DOIs.
- Decide between CRAN and Bioconductor before Stage 2's interface is final, because Bioconductor
  expects its own data structures.

## What stays in this repository

Experimental fitters, diagnostic controls, attempt logs, saved diagnostic results and the progress
book stay in this repository as the record of how the method was developed. They are not cleaned
up or moved. Reducing bloat means moving only what is needed, not reorganizing this repository.

## What to avoid

- Starting the package before Stage 1's trigger, which would package a model that does not
  supervise its loadings.
- Waiting for every result before Stage 2. The later analyses should already run through the
  package.
- Producing paper results from this repository's code after the package exists. That splits
  the method into two implementations that can drift apart.

## Dissertation and journal destinations

The prelim chapter lives in `paper/prelim/`. The final dissertation chapter
will live in `paper/thesis/` and be ported to bios-dissertation.
`SSBMF-paper` builds the journal submission from package outputs.

Reliable survival supervision of the molecular loadings is the top priority.
Start with the S-9 YFB test and its registered recovery, false-active and
convergence criteria. Investigate LB with fold-in prediction if YFB stalls.
Multimodal development is deferred except when a unified fitter run with
M = 1 offers a concrete route to resolving supervision. The single-modality
R package can start once the model meets its acceptance criteria.
