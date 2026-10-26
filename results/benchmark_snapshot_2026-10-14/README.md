# Benchmark snapshot: the approved interim write-up state (2026-10-14)

The current best state of the project, against which every new result is compared. Git tag
`benchmark-interim-writeup-2026-10-14` (commit `a352f17` on `main`) holds the code and the
approved write-up (`paper/writeup/`).

| Model | Configuration | Reported result | Files here |
|---|---|---|---|
| Single modality (D4) | YFB, per-platform z-std, combined-rank top 3,000 genes per cohort, K = 7, α_F = 0, point-exponential L/F | Mean external C 0.627 over 5 cohorts | `single_modality/` (fits, risk scores, gene names, result tables) |
| Multimodal (10/1 framework) | Intercept, point-Laplace F, ebnm priors, ELBO pruning, per-factor EB β, K_init = 7 | ICGC primary C 0.685 (same donors: single modality 0.685) | `multimodal/` (K = 7 fit, expression-only fit, sweep and program tables) |

Rules:
- These files are read-only copies. No runner writes here.
- `MANIFEST.sha256` (tracked in git) detects any change: `shasum -a 256 -c MANIFEST.sha256`.
- The `.rds` fits are git-ignored (they derive from controlled data). Back this folder up outside
  the repository, for example to OneDrive.
- A new configuration replaces the benchmark only if it is at least as good on every project goal
  and better on at least one, with Andrew's approval.
