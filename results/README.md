# results/ — Directory Layout

## Active structure

```
results/
  benchmark_sim/         ← canonical home for all formal benchmark work
    run_LB_benchmark.R   ← Cluster A runner (Lβ linear predictor)
    run_YFB_benchmark.R  ← Cluster B runner (YFβ / Cox-on-YF reformulation)
    run_phase1_diagnostics.R  ← loading heatmaps (merged / tcga_only / cptac_only × prior)
    archive/             ← 9 retired scripts (see archive/README.md for details)
    outputs/             ← benchmark CSV and figure outputs
      LB_benchmark/      ← Cluster A outputs (untracked)
      YFB_benchmark/     ← Cluster B outputs (untracked)

  figures/               ← active cohort-level figure outputs
    CPTAC/, TCGA_PAAD/, Dijk/, Moffitt_GEO_array/,
    PACA_AU_array/, PACA_AU_seq/, PDAC_pooled_rnaseq/, Puleo_array/

  tables/                ← active cohort-level table outputs
    CPTAC/, CPTAC_Keff/, CPTAC_pl/, CPTAC_pl_Keff/,
    Dijk/, Dijk_Keff/, Dijk_pl/, Dijk_pl_Keff/,
    Moffitt_GEO_array/, Moffitt_GEO_array_Keff/, ...
    PDAC_cross_dataset/, PDAC_pooled_rnaseq/, TCGA_PAAD/, ...

  legacy/                ← deprecated simulation generations; not re-run
    full_sim/            ← V1 monolithic runner + report (uses Supervised_Bayesian_MF_V2.R)
    modular_sim_block/   ← block-wise CAVI runner + report (superseded by factor-wise)
    modular_sim_factor/  ← exploratory factor-wise runner + PDAC / synthetic reports
    figures/
      full_sim/, modular_sim/, synthetic/
    tables/
      full_sim/, modular_sim/, synthetic/
```

## Benchmark reports

Formal results reports live in `docs/reports/` with a `_MM_DD_YY` date suffix so versions
are unambiguous as the benchmark evolves.

| Report | Description |
|--------|-------------|
| `docs/reports/ssbmf_summary_report_04_29_26.{qmd,pdf,html}` | Frozen DeSurv benchmark (Cluster A, 2026-04-29). Internal paths are stale post-move — static reference only. |
| `docs/reports/ssbmf_summary_report_MM_DD_YY.qmd` | Future LB/YFB benchmark report will go here. |

## Where new outputs land

- **Cluster A (LB) outputs** → `results/benchmark_sim/outputs/LB_benchmark/`
- **Cluster B (YFB) outputs** → `results/benchmark_sim/outputs/YFB_benchmark/`
- **New benchmark reports** → `docs/reports/ssbmf_summary_report_MM_DD_YY.{qmd,pdf,html}`
- **New figure outputs** → `results/figures/<cohort>/`
- **New table outputs** → `results/tables/<cohort>/`

## What NOT to touch

- `results/legacy/` — retired outputs; no active scripts depend on them
- `results/benchmark_sim/archive/` — retired runner scripts; see `archive/README.md`
