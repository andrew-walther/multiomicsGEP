# Benchmark Archive

Historical diagnostic and one-off scripts retained for reproducibility. These
scripts were used during Cluster A and Cluster B development and are superseded
by the permanent benchmark runners in `results/benchmark_sim/`.

## Superseded by permanent runners

| Script | Replaced by |
|--------|-------------|
| `run_cluster_a_smoke.R` | `run_LB_benchmark.R` |
| `run_cluster_a_external.R` | `run_LB_benchmark.R` |
| `run_cox_on_yf_benchmark.R` | `run_YFB_benchmark.R` |
| `run_ssbmf_benchmark.R` | `run_LB_benchmark.R` + `run_YFB_benchmark.R` (unified runner split into LB/YFB pair) |
| `run_lambda_sweep.R` | N/A — one-off evaluation; lambda ∈ {1, p/n, 2p/n} tested; lambda=1 retained as default |

## One-off diagnostics (β=0 investigation)

| Script | Purpose |
|--------|---------|
| `run_ebmf_diagnostic.R` | EBMF unsupervised fit — Cox-associates EBMF factors with survival; used to identify which factors carry survival signal before CAVI |
| `run_ebmf_warmstart.R` | Warm-start experiments: Exp1 (β-only, EL/EF fixed from EBMF) + Exp2 (full CAVI from EBMF init) |
| `compute_ph_diagnostics.R` | Proportional-hazards assumption checks (Schoenfeld residuals) |
| `sandbox_lambda_test.R` | Lambda multiplier sensitivity: tested lambda ∈ {1, p/n, 2p/n}; lambda=1 best or tied |

These scripts are not actively maintained but can be reproduced from git history.
