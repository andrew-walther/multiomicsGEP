# Rank-5 program variance shares (2026-11-15)

Source for the variance shares quoted in Chapter 4, Training Cross-Validation (committee item C-12).

Fit: `results/benchmark_sim/outputs/k5_working_2026-11-12/s11/fits/size_D4_N2218_w1_f0.rds`
(the full 2,218-gene joint SSBMF fit at w_S = 1, rank 5, 273 patients; the same file used by
`results/benchmark_sim/k5_gene_lists.R`).

Variance share of program k: v_k = ||E[l_k]||^2 ||E[f_k]||^2 / (n p), the definition used for the
survival-share approximation in the derivation (Y is standardized within cohort, so its total sum of
squares is about n p).

```r
s <- readRDS(".../size_D4_N2218_w1_f0.rds")$state
v <- colSums(s$EL^2) * colSums(s$EF^2) / (nrow(s$EL) * nrow(s$EF))
round(100 * v, 2)   # 14.61 10.97 6.48 4.91 3.61  (Programs 1 to 5; total 40.6%)
```

Survival-active Programs 2, 4 and 5: 11.0%, 4.9% and 3.6%. The earlier chapter statement that the real
programs carry "far more" variance than the simulated 0.25% to 5% overstated this; Program 5 is below
the largest simulated share.
