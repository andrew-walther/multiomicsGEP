# Robustness execution status

Exit code: 1, after all final tables were written and printed.

Final console lines:

```text
Warning message:
In run_supervised_robustness(if (length(a)) as.integer(a[1]) else 6L) :
  Simulation warnings are retained in all_starts.csv
Error: unexpected 'else' in "nteger(a[1])else"
Execution halted
```

The aggregate warning is spurious: all 880 saved warning fields are empty/NA.
A diagnostic-only guard was changed in the source while the Rscript process was
running, likely shifting its source read position and causing the trailing error.
The edit changed no fit or metric calculations; no tuning or rerun was attempted.

Independent checks after the error: current runner parses; full suite 550/0;
880 distinct fit rows, 40 successful dataset statuses, 40 selected starts;
30 nonconverged starts; exactly 20 distinct feature folds with finite C and
four full-training fits, all converged. Feature means, SEs and cap counts
recomputed from fold rows agree. Preserved CV prefix CSVs match the branch
baseline byte for byte. All final result tables exist.

The nonzero exit remains part of the result and must not be called a clean run.
