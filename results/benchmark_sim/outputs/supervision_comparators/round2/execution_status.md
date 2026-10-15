# Task4 environment failure before any dataset

Commit dd11923 preregistered the fixed40-dataset/five-arm panel.
The sandboxed launch exited1 while opening the PSOCK listener:
`creation of server socket failed: port 11504 cannot be opened`.
Only design.csv was written; no dataset or fit ran. Preserve this folder.
Retry the exact committed code/design with host socket access and fresh
output `results/benchmark_sim/outputs/supervision_comparators/round2_retry1/`.
No statistical setting, arm, seed, threshold or default changes.
