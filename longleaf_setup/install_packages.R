# =============================================================================
# longleaf_setup/install_packages.R
#
# One-time R package installation for multiomicsGEP on Longleaf.
#
# Usage (in an interactive session):
#   srun -p interact -n 1 --mem=4G -t 00:30:00 --pty bash
#   module add r/4.4.0
#   Rscript longleaf_setup/install_packages.R
#   exit
# =============================================================================

cat("Installing multiomicsGEP dependencies...\n\n")

pkgs <- c("ebnm", "survival")

install.packages(pkgs, repos = "https://cloud.r-project.org")

cat("\n--- Verification ---\n")
for (pkg in pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-12s %s\n", pkg, if (ok) "OK" else "FAILED"))
}
cat("\nDone.\n")
