#!/bin/sh
# ============================================================
# Script: post_commit_hook.sh
# Purpose: multiomicsGEP post-commit logic. If HEAD touched the Chapter 4
#          sources, run sync_to_prelim.sh. Never blocks, fails, or undoes the
#          multiomicsGEP commit (the commit already exists when this runs).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-10-02
# Dependencies: POSIX sh, git
# ============================================================
# Installed as a thin shim by tools/install_hooks.sh; the shim calls this file
# from the current worktree, so the logic stays versioned here.
#
# Part III slides: a commit touching presentation/prelim-slides/ also runs
# presentation/prelim-slides/tools/sync_slides_to_prelim.sh (added 2026-11-03).
#
# Chapter 4 trigger paths (relative to the multiomicsGEP root):
#   paper/prelim/project3-ssbmf.qmd
#   paper/prelim/figures/
#   paper/prelim/tools/
#   paper/ssbmf-refs.bib
# Output is appended to tools/sync.log (gitignored).

root=$(git rev-parse --show-toplevel 2>/dev/null) || exit 0
tools="$root/paper/prelim/tools"
log="$tools/sync.log"
sync="$tools/sync_to_prelim.sh"

# Files changed by HEAD (first parent for merges; --root for a first commit).
changed=$(git diff-tree --no-commit-id --name-only -r --root -m --first-parent HEAD 2>/dev/null)

# Part III prelim slides: mirror presentation/prelim-slides/ into
# bios-dissertation/slides/prelim/part-3/ (added 2026-11-03). Runs before the
# Chapter 4 check below and, like it, never blocks the commit.
slides_sync="$root/presentation/prelim-slides/tools/sync_slides_to_prelim.sh"
if echo "$changed" | grep -q '^presentation/prelim-slides/'; then
  if [ -d "$(git rev-parse --git-path rebase-merge)" ] || [ -d "$(git rev-parse --git-path rebase-apply)" ]; then
    echo "post-commit: rebase in progress; Part III slide sync skipped. Run $slides_sync afterwards." >&2
  elif [ ! -x "$slides_sync" ]; then
    echo "post-commit: WARNING: Part III slide sync script not found or not executable: $slides_sync" >&2
  else
    echo "post-commit: syncing Part III slides to bios-dissertation (log: $log) ..." >&2
    {
      echo "=== $(date '+%Y-%m-%d %H:%M:%S') multiomicsGEP slides $(git rev-parse --short HEAD)"
      "$slides_sync"
      echo "=== exit $?"
    } >>"$log" 2>&1
    case "$(tail -1 "$log")" in
      "=== exit 0") echo "post-commit: $(tail -2 "$log" | head -1)" >&2 ;;
      *) echo "post-commit: WARNING: Part III slide sync FAILED; bios-dissertation was not updated. See $log, then run $slides_sync by hand." >&2 ;;
    esac
  fi
fi

echo "$changed" | grep -q -E '^paper/(prelim/(project3-ssbmf\.qmd$|figures/|tools/)|ssbmf-refs\.bib$)' || exit 0

sha=$(git rev-parse --short HEAD 2>/dev/null)

warn() {
  {
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!! Chapter 4 prelim sync FAILED for multiomicsGEP commit $sha."
    echo "!! $1"
    echo "!! Your multiomicsGEP commit is fine; bios-dissertation was NOT updated."
    echo "!! Fix the problem, then run by hand:"
    echo "!!   $sync"
    echo "!! Log: $log"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""
  } >&2
}

# Mid-rebase every replayed commit would trigger a sync; skip and say so.
if [ -d "$(git rev-parse --git-path rebase-merge)" ] || [ -d "$(git rev-parse --git-path rebase-apply)" ]; then
  echo "post-commit: rebase in progress; Chapter 4 sync skipped. Run $sync afterwards." >&2
  exit 0
fi

if [ ! -x "$sync" ]; then
  warn "sync script not found or not executable: $sync"
  exit 0
fi

echo "post-commit: $sha touched Chapter 4 sources; syncing to bios-dissertation (log: $log) ..." >&2
{
  echo "=== $(date '+%Y-%m-%d %H:%M:%S') multiomicsGEP $sha"
  "$sync"
  rc=$?
  echo "=== exit $rc"
} >>"$log" 2>&1

last=$(tail -2 "$log" | head -1)
case "$(tail -1 "$log")" in
  "=== exit 0") echo "post-commit: $last" >&2 ;;
  *)
    # Quote this run's error lines (from its "===" header onward), since the
    # sync's own stderr went to the log, not the terminal.
    reason=$(awk '/^=== .* multiomicsGEP /{buf=""} {buf=buf $0 "\n"} END{printf "%s", buf}' "$log" |
             grep -E 'ERROR|FAILED|^  - ' | sed 's/^/!!   /')
    warn "Reason (from the log):
$reason" ;;
esac
exit 0
