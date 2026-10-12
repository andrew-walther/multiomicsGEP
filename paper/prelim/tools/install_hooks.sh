#!/usr/bin/env bash
# ============================================================
# Script: install_hooks.sh
# Purpose: Install multiomicsGEP's post-commit hook, which keeps
#          bios-dissertation's Chapter 4 in sync (tools/post_commit_hook.sh).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-10-02
# Dependencies: bash, git
# ============================================================
# Adapted from SpatialCRT's installer. The hook goes in core.hooksPath when
# that is set (this repo sets it to .git/hooks), otherwise in the common hooks
# dir, so every worktree shares it. The installed file is a thin shim that runs
# tools/post_commit_hook.sh from whichever worktree made the commit, so the
# hook's logic stays versioned. Re-running is safe: it replaces only a hook it
# installed itself (recognised by the marker line) and refuses to overwrite
# anyone else's post-commit hook.

set -euo pipefail

MARKER="# multiomicsgep-chapter4-prelim-sync"
cd "$(dirname "${BASH_SOURCE[0]}")"

if hooks="$(git config --get core.hooksPath)"; then
  case "$hooks" in /*) ;; *) hooks="$(git rev-parse --show-toplevel)/$hooks" ;; esac
else
  hooks="$(git rev-parse --path-format=absolute --git-common-dir)/hooks"
fi
hook="$hooks/post-commit"
mkdir -p "$hooks"

if [ -e "$hook" ] && ! grep -q "^$MARKER\$" "$hook"; then
  echo "install_hooks: $hook already exists and was not installed by this script." >&2
  echo "  Not overwriting it. Add this line to it by hand instead:" >&2
  echo "    sh \"\$(git rev-parse --show-toplevel)/paper/prelim/tools/post_commit_hook.sh\"" >&2
  exit 1
fi

cat >"$hook" <<HOOK
#!/bin/sh
$MARKER
# Installed by paper/prelim/tools/install_hooks.sh.
# Runs the versioned hook logic from the worktree that made the commit; it never
# blocks or undoes the commit.
h="\$(git rev-parse --show-toplevel)/paper/prelim/tools/post_commit_hook.sh"
[ -f "\$h" ] && sh "\$h"
exit 0
HOOK
chmod +x "$hook"
echo "install_hooks: installed $hook"
