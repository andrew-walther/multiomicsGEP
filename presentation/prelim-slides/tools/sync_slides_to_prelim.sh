#!/usr/bin/env bash
# ============================================================
# Script: sync_slides_to_prelim.sh
# Purpose: Mirror the Part III prelim deck (presentation/prelim-slides/) into
#          bios-dissertation/slides/prelim/part-3/ and commit the mirrored
#          sources there. Never pushes.
# Author: Andrew Walther
# Created: 2026-11-03
# Dependencies: bash, git, rsync
# ============================================================
#
# Usage:  presentation/prelim-slides/tools/sync_slides_to_prelim.sh
#
# Normally run by the multiomicsGEP post-commit hook
# (paper/prelim/tools/post_commit_hook.sh) when a commit touches
# presentation/prelim-slides/; safe to run by hand at any time.
#
# What is mirrored into part-3/ (part-3/README.md belongs to bios-dissertation
# and is never touched):
#   ssbmf_talk.qmd, OUTLINE.md, theme/ and figures/ (exact mirrors)
#       -> copied and committed
#   ssbmf_talk.html, outline.html
#       -> copied for local viewing only; bios-dissertation ignores *.html
# The commit names only these paths ("Sync Part III slides from multiomicsGEP
# <sha>") and is made only when something changed. Joining Part III into the
# combined prelim deck (dividers, footer prefixes, backup section) is a
# separate manual step.
#
# Environment:
#   BIOS_DISSERTATION  bios-dissertation checkout (default
#                      ~/GithubProjects/bios-dissertation); override for tests.
# ============================================================
set -euo pipefail

fail() { echo "sync_slides_to_prelim: FAILED: $*" >&2; exit 1; }

# Paths ----------------------------------------------------------------------
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DECK_DIR="$(dirname "$TOOLS_DIR")"
ROOT="$(git -C "$DECK_DIR" rev-parse --show-toplevel)"
SHA="$(git -C "$ROOT" rev-parse --short HEAD)"

# Inside a git hook, GIT_DIR / GIT_INDEX_FILE point at multiomicsGEP; clear them
# so every git call below acts only on the repository named with -C.
unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_PREFIX GIT_OBJECT_DIRECTORY

BIOS="${BIOS_DISSERTATION:-$HOME/GithubProjects/bios-dissertation}"
TARGET_REL="slides/prelim/part-3"
TARGET="$BIOS/$TARGET_REL"

[ -d "$BIOS/.git" ] || [ -f "$BIOS/.git" ] || fail "not a git checkout: $BIOS"
[ -d "$TARGET" ] || fail "target directory missing: $TARGET"
for f in ssbmf_talk.qmd OUTLINE.md theme figures; do
  [ -e "$DECK_DIR/$f" ] || fail "missing deck input: $DECK_DIR/$f"
done

dirty="$(git -C "$ROOT" status --porcelain -- "$DECK_DIR")"
if [ -n "$dirty" ]; then
  echo "sync_slides_to_prelim: WARNING: deck has uncommitted changes; the copy is taken" >&2
  echo "  from the working tree, not exactly from $SHA:" >&2
  echo "$dirty" | sed 's/^/    /' >&2
fi

# Copy ---------------------------------------------------------------------------
mkdir -p "$TARGET/theme" "$TARGET/figures"
cp "$DECK_DIR/ssbmf_talk.qmd" "$DECK_DIR/OUTLINE.md" "$TARGET/"
# Exact mirrors, so a file removed from the deck is removed here too. theme/
# holds the shared prelim theme copies and deck.css, so part-3/ renders on its own.
rsync -a --delete "$DECK_DIR/theme/" "$TARGET/theme/"
rsync -a --delete "$DECK_DIR/figures/" "$TARGET/figures/"
# Rendered pages for local viewing (ignored by bios-dissertation's *.html rule)
for h in ssbmf_talk.html outline.html; do
  [ ! -f "$DECK_DIR/$h" ] || cp "$DECK_DIR/$h" "$TARGET/$h"
done

# Commit the mirrored sources only ---------------------------------------------
PATHS=("$TARGET_REL/ssbmf_talk.qmd" "$TARGET_REL/OUTLINE.md" "$TARGET_REL/theme" "$TARGET_REL/figures")
ignored="$(git -C "$BIOS" ls-files --others --ignored --exclude-standard -- "${PATHS[@]}")"
[ -z "$ignored" ] || fail "files are gitignored in bios-dissertation and would not be committed: $ignored"
git -C "$BIOS" add -A -- "${PATHS[@]}"

RESULT="no commit (nothing changed)"
if ! git -C "$BIOS" diff --cached --quiet -- "${PATHS[@]}"; then
  git -C "$BIOS" commit -q -m "Sync Part III slides from multiomicsGEP $SHA" -- "${PATHS[@]}"
  # The commit must touch nothing outside the Part III folder.
  outside="$(git -C "$BIOS" show --name-only --format= HEAD | grep -v "^$TARGET_REL/" || true)"
  [ -z "$outside" ] || fail "commit $(git -C "$BIOS" rev-parse --short HEAD) touched paths outside $TARGET_REL: $outside"
  RESULT="bios-dissertation $(git -C "$BIOS" rev-parse --short HEAD)"
fi
n_fig="$(find "$TARGET/figures" -type f | wc -l | tr -d ' ')"
echo "sync_slides_to_prelim: OK: Part III slides from multiomicsGEP $SHA -> $RESULT; $n_fig figures"
