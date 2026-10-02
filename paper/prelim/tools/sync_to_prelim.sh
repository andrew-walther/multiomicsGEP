#!/usr/bin/env bash
# ============================================================
# Script: sync_to_prelim.sh
# Purpose: Regenerate bios-dissertation's Chapter 4 draft from
#          project3-ssbmf.qmd (the ONLY file anyone edits), render it,
#          and commit the generated copy in bios-dissertation.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-25
# Dependencies: bash, git, python3 (stdlib only), RStudio's bundled Quarto,
#               TinyTeX
# ============================================================
#
# Usage:  tools/sync_to_prelim.sh [--no-render]
#
# Normally run by the multiomicsGEP post-commit hook (tools/install_hooks.sh);
# safe to run by hand at any time.
#
# Steps:
#   1. Bib check + transform (prelim_transform.py). On any citekey or bib-field
#      mismatch it exits non-zero BEFORE writing anything.
#   2. Write draft/project3-ssbmf-draft.qmd, sync draft/figures/, copy
#      manuscript-declarations.md, make sure the bios-prelim.cls symlink exists.
#   3. Render the draft with RStudio's Quarto (skipped with --no-render).
#   4. Stage ONLY the Chapter 4 paths under
#      prelim/project-proposals/project3-ssbmf/, explicitly by path, and
#      commit them with "Sync Chapter 4 from multiomicsGEP <short-sha>" (only if
#      something changed). The rendered draft/project3-ssbmf-draft.pdf is
#      added only when this run rendered it AND some other Chapter 4 path
#      changed (or the PDF is not yet committed); a re-render whose bytes
#      differ only by timestamps is discarded by restoring the committed PDF.
#      `git commit -- <paths>` commits those paths alone, even if something
#      else happens to be staged. Never pushes.
#   5. Print a one-line summary.
#
# Environment:
#   BIOS_DISSERTATION  bios-dissertation checkout (default
#                      ~/GithubProjects/bios-dissertation); override for tests.
# ============================================================

set -euo pipefail

NO_RENDER=0
for arg in "$@"; do
  case "$arg" in
    --no-render) NO_RENDER=1 ;;
    -h|--help) sed -n '2,35p' "$0"; exit 0 ;;
    *) echo "sync_to_prelim: unknown argument: $arg" >&2; exit 64 ;;
  esac
done

fail() { echo "sync_to_prelim: FAILED: $*" >&2; exit 1; }

# Paths -------------------------------------------------------------------
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHAPTER_DIR="$(dirname "$TOOLS_DIR")"
PAPER_DIR="$(dirname "$CHAPTER_DIR")"
CHAPTER_QMD="$CHAPTER_DIR/project3-ssbmf.qmd"
CHAPTER_BIB="$PAPER_DIR/ssbmf-refs.bib"
DECLARATIONS="$CHAPTER_DIR/manuscript-declarations.md"
SPATIAL_ROOT="$(git -C "$CHAPTER_DIR" rev-parse --show-toplevel)"
CHAPTER_REL="${CHAPTER_QMD#"$SPATIAL_ROOT"/}"
SHA="$(git -C "$SPATIAL_ROOT" rev-parse --short HEAD)"

# Inside a git hook, GIT_DIR / GIT_INDEX_FILE etc. point at multiomicsGEP; clear
# them so every git call below acts only on the repo named with -C.
unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_PREFIX GIT_OBJECT_DIRECTORY

BIOS="${BIOS_DISSERTATION:-$HOME/GithubProjects/bios-dissertation}"
TARGET_REL="prelim/project-proposals/project3-ssbmf"
TARGET="$BIOS/$TARGET_REL"
DRAFT_DIR="$TARGET/draft"
MASTER_BIB="$BIOS/prelim/references.bib"
QUARTO="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"

[ -f "$CHAPTER_QMD" ] || fail "chapter not found: $CHAPTER_QMD"
[ -d "$BIOS/.git" ] || [ -f "$BIOS/.git" ] || fail "not a git checkout: $BIOS"
[ -d "$TARGET" ] || fail "target directory missing: $TARGET"
[ -f "$MASTER_BIB" ] || fail "master bib missing: $MASTER_BIB"

dirty="$(git -C "$SPATIAL_ROOT" status --porcelain -- "$CHAPTER_QMD" "$CHAPTER_DIR/figures" \
         "$TOOLS_DIR" "$CHAPTER_BIB" "$DECLARATIONS" | grep -v 'tools/sync\.log$' || true)"
if [ -n "$dirty" ]; then
  echo "sync_to_prelim: WARNING: multiomicsGEP inputs have uncommitted changes; the copy" >&2
  echo "  is generated from the working tree, not exactly from $SHA:" >&2
  echo "$dirty" | sed 's/^/    /' >&2
fi

# 1-2. Bib check, transform, figures ---------------------------------------
RESULT="$(python3 "$TOOLS_DIR/prelim_transform.py" \
  --chapter "$CHAPTER_QMD" --chapter-bib "$CHAPTER_BIB" --master-bib "$MASTER_BIB" \
  --header "$TOOLS_DIR/prelim_header.yml" --out-dir "$DRAFT_DIR" \
  --source-label "multiomicsGEP/$CHAPTER_REL")" || fail "bib check / transform (see message above); nothing written"

# Project 3 has no declarations file yet; copy it once one exists
[ ! -f "$DECLARATIONS" ] || cp "$DECLARATIONS" "$TARGET/manuscript-declarations.md"
if [ ! -L "$DRAFT_DIR/bios-prelim.cls" ]; then
  ln -s ../../../../templates/skeleton/prelim/bios-prelim.cls "$DRAFT_DIR/bios-prelim.cls"
fi
[ -f "$DRAFT_DIR/bios-prelim.cls" ] || fail "bios-prelim.cls symlink does not resolve"

# 3. Render ----------------------------------------------------------------
PAGES="not rendered"
if [ "$NO_RENDER" -eq 0 ]; then
  [ -x "$QUARTO" ] || fail "RStudio Quarto not found at $QUARTO"
  TINYTEX_BIN="$(ls -d "$HOME"/Library/TinyTeX/bin/* 2>/dev/null | head -1)"
  [ -n "$TINYTEX_BIN" ] || fail "TinyTeX not found under ~/Library/TinyTeX/bin"
  export PATH="$TINYTEX_BIN:$PATH"
  RENDER_LOG="$(mktemp -t sync_to_prelim_render)"
  if ! (cd "$DRAFT_DIR" && "$QUARTO" render project3-ssbmf-draft.qmd) >"$RENDER_LOG" 2>&1; then
    tail -40 "$RENDER_LOG" >&2
    fail "quarto render failed (full log: $RENDER_LOG); nothing committed"
  fi
  # Unresolved citations / references exit 0, so check the output explicitly.
  if grep -i -E 'citeproc: reference .* not found|unable to resolve|undefined references|Reference .* undefined|Citation .* undefined' "$RENDER_LOG" >&2; then
    fail "render reported unresolved citations or references (log: $RENDER_LOG); nothing committed"
  fi
  PDF="$DRAFT_DIR/project3-ssbmf-draft.pdf"
  if command -v pdftotext >/dev/null 2>&1; then
    n_qq="$(pdftotext "$PDF" - | grep -c '??' || true)"
    [ "$n_qq" = "0" ] || fail "rendered PDF contains '??' $n_qq time(s) (an unresolved \\ref); nothing committed"
  fi
  if command -v pdfinfo >/dev/null 2>&1; then
    PAGES="$(pdfinfo "$PDF" | awk '/^Pages:/{print $2}') pp"
  else
    PAGES="rendered"
  fi
  rm -f "$RENDER_LOG"
fi

# 4. Commit the Chapter 4 paths only ----------------------------------------
CANDIDATES=(
  "$TARGET_REL/README.md"
  "$TARGET_REL/project3-ssbmf.qmd"
  "$TARGET_REL/manuscript-declarations.md"
  "$TARGET_REL/.gitignore"
  "$TARGET_REL/draft/project3-ssbmf-draft.qmd"
  "$TARGET_REL/draft/bios-prelim.cls"
  "$TARGET_REL/draft/figures"
)
# The rendered PDF is handled separately below: it is committed only when this
# run rendered it (so a --no-render sync never commits a PDF that predates the
# .qmd it just wrote) and only when something else changed (so a re-render
# that differs only in embedded timestamps does not bloat the repo).
PDF_REL="$TARGET_REL/draft/project3-ssbmf-draft.pdf"
[ "$NO_RENDER" -eq 1 ] || CANDIDATES+=("$PDF_REL")
PATHS=()
for p in "${CANDIDATES[@]}"; do
  if [ -e "$BIOS/$p" ] || [ -L "$BIOS/$p" ] || [ -n "$(git -C "$BIOS" ls-files -- "$p")" ]; then
    PATHS+=("$p")
  fi
done

# Checked before `git add`: files inside an ignored directory are skipped
# silently, and an explicitly named ignored file (the rendered PDF) makes
# `git add` exit with a less clear message. Either way a figure or the PDF that
# some bios-dissertation ignore rule matches (e.g. a *.pdf rule coming back)
# would not be committed.
ignored="$(git -C "$BIOS" ls-files --others --ignored --exclude-standard -- "${PATHS[@]}")"
[ -z "$ignored" ] || fail "files are gitignored in bios-dissertation and would not be committed: $ignored"
# Split the PDF off from the other paths.
OTHER=()
HAVE_PDF=0
for p in "${PATHS[@]}"; do
  if [ "$p" = "$PDF_REL" ]; then HAVE_PDF=1; else OTHER+=("$p"); fi
done
git -C "$BIOS" add -A -- "${OTHER[@]}"
PDF_NOTE=""
if [ "$HAVE_PDF" -eq 1 ]; then
  if ! git -C "$BIOS" diff --cached --quiet -- "${OTHER[@]}" ||
     ! git -C "$BIOS" cat-file -e "HEAD:$PDF_REL" 2>/dev/null; then
    git -C "$BIOS" add -A -- "$PDF_REL"
  else
    # Nothing but the PDF would change: restore the committed PDF (index and
    # working tree, this one path only) and commit nothing.
    git -C "$BIOS" checkout HEAD -- "$PDF_REL"
    PATHS=("${OTHER[@]}")
    PDF_NOTE=" (PDF-only re-render discarded; committed PDF restored)"
  fi
fi
BIOS_COMMIT="no commit (nothing changed)$PDF_NOTE"
if ! git -C "$BIOS" diff --cached --quiet -- "${PATHS[@]}"; then
  git -C "$BIOS" commit -q -m "Sync Chapter 4 from multiomicsGEP $SHA" -- "${PATHS[@]}"
  # Belt and braces: the commit must touch nothing outside the target dir.
  outside="$(git -C "$BIOS" show --name-only --format= HEAD | grep -v "^$TARGET_REL/" || true)"
  [ -z "$outside" ] || fail "commit $(git -C "$BIOS" rev-parse --short HEAD) touched paths outside $TARGET_REL: $outside"
  BIOS_COMMIT="bios-dissertation $(git -C "$BIOS" rev-parse --short HEAD)"
fi

# 5. Summary -----------------------------------------------------------------
# RESULT looks like: RESULT qmd=changed figs_copied=N figs_removed=M keys=K
echo "sync_to_prelim: OK: Chapter 4 from multiomicsGEP $SHA -> $BIOS_COMMIT; ${RESULT#RESULT }; $PAGES"
