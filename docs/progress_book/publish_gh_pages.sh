#!/usr/bin/env bash
# ============================================================
# Script: docs/progress_book/publish_gh_pages.sh
# Purpose: Render the progress-book Quarto site and publish it to the
#          gh-pages branch, so the rendered HTML is browsable at
#          https://andrew-walther.github.io/multiomicsGEP/ without anyone
#          needing to clone the repo and render it themselves.
#
#          `quarto publish gh-pages` was not used here because it requires
#          an interactively-created _publish.yml on first use, which isn't
#          practical from a non-interactive session; this script does the
#          equivalent manually: render, then push _book/'s contents to an
#          orphan gh-pages branch via a throwaway clone (never touches this
#          working copy's git state).
#
#          _book/ is intentionally gitignored on main (build artifact) --
#          only gh-pages carries the rendered output.
#
# Usage:   bash docs/progress_book/publish_gh_pages.sh
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# ============================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BOOK_DIR="$REPO_ROOT/docs/progress_book"
REMOTE_URL="$(git -C "$REPO_ROOT" remote get-url origin)"
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT

echo "--- Rendering $BOOK_DIR ---"
(cd "$BOOK_DIR" && quarto render)

echo "--- Cloning into throwaway workspace ---"
git clone --depth 1 --no-checkout "$REMOTE_URL" "$SCRATCH/publish"
cd "$SCRATCH/publish"
git checkout --orphan gh-pages
git reset --hard
cp -r "$BOOK_DIR/_book/." .
git add -A
git commit -m "Publish progress book to GitHub Pages ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
git push origin gh-pages --force

echo ""
echo "Published: https://andrew-walther.github.io/multiomicsGEP/"
