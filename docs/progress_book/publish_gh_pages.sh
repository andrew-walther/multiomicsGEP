#!/usr/bin/env bash
# ============================================================
# Script: docs/progress_book/publish_gh_pages.sh
# Purpose: Render a local preview of the progress-book Quarto site.
#
# Local preview. The GitHub Actions workflow handles publication.
#
# Usage:   bash docs/progress_book/publish_gh_pages.sh
# Author: Andrew Walther
# Created: 2026-09-04
# ============================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BOOK_DIR="$REPO_ROOT/docs/progress_book"
QUARTO="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"
cd "$BOOK_DIR"
"$QUARTO" render
