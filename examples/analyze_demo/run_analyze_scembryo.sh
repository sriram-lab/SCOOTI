#!/usr/bin/env bash
# Archived: moved under examples/_archive_runs/. Use `scooti -A` instead.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ARCH="$REPO_ROOT/examples/_archive_runs/analyze_demo_run_analyze_scembryo.sh"
if [[ -x "$ARCH" ]]; then
  exec "$ARCH" "$@"
else
  echo "Archived script not found: $ARCH" >&2
  exit 1
fi
