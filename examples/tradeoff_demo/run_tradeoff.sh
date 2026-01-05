#!/usr/bin/env bash
# Archived: moved under examples/_archive_runs/. Use `scooti -T ./examples/tradeoff_demo/tradeoff_config.json` instead.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ARCH="$REPO_ROOT/examples/_archive_runs/tradeoff_demo_run_tradeoff_legacy.sh"
if [[ -x "$ARCH" ]]; then
  exec "$ARCH" "$@"
else
  echo "Archived script not found: $ARCH" >&2
  exit 1
fi
