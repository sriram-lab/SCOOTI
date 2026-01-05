#!/usr/bin/env bash
# Archived: use `scooti -A ./examples/analyze_demo/analyze_config.json` instead.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_JSON="${1:-$REPO_ROOT/examples/analyze_demo/analyze_config.json}"
bash "$REPO_ROOT/scooti/run_analyze.sh" "$CONFIG_JSON"

