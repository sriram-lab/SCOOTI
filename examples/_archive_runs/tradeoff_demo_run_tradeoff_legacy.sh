#!/usr/bin/env bash
# Archived: use `scooti -T ./examples/tradeoff_demo/tradeoff_config.legacy.json` instead.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_JSON="${1:-$REPO_ROOT/examples/tradeoff_demo/tradeoff_config.legacy.json}"
bash "$REPO_ROOT/scooti/run_tradeoff.sh" "$CONFIG_JSON"

