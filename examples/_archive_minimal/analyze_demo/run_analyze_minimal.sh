#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONFIG_JSON="${1:-$SCRIPT_DIR/analyze_config.minimal.json}"
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[analyze-demo:minimal] Running analysis with: $CONFIG_JSON"
bash "$REPO_ROOT/scooti/run_analyzer.sh" "$CONFIG_JSON"
echo "[analyze-demo:minimal] Done. Outputs under $(jq -r '.save_root_path' "$CONFIG_JSON")"

