#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

CONFIG_JSON="${1:-$SCRIPT_DIR/tradeoff_config.legacy.json}"
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[tradeoff-demo:legacy] Running analysis with: $CONFIG_JSON"
bash "$REPO_ROOT/scooti/run_analyzer.sh" "$CONFIG_JSON"
echo "[tradeoff-demo:legacy] Done. Outputs under $(jq -r '.save_root_path' "$CONFIG_JSON")"

