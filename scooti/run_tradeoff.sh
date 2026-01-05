#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scooti/run_tradeoff.sh path/to/tradeoff_config.json

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_JSON="${1:-}"
if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/tradeoff_config.json" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[tradeoff] Running analysis with: $CONFIG_JSON"
bash "$SCRIPT_DIR/run_analyzer.sh" "$CONFIG_JSON"
echo "[tradeoff] Done. Outputs under $(command -v jq >/dev/null 2>&1 && jq -r '.save_root_path // empty' "$CONFIG_JSON" || echo '<see config>')"

