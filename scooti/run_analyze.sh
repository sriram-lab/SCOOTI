#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scooti/run_analyze.sh path/to/analyze_config.json

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_JSON="${1:-}"
if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/analyze_config.json" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

bash "$SCRIPT_DIR/run_analyzer.sh" "$CONFIG_JSON"

