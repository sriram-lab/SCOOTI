#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../" && pwd)"
CONFIG_JSON="${1:-$SCRIPT_DIR/analyze_reproduce_config.json}"

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[analyze-reproduce] Running analysis with: $CONFIG_JSON"
# Quick preflight: show resolved coef CSV if directory provided
if command -v jq >/dev/null 2>&1; then
  COEF_DIR=$(jq -r '.coef_paths.exp' "$CONFIG_JSON")
  if [[ -d "$COEF_DIR" ]]; then
    CAND=$(ls "$COEF_DIR"/*.csv 2>/dev/null | head -n 1 || true)
    if [[ -n "$CAND" ]]; then
      echo "[analyze-reproduce] Found coef CSV: $CAND"
    else
      echo "[analyze-reproduce][WARN] No CSV found under $COEF_DIR; analysis_runner will attempt to resolve automatically."
    fi
  fi
fi
bash "$REPO_ROOT/scooti/run_analyzer.sh" "$CONFIG_JSON"
echo "[analyze-reproduce] Done. Outputs under $(jq -r '.save_root_path' "$CONFIG_JSON")"
