#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scooti/run_batch.sh <samples> <config.json>
# Example: bash scooti/run_batch.sh 10 examples/run_flux/constrained_demo_config.json

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if ! command -v jq >/dev/null 2>&1; then
  echo "Error: jq is required for run_batch.sh" >&2
  exit 1
fi

SAMPLES="${1:-}"
CONFIG_JSON="${2:-}"

if [[ -z "$SAMPLES" || -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 <samples> <config.json>" >&2
  exit 1
fi
if [[ ! "$SAMPLES" =~ ^[0-9]+$ ]] || [[ "$SAMPLES" -lt 1 ]]; then
  echo "Error: samples must be a positive integer" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: Config file not found: $CONFIG_JSON" >&2
  exit 1
fi

BASE_DIR="$(cd "$(dirname "$CONFIG_JSON")" && pwd)"
BASENAME="$(basename "$CONFIG_JSON")"
NAME_NO_EXT="${BASENAME%.json}"

echo "[run_batch] Running $SAMPLES jobs from $CONFIG_JSON"

for (( jj=1; jj<=SAMPLES; jj++ )); do
  OUT_CFG="$BASE_DIR/${NAME_NO_EXT}_jj${jj}.json"
  # Set .jj; set .paraLen if missing to total samples
  if jq 'has("paraLen") | not' "$CONFIG_JSON" >/dev/null; then
    jq --argjson jj "$jj" --argjson para "$SAMPLES" '.jj = $jj | .paraLen = $para' \
      "$CONFIG_JSON" > "$OUT_CFG"
  else
    jq --argjson jj "$jj" '.jj = $jj' "$CONFIG_JSON" > "$OUT_CFG"
  fi
  echo "[run_batch] jj=$jj -> $OUT_CFG"
  bash "$SCRIPT_DIR/run_flux.sh" "$OUT_CFG"
done

echo "[run_batch] Completed $SAMPLES jobs."

