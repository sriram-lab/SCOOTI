#!/usr/bin/env bash
set -euo pipefail

# Demo: generate omics-constrained flux models via SCOOTI (CFR/DFA framework)
# Uses MATLAB via scooti/run_flux.sh with the provided config JSON.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_JSON="${1:-$SCRIPT_DIR/constrained_demo_config.json}"

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[constrained-demo] Running flux generation with: $CONFIG_JSON"
bash "$REPO_ROOT/scooti/run_flux.sh" "$CONFIG_JSON"
OUT_DIR="$REPO_ROOT/examples/constrained_demo/out/constrained_models"
echo "[constrained-demo] Done. Outputs under examples/constrained_demo/out/constrained_models/"
if [[ -d "$OUT_DIR" ]]; then
  echo "[constrained-demo] Showing a few generated files:"
  find "$OUT_DIR" -maxdepth 2 -type f -name '*_fluxes.csv*' | sort | head -n 10
fi

