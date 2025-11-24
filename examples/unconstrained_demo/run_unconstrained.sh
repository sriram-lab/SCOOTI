#!/usr/bin/env bash
set -euo pipefail

# Demo: generate unconstrained flux models optimizing ideal objectives
# Uses MATLAB via scooti/run_flux.sh with the provided config JSON.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_JSON="${1:-$SCRIPT_DIR/unconstrained_demo_config.json}"

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[unconstrained-demo] Running flux generation with: $CONFIG_JSON"
bash "$REPO_ROOT/scooti/run_flux.sh" "$CONFIG_JSON"
OUT_DIR="$REPO_ROOT/examples/unconstrained_demo/out/unconstrained_models"
echo "[unconstrained-demo] Done. Outputs under examples/unconstrained_demo/out/unconstrained_models/"
if [[ -d "$OUT_DIR" ]]; then
  echo "[unconstrained-demo] Showing a few generated files:"
  find "$OUT_DIR" -maxdepth 2 -type f -name '*_fluxes.csv*' | sort | head -n 10
fi
