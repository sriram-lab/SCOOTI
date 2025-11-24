#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

echo "[quickstart] Generating demo flux predictions (unconstrained)"
bash "$ROOT_DIR/scooti/run_flux.sh" "$ROOT_DIR/examples/quickstart/configs/unconstrained.json"

echo "[quickstart] Generating demo flux predictions (constrained)"
bash "$ROOT_DIR/scooti/run_flux.sh" "$ROOT_DIR/examples/quickstart/configs/constrained.json"

echo "[quickstart] Running objective inference"
bash "$ROOT_DIR/scooti/run_trainer.sh" "$ROOT_DIR/examples/quickstart/configs/inference.json"

echo "[quickstart] Done. Outputs in examples/quickstart/out/"

