#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"
python3 "$SCRIPT_DIR/analyze_scembryo.py"
echo "[analyze-scembryo] Done. Outputs under ./examples/analyze_demo/out_scEmbryo/"
