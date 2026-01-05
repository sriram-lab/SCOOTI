#!/usr/bin/env bash
set -euo pipefail
# Default config now targets single-cell embryogenesis transitions.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CFG="${1:-$SCRIPT_DIR/identify_siggenes_config.json}"

if [[ ! -f "$CFG" ]]; then
  echo "Config not found: $CFG" >&2
  exit 1
fi

# Activate environment if needed (optional)
# source "$REPO_ROOT/.venv/bin/activate" || true

echo "[identifySigGenes] Using config: $CFG"
python3 "$SCRIPT_DIR/demo_identify_siggenes.py" "$CFG"
