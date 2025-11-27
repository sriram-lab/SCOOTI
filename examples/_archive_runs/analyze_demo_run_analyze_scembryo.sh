#!/usr/bin/env bash
# Archived: use the unified `scooti` CLI instead.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
python3 "$REPO_ROOT/examples/analyze_demo/analyze_scembryo.py"

