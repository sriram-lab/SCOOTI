#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scooti/run_analyzer.sh path/to/analyze_config.json

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

CONFIG_JSON="${1:-}"
if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/analyze_config.json" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: Config file not found: $CONFIG_JSON" >&2
  exit 1
fi

# Prefer 'scooti' conda env's python if present; also preload libittnotify if available
PYTHON="python3"
if command -v conda >/dev/null 2>&1; then
  ENV_PATH=$(conda env list | awk '$1=="scooti"{print $NF}')
  if [[ -n "$ENV_PATH" && -x "$ENV_PATH/bin/python" ]]; then
    PYTHON="$ENV_PATH/bin/python"
  fi
fi

PRELOAD_ARG=()
if [[ -n "${CONDA_PREFIX:-}" ]]; then
  for cand in "$CONDA_PREFIX"/lib/libittnotify.so "$CONDA_PREFIX"/lib/libittnotify.so.*; do
    if [[ -f "$cand" ]]; then
      echo "Using libittnotify at $cand"
      PRELOAD_ARG=("LD_PRELOAD=$cand")
      break
    fi
  done
fi

env -u LD_LIBRARY_PATH "${PRELOAD_ARG[@]}" "$PYTHON" -m scooti.analysis_runner --config "$CONFIG_JSON"
