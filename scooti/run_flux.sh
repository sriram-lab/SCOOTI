#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scooti/run_flux.sh path/to/config.json

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_JSON="${1:-}"
if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/config.json" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: Config file not found: $CONFIG_JSON" >&2
  exit 1
fi

# Resolve absolute path to config
if command -v realpath >/dev/null 2>&1; then
  CONFIG_ABS="$(realpath "$CONFIG_JSON")"
else
  CONFIG_ABS="$(python3 - "$CONFIG_JSON" <<'PY'
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
)"
fi

# Default COBRA toolbox path (can be overridden in JSON config)
COBRA_DEFAULT="${COBRA_DEFAULT:-$HOME/cobratoolbox}"

# Run MATLAB with repo-anchored path to metabolicModel and default COBRA_path fallback
matlab -nodisplay -nosplash -r "addpath(genpath(fullfile('$REPO_ROOT','scooti','metabolicModel'))); \
  config = jsondecode(fileread('$CONFIG_ABS')); \
  if ~isfield(config,'COBRA_path') || isempty(config.COBRA_path), config.COBRA_path = '$COBRA_DEFAULT'; end; \
  if contains(config.COBRA_path, '~'), config.COBRA_path = strrep(config.COBRA_path, '~', getenv('HOME')); end; \
  validate_config(config); \
  multiObj_CBM(config); \
  exit;"

# Post-process: gzip flux CSVs to .csv.gz for inference demo
# Read save_root_path from config and resolve to absolute path
if command -v jq >/dev/null 2>&1; then
  OUT_DIR=$(jq -r '.save_root_path // .save_root_path // empty' "$CONFIG_ABS")
else
  OUT_DIR=$(python3 - "$CONFIG_ABS" <<'PY'
import json, os, sys
cfg=json.load(open(sys.argv[1]))
print(cfg.get('save_root_path',''))
PY
  )
fi

if [[ -n "$OUT_DIR" ]]; then
  if [[ "$OUT_DIR" != /* ]]; then
    OUT_DIR_ABS="$REPO_ROOT/${OUT_DIR#./}"
  else
    OUT_DIR_ABS="$OUT_DIR"
  fi
  if [[ -d "$OUT_DIR_ABS" ]]; then
    echo "[run_flux] Compressing flux CSVs under: $OUT_DIR_ABS"
    # Compress any *_fluxes.csv within two levels
    find "$OUT_DIR_ABS" -maxdepth 2 -type f -name '*_fluxes.csv' -print0 | xargs -0 -r gzip -f || true
  fi
fi
