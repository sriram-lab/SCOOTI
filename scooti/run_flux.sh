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

run_matlab_with_config() {
  local cfg_abs="$1"
  matlab -nodisplay -nosplash -r "addpath(genpath(fullfile('$REPO_ROOT','scooti','metabolicModel'))); \
    config = jsondecode(fileread('$cfg_abs')); \
    if ~isfield(config,'COBRA_path') || isempty(config.COBRA_path), config.COBRA_path = '$COBRA_DEFAULT'; end; \
    if contains(config.COBRA_path, '~'), config.COBRA_path = strrep(config.COBRA_path, '~', getenv('HOME')); end; \
    validate_config(config); \
    multiObj_CBM(config); \
    exit;"
}

# Auto-batch: if jj is absent and parameter lists/logspace scanning are detected, run jj=1..N
if command -v jq >/dev/null 2>&1; then
  HAS_JJ=$(jq 'has("jj")' "$CONFIG_ABS")
  PARA_LEN=$(jq -r '.paraLen // 1' "$CONFIG_ABS")
  # Derive list lengths from kappaArr/rhoArr or CFR_kappa/CFR_rho strings with commas
  LEN_K=$(jq -r '((.kappaArr // .CFR_kappa // "" ) | tostring) as $s | (if ( ($s|type)=="string" and ($s|contains(",")) ) then ($s|split(",")|length) else 1 end)' "$CONFIG_ABS")
  LEN_R=$(jq -r '((.rhoArr   // .CFR_rho   // "" ) | tostring) as $s | (if ( ($s|type)=="string" and ($s|contains(",")) ) then ($s|split(",")|length) else 1 end)' "$CONFIG_ABS")
  TOTAL=$(( LEN_K * LEN_R ))
  CFR_NEG1=$(jq -r '(.CFR_kappa // 0) == -1' "$CONFIG_ABS")

  if [[ "$HAS_JJ" != true ]]; then
    if [[ $TOTAL -gt 1 ]]; then
      SAMPLES=$TOTAL
    elif [[ "$CFR_NEG1" == true && $PARA_LEN -gt 1 ]]; then
      SAMPLES=$PARA_LEN
    else
      SAMPLES=1
    fi
  else
    SAMPLES=1
  fi
else
  SAMPLES=1
fi

if [[ ${SAMPLES:-1} -gt 1 ]]; then
  echo "[run_flux] Auto-batch detected: running jj=1..$SAMPLES"
  BASE_DIR="$(dirname "$CONFIG_ABS")"
  BASENAME="$(basename "$CONFIG_ABS")"
  NAME_NO_EXT="${BASENAME%.json}"
  for (( jj=1; jj<=SAMPLES; jj++ )); do
    OUT_CFG="$BASE_DIR/${NAME_NO_EXT}_jj${jj}.json"
    if command -v jq >/dev/null 2>&1; then
      jq --argjson jj "$jj" --argjson para "$SAMPLES" '.jj = $jj | (.paraLen = (.paraLen // $para))' "$CONFIG_ABS" > "$OUT_CFG"
    else
      cp "$CONFIG_ABS" "$OUT_CFG"
    fi
    # Resolve absolute path to generated config
    if command -v realpath >/dev/null 2>&1; then
      OUT_ABS="$(realpath "$OUT_CFG")"
    else
      OUT_ABS="$(python3 - "$OUT_CFG" <<'PY'
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
)"
    fi
    echo "[run_flux] jj=$jj -> $OUT_ABS"
    run_matlab_with_config "$OUT_ABS"
  done
else
  run_matlab_with_config "$CONFIG_ABS"
fi

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
