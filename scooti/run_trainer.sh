#!/usr/bin/env bash
set -euo pipefail

# Example: bash scooti/run_trainer.sh ./examples/run_inference/demo_inference_config.json

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

CONFIG_JSON="${1:-}"
if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/config.json" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: Config file not found: $CONFIG_JSON" >&2
  exit 1
fi

# Extract fields from JSON (prefer jq; fallback to Python)
if command -v jq >/dev/null 2>&1; then
  JQ() { jq -r "$1" "$CONFIG_JSON"; }
  unconModel=$(JQ '.unconModel')
  conModel=$(JQ '.conModel')
  savePath=$(JQ '.savePath')
  kappaArr=$(JQ '.kappaArr')
  rhoArr=$(JQ '.rhoArr')
  dkappaArr=$(JQ '.dkappaArr')
  expName=$(JQ '.expName')
  unconNorm=$(JQ '.unconNorm')
  conNorm=$(JQ '.conNorm')
  medium=$(JQ '.medium')
  method=$(JQ '.method')
  model=$(JQ '.model')
  inputType=$(JQ '.inputType')
  clusterPath=$(JQ '.clusterPath')
  objListPath=$(JQ '.objListPath')
  rank=$(JQ '.rank')
  stackModel=$(JQ '.stackModel')
  sampling=$(JQ '.sampling')
  learner=$(JQ '.learner')
  geneKO=$(JQ '.geneKO')
  geneListPath=$(JQ '.geneListPath')
  learningRate=$(JQ '.learningRate')
  epo=$(JQ '.epo')
  fileSuffix=$(JQ '.fileSuffix // "_fluxes.csv.gz"')
else
  # Python fallback
  mapfile -t KV < <(python3 - "$CONFIG_JSON" <<'PY'
import json, sys
cfg = json.load(open(sys.argv[1]))
  keys = [
    "unconModel","conModel","savePath","kappaArr","rhoArr","dkappaArr","expName",
    "unconNorm","conNorm","medium","method","model","inputType","clusterPath",
    "objListPath","rank","stackModel","sampling","learner","geneKO","geneListPath",
    "learningRate","epo","fileSuffix"
  ]
for k in keys:
  v = cfg.get(k, "")
  # Keep arrays/objects as compact JSON strings to match prior usage
  import json as _json
  if isinstance(v, (list, dict)):
    v = _json.dumps(v)
  print(f"{k}={v}")
PY
  )
  for line in "${KV[@]}"; do eval "$line"; done
fi

# Ensure output directory exists (savePath)
if [[ -n "$savePath" ]]; then
  mkdir -p "$savePath"
fi

# Prefer conda's libittnotify to satisfy iJIT_NotifyEvent if present (match any soname)
PRELOAD_ARG=()
if [[ -n "${CONDA_PREFIX:-}" ]]; then
  for cand in "$CONDA_PREFIX"/lib/libittnotify.so "$CONDA_PREFIX"/lib/libittnotify.so.*; do
    if [[ -f "$cand" ]]; then
      echo "Using libittnotify at $cand to satisfy iJIT_NotifyEvent"
      PRELOAD_ARG=("LD_PRELOAD=$cand")
      break
    fi
  done
fi

# Run Python trainer with a clean library path to avoid conflicts, but preload libittnotify if available
env -u LD_LIBRARY_PATH "${PRELOAD_ARG[@]}" python3 -m scooti.SCOOTI_trainer \
  --unconModel "$unconModel" \
  --conModel "$conModel" \
  --savePath "$savePath" \
  --kappaArr "$kappaArr" \
  --rhoArr "$rhoArr" \
  --dkappaArr "$dkappaArr" \
  --expName "$expName" \
  --unconNorm "$unconNorm" \
  --conNorm "$conNorm" \
  --medium "$medium" \
  --method "$method" \
  --model "$model" \
  --inputType "$inputType" \
  --clusterPath "$clusterPath" \
  --objListPath "$objListPath" \
  --rank "$rank" \
  --stackModel "$stackModel" \
  --sampling "$sampling" \
  --learner "$learner" \
  --geneKO "$geneKO" \
  --geneListPath "$geneListPath" \
  --learningRate "$learningRate" \
  --epo "$epo" \
  --fileSuffix "$fileSuffix"
