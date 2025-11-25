#!/usr/bin/env bash
set -euo pipefail
# All-in-one: run inference then analysis under reproduce_demo
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INF_CFG_DEFAULT="$SCRIPT_DIR/inference_reproduce/reproduce_inference_config.json"
AN_CFG_DEFAULT="$SCRIPT_DIR/analyze_reproduce/analyze_reproduce_config.json"
ENGINE="legacy" # or minimal

# Parse simple flags: --engine=legacy|minimal; --inf=<cfg>; --an=<cfg>
for arg in "$@"; do
  case "$arg" in
    --engine=*) ENGINE="${arg#*=}" ;;
    --inf=*) INF_CFG_DEFAULT="${arg#*=}" ;;
    --an=*) AN_CFG_DEFAULT="${arg#*=}" ;;
  esac
done

INF_RUN="$SCRIPT_DIR/inference_reproduce/run_reproduce_inference.sh"
AN_RUN_LEG="$SCRIPT_DIR/analyze_reproduce/run_analyze_reproduce.sh"
AN_RUN_MIN="$SCRIPT_DIR/analyze_reproduce/run_analyze_reproduce_minimal.sh"

if [[ ! -f "$INF_CFG_DEFAULT" ]]; then
  echo "Inference config not found: $INF_CFG_DEFAULT" >&2; exit 1
fi

echo "[reproduce-all] 1/2 Inference ..."
bash "$INF_RUN" "$INF_CFG_DEFAULT"

if [[ "$ENGINE" == "minimal" ]]; then
  if [[ -f "$AN_RUN_MIN" ]]; then
    echo "[reproduce-all] 2/2 Analysis (minimal) ..."
    bash "$AN_RUN_MIN" "$AN_CFG_DEFAULT"
  else
    echo "[reproduce-all][WARN] minimal runner missing; falling back to legacy."
    bash "$AN_RUN_LEG" "$AN_CFG_DEFAULT"
  fi
else
  echo "[reproduce-all] 2/2 Analysis (legacy) ..."
  bash "$AN_RUN_LEG" "$AN_CFG_DEFAULT"
fi

echo "[reproduce-all] Done."
