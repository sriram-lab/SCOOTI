#!/usr/bin/env bash
set -euo pipefail

# scooti.sh — Mother wrapper for SCOOTI workflows
#
# Quick usage examples:
#   scooti quickstart
#   scooti -U ./examples/unconstrained_demo/unconstrained_demo_config.json
#   scooti -C ./examples/constrained_demo/constrained_demo_config.json
#   scooti -I ./examples/inference_demo/infer_config.json
#   scooti -A ./examples/analyze_demo/analyze_config.json
#   scooti -T ./examples/tradeoff_demo/tradeoff_config.json
#   scooti -U uncon.json -C con.json -I infer.json -A analyze.json
#
# Flags (order determines execution):
#   -U, --uncon   PATH   Unconstrained/constrained modeling (MATLAB) via run_flux.sh
#   -C, --con     PATH   Constrained modeling (same engine; config controls constraint)
#   -I, --infer   PATH   Training/inference via run_trainer.sh
#   -A, --analyze PATH   Post-analysis via run_analyze.sh
#   -T, --tradeoff PATH  Tradeoff analysis via run_tradeoff.sh
#   -h, --help           Show this help message and exit
#
# Special:
#   quickstart           Run the quickstart reproduction script
#   quickstart-scembryo  Run scEmbryo quickstart (inference + analysis) using archived models

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

print_banner() {
  cat <<'BANNER'

 ░▒▓███████▓▒░░▒▓██████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░▒▓████████▓▒░▒▓█▓▒░ 
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░ 
░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░ 
 ░▒▓██████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░ 
       ░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░ 
       ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░ 
░▒▓███████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░ 
     __o      
   _ \<_      SCOOTI: Single-Cell Optimization Objective and Tradeoff Inference
  (_)/(_)     

BANNER
}

usage() {
  cat <<EOF
Usage:
  scooti quickstart
  scooti quickstart-scembryo
  scooti [OPTIONS]

Options (each expects a JSON config; executed in given order):
  -U, --uncon   PATH   Un/constrained modeling via MATLAB engine (run_flux.sh)
  -C, --con     PATH   Constrained modeling (same engine; config sets constraint)
  -I, --infer   PATH   Training/inference (run_trainer.sh)
  -A, --analyze PATH   Post-analysis (run_analyze.sh)
  -T, --tradeoff PATH  Tradeoff analysis (run_tradeoff.sh)
  -h, --help           Show this help and exit

Examples:
  scooti -U ./examples/unconstrained_demo/unconstrained_demo_config.json
  scooti -U uncon.json -C con.json -I infer.json -A analyze.json
  scooti -T ./examples/tradeoff_demo/tradeoff_config.json
  scooti quickstart-scembryo
EOF
}

declare -a QUEUE

if [[ ${1:-} == "quickstart" ]]; then
  print_banner
  LOG_DIR="${SCOOTI_LOG_DIR:-$REPO_ROOT/logs}"
  mkdir -p "$LOG_DIR"
  LOG_FILE="$LOG_DIR/scooti-$(date +%Y%m%d-%H%M%S)-$$.log"
  echo "[scooti] Logging to: $LOG_FILE"
  exec > >(tee -a "$LOG_FILE") 2>&1
  QS_SCRIPT="$REPO_ROOT/examples/quickstart/run_reproduce_all.sh"
  if [[ ! -x "$QS_SCRIPT" ]]; then
    echo "Quickstart script not found or not executable: $QS_SCRIPT" >&2
    exit 1
  fi
  bash "$QS_SCRIPT"
  exit 0
fi

if [[ ${1:-} == "quickstart-scembryo" ]]; then
  shift
  print_banner
  LOG_DIR="${SCOOTI_LOG_DIR:-$REPO_ROOT/logs}"
  mkdir -p "$LOG_DIR"
  LOG_FILE="$LOG_DIR/scooti-$(date +%Y%m%d-%H%M%S)-$$.log"
  echo "[scooti] Logging to: $LOG_FILE"
  exec > >(tee -a "$LOG_FILE") 2>&1

  QSDIR="$REPO_ROOT/examples/quickstart/_archive_scEmbryo"
  UDIR="$QSDIR/unconstrained_models"
  CDIR="$QSDIR/constrained_models"
  OUTI="$QSDIR/out/regression_models"
  OUTA="$QSDIR/out/analyze"
  mkdir -p "$OUTI" "$OUTA"

  INF_CFG="$QSDIR/quickstart_scembryo_infer.json"
  AN_CFG="$QSDIR/quickstart_scembryo_analyze.json"

  echo "[scooti] Preparing scEmbryo quickstart configs under $QSDIR"
  cat > "$INF_CFG" <<EOF
{
  "unconModel": "$UDIR/",
  "conModel": "$CDIR/",
  "savePath": "$OUTI/",
  "kappaArr": "0.1",
  "rhoArr": "0.01",
  "dkappaArr": "0.1",
  "expName": "scembryo_quickstart",
  "unconNorm": "T",
  "conNorm": "F",
  "medium": "KSOM",
  "method": "cfr",
  "model": "recon1",
  "inputType": "flux",
  "clusterPath": "",
  "objListPath": "",
  "rank": "F",
  "stackModel": "F",
  "sampling": "F",
  "learner": "L",
  "geneKO": "F",
  "geneListPath": "",
  "learningRate": 0.001,
  "epo": 5000,
  "fileSuffix": "_fluxes.csv.gz"
}
EOF

  cat > "$AN_CFG" <<EOF
{
  "flux_paths": { "exp": "$CDIR/" },
  "coef_paths": { "exp": "$OUTI/" },
  "save_root_path": "$OUTA/",
  "engine": "legacy",
  "reduction": "auto",
  "GEM_path": "./scooti/metabolicModel/GEMs/Shen2019.mat",
  "uncon_model_path": "$UDIR/",
  "col_map": {},
  "samplingFlux_path": "",
  "sel_para": "k0.1_r0.01",
  "prefix": "scembryo_quickstart",
  "medium": "KSOM",
  "labels": {
    "mode": "regex",
    "regex": "(1C_to_2cell|2C_to_32cell)",
    "regex_group": 1,
    "group_map": {
      "1C_to_2cell": "1C2C",
      "2C_to_32cell": "2CBC"
    }
  },
  "get_coef": { "metType_cluster": false },
  "coef_analysis": {
    "unknown_clustering": false,
    "clustering": true,
    "entropy": true,
    "distance": true,
    "compare": true,
    "umap_para": [5, 50],
    "method": "average",
    "ref_col": null
  }
}
EOF

  echo "[scooti] 1/2 Inference (scEmbryo)"
  bash "$SCRIPT_DIR/run_trainer.sh" "$INF_CFG"
  echo "[scooti] 2/2 Analysis (scEmbryo)"
  bash "$SCRIPT_DIR/run_analyze.sh" "$AN_CFG"
  echo "[scooti] scEmbryo quickstart completed. Outputs under $QSDIR/out/"
  exit 0
fi

if [[ $# -eq 0 ]]; then
  usage
  exit 0
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
  -U | --uncon)
    [[ $# -lt 2 ]] && {
      echo "Missing path for $1"
      exit 1
    }
    QUEUE+=("$SCRIPT_DIR/run_flux.sh|$2")
    shift 2
    ;;
  -C | --con)
    [[ $# -lt 2 ]] && {
      echo "Missing path for $1"
      exit 1
    }
    QUEUE+=("$SCRIPT_DIR/run_flux.sh|$2")
    shift 2
    ;;
  -I | --infer)
    [[ $# -lt 2 ]] && {
      echo "Missing path for $1"
      exit 1
    }
    QUEUE+=("$SCRIPT_DIR/run_trainer.sh|$2")
    shift 2
    ;;
  -A | --analyze)
    [[ $# -lt 2 ]] && {
      echo "Missing path for $1"
      exit 1
    }
    QUEUE+=("$SCRIPT_DIR/run_analyze.sh|$2")
    shift 2
    ;;
  -T | --tradeoff)
    [[ $# -lt 2 ]] && {
      echo "Missing path for $1"
      exit 1
    }
    QUEUE+=("$SCRIPT_DIR/run_tradeoff.sh|$2")
    shift 2
    ;;
  -h | --help)
    usage
    exit 0
    ;;
  *)
    echo "Unknown arg: $1" >&2
    usage
    exit 1
    ;;
  esac
done

print_banner

# Initialize logging
LOG_DIR="${SCOOTI_LOG_DIR:-$REPO_ROOT/logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/scooti-$(date +%Y%m%d-%H%M%S)-$$.log"
echo "[scooti] Logging to: $LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

run_with_config() {
  local script="$1"
  shift
  local cfg="$1"
  shift || true
  if [[ -z "$cfg" ]]; then
    echo "Error: missing --config for $script" >&2
    exit 1
  fi
  if [[ ! -f "$cfg" ]]; then
    echo "Error: config not found: $cfg" >&2
    exit 1
  fi
  bash "$script" "$cfg"
}

step=0
for item in "${QUEUE[@]}"; do
  step=$((step + 1))
  IFS='|' read -r script cfg <<<"$item"
  echo "[scooti] Step $step/${#QUEUE[@]}: $(basename "$script") with $cfg"
  run_with_config "$script" "$cfg"
done

echo "[scooti] All requested stages completed."
