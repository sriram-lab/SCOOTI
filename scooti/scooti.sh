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
#   scooti -T ./examples/tradeoff_demo/tradeoff_config.legacy.json
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
  scooti -T ./examples/tradeoff_demo/tradeoff_config.legacy.json
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
