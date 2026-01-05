#!/usr/bin/env bash
set -euo pipefail

# Reproduce: run meta-regression with user-provided unconstrained models
# - unconstrained: ./examples/quickstart/unconstrained_models/
# - constrained:   ./examples/constrained_demo/out/constrained_models/
# Uses scooti/run_trainer.sh with the provided config.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../" && pwd)"
CONFIG_JSON="${1:-$SCRIPT_DIR/reproduce_inference_config.json}"

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[inference-reproduce] Running meta-regression with: $CONFIG_JSON"

UM=$(jq -r '.unconModel' "$CONFIG_JSON")
CM=$(jq -r '.conModel' "$CONFIG_JSON")
MED=$(jq -r '.medium' "$CONFIG_JSON")
SFX=$(jq -r '.fileSuffix // "_fluxes.csv.gz"' "$CONFIG_JSON")

UM_ABS="$UM"
CM_ABS="$CM"
if [[ -n "$UM_ABS" && "$UM_ABS" != /* ]]; then UM_ABS="$REPO_ROOT/${UM_ABS#./}"; fi
if [[ -n "$CM_ABS" && "$CM_ABS" != /* ]]; then CM_ABS="$REPO_ROOT/${CM_ABS#./}"; fi

echo "[inference-reproduce] unconstrained dir: $UM_ABS"
echo "[inference-reproduce] constrained   dir: $CM_ABS"
echo "[inference-reproduce] medium='$MED' fileSuffix='$SFX'"

# Preflight: check unconstrained path
if [[ ! -d "$UM_ABS" ]]; then
  echo "[inference-reproduce][WARN] Unconstrained path not found: $UM_ABS"
  echo "[inference-reproduce][HINT] Create it and place *_metadata.json + *${SFX} files."
else
  echo "[inference-reproduce] Scanning unconstrained metadata in: $UM_ABS"
  if command -v jq >/dev/null 2>&1; then
    mapfile -t MEDS < <(find "$UM_ABS" -maxdepth 1 -type f -name '*_metadata.json' -print0 | xargs -0 -r jq -r '.medium' 2>/dev/null | sort | uniq -c)
    if [[ ${#MEDS[@]} -gt 0 ]]; then
      printf '[inference-reproduce] Found medium tags in unconstrained: %s\n' "${MEDS[@]}"
      if ! printf '%s\n' "${MEDS[@]}" | grep -q "$MED"; then
        echo "[inference-reproduce][WARN] Config medium='$MED' not present in unconstrained metadata above."
      fi
    else
      echo "[inference-reproduce][WARN] No *_metadata.json found under $UM_ABS"
    fi
  fi
  # Suffix check (unconstrained)
  if ! find "$UM_ABS" -type f -name "*${SFX}" | grep -q .; then
    if [[ "$SFX" == *_fluxes.csv.gz ]]; then
      if find "$UM_ABS" -type f -name "*_fluxes.csv" | grep -q .; then
        echo "[inference-reproduce] No ${SFX} found, but *_fluxes.csv exist. Loader should fall back automatically."
      else
        echo "[inference-reproduce][WARN] No *_fluxes.csv.gz or *_fluxes.csv found in $UM_ABS"
        echo "[inference-reproduce][INFO] Listing a few files in $UM_ABS for diagnosis:"
        find "$UM_ABS" -maxdepth 1 -type f | head -n 10 || true
      fi
    fi
  fi
fi

# Preflight: check constrained path exists
if [[ ! -d "$CM_ABS" ]]; then
  echo "[inference-reproduce][ERROR] Constrained path not found: $CM_ABS" >&2
  exit 1
fi

bash "$REPO_ROOT/scooti/run_trainer.sh" "$CONFIG_JSON"

OUT_DIR=$(jq -r '.savePath' "$CONFIG_JSON")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != /* ]]; then
  OUT_DIR="$REPO_ROOT/${OUT_DIR#./}"
fi
echo "[inference-reproduce] Done. Outputs in: ${OUT_DIR:-<unknown>}"
if [[ -n "${OUT_DIR:-}" && -d "$OUT_DIR" ]]; then
  echo "[inference-reproduce] Showing result files:"
  find "$OUT_DIR" -maxdepth 1 -type f -name '*.csv' | sort | head -n 10
fi
