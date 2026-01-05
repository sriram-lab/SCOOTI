#!/usr/bin/env bash
set -euo pipefail

# Demo: Inference of Metabolic Objectives (meta-regression)
# Uses scooti/run_trainer.sh with the provided demo config.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_JSON="${1:-$SCRIPT_DIR/demo_inference_config.json}"

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Config not found: $CONFIG_JSON" >&2
  exit 1
fi

echo "[inference-demo] Running meta-regression with: $CONFIG_JSON"

# Preflight: check medium + suffix alignment
UM=$(jq -r '.unconModel' "$CONFIG_JSON")
CM=$(jq -r '.conModel' "$CONFIG_JSON")
MED=$(jq -r '.medium' "$CONFIG_JSON")
SFX=$(jq -r '.fileSuffix' "$CONFIG_JSON")

echo "[inference-demo] Config: medium='$MED' fileSuffix='$SFX'"
UM_ABS="$UM"
CM_ABS="$CM"
if [[ -n "$UM_ABS" && "$UM_ABS" != /* ]]; then UM_ABS="$REPO_ROOT/${UM_ABS#./}"; fi
if [[ -n "$CM_ABS" && "$CM_ABS" != /* ]]; then CM_ABS="$REPO_ROOT/${CM_ABS#./}"; fi

if [[ -n "$UM_ABS" && -d "$UM_ABS" ]] && command -v jq >/dev/null 2>&1; then
  echo "[inference-demo] Scanning unconstrained metadata in: $UM_ABS"
  mapfile -t MEDS < <(find "$UM_ABS" -maxdepth 1 -type f -name '*_metadata.json' -print0 | xargs -0 -r jq -r '.medium' 2>/dev/null | sort | uniq -c)
  if [[ ${#MEDS[@]} -gt 0 ]]; then
    printf '[inference-demo] Found medium tags in unconstrained: %s\n' "${MEDS[@]}"
    if ! printf '%s\n' "${MEDS[@]}" | grep -q "$MED"; then
      echo "[inference-demo][WARN] Config medium='$MED' not found in unconstrained metadata above. Update the config or regenerate fluxes with matching medium."
    fi
  else
    echo "[inference-demo][WARN] No *_metadata.json found under $UM_ABS"
  fi
  # Suffix check (unconstrained)
  if ! find "$UM_ABS" -type f -name "*${SFX}" | grep -q .; then
    if [[ "$SFX" == *_fluxes.csv.gz ]]; then
      if find "$UM_ABS" -type f -name "*_fluxes.csv" | grep -q .; then
        echo "[inference-demo] No ${SFX} found, but *_fluxes.csv exist. Loader will fall back automatically."
      else
        echo "[inference-demo][WARN] No *_fluxes.csv.gz or *_fluxes.csv found in $UM_ABS"
        echo "[inference-demo][INFO] Listing a few files in $UM_ABS for diagnosis:"
        find "$UM_ABS" -maxdepth 1 -type f | head -n 10 || true
      fi
    fi
  fi
fi

bash "$REPO_ROOT/scooti/run_trainer.sh" "$CONFIG_JSON"

OUT_DIR=$(jq -r '.savePath' "$CONFIG_JSON")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != /* ]]; then
  OUT_DIR="$REPO_ROOT/${OUT_DIR#./}"
fi
echo "[inference-demo] Done. Outputs in: ${OUT_DIR:-<unknown>}"
if [[ -n "${OUT_DIR:-}" && -d "$OUT_DIR" ]]; then
  echo "[inference-demo] Showing result files:"
  find "$OUT_DIR" -maxdepth 1 -type f -name '*.csv' | sort | head -n 10
fi
