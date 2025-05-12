#!/bin/bash

# Example: bash run_trainer.sh ./example/demo_trainer_config.json

# ======== Set up ==========

CONFIG_JSON="$1"

if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/config.json"
  exit 1
fi

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: Config file not found: $CONFIG_JSON"
  exit 1
fi

# ======== Activate environment (optional) ==========
# If you have a virtualenv or conda environment, activate it here
# source ~/your_env/bin/activate
# conda activate your_env

# ======== Parse JSON and Run ==========

# Parse JSON fields
# You need to install jq if you don't have it: sudo apt install jq
unconModel=$(jq -r '.unconModel' "$CONFIG_JSON")
conModel=$(jq -r '.conModel' "$CONFIG_JSON")
savePath=$(jq -r '.savePath' "$CONFIG_JSON")
kappaArr=$(jq -r '.kappaArr' "$CONFIG_JSON")
rhoArr=$(jq -r '.rhoArr' "$CONFIG_JSON")
dkappaArr=$(jq -r '.dkappaArr' "$CONFIG_JSON")
expName=$(jq -r '.expName' "$CONFIG_JSON")
unconNorm=$(jq -r '.unconNorm' "$CONFIG_JSON")
conNorm=$(jq -r '.conNorm' "$CONFIG_JSON")
medium=$(jq -r '.medium' "$CONFIG_JSON")
method=$(jq -r '.method' "$CONFIG_JSON")
model=$(jq -r '.model' "$CONFIG_JSON")
inputType=$(jq -r '.inputType' "$CONFIG_JSON")
clusterPath=$(jq -r '.clusterPath' "$CONFIG_JSON")
objListPath=$(jq -r '.objListPath' "$CONFIG_JSON")
rank=$(jq -r '.rank' "$CONFIG_JSON")
stackModel=$(jq -r '.stackModel' "$CONFIG_JSON")
sampling=$(jq -r '.sampling' "$CONFIG_JSON")
learner=$(jq -r '.learner' "$CONFIG_JSON")
geneKO=$(jq -r '.geneKO' "$CONFIG_JSON")
geneListPath=$(jq -r '.geneListPath' "$CONFIG_JSON")
learningRate=$(jq -r '.learningRate' "$CONFIG_JSON")
epo=$(jq -r '.epo' "$CONFIG_JSON")

# ======== Run Python Training ==========
python3 SCOOTI_trainer.py \
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
  --epo "$epo"


