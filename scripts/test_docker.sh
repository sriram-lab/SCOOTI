#!/usr/bin/env bash
set -euo pipefail

if ! command -v docker >/dev/null; then
  echo "docker not found; please install Docker first." >&2
  exit 1
fi

image_name="${1:-scooti-test:ubuntu24}"

docker build -f docker/Dockerfile -t "${image_name}" .

docker run --rm "${image_name}" /bin/bash -lc "conda run -n scooti scooti -h"
docker run --rm "${image_name}" /bin/bash -lc "NUMBA_CACHE_DIR=/tmp/numba_cache conda run -n scooti python -c \"import scanpy as sc; print(sc.__version__)\""
docker run --rm "${image_name}" /bin/bash -lc "NUMBA_CACHE_DIR=/tmp/numba_cache conda run -n scooti python -c \"import scooti; scooti.load(); print('scooti load ok')\""
