#!/usr/bin/env bash
set -Eeuo pipefail

IMAGE_NAME="mwdrop"
HOST_OUTPUT_DIR="${HOST_OUTPUT_DIR:-$PWD/test_output}"
CONTAINER_OUTPUT_DIR="/data/test_runs"

mkdir -p "$HOST_OUTPUT_DIR"

docker build --no-cache -t "$IMAGE_NAME" ../../

docker run --rm -it \
-v "${HOST_OUTPUT_DIR}:${CONTAINER_OUTPUT_DIR}" \
"$IMAGE_NAME" \
/bin/bash -lc "cd /opt/mwdrop/testing/load-test && ./test.sh -p ${CONTAINER_OUTPUT_DIR} $*"