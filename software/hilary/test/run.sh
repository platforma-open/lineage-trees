#!/usr/bin/env bash
# Run the cluster.py behaviour checks inside the block's HILARy image.
#
#   software/hilary/test/run.sh [image]
#
# Without an argument the image is the one the last build recorded in
# dist/artifacts/main/docker_x64.json: used if present locally, pulled if not,
# built from the generated Dockerfile as a last resort. The working tree's
# cluster.py is mounted over the image's copy.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pkg="$(dirname "$here")"

image="${1:-}"
if [ -z "$image" ]; then
  info="$pkg/dist/artifacts/main/docker_x64.json"
  if [ -f "$info" ]; then
    image="$(sed -n 's/.*"remoteArtifactLocation":"\([^"]*\)".*/\1/p' "$info")"
  fi
  if [ -n "$image" ] && docker image inspect "$image" >/dev/null 2>&1; then
    :
  elif [ -n "$image" ] && docker pull -q "$image" >/dev/null 2>&1; then
    echo "pulled $image"
  elif [ -f "$pkg/dist/docker/Dockerfile-main" ]; then
    image="lineage-trees-hilary:test"
    echo "no built image available, building $image"
    docker build -q -t "$image" -f "$pkg/dist/docker/Dockerfile-main" "$pkg/src" >/dev/null
  else
    echo "no HILARy image: build the package first" >&2
    exit 1
  fi
fi
echo "image: $image"

start=$SECONDS
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$pkg/src:/sw/src:ro" -v "$here:/sw/test:ro" \
  "$image" python /sw/test/test_cluster.py
echo "HILARy checks finished in $((SECONDS - start)) s"
