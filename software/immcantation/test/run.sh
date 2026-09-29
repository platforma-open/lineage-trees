#!/usr/bin/env bash
# Run the trees.R and alleles.R behaviour checks inside the block's immcantation
# image, both suites at once.
#
#   software/immcantation/test/run.sh [image]
#
# Without an argument the image is the one the last build recorded in
# dist/artifacts/main/docker_x64.json: used if present locally, pulled if not,
# built from the Dockerfile as a last resort. The scripts under test are mounted
# over the image's copies, so the working tree is what gets checked.
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
  else
    image="lineage-trees-immcantation:test"
    echo "no built image available, building $image"
    docker build -q -t "$image" -f "$pkg/Dockerfile" "$pkg/context" >/dev/null
  fi
fi
echo "image: $image"

# Under the package rather than /tmp: on a docker-in-docker runner only the
# workspace is shared with the daemon. Ignored by the repo's `work/` rule.
work="$here/work"
rm -rf "$work"
mkdir -p "$work/trees" "$work/alleles"
trap 'rm -rf "$work"' EXIT

# Run as the invoking user so the fixtures stay removable; R then wants a
# writable HOME.
run() {
  local entrypoint="$1"; shift
  docker run --rm \
    --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$here:/test:ro" \
    -v "$pkg/context/trees.R:/app/trees.R:ro" \
    -v "$pkg/context/alleles.R:/app/alleles.R:ro" \
    -v "$pkg/context/common.R:/app/common.R:ro" \
    -v "$work:/work" \
    --entrypoint "$entrypoint" "$image" "$@" 2>&1 | grep -v "replacing previous import"
}

start=$SECONDS
run bash -c "Rscript /test/fixture.R /work/trees && Rscript /test/test_trees.R /work/trees" \
  > "$work/trees.log" &
trees=$!
run Rscript /test/test_alleles.R /work/alleles > "$work/alleles.log" &
alleles=$!

status=0
wait "$trees" || status=1
wait "$alleles" || status=1
echo "== trees.R =="; cat "$work/trees.log"
echo "== alleles.R =="; cat "$work/alleles.log"
echo "immcantation checks finished in $((SECONDS - start)) s"
exit "$status"
