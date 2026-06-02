#!/usr/bin/env bash
# build_llmap_container.sh — build llmap.sif locally with apptainer (or docker).
#
# Two paths:
#   1. Native apptainer present -> build directly from containers/llmap.def
#      via `apptainer build`.
#   2. No apptainer, but docker available -> build the rootfs with a
#      Dockerfile (which has working DNS for apt), then convert the docker
#      image to a .sif via the quay.io/singularity/singularity:v3.11.4
#      container. This sidesteps singularity 3.11's %post DNS isolation.
#
# If neither is available, fail.
#
# Outputs:
#   /home/<user>/llmap-local/containers/llmap.sif
#
# Usage:
#   bash scripts/build_llmap_container.sh            # plain build
#   bash scripts/build_llmap_container.sh --force    # rebuild even if .sif exists
#   bash scripts/build_llmap_container.sh --fakeroot # apptainer-only: fakeroot

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLMAP_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEF="$LLMAP_ROOT/containers/llmap.def"
DOCKERFILE="$LLMAP_ROOT/containers/Dockerfile"
SIF="$LLMAP_ROOT/containers/llmap.sif"

FORCE=false
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)    FORCE=true; shift ;;
    --fakeroot) EXTRA_ARGS+=(--fakeroot); shift ;;
    -h|--help)
      sed -n '1,22p' "$0"; exit 0 ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

say() { echo "[build_llmap_container] $*" >&2; }

# Preflight ------------------------------------------------------------------
if [[ ! -f "$DEF" ]]; then
  echo "ERROR: definition file missing: $DEF" >&2
  exit 1
fi

USE_DOCKER=false
if ! command -v apptainer >/dev/null 2>&1; then
  if command -v docker >/dev/null 2>&1; then
    say "no native apptainer — falling back to docker + singularity"
    USE_DOCKER=true
  else
    echo "ERROR: neither apptainer nor docker found on this host." >&2
    echo "       Install one, or build the .sif elsewhere and copy it in." >&2
    exit 1
  fi
fi

LLMAP_BIN_LOCAL="$LLMAP_ROOT/build/src/llmap"
if [[ ! -x "$LLMAP_BIN_LOCAL" ]]; then
  echo "ERROR: llmap binary not found: $LLMAP_BIN_LOCAL" >&2
  echo "       Build it first: cmake --build build -j --target llmap" >&2
  exit 1
fi

ONNX_LIB="$LLMAP_ROOT/third_party/onnxruntime/lib/libonnxruntime.so.1"
if [[ ! -f "$ONNX_LIB" ]]; then
  echo "ERROR: onnxruntime not found: $ONNX_LIB" >&2
  exit 1
fi

if [[ -f "$SIF" && "$FORCE" != true ]]; then
  say "$SIF already exists — pass --force to rebuild"
  ls -lh "$SIF"
  exit 0
fi

mkdir -p "$LLMAP_ROOT/containers"

# Build ----------------------------------------------------------------------
if $USE_DOCKER; then
  if [[ ! -f "$DOCKERFILE" ]]; then
    echo "ERROR: Dockerfile not found: $DOCKERFILE" >&2
    exit 1
  fi

  SIF_IMG="quay.io/singularity/singularity:v3.11.4"
  DOCKER_TAG="localhost/llmap:1.0.0"

  say "stage 1: docker build $DOCKER_TAG (apt has working DNS here)"
  docker build \
    --network host \
    --platform linux/amd64 \
    -t "$DOCKER_TAG" \
    -f "$DOCKERFILE" \
    "$LLMAP_ROOT"

  say "stage 2: docker pull $SIF_IMG (if not cached)"
  docker pull "$SIF_IMG" >/dev/null 2>&1 || true

  DOCKER_TAR="$(mktemp /tmp/llmap-docker.XXXXXX.tar)"
  trap "rm -f '$DOCKER_TAR'" EXIT
  say "stage 2: docker save $DOCKER_TAG -> $DOCKER_TAR"
  docker save "$DOCKER_TAG" -o "$DOCKER_TAR"

  say "stage 2: singularity build $SIF from docker-archive (no %post DNS)"
  docker run --rm --privileged \
    --network host \
    -v "$LLMAP_ROOT:/work" \
    -v "$DOCKER_TAR:/tmp/llmap-docker.tar:ro" \
    -w /work \
    "$SIF_IMG" \
    build --force /work/containers/llmap.sif docker-archive:///tmp/llmap-docker.tar

else
  say "apptainer build: $SIF <- $DEF"
  say "apptainer version: $(apptainer --version)"
  if apptainer build "${EXTRA_ARGS[@]}" "$SIF" "$DEF"; then
    :
  elif [[ " ${EXTRA_ARGS[*]} " != *" --fakeroot "* ]]; then
    say "plain build failed — retrying with --fakeroot"
    apptainer build --fakeroot "$SIF" "$DEF"
  else
    echo "ERROR: apptainer build failed" >&2
    exit 1
  fi
fi

say "built: $(ls -lh "$SIF")"

# Smoke test -----------------------------------------------------------------
if $USE_DOCKER; then
  SIF_IMG="quay.io/singularity/singularity:v3.11.4"
  say "self-test (via docker): singularity exec llmap.sif llmap --version"
  docker run --rm --privileged \
    -v "$LLMAP_ROOT:/work:ro" \
    -w /work \
    "$SIF_IMG" \
    exec /work/containers/llmap.sif /opt/llmap/bin/llmap --version
else
  say "self-test: apptainer exec $SIF llmap --version"
  apptainer exec "$SIF" llmap --version
fi

say "OK"
