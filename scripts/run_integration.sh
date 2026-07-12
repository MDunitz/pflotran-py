#!/usr/bin/env bash
# Build the PFLOTRAN test image and run integration tests.
#
# Usage:
#   ./scripts/run_integration.sh              # auto-detect docker or podman
#   ./scripts/run_integration.sh docker
#   ./scripts/run_integration.sh podman
#   RUNTIME=podman ./scripts/run_integration.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

pick_runtime() {
  if [[ -n "${1:-}" ]]; then
    echo "$1"
    return
  fi
  if [[ -n "${RUNTIME:-}" ]]; then
    echo "$RUNTIME"
    return
  fi
  if command -v docker >/dev/null 2>&1; then
    echo docker
    return
  fi
  if command -v podman >/dev/null 2>&1; then
    echo podman
    return
  fi
  echo "ERROR: neither docker nor podman found in PATH" >&2
  exit 1
}

RUNTIME="$(pick_runtime "${1:-}")"
IMAGE="${IMAGE:-pflotran-py-test}"

if [[ "$RUNTIME" == "podman" ]]; then
  # :Z is required on SELinux hosts (Fedora/RHEL); harmless elsewhere.
  VOLUME_SPEC="${ROOT}:/work:Z"
else
  VOLUME_SPEC="${ROOT}:/work"
fi

echo "==> Building with ${RUNTIME}: ${IMAGE}"
"$RUNTIME" build -t "$IMAGE" -f Containerfile .

echo "==> Running integration tests with ${RUNTIME}"
"$RUNTIME" run --rm \
  -v "$VOLUME_SPEC" \
  -w /work \
  "$IMAGE" \
  pytest tests/ -v --tb=short -m integration
