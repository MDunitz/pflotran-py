#!/usr/bin/env bash
# Build PFLOTRAN with pflotran-py custom AWINHIBIT reaction sandboxes.
#
# Requires PETSC_DIR / PETSC_ARCH (set automatically in the base container).
# Writes the binary to $OUTPUT_DIR/pflotran by default.
#
# Usage:
#   ./scripts/build_pflotran_custom.sh
#   PFLOTRAN_SRC=/path/to/pflotran/src/pflotran OUTPUT_DIR=/work/build ./scripts/build_pflotran_custom.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PFLOTRAN_SRC="${PFLOTRAN_SRC:-/scratch/pflotran/src/pflotran}"
SANDBOX_DIR="${SANDBOX_DIR:-$REPO_ROOT/sandbox}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/build}"
PETSC_DIR="${PETSC_DIR:-/scratch/petsc}"
PETSC_ARCH="${PETSC_ARCH:-petsc-arch}"

echo "Patching PFLOTRAN source at $PFLOTRAN_SRC ..."
python3 "$REPO_ROOT/scripts/patch_pflotran_sandboxes.py" \
    --pflotran-src "$PFLOTRAN_SRC" \
    --sandbox-dir "$SANDBOX_DIR"

echo "Building PFLOTRAN (PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH) ..."
export PETSC_DIR PETSC_ARCH
cd "$PFLOTRAN_SRC"
make clean >/dev/null 2>&1 || true
make -j"$(nproc 2>/dev/null || echo 2)" pflotran

mkdir -p "$OUTPUT_DIR"
cp "$PFLOTRAN_SRC/pflotran" "$OUTPUT_DIR/pflotran"
echo "Custom PFLOTRAN binary: $OUTPUT_DIR/pflotran"
