#!/usr/bin/env bash
# Finish the bottle-batch comparison figures end-to-end.
# Run from the repo root in a terminal that can talk to Docker Desktop:
#
#   bash scripts/finish_comparison_figures.sh
#
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=src

echo "==> Checking Docker"
docker info >/dev/null
if ! docker image inspect pflotran-py-test >/dev/null 2>&1; then
  echo "==> Building pflotran-py-test image (first time can take a while)"
  docker build -t pflotran-py-test -f Containerfile .
fi

# Monod-only salt story: the fitted Cl- smoothstep (threshold 0.75) is not
# applied. The inhibition diagnostic showed it was responsible for the ~850x
# cliff onto a methane floor; the network's own 0.2 M Cl- Monod remains.
if [[ ! -d decks ]] || [[ -z "$(ls decks/*.in 2>/dev/null || true)" ]]; then
  echo "==> Generating decks (cellulose hydrolysis; no fitted Cl- smoothstep)"
  python -m pflotran_py.comparison.decks --output-dir decks \
    --cellulose-hydrolysis
else
  echo "==> Decks already present ($(ls decks/*.in | wc -l | tr -d ' ') files)"
fi

echo "==> Running 15 closed-batch decks in Docker"
python -m pflotran_py.comparison.run_decks --run-root runs --clean

echo "==> Building absolute + per-starting-C figures"
python -m pflotran_py.comparison.figures

echo
echo "Done. Figures in output/comparison/:"
ls -1 output/comparison/*.png
