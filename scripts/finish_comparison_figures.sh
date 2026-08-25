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
# Always rebuild: AWINHIBIT Fortran carries network Monod keywords that an
# older image will not recognise.
echo "==> Building pflotran-py-test image with patched AWINHIBIT sandboxes"
docker build -t pflotran-py-test -f Containerfile .

echo "==> Generating decks (cellulose; a_w sandboxes; no Cl- Monod)"
# Meter-read a_w is FIXED_WATER_ACTIVITY by default. Add --use-computed-aw
# only to feed PHREEQC/pitzer.dat a_w into the sandboxes instead.
rm -rf decks
python -m pflotran_py.comparison.decks --output-dir decks \
  --cellulose-hydrolysis --no-cl-inhibition --aw-threshold 0.95

echo "==> Running 15 closed-batch decks in Docker"
python -m pflotran_py.comparison.run_decks --run-root runs --clean

echo "==> Building absolute + per-starting-C figures"
python -m pflotran_py.comparison.figures

echo
echo "Done. Figures in output/comparison/:"
ls -1 output/comparison/*.png
