#!/usr/bin/env bash
# Decompose the modelled salt suppression into the terms that cause it.
# Run from the repo root in a terminal that can talk to Docker Desktop:
#
#   bash scripts/run_inhibition_diagnostic.sh
#
# Five deck variants x 15 batches = 75 runs, sequential. Budget an hour or two.
# Nothing is fitted: each variant is the existing deck with one term switched
# off. Results land in output/comparison/diagnostic/.
#
# To re-read finished runs without re-running them:
#
#   PYTHONPATH=src python -m pflotran_py.comparison.inhibition_diagnostic \
#     --stage summarise
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

python -m pflotran_py.comparison.inhibition_diagnostic --stage decks >/dev/null
echo "==> Decks: $(find decks/diagnostic -name '*.in' | wc -l | tr -d ' ') across 5 variants"

python -m pflotran_py.comparison.inhibition_diagnostic --stage runs
python -m pflotran_py.comparison.inhibition_diagnostic --stage summarise

echo
echo "Done. Tables and figure in output/comparison/diagnostic/:"
ls -1 output/comparison/diagnostic/
