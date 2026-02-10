#!/usr/bin/env bash
# Sequential PFLOTRAN batch runner (safe memory usage)
##### future users: you CAN make it run in parallel, but some of them are computationally expensive
##### so you're better off not crashing your laptop/venv
#############################################################################
# the script should run all .in files in the specified directory


##### instructions
"""
This assumes you have already set up your PFLOTRAN environment...

1) Copy in the file, e.g:
nano run_pflotran_batch.sh


2) Make it executable:
chmod +x run_pflotran_batch.sh

3) Run it, from the directory where your .in files are, or specify the directory as shown in the second line:
./run_pflotran_batch.sh
  OR if you have a specific directory: ./run_pflotran_batch.sh /path/to/inputs



"""


###########



set -euo pipefail

### User-adjustable defaults
SIM_DIR="${1:-.}"                 # Directory to search for *.in files (arg1 or current dir)
NPROCS="${NPROCS:-1}"             # Override by environment: NPROCS=4 ./run_pflotran_batch.sh
USE_INPUT_PREFIX="${USE_INPUT_PREFIX:-1}"  # 1 = use -input_prefix; 0 = use -pflotranin
PFLOTRAN_BIN="${PFLOTRAN_BIN:-}"  # Optional: full path to pflotran binary
                                   # If empty, will default to $PFLOTRAN_DIR/src/pflotran/pflotran

### make sure pflotran is executable + other stuff
if [[ -z "${PFLOTRAN_BIN}" ]]; then
  if [[ -z "${PFLOTRAN_DIR:-}" ]]; then
    echo "ERROR: PFLOTRAN_BIN not set and PFLOTRAN_DIR not set."
    echo "Set PFLOTRAN_BIN to the pflotran executable, or export PFLOTRAN_DIR."
    exit 1
  fi
  PFLOTRAN_BIN="${PFLOTRAN_DIR}/src/pflotran/pflotran"
fi

if [[ ! -x "${PFLOTRAN_BIN}" ]]; then
  echo "ERROR: PFLOTRAN binary not executable: ${PFLOTRAN_BIN}"
  exit 1
fi

### get the input files
shopt -s nullglob
mapfile -t INFILES < <(find "${SIM_DIR}" -maxdepth 1 -type f -name "*.in" | sort)
shopt -u nullglob

if (( ${#INFILES[@]} == 0 )); then
  echo "No .in files found in: ${SIM_DIR}"
  exit 0
fi

### prepare the logs! for troubleshooting if need be
LOG_DIR="${SIM_DIR%/}/logs"
mkdir -p "${LOG_DIR}"

echo "Found ${#INFILES[@]} .in files in ${SIM_DIR}"
echo "Running sequentially with NPROCS=${NPROCS}"
echo "Logging to ${LOG_DIR}"
echo

FAILED=()

for infile in "${INFILES[@]}"; do
  # get a prefix (basename without .in)
  base="$(basename "${infile}")"
  prefix="${base%.in}"

  # Build the command (sequential, safe memory)
  if [[ "${USE_INPUT_PREFIX}" == "1" ]]; then
    CMD=(mpirun -n "${NPROCS}" "${PFLOTRAN_BIN}" -input_prefix "${prefix}")
    pretty_target="${prefix} (via -input_prefix)"
  else
    # Alternatively, run explicitly with -pflotranin <file.in>
    CMD=(mpirun -n "${NPROCS}" "${PFLOTRAN_BIN}" -pflotranin "${infile}")
    pretty_target="${infile} (via -pflotranin)"
  fi

  LOG_FILE="${LOG_DIR}/${prefix}.log"
  echo "=== Running: ${pretty_target}"
  echo "    Log: ${LOG_FILE}"
  echo "    Command: ${CMD[*]}"
  echo

  # Run and send output to log; capture failures but continue to next
  # this allows us to run all simulations even if some fail
  # and report all failures at the end
  if "${CMD[@]}" 2>&1 | tee "${LOG_FILE}"; then
    echo "--- Completed: ${pretty_target}"
  else
    echo "XXX FAILED: ${pretty_target}"
    FAILED+=("${pretty_target}")
  fi
  echo
done

if (( ${#FAILED[@]} > 0 )); then
  echo "Some runs failed:"
  for f in "${FAILED[@]}"; do
    echo " - ${f}"
  done
  exit 2
else
  echo "All simulations completed successfully."
fi
