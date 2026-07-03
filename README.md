# PFLOTRAN Reactive Transport Modeling

Simulation of microbial redox networks (methanogenesis, sulfate reduction, iron cycling) in saline sediments using [PFLOTRAN](https://www.pflotran.org/).

## Directory Structure

```
├── README.md
├── requirements.txt          ← Python dependencies
├── Containerfile             ← Docker/Podman image with custom-sandbox PFLOTRAN
├── scripts/                  ← PFLOTRAN sandbox patch + build scripts
├── generator/                ← Python code to produce PFLOTRAN .in files
├── sandbox/                  ← Custom Fortran 90 reaction modules (water activity inhibition)
├── visualization/            ← Post-processing pipeline (extract → plot → gradient/flux)
├── batch/                    ← Batch file generation + inhibition diagnostics
├── sample_data/              ← Example PFLOTRAN Tecplot output files
├── reference/                ← Historical/reference input decks
├── compare/                  ← Comparison notebooks
└── exploratory/              ← Exploratory analyses (Sanskriti)
    ├── pflotran/             ← PFLOTRAN testing, sandbox builds, visualization
    │   ├── testing/          ← Iterative .in file development (1–10)
    │   ├── sandbox/          ← Fortran reaction sandbox source files
    │   └── visualization/    ← Tecplot post-processing and Plotly 3D viz
    ├── long_term_isq/        ← ISQ gas concentration dashboard (Streamlit)
    ├── CCS_LT.in             ← Long-term CO2 sequestration input deck
    ├── pflotran_vars.py      ← Variable definitions and constants
    ├── constants.py          ← Unit conversion constants
    └── *.ipynb               ← Notebooks (water activity curves, 100-year projections)
```

## Quick Start

### 1. Generate input files at varying salinity

```bash
cd batch/
python create_modified_files.py --max-multiplier 10  # generates 1×–10× seawater .in files
```

### 2. Run the three-variant inhibition diagnostic

```bash
cd batch/
python run_inhibition_test.py --multiplier 5  # A: Cl⁻ only, B: aₓ only, C: both
```

### 3. Run PFLOTRAN simulations

```bash
cd batch/
./run_pflotran_batch.sh /path/to/input/files
```

### 4. Post-process and visualize results

```bash
cd visualization/
python step_orchestra.py  # runs step1 → step2 → step3 → step4
```

---

## Running PFLOTRAN

Most workflows in this repo need a PFLOTRAN binary. Choose the option that matches
your input decks:

Also, note that if you are using Docker, you don't NEED to install PFLOTRAN separately.

| Decks | Custom sandboxes? | Binary |
|-------|-------------------|--------|
| `exploratory/pflotran/testing/1`–`6` | No | Stock PFLOTRAN |
| `exploratory/pflotran/testing/7`–`10`, `batch/9_addnitrogen.in` | Yes (AWINHIBIT) | Custom build below |
| `generator/`-produced decks with `enable_aw_sandbox=True` | Yes | Custom build below |

All decks should point `DATABASE` at `sandbox/hanford.dat` (the committed decks
use developer-specific absolute paths; the integration tests rewrite this
automatically).

### Option A: Docker / Podman (recommended)

Build the test image once. It includes PFLOTRAN v6 with the AWINHIBIT sandboxes
from `sandbox/` pre-compiled at `/opt/pflotran-py/pflotran`:

```bash
# Docker
docker build -t pflotran-py-test -f Containerfile .

# Podman (Linux)
podman build -t pflotran-py-test -f Containerfile .
```

**Stock PFLOTRAN tests** (no custom sandboxes, e.g. deck `3_smaller_grid.in`):

```bash
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  pytest tests/test_docker_e2e.py -v --tb=short
```

**Custom sandbox tests** (AWINHIBIT decks 7–10):

```bash
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  pytest tests/test_custom_docker_e2e.py -v --tb=short
```

**All integration tests:**

```bash
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  pytest tests/ -v --tb=short -m integration
```

**Run the full custom-sandbox pipeline** (prep deck, simulate, extract, plot):

```bash
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  python3 tests/test_custom_docker_e2e.py
```

Output lands in `tests/custom_e2e_output/` (CSV + PNGs).

On Linux with Podman, add `:Z` to the volume mount for SELinux:
`-v "$(pwd)":/work:Z`.

### Option B: Custom sandbox build (manual, inside base image)

If you already have the base PFLOTRAN image but haven't built `Containerfile`:

```bash
docker run --rm -v "$(pwd)":/work -w /work \
  pshuai/jupyter-pflotran-multiplatform:base_v6 \
  bash -c "./scripts/build_pflotran_custom.sh && python3 tests/test_custom_docker_e2e.py"
```

This patches PFLOTRAN with `scripts/patch_pflotran_sandboxes.py` and writes
`build/pflotran`. Set `PFLOTRAN_CUSTOM_EXE` to point at that binary.

### Option C: Stock PFLOTRAN only (no Docker)

For decks without AWINHIBIT sandboxes, use a standard PFLOTRAN install. See
[Installing PFLOTRAN](#installing-pflotran) below for building from source on
Linux/HPC. Point `DATABASE` at `sandbox/hanford.dat` and run:

```bash
mpirun -n 1 /path/to/pflotran -input_prefix my_simulation
```

To patch an existing PFLOTRAN source tree with the custom sandboxes by hand:

```bash
python3 scripts/patch_pflotran_sandboxes.py \
  --pflotran-src /path/to/pflotran/src/pflotran \
  --sandbox-dir sandbox/
cd /path/to/pflotran/src/pflotran && make clean && make pflotran
```

---

## Installing PFLOTRAN
Remember, this is only if you don't want to use the docker stuff.
### Linux

See the [official docs](https://documentation.pflotran.org/user_guide/how_to/installation/linux.html#linux-install) for details.

#### PETSc

```bash
git clone https://gitlab.com/petsc/petsc petsc
cd petsc
git checkout v3.21.5
```

If running Python 3.13, create a 3.11 environment first:

```bash
conda create -n petscpy311 python=3.11 -y
conda activate petscpy311
cd ~/petsc
```

Configure PETSc:

```bash
./configure \
  --with-python-exec=$(which python) \
  --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
  --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3' \
  --with-debugging=0 \
  --download-mpich \
  --download-hdf5 --download-hdf5-fortran-bindings \
  --download-fblaslapack \
  --download-metis --download-parmetis \
  --download-hdf5-configure-arguments="--with-zlib=yes"
```

Build and check:

```bash
make PETSC_DIR=/home/mdunitz/2025/pflotran/petsc PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=/home/mdunitz/2025/pflotran/petsc PETSC_ARCH=arch-linux-c-opt check
```

Set environment variables:

```bash
export PETSC_DIR=$PWD
export PETSC_ARCH=$(ls -d arch-* | head -n1)
```

#### PFLOTRAN

```bash
cd ..
git clone https://bitbucket.org/pflotran/pflotran
cd pflotran/src/pflotran
make pflotran
```

#### Running a simulation

```bash
export PFLOTRAN_DIR=/home/mdunitz/2025/pflotran
module load mpi/openmpi-x86_64
mpirun -n 1 $PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix filenamehere
```

Use one of the test files in `reference/` to verify installation. You will also need `sandbox/hanford.dat` — it includes species not in PFLOTRAN's default database. Make sure the `.in` file's `DATABASE` line points to the correct path.

---

## Module Details

### generator/

`PFLOTRANGenerator` produces complete `.in` files via a kwargs interface:

```python
from pflotran_generator import PFLOTRANGenerator

gen = PFLOTRANGenerator(
    concentrations={'Cl-': '2.68 Z', 'Na+': '2.295 T'},
    aw_threshold=0.6,
    dimensions='1d',        # '1d' | '2d' | '3d'
    temperature=8.0,
    enable_cl_inhibition=True,
    enable_aw_sandbox=True,
)
gen.generate('simulation.in')
```

Key files:
- `pflotran_generator.py` — Generator class with configurable rate constants, half-saturations, inhibition thresholds, grid presets (1D/2D/3D)
- `pflotran_templates.py` — Static PFLOTRAN input blocks (17 primary species, 62 secondary, 5 gas, 3 mineral)
- `REFERENCES.md` — Full citations for all rate constants and parameters
- `9_addnitrogen_example.in` — Example generated output

### sandbox/

Fortran 90 modules that extend PFLOTRAN with water-activity-dependent inhibition of methanogenesis. Must be compiled into PFLOTRAN (see below and `step_by_step_instructions.md`).

| Module | Pathway | Reaction |
|--------|---------|----------|
| `reaction_sandbox_awinhibit.F90` | Hydrogenotrophic | CO₂ + 4H₂ → CH₄ + 2H₂O |
| `reaction_sandbox_awinhibitacetate.F90` | Acetoclastic | CH₃COO⁻ + H⁺ → CH₄ + CO₂ |
| `reaction_sandbox_awinhibitmethyl.F90` | Methylotrophic | CH₃OH → ¾CH₄ + ¼CO₂ + ½H₂O |

Each implements threshold and smoothstep inhibition modes. `hanford.dat` is the thermodynamic database.

#### Compiling reaction sandboxes into PFLOTRAN

To add a sandbox (e.g. `awinhibit`) to your PFLOTRAN build:

1. Copy the `.F90` file into PFLOTRAN's `src/pflotran/` directory alongside the other reaction sandboxes.

2. Add to `reaction_sandbox.F90`:
   ```fortran
   use Reaction_Sandbox_AWInhibit_class
   ```
   and in the select-case block:
   ```fortran
   case('AWINHIBIT')
       new_sandbox => AWInhibitCreate()
   ```

3. Add to `pflotran_object_files.txt` (in the `chem_obj` section):
   ```
   ${common_src}reaction_sandbox_awinhibit.o \
   ```

4. Add to `pflotran_dependencies.txt`:
   ```
   reaction_sandbox_awinhibit.o : \
     reaction_sandbox_base.o \
     reactive_transport_aux.o \
     global_aux.o \
     reaction_aux.o
   ```
   Also add all three to the dependency list of `reaction_sandbox.o`:
   ```
   reaction_sandbox.o : \
     global_aux.o \
     input_aux.o \
     material_aux.o \
     option.o \
     output_aux.o \
     pflotran_constants.o \
     reaction_aux.o \
     reaction_sandbox_base.o \
     reaction_sandbox_bioTH.o \
     reaction_sandbox_biohill.o \
     reaction_sandbox_calcite.o \
     reaction_sandbox_chromium.o \
     reaction_sandbox_clm_cn.o \
     reaction_sandbox_equilibrate.o \
     reaction_sandbox_example.o \
     reaction_sandbox_flexbiohill.o \
     reaction_sandbox_gas.o \
     reaction_sandbox_pnnl_cyber.o \
     reaction_sandbox_pnnl_lambda.o \
     reaction_sandbox_radon.o \
     reaction_sandbox_simple.o \
     reaction_sandbox_ufd_wp.o \
     reaction_sandbox_awinhibit.o\
     reaction_sandbox_awinhibitacetate.o\
     reaction_sandbox_awinhibitmethyl.o\
     reactive_transport_aux.o \
     string.o \
     utility.o
   ```

5. Repeat for all three sandbox files (`awinhibit`, `awinhibitacetate`, `awinhibitmethyl`).

6. Recompile: `cd pflotran/src/pflotran && make pflotran`

Constants in the reaction sandboxes can be changed without recompiling PFLOTRAN — they're read from the `.in` file at runtime.

### visualization/

Post-processing pipeline for PFLOTRAN Tecplot output:

| Step | Script | Output |
|------|--------|--------|
| 1 | `step1_extract.py` | `.tec` → DataFrame → `pflotran_data.csv` |
| 2 | `step2_plot.py` | 3D Plotly scatter with time animation |
| 3 | `step3_flux.py` | Bokeh gradient + diffusive flux surface maps |
| 4 | `step4_plotflux.py` | Bokeh gradient + flux time series |

`shared_utils.py` is the single source of truth for:
- Concentration gradient computation (∇C)
- Diffusive flux via Fick's First Law: J = −D·∇C (with Stokes-Einstein temperature correction)
- Species-specific diffusion coefficients from Boudreau (1997)
- Consistent column naming and Bokeh tooltip generation

`step_orchestra.py` runs the full pipeline with configurable temperature (default 8°C) and flux computation. By default it reads sample data from `../sample_data/`.

### sample_data/

Example PFLOTRAN Tecplot output files (`test29-000.tec` through `test29-005.tec`) and extracted CSV (`pflotran_data.csv`). Used by the visualization pipeline for testing. Replace with your own simulation output for real analysis.

### batch/

- `create_modified_files.py` — Generates `.in` files at 1×–N× seawater concentration. Seawater baseline from Millero (2013). Multiplier math stays here, not in generator.
- `run_inhibition_test.py` — Three-variant diagnostic (Cl⁻ only / aₓ only / both) for detecting double-counting of salinity inhibition.
- `run_pflotran_batch.sh` — Sequential runner: finds `*.in`, runs each via `mpirun`, logs to `logs/`. Runs in sequence (not parallel) to avoid memory issues. Skips failed runs and reports them at the end.

### reference/

Historical input decks from prior PFLOTRAN runs (not generated by the Python tooling):
- `Erin_UC_increasedT_*.in` — Original and updated-syntax versions
- `test9_w_nitrogen.in` — Nitrogen-inclusive variant

### compare/

- `comparing_aw.ipynb` — Plot CH₄/CO₂ flux as a function of water activity across multiple simulation runs. Update `working_directory` in the notebook to point to your `.tec` output directory.

### exploratory/

Sanskriti's exploratory PFLOTRAN work, covering iterative input deck development, custom reaction sandbox testing, and long-term CO2 sequestration modeling.

- `pflotran/testing/` -- Sequence of `.in` files (1-10) building up from a simple rain-grid methanol model to full nitrogen + sulfate inhibition. Each file adds one new feature to the previous.
- `pflotran/sandbox/` -- Source copies of the Fortran reaction sandbox modules used during development.
- `pflotran/visualization/` -- Tecplot extraction and Plotly 3D visualization scripts, plus `.tec` output files from test runs (test29, test30 series).
- `long_term_isq/` -- Streamlit dashboard for ISQ gas concentration data (CH4, CO2, H2S) across experimental conditions (methanogen, spirulina, mix). Includes CSV datasets and interactive plotting.
- `CCS_LT.in` + `ccs_lt.grdecl` -- Long-term (100-year) CO2 capture and sequestration input deck with Eclipse grid.
- `pflotran_vars.py` -- PFLOTRAN variable definitions and parameter sets.
- `constants.py` -- Unit conversion constants (molar masses, time conversions).
- Notebooks: water activity curves, 100-year linear/exponential projection models.

---

## Development

### CI

This repo uses GitHub Actions on all PRs against `main`:

| Workflow | What it runs |
|----------|--------------|
| `ci.yml` | `black`, `flake8`, `pytest -m "not integration"` (no PFLOTRAN needed) |
| `integration.yml` | Builds `Containerfile`, runs `pytest -m integration` (real PFLOTRAN simulations) |

### Running checks locally

```bash
pip install -r requirements.txt
black --check .
flake8
pytest -m "not integration"    # fast: no PFLOTRAN binary required
```

Integration tests (require Docker — see [Running PFLOTRAN](#running-pflotran)):

```bash
docker build -t pflotran-py-test -f Containerfile .
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  pytest tests/ -v --tb=short -m integration
```

---

## Key Physical Parameters

**Default temperature:** 8°C (isothermal, coastal sediment)

**Seawater 1× baseline** (Millero 2013):

| Ion | Concentration [mol/L] |
|-----|----------------------|
| Cl⁻ | 0.536 |
| Na⁺ | 0.459 |
| Mg²⁺ | 0.0523 |
| SO₄²⁻ | 0.0276 |
| Ca²⁺ | 0.0100 |
| K⁺ | 0.00972 |

**Inhibition mechanisms:**
- **Cl⁻ Monod inhibition** — empirical, half-saturation at 0.2 M, applied to all methanogenesis reactions
- **Water activity (aₓ) sandbox** — thermodynamic, threshold-based, captures effect of ALL dissolved ions

See `generator/REFERENCES.md` for complete parameter sourcing.

## Workflows

Adapt paths below if your clone location differs from the defaults.

### Regenerating the main input file with different constants

Navigate to `generator/`. Change the constants you want using the dictionaries at the top of `pflotran_generator.py`. Default output is `.tec` files. Follow the comments in the generator if you want HDF5 output (recommended for multi-salinity runs).

### Generating files with different salt concentrations

Navigate to `batch/`. Copy in the `9_addnitrogen.in` generated earlier (or use the one already there if you haven't changed constants). Run:

```bash
python create_modified_files.py
```

By default it creates 20 files from 1x to 20x seawater. The seawater ion concentrations were calculated from: https://docs.google.com/spreadsheets/d/1iNVlg_OOcvQkkKXAuV_2iWS9l-c619CPPVjG7pkcaQE/edit?pli=1&gid=0#gid=0

### Running multiple files

```bash
cd batch/
chmod +x run_pflotran_batch.sh
./run_pflotran_batch.sh
```

Runs all `.in` files in the directory sequentially (not parallel). Skips failures and reports them at the end.

### Visualizing single PFLOTRAN runs

```bash
cd visualization/
python step_orchestra.py
```

Edit the config at the top of `step_orchestra.py` to point to your data directory, set `tec` or `hdf5` format, and choose species. Set `verbose=True` in individual step functions for debugging. Opens three HTML files.

### Visualizing multiple conditions (varying water activity)

Run the notebook at `compare/comparing_aw.ipynb`. Update `working_directory` to point to your `.tec` output directory:

```python
working_directory = "batch/modified_pflotran_files"
```

## General PFLOTRAN Tips

- Any new reactions must be stoichiometrically balanced.
- Any new species must exist in `hanford.dat`. If not, add them manually.
- Default output is `.tec` (Tecplot). For multi-condition runs, HDF5 output is recommended — set it in the `.in` file's `OUTPUT` block.
- Check the [PFLOTRAN documentation](https://documentation.pflotran.org/) before adding new chemistry.

