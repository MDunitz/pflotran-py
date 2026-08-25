# PFLOTRAN Reactive Transport Modeling

Simulation of microbial redox networks (methanogenesis, sulfate reduction, iron cycling) in saline sediments using [PFLOTRAN](https://www.pflotran.org/).

## Directory Structure

```
├── README.md
├── requirements.txt          ← Python dependencies
├── Containerfile             ← Docker/Podman image with custom-sandbox PFLOTRAN
├── scripts/                  ← PFLOTRAN sandbox patch + build scripts
├── sandbox/                  ← Custom Fortran 90 reaction modules (water activity inhibition)
├── src/pflotran_py/          ← Installable package
│   ├── generator/            ← Python code to produce PFLOTRAN .in files
│   │                           (including bottle_generator.py, closed-batch decks)
│   │                           constants.py: AWINHIBIT a_crit / inhibition type
│   ├── analysis/             ← Post-PFLOTRAN compute (extract, gradients, transforms)
│   │                           constants.py: diffusion coefficients, water viscosity
│   ├── comparison/           ← Model vs. measured laboratory incubations
│   ├── visualization/        ← Presentation (Bokeh / Plotly rendering)
│   └── config.py             ← Pipeline knobs (species map, temperature, output paths)
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

For the bottle decks there is a wider five-variant version that also separates
the two chloride terms from sulfate competition and from the carbon-supply
ceiling; see [Attributing the modelled salt
suppression](#attributing-the-modelled-salt-suppression).

### 3. Run PFLOTRAN simulations

```bash
cd batch/
./run_pflotran_batch.sh /path/to/input/files
```

### 4. Post-process and visualize results

```bash
python -m pflotran_py.visualization.step_orchestra  # runs step1 → step2 → step3 → step4
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
# Docker
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  pytest tests/test_docker_e2e.py -v --tb=short

# Podman (add :Z on SELinux hosts)
podman run --rm -v "$(pwd)":/work:Z -w /work pflotran-py-test \
  pytest tests/test_docker_e2e.py -v --tb=short
```

**Custom sandbox tests** (AWINHIBIT decks 7–10):

```bash
# Docker
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  pytest tests/test_custom_docker_e2e.py -v --tb=short

# Podman
podman run --rm -v "$(pwd)":/work:Z -w /work pflotran-py-test \
  pytest tests/test_custom_docker_e2e.py -v --tb=short
```

**All integration tests:**

```bash
# Auto-detect docker or podman (recommended)
./scripts/run_integration.sh

# Or pick explicitly
./scripts/run_integration.sh docker
./scripts/run_integration.sh podman
```

**Run the full custom-sandbox pipeline** (prep deck, simulate, extract, plot):

```bash
# Docker
docker run --rm -v "$(pwd)":/work -w /work pflotran-py-test \
  python3 tests/test_custom_docker_e2e.py

# Podman
podman run --rm -v "$(pwd)":/work:Z -w /work pflotran-py-test \
  python3 tests/test_custom_docker_e2e.py
```

Output lands in `tests/custom_e2e_output/` (CSV + PNGs).

On Linux with Podman, `:Z` on the volume mount is required for SELinux
(Fedora/RHEL). It is harmless on hosts without SELinux.

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
- `constants.py` — Deck defaults with citations: AWINHIBIT `a_crit` per pathway (`AW_CRIT_*`) and `AW_INHIBITION_TYPE`
- `pflotran_generator.py` — Generator class with configurable rate constants, half-saturations, inhibition thresholds, grid presets (1D/2D/3D)
- `pflotran_templates.py` — Static PFLOTRAN input blocks (17 primary species, 62 secondary, 5 gas, 3 mineral)
- `REFERENCES.md` — Full citations for all rate constants and parameters
- `9_addnitrogen_example.in` — Example generated output

### Where constants live

Three files named `constants.py`; they are not interchangeable:

| File | What it holds |
|------|----------------|
| `src/pflotran_py/generator/constants.py` | Deck-generator defaults written into `.in` files: sealed-bottle geometry/run setup (`VIAL_VOLUME_L`, `BOTTLE_*`) and pathway `a_crit` / `ONE_MINUS_AW` (Oren 1999/2011). CLI flags and generators import from here. |
| `src/pflotran_py/analysis/constants.py` | Post-processing physics: Boudreau 25 °C diffusion coefficients, Vogel water viscosity, unit-conversion factors. Used by the visualization physics layer, not by deck generation. |
| `exploratory/constants.py` | Notebook/unit-conversion leftovers (molar mass of C, sample volume, time factors). Exploratory only; the package does not import it. |

Pipeline settings that are not literature constants (species map, default temperature, output paths) live in `src/pflotran_py/config.py`.

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

### src/pflotran_py/visualization/

Post-processing pipeline for PFLOTRAN Tecplot output (importable as `pflotran_py.visualization`):

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

`step_orchestra.py` runs the full pipeline with configurable temperature (default 8°C) and flux computation; run it as `python -m pflotran_py.visualization.step_orchestra`. By default it reads sample data from `sample_data/`.

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
- `constants.py` -- Exploratory unit conversions (molar mass of C, sample volume, time factors). Not imported by the installable package; see `src/pflotran_py/generator/constants.py` and `analysis/constants.py` for those.
- Notebooks: water activity curves, 100-year linear/exponential projection models.

---

## Comparing the model against the laboratory incubations

Everything above models an open sediment column. This section is about a
different physical system: the sealed 125 mL vials the laboratory actually
incubates, and comparing what the model predicts for them against what the gas
chromatograph measured.

### Why the column deck cannot be used for this

The column decks put a one-cubic-metre domain under an atmospheric boundary, so
water and solutes move through it. Run against a sealed bottle that is not a
small inaccuracy, it is the wrong system. In the five-day test output committed
with the repository, chloride falls from 6.2 M to 7.6e-14 M within a day as the
boundary flushes the brine out of the domain, taking with it the salt stress the
experiment exists to study, and methane never leaves its 1e-15 numerical floor.

`BottleGenerator` (in `generator/bottle_generator.py`) subclasses
`PFLOTRANGenerator` and changes only the physical setup, inheriting the reaction
network, the water-activity sandboxes and the species lists unchanged:

| | column deck | bottle deck |
|---|---|---|
| Domain | 1 m³, 10–64 cells | one cell of 125 mL |
| Boundaries | atmospheric top | **none** — sealed |
| Gas / liquid | pore space | 100 mL headspace, 25 mL liquid |
| Duration | 31 days | 130 days |
| Temperature | 8 °C | 18 °C |

Removing every boundary condition is what closes the bottle: PFLOTRAN treats a
face with no boundary condition as zero-flux. With that change chloride is
retained and methane is produced.

A single-cell domain has no spatial extent, so all concentration gradients in it
are zero and the flux stages of the post-processing pipeline produce nothing
meaningful for these runs. That is correct for a bottle. What a bottle gives you
is a concentration time series.

### Two chemistry changes the bottle decks make

Both are off by default in `PFLOTRANGenerator`, so column decks are unaffected,
and on by default in `BottleGenerator`.

**`couple_carbonate`** moves dissolved carbon dioxide from the primary species
list to the secondary list, so PFLOTRAN computes it from bicarbonate and pH
through the database reaction `CO2(aq) = HCO3- + H+ - H2O`. Without this the
model cannot predict headspace carbon dioxide at all: no reaction in the network
produces `CO2(aq)` — every carbon-oxidising step yields bicarbonate — and as a
decoupled primary species it therefore never moves from its initial value
however much carbon is respired. Coupling it also repairs an inconsistent
initial state, where the deck began with dissolved carbon dioxide that
contradicted its own bicarbonate and pH.

Note that the same argument does **not** apply to `CH4(aq)`, whose database
entry is a redox reaction rather than an acid-base one. Decoupling methane is
what allows the kinetic network to produce it and must stay.

**`methane_gas_phase`** declares `CH4(g)` as an active gas species and removes
the ebullition proxy reaction. The proxy, `CH4(aq) -> Tracer2` gated at
2.5e-3 M, models methane leaving as bubbles once the liquid is supersaturated
enough for them to nucleate. That is a reasonable picture of submerged sediment
and the wrong picture of a bottle, where the headspace already exists and
methane partitions into it from the first molecule produced, with no threshold.

This needs a small database change, because `hanford.dat` defines `CH4(g)`
against `Methane(aq)`, an organic species on an ethane basis that this network
does not carry, while the deck's methane is `CH4(aq)`, a redox species on a
bicarbonate basis. The two are the same molecule — both are recorded at
16.0428 g/mol — written against different basis species. What transfers is the
equilibrium constant: the tabulated `CH4(g)` log K series *is* the methane
Henry constant, giving 1.41e-3 mol/(L·atm) at 25 °C against Sander (2023)'s
1.4e-3, and agreeing to two percent at 0 °C. Volatility belongs to the molecule,
not to the basis chosen to describe its dissolved form.

`write_bottle_database()` therefore writes `sandbox/hanford_bottle.dat`, a copy
of the database with that one line's aqueous partner rewritten and its constants
untouched. The shared database is never modified, and the change is visible as a
diff between two files.

As a check on the result, the partition PFLOTRAN then performs internally agrees
with an independent Henry's-law calculation to within 16 percent at 18 °C.

### The comparison package

`comparison/` holds everything that touches measured data.

| Module | Role |
|---|---|
| `brines.py` | Reads the measured salt recipes from the batch workbook and derives per-ion concentrations for each batch |
| `corrections.py` | Documented corrections to the recorded data, each with its evidence |
| `decks.py` | One closed-batch deck per measured batch, built from the weighed salt |
| `run_decks.py` | Runs the decks through PFLOTRAN in the container |
| `headspace.py` | Converts between model concentrations and headspace moles |
| `figures.py` | The comparison figures |
| `palette.json` | Colour roles, anchored on the colourblind-safe palette used in the measurement repository |

Decks are built from the salt that was weighed out, not from a target water
activity inverted through an idealised salt. Water activity is what the model is
being asked to predict — PFLOTRAN computes it from composition at each timestep,
and the inhibition sandboxes act on the value it computes. Feeding in a
composition reverse-engineered from the measured water activity would hand the
model part of its own answer. Building from the recipe keeps the measured water
activity as an independent check.

### Two further options for the carbon inventory and the inhibition

Both default off, so nothing changes unless asked for.

**`cellulose_hydrolysis`** holds the substrate carbon in a solid pool that
dissolves into solution, rather than as 5 mol/L of dissolved DOM1. This matters
because DOM1 is glucose -- the database gives its molar mass as 180.1566 and
labels the related solid pool "TAO-glucose" -- so 5 mol/L is 901 g/L, which is
glucose's solubility limit. PFLOTRAN computes water activity as
`1 - 0.017 * sum(molality)`, counting every solute, so at that concentration the
organic pool contributed about ninety percent of the osmolality in an unsalted
bottle: the model's control sat at a water activity of 0.909 while the meter
read 1.000. With hydrolysis on, dissolved DOM1 is millimolar and water activity
is salt-driven; the solid pool's volume fraction is set so total starting C
matches the incubations' recipe-derived biomass carbon (~0.0565 mol; VF
0.012182). See "Starting carbon" below.

This uses the `Cellulose_min` record already in the database, which needed one
repair: it declares two species but carries two surplus fields, and its second
species has a coefficient of zero that PFLOTRAN drops, so the deck is refused
with a species-count mismatch. `write_bottle_database()` fixes the record
alongside the methane one, leaving the chemistry unchanged.

**`salinity_inhibition`** (optional) adds a sigmoidal Cl⁻ term on top of the
network's own methanogenesis reactions. The comparison pipeline **does not**
pass it any more. The inhibition diagnostic showed that fitted smoothstep
(threshold 0.75 mol/L, interval 1.0) was responsible for an ~850× collapse onto
a methane floor at mid/high salt, while the network's legacy 0.2 M Cl⁻ Monod
alone gives a gradual decline closer to the measured shape. Keeping both was
double-counting chloride. The flag remains available for attribution runs:

```python
salinity_inhibition={"species": "Cl-", "threshold": 0.75, "interval": 1.0}
```

The AWINHIBIT sandboxes **own** methanogenesis: each carries the matching
network Monod rate law (same rate constant, half-saturations, and O₂/Fe/H⁺
inhibition) and multiplies by an a_w smoothstep. The network's three
methane-producing `MICROBIAL_REACTION` blocks are omitted so the two do not
double-produce methane. Comparison decks therefore use `--no-cl-inhibition`
and `--aw-inhibition-type ONE_MINUS_AW` with pathway-specific
`a_crit` from `pflotran_py.generator.constants` (hydrogenotrophic 0.80,
methylotrophic 0.85, acetoclastic 0.90 — acetoclasts fail first under salt;
Oren 1999/2011, not a methane fit),
and pass each batch's **meter-read** water activity as
`FIXED_WATER_ACTIVITY`. PHREEQC/`pitzer.dat` a_w computed from the weighed
recipe is an independent oracle, not the default inhibition input: it is
near-exact for 1:1 NaCl and ~0.02 high for the Mg brines, so feeding it in
would bias the Na-vs-Mg contrast. Pass `--use-computed-aw` only for
sensitivity checks. Rebuild the container after pulling sandbox Fortran
changes.

### Running it

```bash
# 1. Pull the measured recipes and derive per-batch ion concentrations.
python -m pflotran_py.comparison.brines --output data/incubation_batch_composition.csv

# 2. Build one closed-batch deck per measured batch.
#    AWINHIBIT sandboxes own methanogenesis (network Monod rates + a_w).
#    Cl- Monod is off so salt is not double-counted.
python -m pflotran_py.comparison.decks --output-dir decks \
    --cellulose-hydrolysis --no-cl-inhibition \
    --aw-inhibition-type ONE_MINUS_AW \
    --aw-threshold 0.80 --aw-threshold-methyl 0.85 --aw-threshold-acetate 0.90

# 3. Run them (needs the container image built; see Running PFLOTRAN above).
python -m pflotran_py.comparison.run_decks --run-root runs --clean

# 4. Build the figures (absolute moles and per-starting-C companions).
python -m pflotran_py.comparison.figures
```

Step 1 needs network access; the workbook is published as CSV and needs no
credentials. Step 4 additionally needs the measurement pipeline's `.ecsv`
output, produced in the saltyBiomass repository by
`make run-mzml-pipeline FULL=1`; point at it with `--ecsv-glob`.

Each deck runs in about a second, because the domain is a single cell.

### Converting between the two sides

The model reports aqueous concentrations; the chromatograph reports headspace
moles. `headspace.py` converts the first into the second.

Where a deck carries a gas phase, the gas-phase concentration is read directly
and multiplied by the headspace volume — PFLOTRAN has already done the
partition, and reporting it is the whole job. Where a deck does not, the model's
entire inventory of the gas is summed and then **partitioned** between headspace
and liquid using Henry's law with a Setschenow salting-out correction.

That distinction is worth stating plainly, because getting it wrong is easy and
quiet. Multiplying a dissolved concentration by the partition coefficient and
the headspace volume answers the question "if this liquid were equilibrated
against a headspace, what would the headspace hold?" For a deck whose dissolved
concentration is not in equilibrium with any headspace — because the deck has
none — that answer exceeds the total gas the model ever made, by roughly
`(K·V_gas + V_liquid) / V_liquid`, about a hundredfold for methane in this vial.
An earlier version of this comparison did exactly that and overstated the model
by a factor of about thirty. `partition_total_moles` is mass-conserving by
construction and the test suite asserts that the prediction can never exceed the
gas the model produced.

At 18 °C in this vial, 99 percent of methane sits in the headspace and about
80 percent of carbon dioxide; salt raises both, because dissolved salt reduces
gas solubility.

The reverse conversion — measured headspace carbon dioxide back to total carbon
produced — is deliberately **not** offered. It needs carbonate speciation, and
the available equilibrium constants are calibrated for seawater to an ionic
strength near 0.7 mol/L while these brines run from 1.2 to 5.8 mol/L. That is an
extrapolation of up to eightfold, and the result would look like a measurement
while being closer to a guess.

### Starting carbon: what we have, what we plot, what we claim

**Lab TOC / volatile solids were not assayed** for Exp003 / Exp004. What the
measurement pipeline *does* carry is a recipe-derived starting carbon inventory:

```text
dry biomass (g) = incubation mass
                × (sludge fraction × (1 − 0.732)
                   + spirulina fraction × (1 − 0.05))
moles starting C = dry biomass × 0.5 / 12.011 g mol⁻¹
```

(`calculate_dry_biomass_in_incubation` / `convert_dry_biomass_to_moles_C` in
saltyBiomass). That is dry sludge (+ spirulina when present) assumed 50% C by
mass — the same denominator the GC dashboard uses for “percent of starting
biomass C released.” Across Exp003/Exp004 it is nearly constant at about
**0.0565 mol C per bottle** (batch-to-batch spread only a few percent); see
`MEASURED_INCUBATION_STARTING_CARBON_MOLES` in `comparison/carbon_inventory.py`.

**The closed-batch decks match that inventory.** With `--cellulose-hydrolysis`,
solid `Cellulose_min` uses volume fraction **0.012182** so model starting C is
**≈ 0.0565 mol** (six carbons per glucose unit, plus 1 mM dissolved DOM1). That
is inventory normalisation — one geometric factor, no rate constants touched —
not a kinetics fit. An earlier default of VF 0.2 held ~0.925 mol C (~16× too
much) and made absolute mole overlays agree with the chromatograph only by luck.

`cellulose_volume_fraction_for_starting_carbon()` converts a target moles-C
into the VF; the default is the value for 0.0565 mol. After changing VF,
re-check control water activity and dissolved DOC under the existing
hydrolysis-rate constraints (the rate was tuned at the old inventory).

**Figures.** Absolute-mole overlays and per-starting-C companions now share the
same matched denominator (~0.0565 mol C). Prefer the companions for yield /
Madison-style questions; use the absolute overlays for trajectory shape and
salt ranking.

| Figure | Role |
|---|---|
| `*_per_starting_c.png` | Yield (mol gas / mol starting C) |
| `methane_over_time.png`, `carbon_dioxide.png`, … | Trajectory shape and salt ranking |

### Which numbers were fitted, and which were not

Read this before quoting any agreement statistic from this comparison.

**One parameter is still fitted against the measurements for the carbon
inventory. The fitted Cl⁻ smoothstep is no longer used in the comparison.**

| Parameter | Value | Fitted against | How | Status |
|---|---|---|---|---|
| Cellulose hydrolysis rate | `2.d-8` mol/m²/s | Exp004 **water activity**, not methane | Sweep of four values | In use |
| Salinity inhibition threshold | 0.75 mol/L Cl⁻ | Exp004 methane | Grid search, 24 combinations | **Retired** from comparison default |
| Salinity inhibition interval | 1.0 decades | Exp004 methane | Same grid search | **Retired** from comparison default |

The smoothstep pair was selected by `pflotran_py.comparison.calibrate` on
Exp004 alone. The inhibition diagnostic later showed it produced the modelled
methane cliff (~850× at mid/high salt) and that the network's own 0.2 M Cl⁻
Monod, used alone, recovers a gradual decline. Those two fitted numbers remain
in `calibrate` / `forecast` for reproducibility of the old protocol; the
comparison decks omit `--salinity-threshold`.

The hydrolysis rate was chosen so that the unsalted bottle's modelled water
activity matched the meter reading of 1.000, and so that the dissolved organic
pool landed at a concentration an active sludge porewater plausibly holds.
Methane was not consulted. A faster rate leaves the organic pool depressing
water activity; a slower one makes carbon supply itself the limiting factor.

A fourth number, the solid carbon volume fraction (**0.012182**), is
**not** a methane fit: it is set so model starting C matches the
recipe-derived incubation inventory (~0.0565 mol C). See "Starting carbon"
above.

**Not fitted, and not adjusted at any point:** the sixteen-reaction network and
every rate constant and half-saturation in it (see `generator/REFERENCES.md`),
the Henry solubilities and Setschenow coefficients, all fifteen brine
compositions (derived from the weighed recipes), the bottle geometry, the
temperature, and the run duration. Nor are the four structural changes fits --
sealing the domain, giving methane a gas phase, coupling carbonate, and matching
the solid carbon pool to measured starting C each correct a specific defect and
introduce no free kinetic parameter. Each was checked against something
independent: chloride retention, agreement between PFLOTRAN's internal methane
partition and Henry's law to within 16 percent, mass conservation across deck
versions, the measured water activity of the control, and (for the carbon pool)
equality with recipe-derived dry-biomass moles C.

**Carbon dioxide has no fitted parameters at all.** Nothing was ever tuned
against it. Its panel in the parity figure is a parameter-free prediction,
though an indirect one: the hydrolysis rate and the methane inhibition both
change how carbon is routed.

### Fitting on one experiment and testing on the other

```bash
python -m pflotran_py.comparison.calibrate
```

Exp004 uses sodium chloride and magnesium chloride brines; Exp003 uses
artificial sea salt, which carries sulfate and a different divalent balance.
Parameters fitted on the first and applied unchanged to the second are being
asked to transfer across salt chemistry, not merely across replicates.

| | batches | score | median miss | within a factor of ten |
|---|---|---|---|---|
| Exp004, fitted on | 7 | 0.79 | — | — |
| Exp003, held out | 8 | **1.01** | factor of 7.6 | 50% |

The score degrades from 0.79 to 1.01 on transfer, which is the expected
direction and a modest amount: a typical miss grows from about sixfold to about
tenfold.

Two things temper this. The objective surface is flat -- seven of the
twenty-four grid points score between 0.79 and 0.80 -- so the two parameters are
only weakly identified, and a different tie-break would have chosen differently
without changing the fit quality. And the held-out misses are not random: the
model under-predicts almost every sea-salt condition, worst at the mid strength
by a factor of about forty. Sea salt carries sulfate, so sulfate reduction is
already competing with methanogenesis there, and a chloride-keyed inhibition
term penalises those bottles a second time for the same salt.

**This is a pre-registered protocol, not a blind prediction.** Exp003 was
examined during the work that produced this design, and that knowledge cannot be
unlearned; it can leak into choices as ordinary as which grid to sweep. Read the
held-out number as evidence about whether two parameters transfer across salt
systems, which it can genuinely answer, and not as a forecast of unseen data,
which it cannot.

### What the comparison currently shows

With the domain sealed, a real methane gas phase, coupled carbonate, the carbon
pool matched to recipe-derived starting C (~0.0565 mol), and **AWINHIBIT
sandboxes owning methanogenesis** (network Monod kinetics × a_w smoothstep at
threshold 0.95; network methanogenesis and Cl⁻ Monod off), salt stress is keyed
on water activity rather than double-counted chloride.

The fitted Cl⁻ smoothstep remains retired (see the inhibition diagnostic). The
a_w threshold is placed at the top of the measured salted range so inhibition
engages across Exp003/Exp004; it is not fitted to methane yields. Re-run decks
after rebuilding the container image that patches the sandbox Fortran.

One limit is not addressable from a deck. PFLOTRAN computes water activity as
`1 - 0.017 * sum(molality)`, ideal Raoult with no osmotic coefficient, so it
cannot capture the non-ideality of concentrated brine: the strongest magnesium
chloride bottle computes 0.884 against a measured 0.824. Water activity is now a
salt-driven quantity and a usable diagnostic, but it is systematically too high
at the salty end, and correcting that is a Fortran change.

### Attributing the modelled salt suppression

Modelled methane falls off a cliff between water activity 0.94 and 0.93 and then
sits on a floor near 10⁻⁶ mol per mol starting C, while the measurements decline
gradually across the whole range. Four things could produce that, and all four
are currently switched on at once, so the comparison figures cannot say which:

1. **A chloride Monod inhibition baked into the reaction network** at 0.2 mol/L,
   attached to the three methanogenesis pathways *and to fermentation*. Every
   salted batch carries at least 1.2 mol/L chloride, so `K / (K + C)` runs from
   about a seventh to a thirtieth — and because fermentation is throttled too,
   that suppression compounds through the carbon chain. `REFERENCES.md` records
   the threshold as empirical, with no citation. A generated deck contains four
   copies of this term and three of the smoothstep, which is the double count in
   plain sight.
2. **The fitted chloride smoothstep**, on the three methanogenesis reactions and
   deliberately not on fermentation. It was fitted on Exp004, whose brines carry
   no sulfate, so it had to absorb the entire salt effect of a sulfate-free
   solution.
3. **Sulfate competition.** Sulfate reducers outcompete methanogens for acetate
   and hydrogen, and sulfate-dependent anaerobic methane oxidation destroys
   methane after it is made. Both are real chemistry, and both act only on the
   Exp003 sea-salt brines — which is how a term fitted on sulfate-free brines
   ends up over-suppressing sea salt when carried across.
4. **Carbon supply.** With the pool now matched to ~0.0565 mol C, hydrolysis may
   be the binding constraint, which would flatten yield across brines for
   reasons having nothing to do with salt.

`comparison/inhibition_diagnostic.py` runs the five decks that separate these.
It fits nothing and adds no parameter: each variant is the existing deck with
one term switched off.

| Variant | Legacy Cl⁻ Monod | Fitted smoothstep | Sulfate pathways |
| --- | --- | --- | --- |
| `baseline` | on | on | on |
| `no_cl_monod` | **off** | on | on |
| `no_smoothstep` | on | **off** | on |
| `no_salt_terms` | off | off | on |
| `ceiling` | off | off | **off** |

Because inhibition factors multiply, each pair that differs in exactly one
switch gives the fold change that term is responsible for, per batch. The point
is the *shape* across water activity rather than the size anywhere: if removing
one term turns the cliff into a gradual decline, that term produced the cliff.

```bash
bash scripts/run_inhibition_diagnostic.sh          # 75 runs, budget an hour or two
```

Or a stage at a time, which is useful because only the middle one needs Docker:

```bash
python -m pflotran_py.comparison.inhibition_diagnostic --stage decks
python -m pflotran_py.comparison.inhibition_diagnostic --stage runs
python -m pflotran_py.comparison.inhibition_diagnostic --stage summarise
```

Results land in `output/comparison/diagnostic/`: `per_variant.csv` (final
methane per variant per batch), `attribution.csv` (the fold changes), and
`inhibition_decomposition.png` (yield against water activity per variant, with
the measurements overlaid, beside a per-batch attribution panel).

The last two variants omit sulfate reduction and anaerobic methane oxidation,
which these incubations certainly perform. They exist to bound the carbon-supply
ceiling and are diagnostics, not predictions; nothing about them should be
carried into a deck used for prediction. The same caution applies to the
`--disable-reactions` flag on `comparison.decks` that makes them possible.

The likely honest outcome is that one of the two chloride terms is retired and
the survivor refitted, with Exp003 still held out — which lowers the parameter
count rather than raising it.

**Update after the diagnostic ran:** that is what happened. The fitted smoothstep
was the cliff (~850×); sulfate competition was ~1×; carbon supply was not
limiting. The comparison pipeline now matches the `no_smoothstep` variant
(legacy Cl⁻ Monod only). The smoothstep remains available for attribution and
for reproducing the old calibrate/forecast protocol.

---

## Development

### CI

This repo uses GitHub Actions on all PRs against `main`:

| Workflow | What it runs |
|----------|--------------|
| `ci.yml` | `black`, `flake8`, `pytest -m "not integration"` (no PFLOTRAN needed) |
| `integration.yml` | Builds `Containerfile` and runs `pytest -m integration` with **both Docker and Podman** (matrix) |

### Running checks locally

```bash
pip install -r requirements.txt
black --check .
flake8
pytest -m "not integration"    # fast: no PFLOTRAN binary required
```

Integration tests (require Docker or Podman — see [Running PFLOTRAN](#running-pflotran)):

```bash
./scripts/run_integration.sh          # auto-detects docker or podman
./scripts/run_integration.sh docker
./scripts/run_integration.sh podman
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

Navigate to `src/pflotran_py/generator/`. Pathway a_w defaults are in `constants.py`. Change other kinetic values using the dictionaries at the top of `pflotran_generator.py`. Default output is `.tec` files. Follow the comments in the generator if you want HDF5 output (recommended for multi-salinity runs).

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
python -m pflotran_py.visualization.step_orchestra
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

