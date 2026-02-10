# PFLOTRAN Reactive Transport Modeling

Simulation of microbial redox networks (methanogenesis, sulfate reduction, iron cycling) in saline sediments using [PFLOTRAN](https://www.pflotran.org/).

## Directory Structure

```
pflotran/
├── README.md
├── requirements.txt       ← Python dependencies for this module
├── generator/             ← Python code to produce PFLOTRAN .in files
├── sandbox/               ← Custom Fortran 90 reaction modules (water activity inhibition)
├── visualization/         ← Post-processing pipeline (extract → plot → gradient/flux)
├── batch/                 ← Batch file generation + inhibition diagnostics
├── sample_data/           ← Example PFLOTRAN Tecplot output files
├── reference/             ← Historical/reference input decks
└── compare/               ← Comparison notebooks
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

## Installing PFLOTRAN

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

Build and check (use the exact paths your configure step prints):

```bash
make PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-linux-c-opt check
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
export PFLOTRAN_DIR=/path/to/pflotran
mpirun -n 1 $PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix filenamehere
```

Use one of the test files in `reference/` to verify installation. You will also need `sandbox/hanford.dat` — it includes species not in PFLOTRAN's default database. Make sure the `.in` file's `DATABASE` line points to the correct path.

### HPC (Slurm)

```bash
srun -A ccsi -p burst -N 1 -n 1 -c 1 --mem=8g -t 00:30:00 --pty bash -i
module load mpi/openmpi-x86_64
$PFLOTRAN_DIR/src/pflotran/pflotran -pflotranin your_file.in
```

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
   Also add `reaction_sandbox_awinhibit.o\` to the dependency list of `reaction_sandbox.o`.

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

## General PFLOTRAN Tips

- Any new reactions must be stoichiometrically balanced.
- Any new species must exist in `hanford.dat`. If not, add them manually.
- Default output is `.tec` (Tecplot). For multi-condition runs, HDF5 output is recommended — set it in the `.in` file's `OUTPUT` block.
- Check the [PFLOTRAN documentation](https://documentation.pflotran.org/) before adding new chemistry.
