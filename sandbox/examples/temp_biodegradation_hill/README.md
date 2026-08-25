# Temperature-Dependent Biodegradation Hill

Literature use case based on **Flexible Biodegradation Hill** from:

> Hammond, G. E. (2022). The PFLOTRAN Reaction Sandbox. *Geosci. Model Dev.*, 15, 1659–1676. https://doi.org/10.5194/gmd-15-1659-2022

## Conceptual model

```text
Aaq + 0.25 Baq  →  0.33 Caq + Daq
```

mediated by immobile biomass `Xim`, with:

- Hill kinetics on electron donor `Aaq`
- Monod on acceptor `Baq`
- product inhibition by `Caq`
- biomass growth (yield) and first-order decay

**Custom extension:** temperature-dependent activity via Q₁₀:

```text
r(T) = r_ref * Q10 ** ((T - T_ref) / 10)
```

Input keyword: `TEMPBIOHILL`  
Fortran module: [`reaction_sandbox_tempbiohill.F90`](reaction_sandbox_tempbiohill.F90)

Setting `Q10 1.d0` recovers stock `FLEXIBLE_BIODEGRADATION_HILL` exactly.

## Layout

```text
temp_biodegradation_hill/
  reaction_sandbox_tempbiohill.F90   ← custom sandbox (extends BioHill)
  reaction_sandbox_block.in          ← pasteable REACTION_SANDBOX block
  reaction_def.json
  reference_biohill.F90              ← stock PFLOTRAN (for reading)
  reference_flexbiohill.F90          ← stock PFLOTRAN (for reading)
  gold/flexible_biodegradation_hill.regression.gold
  decks/
    01_batch_flexbiohill_reference.in   ← stock FlexBioHill (paper)
    02_batch_tempbiohill_Q10_1.in       ← TempBioHill, Q10=1 (match 01)
    03_batch_tempbiohill_warm.in        ← Q10=2, T=35 C
    04_column_1d_tempbiohill.in         ← 1D flow column, Q10=2
    05_column_1d_Q10_1.in               ← 1D flow column, Q10=1
  scripts/
    compare_mass_balance.py
    run_example.sh
  compare_and_plot.ipynb
```

## Install into PFLOTRAN

From the repo root (after a normal AWINHIBIT patch):

```bash
python3 scripts/patch_pflotran_sandboxes.py \
  --pflotran-src /path/to/pflotran/src/pflotran \
  --sandbox-dir sandbox/ \
  --extra-sandbox-dir sandbox/examples/temp_biodegradation_hill

cd /path/to/pflotran/src/pflotran && make clean && make pflotran
```

`TEMPBIOHILL` depends on stock `reaction_sandbox_biohill.F90` (already in PFLOTRAN).

## Run

```bash
# Stock paper reproduction (any PFLOTRAN with FlexBioHill — including this repo's Docker image)
pflotran -input_prefix decks/01_batch_flexbiohill_reference

# After patching TEMPBIOHILL:
pflotran -input_prefix decks/02_batch_tempbiohill_Q10_1
pflotran -input_prefix decks/03_batch_tempbiohill_warm
pflotran -input_prefix decks/04_column_1d_tempbiohill
pflotran -input_prefix decks/05_column_1d_Q10_1
```

Or use `scripts/run_example.sh` (sets `PFLOTRAN` / Docker).

## Checks

1. **Gold regression (batch):** final free ions from `01` (and `02` with Q10=1) should match [`gold/flexible_biodegradation_hill.regression.gold`](gold/flexible_biodegradation_hill.regression.gold).
2. **Stoichiometric mass balance (batch):**  
   `python scripts/compare_mass_balance.py runs/01_batch_flexbiohill_reference`
3. **Q10 effect:** warm batch (`03`) should deplete `Aaq` faster than `01`/`02`.
4. **1D column:** use [`compare_and_plot.ipynb`](compare_and_plot.ipynb) for profiles vs distance and mid-column time series; overlay `04` vs `05`.

## Note on the notebook generator

[`sandbox/create_custom_sandbox.ipynb`](../../create_custom_sandbox.ipynb) targets AWINHIBIT-style power-law sandboxes. FlexBioHill / TempBioHill need Hill + biomass + analytical Jacobian, so this example is hand-authored in the same `custom_*`-style layout rather than emitted by the v1 generator.
