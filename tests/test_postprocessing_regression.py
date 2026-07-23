"""Numerical regression pins for the post-PFLOTRAN compute chain.

`test_e2e.py` pins extraction against `sample_data/pflotran_data.csv`, but the
gradient and flux stages were only asserted to be finite and correctly named.
That leaves the arithmetic unpinned: a refactor can move or rewrite the compute
path and change every flux value without failing a test.

These pins close that gap. Each column is summarised by two order-sensitive
aggregates -- sum(|v|) and max(|v|) -- captured from the committed sample_data
at TEMPERATURE_C. Any change to the gradient stencil, the Fick conversion, the
Stokes-Einstein correction, the diffusion table, or the M/m -> SI factor moves
at least one of them.

Golden values captured 2026-07-23 on `main` @ 229b2a2 (pre PR-4 layout split).
They are characterisation values, not independently derived ground truth: they
record what the code *does*, so a deliberate physics change is expected to fail
here and should be re-pinned in the same commit that makes the change.

See issue #9 (regression tests against frozen sample_data output).
"""

import os

import numpy as np
import pytest

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")

from pflotran_py.visualization import extract  # noqa: E402
from pflotran_py.visualization import physics  # noqa: E402

SAMPLE_DATA_DIR = os.path.join(REPO_ROOT, "sample_data")
FILE_TEMPLATE = "test29-{:03d}.tec"
N_FILES = 6
TEMPERATURE_C = 8.0
SPECIES_MAP = {
    "CO2": "CO2(aq) [M]",
    "CH4": "Free CH4(aq) [M]",
}

# column -> (sum(|v|), max(|v|))
GOLDEN = {
    "CO2_grad_x": (2.982049204492e-04, 3.579489932886e-06),
    "CO2_grad_y": (2.982049204492e-04, 3.579489932886e-06),
    "CO2_grad_z": (7.429678111297e-03, 3.279311698113e-05),
    "CO2_grad_magnitude": (7.458127907631e-03, 3.279332484505e-05),
    "CO2_flux_x": (3.478076851758e-10, 4.174894585212e-12),
    "CO2_flux_y": (3.478076851758e-10, 4.174894585212e-12),
    "CO2_flux_z": (8.665514779563e-09, 3.824785348855e-11),
    "CO2_flux_magnitude": (8.698696853796e-09, 3.824809592811e-11),
    "CH4_grad_x": (2.131750485413e-19, 1.267114093957e-20),
    "CH4_grad_y": (2.131750485413e-19, 1.267114093957e-20),
    "CH4_grad_z": (6.564053653774e-18, 1.300830188680e-19),
    "CH4_grad_magnitude": (6.580571701363e-18, 1.300834533462e-19),
    "CH4_flux_x": (1.939606543281e-25, 1.152903590049e-26),
    "CH4_flux_y": (1.939606543281e-25, 1.152903590049e-26),
    "CH4_flux_z": (5.972406951199e-24, 1.183580706525e-25),
    "CH4_flux_magnitude": (5.987436155323e-24, 1.183584659693e-25),
}

# Stokes-Einstein correction factors, pinned independently of the DataFrame path.
GOLDEN_STOKES_EINSTEIN = {
    8.0: 0.610648087951,
    18.0: 0.826262134310,
    25.0: 1.000000000000,
}


@pytest.fixture(scope="module")
def flux_df():
    """Full compute chain: .tec -> DataFrame -> gradients -> fluxes."""
    df = extract.extract_pflotran_data_tec(
        data_dir=SAMPLE_DATA_DIR,
        file_name_template=FILE_TEMPLATE,
        n_files=N_FILES,
    )
    df = physics.calculate_gradients(df, SPECIES_MAP)
    return physics.convert_to_flux(df, list(SPECIES_MAP), temperature_c=TEMPERATURE_C)


@pytest.mark.parametrize("column", sorted(GOLDEN))
def test_gradient_and_flux_values_are_pinned(flux_df, column):
    expected_sum, expected_max = GOLDEN[column]
    values = np.abs(flux_df[column].values)

    np.testing.assert_allclose(
        values.sum(),
        expected_sum,
        rtol=1e-10,
        err_msg=f"sum(|{column}|) changed from the pinned value",
    )
    np.testing.assert_allclose(
        values.max(),
        expected_max,
        rtol=1e-10,
        err_msg=f"max(|{column}|) changed from the pinned value",
    )


@pytest.mark.parametrize("temperature_c", sorted(GOLDEN_STOKES_EINSTEIN))
def test_stokes_einstein_factors_are_pinned(temperature_c):
    np.testing.assert_allclose(
        physics.stokes_einstein_correction(temperature_c),
        GOLDEN_STOKES_EINSTEIN[temperature_c],
        rtol=1e-10,
    )
