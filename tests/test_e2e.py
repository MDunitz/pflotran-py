"""End-to-end test: .tec extraction -> gradient -> flux computation.

Runs the full visualization pipeline (steps 1, 3) against sample_data/
and verifies the outputs are structurally correct and physically
reasonable. Does NOT require a PFLOTRAN binary -- uses the committed
.tec output files as input.
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

# Repo root so imports resolve
REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "visualization"))

import step1_extract  # noqa: E402
import shared_utils  # noqa: E402

# ── Fixtures ──────────────────────────────────────────────────────

SAMPLE_DATA_DIR = os.path.join(REPO_ROOT, "sample_data")
FILE_TEMPLATE = "test29-{:03d}.tec"
N_FILES = 6
TEMPERATURE_C = 8.0
SPECIES_MAP = {
    "CO2": "CO2(aq) [M]",
    "CH4": "Free CH4(aq) [M]",
}


@pytest.fixture(scope="module")
def extracted_df():
    """Extract all test29 .tec files into a single DataFrame."""
    return step1_extract.extract_pflotran_data_tec(
        data_dir=SAMPLE_DATA_DIR,
        file_name_template=FILE_TEMPLATE,
        n_files=N_FILES,
    )


@pytest.fixture(scope="module")
def gradient_df(extracted_df):
    """Compute concentration gradients for CO2 and CH4."""
    return shared_utils.calculate_gradients(extracted_df.copy(), SPECIES_MAP)


@pytest.fixture(scope="module")
def flux_df(gradient_df):
    """Convert gradients to diffusive fluxes with Stokes-Einstein correction."""
    species_list = list(SPECIES_MAP.keys())
    return shared_utils.convert_to_flux(
        gradient_df.copy(), species_list, temperature_c=TEMPERATURE_C
    )


# ── Step 1: Extraction ───────────────────────────────────────────


def test_extraction_produces_dataframe(extracted_df):
    assert isinstance(extracted_df, pd.DataFrame)
    assert len(extracted_df) > 0


def test_extraction_has_expected_shape(extracted_df):
    # 4x4x4 grid = 64 cells, 6 timesteps = 384 rows
    assert len(extracted_df) == 64 * N_FILES


def test_extraction_has_time_index(extracted_df):
    assert "Time Index" in extracted_df.columns
    assert set(extracted_df["Time Index"].unique()) == set(range(N_FILES))


def test_extraction_has_spatial_columns(extracted_df):
    for col in ["X [m]", "Y [m]", "Z [m]"]:
        assert col in extracted_df.columns


def test_extraction_has_species_columns(extracted_df):
    for species_col in SPECIES_MAP.values():
        assert species_col in extracted_df.columns


def test_extraction_matches_reference_csv(extracted_df):
    """Regression: extracted DataFrame matches committed pflotran_data.csv."""
    ref_path = os.path.join(SAMPLE_DATA_DIR, "pflotran_data.csv")
    ref_df = pd.read_csv(ref_path)

    # Same shape
    assert (
        extracted_df.shape == ref_df.shape
    ), f"Shape mismatch: extracted {extracted_df.shape} vs reference {ref_df.shape}"

    # Same columns (order may differ)
    assert set(extracted_df.columns) == set(ref_df.columns)

    # Numerical values match within floating-point tolerance
    for col in ref_df.columns:
        if ref_df[col].dtype in [np.float64, np.float32]:
            np.testing.assert_allclose(
                extracted_df[col].values,
                ref_df[col].values,
                rtol=1e-6,
                err_msg=f"Column {col} values differ from reference",
            )


# ── Step 3: Gradients ────────────────────────────────────────────


def test_gradient_columns_created(gradient_df):
    for species in SPECIES_MAP:
        for component in ["x", "y", "z"]:
            col = shared_utils.gradient_col(species, component)
            assert col in gradient_df.columns, f"Missing gradient column: {col}"


def test_gradients_are_finite(gradient_df):
    for species in SPECIES_MAP:
        for component in ["x", "y", "z"]:
            col = shared_utils.gradient_col(species, component)
            vals = gradient_df[col].dropna()
            assert np.all(np.isfinite(vals)), f"Non-finite values in {col}"


def test_initial_gradients_near_zero(gradient_df):
    """At t=0, uniform initial conditions should give near-zero gradients."""
    t0 = gradient_df[gradient_df["Time Index"] == 0]
    for species in SPECIES_MAP:
        for component in ["x", "y", "z"]:
            col = shared_utils.gradient_col(species, component)
            vals = t0[col].dropna()
            if len(vals) > 0:
                assert np.max(np.abs(vals)) < 1e-3, (
                    f"{col} at t=0 has large gradient ({np.max(np.abs(vals)):.2e}), "
                    "expected near-zero for uniform initial conditions"
                )


# ── Step 3→4: Fluxes ─────────────────────────────────────────────


def test_flux_columns_created(flux_df):
    for species in SPECIES_MAP:
        for component in ["x", "y", "z"]:
            col = shared_utils.flux_col(species, component)
            assert col in flux_df.columns, f"Missing flux column: {col}"


def test_fluxes_are_finite(flux_df):
    for species in SPECIES_MAP:
        for component in ["x", "y", "z"]:
            col = shared_utils.flux_col(species, component)
            vals = flux_df[col].dropna()
            assert np.all(np.isfinite(vals)), f"Non-finite values in {col}"


def test_stokes_einstein_at_reference_is_unity():
    """Stokes-Einstein correction at reference temperature should be 1.0.

    Stokes-Einstein equation:
        D(T) / D(T_ref) = (mu(T_ref) / mu(T)) * (T / T_ref)

    At T = T_ref, this must equal 1.0 by definition.
    """
    correction = shared_utils.stokes_einstein_correction(25.0, reference_c=25.0)
    assert abs(correction - 1.0) < 1e-10


def test_stokes_einstein_at_8c_less_than_unity():
    """At 8C (below 25C reference), diffusivity should be reduced.

    Stokes-Einstein equation:
        D(T) / D(T_ref) = (mu(T_ref) / mu(T)) * (T / T_ref)

    Water viscosity increases at lower temperatures, so D(8C) < D(25C)
    and the correction factor should be < 1.
    """
    correction = shared_utils.stokes_einstein_correction(8.0, reference_c=25.0)
    assert 0 < correction < 1.0, (
        f"Stokes-Einstein correction at 8C = {correction:.4f}, "
        "expected < 1.0 (higher viscosity reduces diffusivity)"
    )
