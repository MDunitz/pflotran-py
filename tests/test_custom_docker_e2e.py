"""End-to-end test for decks with custom AWINHIBIT reaction sandboxes.

Requires a PFLOTRAN binary compiled with the sandbox modules from sandbox/
(see scripts/build_pflotran_custom.sh). Uses the repo's sandbox/hanford.dat.

Default deck: exploratory/pflotran/testing/7_sandbox_try.in
  - 4x4x4 grid, 5 days, AWINHIBIT + AWINHIBITACETATE + AWINHIBITMETHYL
  - ACTIVITY_WATER inhibition with repo hanford.dat

Build the custom binary (once) inside the container:

    ./scripts/build_pflotran_custom.sh
    # -> build/pflotran

Run standalone:

    python3 tests/test_custom_docker_e2e.py
    python3 tests/test_custom_docker_e2e.py exploratory/pflotran/testing/10_addsulfateinhibit.in
"""

import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

import pflotran_e2e_common as e2e  # noqa: E402
from pflotran_py.visualization import shared_utils  # noqa: E402

# 4x4x4, 5 days — first deck with active (non-skipped) AWINHIBIT sandboxes.
# Emits Tecplot output.
DEFAULT_INPUT = os.path.join(
    REPO_ROOT, "exploratory", "pflotran", "testing", "7_sandbox_try.in"
)
# 4x4x4, 5 days — AWINHIBIT sandboxes with FORMAT HDF5 output. Exercises the
# HDF5 extraction path (vs the Tecplot path used by DEFAULT_INPUT).
HDF5_INPUT = os.path.join(REPO_ROOT, "batch", "9_addnitrogen.in")
OUTPUT_DIR = os.path.join(REPO_ROOT, "tests", "custom_e2e_output")
HDF5_OUTPUT_DIR = os.path.join(REPO_ROOT, "tests", "custom_e2e_hdf5_output")
SPECIES_MAP = e2e.SPECIES_MAP

pytestmark = pytest.mark.integration


def _require_custom_pflotran():
    exe = e2e.find_custom_pflotran()
    if exe is None:
        pytest.skip(
            "Custom PFLOTRAN binary not found. Run ./scripts/build_pflotran_custom.sh "
            "or set PFLOTRAN_CUSTOM_EXE."
        )
    return exe


@pytest.fixture(scope="module")
def pipeline():
    exe = _require_custom_pflotran()
    return e2e.run_full_pipeline(
        src_in=DEFAULT_INPUT,
        output_dir=OUTPUT_DIR,
        exe=exe,
        custom_sandboxes=True,
    )


def test_custom_pflotran_produced_snapshots(pipeline):
    expected = e2e.parse_expected_snapshots(DEFAULT_INPUT)
    if expected is None:
        assert pipeline["n_snapshots"] > 0
    else:
        assert pipeline["n_snapshots"] == expected


def test_extracted_grid_shape(pipeline):
    df = pipeline["dataframe"]
    cells = e2e.parse_grid_cells(DEFAULT_INPUT)
    assert cells is not None, "Could not parse NXYZ from deck"
    assert len(df) == cells * pipeline["n_snapshots"]


def test_species_columns_present(pipeline):
    df = pipeline["dataframe"]
    for col in pipeline["species_map"].values():
        assert col in df.columns


def test_concentrations_physical(pipeline):
    df = pipeline["dataframe"]
    for col in pipeline["species_map"].values():
        vals = df[col].values
        assert np.all(np.isfinite(vals)), f"Non-finite values in {col}"
        assert np.all(vals >= 0), f"Negative concentration in {col}"


def test_flux_columns_finite(pipeline):
    df = pipeline["dataframe"]
    for species in pipeline["species_map"]:
        col = shared_utils.flux_col(species, "magnitude")
        assert col in df.columns
        assert np.all(np.isfinite(df[col].dropna().values))


def test_salinity_chemistry_evolved(pipeline):
    """Deck 7 starts at high Cl-; reactive transport should change salinity."""
    df = pipeline["dataframe"]
    cl = "Free Cl- [M]"
    assert cl in df.columns
    t0 = df.loc[df["Time Index"] == 0, cl].mean()
    tf = df.loc[df["Time Index"] == df["Time Index"].max(), cl].mean()
    assert tf < t0 / 100, f"Cl- did not decrease (t0={t0:.3e}, tf={tf:.3e})"


def test_sulfate_was_consumed(pipeline):
    """Sulfate reduction should draw down SO4-- from its initial value."""
    df = pipeline["dataframe"]
    assert "Total SO4-- [M]" in df.columns
    so4 = "Total SO4-- [M]"
    t0 = df.loc[df["Time Index"] == 0, so4].mean()
    tf = df.loc[df["Time Index"] == df["Time Index"].max(), so4].mean()
    assert tf < t0 / 10, f"SO4-- was not consumed (t0={t0:.3e}, tf={tf:.3e})"


def test_images_written(pipeline):
    for name, path in pipeline["images"].items():
        assert os.path.isfile(path), f"Missing image {name}: {path}"
        assert os.path.getsize(path) > 0, f"Empty image {name}: {path}"


def test_output_format_is_tecplot(pipeline):
    """DEFAULT_INPUT (deck 7) uses FORMAT TECPLOT."""
    assert pipeline["output_format"] == "tec"


# ── HDF5 output path ──────────────────────────────────────────────
# batch/9_addnitrogen.in uses FORMAT HDF5, so this exercises the HDF5 reader
# end-to-end through the exact same pipeline (extract -> gradient -> flux ->
# images) as the Tecplot decks.


@pytest.fixture(scope="module")
def hdf5_pipeline():
    exe = _require_custom_pflotran()
    return e2e.run_full_pipeline(
        src_in=HDF5_INPUT,
        output_dir=HDF5_OUTPUT_DIR,
        exe=exe,
        custom_sandboxes=True,
    )


def test_hdf5_output_format_detected(hdf5_pipeline):
    assert hdf5_pipeline["output_format"] == "hdf5"


def test_hdf5_snapshots_match_deck(hdf5_pipeline):
    expected = e2e.parse_expected_snapshots(HDF5_INPUT)
    if expected is None:
        assert hdf5_pipeline["n_snapshots"] > 0
    else:
        assert hdf5_pipeline["n_snapshots"] == expected


def test_hdf5_grid_shape(hdf5_pipeline):
    df = hdf5_pipeline["dataframe"]
    cells = e2e.parse_grid_cells(HDF5_INPUT)
    assert cells is not None, "Could not parse NXYZ from deck"
    assert len(df) == cells * hdf5_pipeline["n_snapshots"]


def test_hdf5_has_spatial_and_species_columns(hdf5_pipeline):
    df = hdf5_pipeline["dataframe"]
    for col in ["X [m]", "Y [m]", "Z [m]"]:
        assert col in df.columns, f"Missing coordinate column {col}"
    for col in hdf5_pipeline["species_map"].values():
        assert col in df.columns, f"Missing species column {col}"


def test_hdf5_concentrations_physical(hdf5_pipeline):
    df = hdf5_pipeline["dataframe"]
    for col in hdf5_pipeline["species_map"].values():
        vals = df[col].values
        assert np.all(np.isfinite(vals)), f"Non-finite values in {col}"
        assert np.all(vals >= 0), f"Negative concentration in {col}"


def test_hdf5_images_written(hdf5_pipeline):
    for name, path in hdf5_pipeline["images"].items():
        assert os.path.isfile(path), f"Missing image {name}: {path}"
        assert os.path.getsize(path) > 0, f"Empty image {name}: {path}"


def main(argv):
    exe = e2e.find_custom_pflotran()
    if exe is None:
        print(
            "ERROR: custom PFLOTRAN not found.\n"
            "  Run: ./scripts/build_pflotran_custom.sh\n"
            "  Or set PFLOTRAN_CUSTOM_EXE to a binary with AWINHIBIT sandboxes."
        )
        return 1

    src_in = argv[1] if len(argv) > 1 else DEFAULT_INPUT
    print(f"Input deck : {src_in}")
    print(f"Database   : {e2e.HANFORD_DB}")
    result = e2e.run_full_pipeline(
        src_in=src_in,
        output_dir=OUTPUT_DIR,
        exe=exe,
        custom_sandboxes=True,
    )
    print(f"Executable : {result['executable']}")
    print(f"Format     : {result['output_format']}")
    print(f"Snapshots  : {result['n_snapshots']}")
    print(f"Rows       : {len(result['dataframe'])}")
    print(f"CSV        : {result['csv']}")
    print("Images     :")
    for name, path in result["images"].items():
        print(f"   - {name}: {path}")
    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
