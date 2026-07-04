"""Full-stack end-to-end test: .in file -> PFLOTRAN -> calculations -> images.

Runs on the stock PFLOTRAN binary (no custom AWINHIBIT sandboxes required).
For decks with custom reaction sandboxes, see test_custom_docker_e2e.py.

Standalone use (also produces images/CSV in tests/e2e_output/):

    python3 tests/test_docker_e2e.py
    python3 tests/test_docker_e2e.py exploratory/pflotran/testing/3_smaller_grid.in
"""

import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "visualization"))

import pflotran_e2e_common as e2e  # noqa: E402
import shared_utils  # noqa: E402

DEFAULT_INPUT = os.path.join(
    REPO_ROOT, "exploratory", "pflotran", "testing", "3_smaller_grid.in"
)
OUTPUT_DIR = os.path.join(REPO_ROOT, "tests", "e2e_output")
SPECIES_MAP = e2e.SPECIES_MAP

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def pipeline():
    exe = e2e.find_pflotran()
    if exe is None:
        pytest.skip("PFLOTRAN executable not available on this machine")
    return e2e.run_full_pipeline(
        src_in=DEFAULT_INPUT,
        output_dir=OUTPUT_DIR,
        exe=exe,
        custom_sandboxes=False,
    )


def test_pflotran_produced_snapshots(pipeline):
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


def test_co2_is_produced(pipeline):
    df = pipeline["dataframe"]
    co2 = pipeline["species_map"]["CO2"]
    t0 = df.loc[df["Time Index"] == 0, co2].mean()
    tf = df.loc[df["Time Index"] == df["Time Index"].max(), co2].mean()
    assert tf > 10 * t0, f"CO2 did not accumulate (t0={t0:.3e}, tf={tf:.3e})"


def test_flux_develops_over_time(pipeline):
    df = pipeline["dataframe"]
    mag = shared_utils.flux_col("CO2", "magnitude")
    tf_max = df.loc[df["Time Index"] == df["Time Index"].max(), mag].max()
    assert tf_max > 0, f"No diffusive flux developed (max |J|={tf_max:.3e})"


def test_images_written(pipeline):
    for name, path in pipeline["images"].items():
        assert os.path.isfile(path), f"Missing image {name}: {path}"
        assert os.path.getsize(path) > 0, f"Empty image {name}: {path}"


def main(argv):
    exe = e2e.find_pflotran()
    if exe is None:
        print(
            "ERROR: no PFLOTRAN executable found. Set PFLOTRAN_EXE or run "
            "inside the PFLOTRAN container."
        )
        return 1

    src_in = argv[1] if len(argv) > 1 else DEFAULT_INPUT
    print(f"Input deck : {src_in}")
    result = e2e.run_full_pipeline(
        src_in=src_in,
        output_dir=OUTPUT_DIR,
        exe=exe,
        custom_sandboxes=False,
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
