"""Unit tests for PFLOTRAN HDF5 extraction (no PFLOTRAN binary required)."""

import os
import sys

import h5py
import numpy as np
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "visualization"))

import step1_extract_hdf5 as hdf5_extract  # noqa: E402


@pytest.fixture
def sample_h5(tmp_path):
    """Write a minimal PFLOTRAN-like HDF5 file (2x2x2 grid, 2 snapshots)."""
    path = tmp_path / "sim.h5"
    nx, ny, nz = 2, 2, 2

    with h5py.File(path, "w") as h5:
        coords = h5.create_group("Coordinates")
        coords.create_dataset("X [m]", data=np.linspace(0.0, 2.0, nx + 1))
        coords.create_dataset("Y [m]", data=np.linspace(0.0, 2.0, ny + 1))
        coords.create_dataset("Z [m]", data=np.linspace(0.0, 2.0, nz + 1))

        for time_idx, gname in enumerate(
            ["Time:  0.00000E+00 d", "Time:  1.00000E+00 d"]
        ):
            group = h5.create_group(gname)
            group.create_dataset(
                "Free_CO2(aq) [M]",
                data=np.full((nx, ny, nz), float(time_idx + 1)),
            )
            group.create_dataset(
                "Free_CH4(aq) [M]",
                data=np.full((nx, ny, nz), 0.1 * (time_idx + 1)),
            )

    return path


def test_find_hdf5_output_prefers_prefix(tmp_path):
    (tmp_path / "sim.h5").write_bytes(b"")
    (tmp_path / "sim-obs-0.h5").write_bytes(b"")
    (tmp_path / "other.h5").write_bytes(b"")

    found = hdf5_extract.find_hdf5_output(str(tmp_path), prefix="sim")
    assert os.path.basename(found) == "sim.h5"


def test_find_hdf5_output_ignores_observation_files(tmp_path):
    (tmp_path / "sim-obs-0.h5").write_bytes(b"")

    assert hdf5_extract.find_hdf5_output(str(tmp_path), prefix="sim") is None


def test_extract_hdf5_grid_shape(sample_h5):
    df = hdf5_extract.extract_pflotran_data_hdf5(str(sample_h5))

    # 2x2x2 cells, 2 snapshots -> 16 rows
    assert len(df) == 16
    assert df["Time Index"].nunique() == 2


def test_extract_hdf5_spatial_columns(sample_h5):
    df = hdf5_extract.extract_pflotran_data_hdf5(str(sample_h5))

    for col in ["X [m]", "Y [m]", "Z [m]"]:
        assert col in df.columns
    # Cell centers from edges [0,1,2] -> [0.5, 1.5]
    assert set(df["X [m]"].unique()) == {0.5, 1.5}


def test_extract_hdf5_normalizes_variable_names(sample_h5):
    df = hdf5_extract.extract_pflotran_data_hdf5(str(sample_h5))

    assert "Free CO2(aq) [M]" in df.columns
    assert "Free CH4(aq) [M]" in df.columns
    assert "Free_CO2(aq) [M]" not in df.columns


def test_extract_hdf5_values_by_snapshot(sample_h5):
    df = hdf5_extract.extract_pflotran_data_hdf5(str(sample_h5))

    t0 = df.loc[df["Time Index"] == 0, "Free CO2(aq) [M]"].iloc[0]
    t1 = df.loc[df["Time Index"] == 1, "Free CO2(aq) [M]"].iloc[0]
    assert t0 == pytest.approx(1.0)
    assert t1 == pytest.approx(2.0)


def test_extract_hdf5_no_snapshots_raises(tmp_path):
    path = tmp_path / "no_times.h5"
    with h5py.File(path, "w") as h5:
        coords = h5.create_group("Coordinates")
        coords.create_dataset("X [m]", data=np.array([0.0, 1.0]))
        coords.create_dataset("Y [m]", data=np.array([0.0, 1.0]))
        coords.create_dataset("Z [m]", data=np.array([0.0, 1.0]))

    with pytest.raises(ValueError, match="No time-series data"):
        hdf5_extract.extract_pflotran_data_hdf5(str(path))


def test_extract_hdf5_warns_on_shape_mismatch(tmp_path):
    """Non-grid datasets are skipped with a UserWarning (not silently)."""
    path = tmp_path / "mixed.h5"
    nx, ny, nz = 2, 2, 2

    with h5py.File(path, "w") as h5:
        coords = h5.create_group("Coordinates")
        coords.create_dataset("X [m]", data=np.linspace(0.0, 2.0, nx + 1))
        coords.create_dataset("Y [m]", data=np.linspace(0.0, 2.0, ny + 1))
        coords.create_dataset("Z [m]", data=np.linspace(0.0, 2.0, nz + 1))

        for gname in ["Time:  0.00000E+00 d", "Time:  1.00000E+00 d"]:
            group = h5.create_group(gname)
            group.create_dataset(
                "Free_CO2(aq) [M]",
                data=np.ones((nx, ny, nz)),
            )
            # Non-grid dataset that should be skipped with a warning
            group.create_dataset("Some_Scalar", data=np.array([42.0]))

    with pytest.warns(UserWarning, match="Skipping 'Some_Scalar'.*does not match"):
        df = hdf5_extract.extract_pflotran_data_hdf5(str(path))

    assert "Free CO2(aq) [M]" in df.columns
    assert "Some_Scalar" not in df.columns
    # One warning per variable name, not one per snapshot
    assert len(df) == nx * ny * nz * 2
