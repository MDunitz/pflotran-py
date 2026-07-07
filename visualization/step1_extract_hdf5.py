"""Extract PFLOTRAN simulation output from HDF5 (.h5) files."""

import glob
import os

import h5py
import numpy as np
import pandas as pd

# PFLOTRAN writes a single .h5 file containing every snapshot as a
# "Time:  <value> <unit>" group. Each group holds (nx, ny, nz) arrays named
# like "Free_CH4(aq) [M]" / "Total_CO2(aq) [M]". The Coordinates group stores
# grid *edges* (length n+1 per axis); cell centers are the edge midpoints.
#
# To stay drop-in compatible with the Tecplot extractor, variable names are
# normalized to the Tecplot spelling: the "Free_" / "Total_" prefix underscore
# becomes a space (e.g. "Free_CH4(aq) [M]" -> "Free CH4(aq) [M]").


def _h5_cell_centers(edges):
    """Convert grid edge coordinates (length n+1) to cell centers (length n)."""
    edges = np.asarray(edges, dtype=float)
    if len(edges) >= 2:
        return 0.5 * (edges[:-1] + edges[1:])
    return edges


def _normalize_h5_var_name(name):
    """Map an HDF5 dataset name to the Tecplot column spelling."""
    if name.startswith("Free_"):
        return "Free " + name[len("Free_") :]
    if name.startswith("Total_"):
        return "Total " + name[len("Total_") :]
    return name


def _h5_time_groups(h5file):
    """Return time-group keys sorted by their numeric time value."""

    def time_value(key):
        # key looks like "Time:  1.00000E+00 d"
        try:
            return float(key.split(":", 1)[1].split()[0])
        except (IndexError, ValueError):
            return float("inf")

    groups = [k for k in h5file.keys() if k.startswith("Time:")]
    return sorted(groups, key=time_value)


def extract_pflotran_data_hdf5(filepath, verbose=False):
    """Extract every snapshot from a PFLOTRAN HDF5 output file.

    Returns a DataFrame with X/Y/Z [m] cell-center columns, a "Time Index"
    column (0-based, ordered by simulation time), and one column per output
    variable using Tecplot-compatible names.
    """
    all_data = []
    with h5py.File(filepath, "r") as h5:
        coords = h5["Coordinates"]
        xc = _h5_cell_centers(coords["X [m]"][:])
        yc = _h5_cell_centers(coords["Y [m]"][:])
        zc = _h5_cell_centers(coords["Z [m]"][:])

        # Cell-center grid, indexing='ij' so [i,j,k] matches the (nx,ny,nz) data
        xx, yy, zz = np.meshgrid(xc, yc, zc, indexing="ij")

        time_groups = _h5_time_groups(h5)
        if verbose:
            print(f"\n HDF5 file: {filepath}")
            print(f"   grid: {len(xc)} x {len(yc)} x {len(zc)}")
            print(f"   snapshots: {len(time_groups)}")

        for time_idx, gname in enumerate(time_groups):
            group = h5[gname]
            frame = {
                "X [m]": xx.ravel(),
                "Y [m]": yy.ravel(),
                "Z [m]": zz.ravel(),
            }
            for var_name in group.keys():
                data = group[var_name][:]
                if getattr(data, "shape", None) != xx.shape:
                    # Skip scalars / non-grid datasets
                    continue
                frame[_normalize_h5_var_name(var_name)] = data.ravel()

            df = pd.DataFrame(frame)
            df["Time Index"] = time_idx
            all_data.append(df)

    if not all_data:
        raise ValueError(f"No time-series data found in {filepath}")

    return pd.concat(all_data, ignore_index=True)


def find_hdf5_output(data_dir, prefix=None):
    """Return the path to a PFLOTRAN snapshot .h5 file, or None.

    Ignores observation/restart files (``*-obs*.h5``, ``*restart*.h5``).
    """
    candidates = sorted(glob.glob(os.path.join(data_dir, "*.h5")))
    candidates = [
        c
        for c in candidates
        if "-obs" not in os.path.basename(c) and "restart" not in os.path.basename(c)
    ]
    if prefix:
        preferred = [c for c in candidates if os.path.basename(c).startswith(prefix)]
        if preferred:
            return preferred[0]
    return candidates[0] if candidates else None
