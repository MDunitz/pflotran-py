#!/usr/bin/env python3
# PFLOTRAN Data Extraction Script

# Should install + import dependencies
# Will show you what variables are available in the tecplot files
# Save the combined data for visualization in step2_plot.py, as a csv

# ##################################################################
#  1) Install and import dependencies
# ##################################################################
import os
import glob

import h5py
import numpy as np
import pandas as pd


# ##################################################################
#  2) Get variable names from a sample tecplot file
# ##################################################################
# get the variable names from a sample tecplot file (there will only be one row)
# useful if you're working with a new PFLOTRAN output file
def extract_plfotran_tecplot_variable_names(filepath):
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("VARIABLES"):
                # Clean and split variable names
                variables = [
                    v.strip().strip('"') for v in line.split("=")[1].split(",")
                ]
                return variables
    return []


# ##################################################################
#  3) Read tecplot files, return a dataframe to give structure
# ##################################################################


def read_tec_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Extract variables
    var_line = next(ln for ln in lines if ln.startswith("VARIABLES"))
    variables = [v.strip().strip('"') for v in var_line.split("=")[1].split(",")]

    # Extract data lines (after the ZONE line)
    data_start_idx = next(i for i, ln in enumerate(lines) if ln.startswith("ZONE"))
    data_lines = lines[data_start_idx + 1 :]
    data = [list(map(float, ln.strip().split())) for ln in data_lines]

    df = pd.DataFrame(data, columns=variables)
    return df


# ##################################################################
#  4) Extract pflotran data from all tecplot files
# ##################################################################
# Will also call functions 2 and 3 (note: will only call 2 if you're doing troubleshooting, and you need to comment out stuff below)
def extract_pflotran_data_tec(
    data_dir=".", file_name_template="test29-{:03d}.tec", n_files=100, verbose=False
):

    # sample file exists to get variable names
    sample_filename = file_name_template.format(0)
    sample_filepath = os.path.join(data_dir, sample_filename)

    if verbose:
        # If you want to see the variable names, uncomment this block, and add it as something you return
        # Get variable names from sample file
        variable_names = extract_plfotran_tecplot_variable_names(sample_filepath)
        print("\n Variables found in tecplot files:")
        for i, var in enumerate(variable_names):
            print(f"   {i:2d}: {var}")

    # Storage for  data
    all_data = []

    # Read all files
    files_read = 0
    for i in range(n_files):
        filename = file_name_template.format(i)
        filepath = os.path.join(data_dir, filename)

        df = read_tec_file(filepath)
        df["Time Index"] = i
        all_data.append(df)
        files_read += 1

    # Combine all dataframes
    full_df = pd.concat(all_data, ignore_index=True)

    return full_df


# ##################################################################
#  4b) Extract pflotran data from an HDF5 file
# ##################################################################
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


# ##################################################################
#  4c) Auto-detecting extractor (TEC or HDF5)
# ##################################################################
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


def extract_pflotran_data(
    data_dir=".",
    file_name_template="{prefix}-{:03d}.tec",
    prefix="sim",
    data_format="auto",
    verbose=False,
):
    """Extract PFLOTRAN output, auto-selecting the Tecplot or HDF5 reader.

    Parameters
    ----------
    data_dir : str
        Directory containing PFLOTRAN output.
    file_name_template : str
        Template for Tecplot snapshot files. If it contains ``{prefix}`` it is
        formatted with ``prefix`` first.
    prefix : str
        Simulation prefix (``-input_prefix``), used to locate files.
    data_format : {"auto", "tec", "hdf5"}
        Output format. "auto" detects by looking for files on disk.
    """
    tec_template = file_name_template
    if "{prefix}" in tec_template:
        tec_template = tec_template.replace("{prefix}", prefix)

    resolved = data_format
    if resolved == "auto":
        tec_files = sorted(
            glob.glob(os.path.join(data_dir, f"{prefix}-[0-9][0-9][0-9].tec"))
        )
        if tec_files:
            resolved = "tec"
        elif find_hdf5_output(data_dir, prefix=prefix) is not None:
            resolved = "hdf5"
        else:
            raise FileNotFoundError(
                f"No PFLOTRAN .tec or .h5 output found in {data_dir}"
            )

    if resolved == "tec":
        n_files = len(
            glob.glob(os.path.join(data_dir, f"{prefix}-[0-9][0-9][0-9].tec"))
        )
        return extract_pflotran_data_tec(
            data_dir=data_dir,
            file_name_template=tec_template,
            n_files=n_files,
            verbose=verbose,
        )
    elif resolved == "hdf5":
        h5_path = find_hdf5_output(data_dir, prefix=prefix)
        if h5_path is None:
            raise FileNotFoundError(f"No PFLOTRAN .h5 output found in {data_dir}")
        return extract_pflotran_data_hdf5(h5_path, verbose=verbose)
    else:
        raise ValueError(f"Unknown data_format: {data_format!r}")


# ##################################################################
#  5) Save the combined dataframe
# ##################################################################
def save_data(df, filename="pflotran_data.pkl"):
    df.to_pickle(filename)

    # Also save as CSV for convenience (though it will be larger)
    csv_filename = filename.replace(".pkl", ".csv")
    df.to_csv(csv_filename, index=False)


# ##################################################################
#  0) Main
# ##################################################################
def main(data_format, data_dir, file_name_template, n_files):
    if data_format == "tec":
        print("\nExtracting data from tecplot files...")
        full_df = extract_pflotran_data_tec(
            data_dir=data_dir, file_name_template=file_name_template, n_files=n_files
        )
    elif data_format == "hdf5":
        print("\nExtracting data from hdf5 file...")
        h5_path = find_hdf5_output(data_dir)
        if h5_path is None:
            raise FileNotFoundError(f"No PFLOTRAN .h5 output found in {data_dir}")
        full_df = extract_pflotran_data_hdf5(h5_path)
    else:
        raise ValueError(f"Unknown data_format: {data_format!r}")

    print("\nSaving data...")
    save_data(full_df, "pflotran_data.pkl")


if __name__ == "__main__":
    main()
