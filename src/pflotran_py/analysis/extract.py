"""Extract PFLOTRAN simulation output from Tecplot (.tec) files.

Auto-detects Tecplot vs HDF5 and dispatches to the matching reader. The HDF5
helpers are re-exported so callers can import either format from here.
"""

import glob
import logging
import os
import re
import warnings

import pandas as pd

from . import extract_hdf5
from .columns import TIME_COL, time_to_days

logger = logging.getLogger(__name__)

# Re-export HDF5 helpers so callers can import either format from `extract`.
extract_pflotran_data_hdf5 = extract_hdf5.extract_pflotran_data_hdf5
find_hdf5_output = extract_hdf5.find_hdf5_output


def extract_plfotran_tecplot_variable_names(filepath):
    """Return the VARIABLES list from a Tecplot header (one row expected)."""
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("VARIABLES"):
                return [v.strip().strip('"') for v in line.split("=")[1].split(",")]
    return []


def parse_tecplot_time_days(filepath):
    """Parse simulation time in days from a Tecplot snapshot header.

    Prefers ``TITLE = "  1.00000E+00 [d]"`` (value + unit). Falls back to
    ``SOLUTIONTIME=...`` treated as days if no unit is present.
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("TITLE"):
            match = re.search(
                r"([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?\d+)?)\s*\[([A-Za-z]+)\]",
                line,
            )
            if match:
                return time_to_days(match.group(1), match.group(2))

    for line in lines:
        if "SOLUTIONTIME" in line:
            match = re.search(r"SOLUTIONTIME\s*=\s*([0-9.eE+-]+)", line)
            if match:
                warnings.warn(
                    f"{filepath}: no time unit in TITLE; "
                    f"treating SOLUTIONTIME={match.group(1)} as days",
                    UserWarning,
                    stacklevel=2,
                )
                return float(match.group(1))

    raise ValueError(f"Could not parse simulation time from {filepath}")


def read_tec_file(filepath):
    """Read one Tecplot snapshot into a DataFrame using its VARIABLES header."""
    with open(filepath, "r") as f:
        lines = f.readlines()

    var_line = next(ln for ln in lines if ln.startswith("VARIABLES"))
    variables = [v.strip().strip('"') for v in var_line.split("=")[1].split(",")]

    data_start_idx = next(i for i, ln in enumerate(lines) if ln.startswith("ZONE"))
    data_lines = lines[data_start_idx + 1 :]
    data = [list(map(float, ln.strip().split())) for ln in data_lines]

    return pd.DataFrame(data, columns=variables)


def extract_pflotran_data_tec(
    data_dir=".", file_name_template="test29-{:03d}.tec", n_files=100, verbose=False
):
    """Read n_files Tecplot snapshots into a single DataFrame.

    Adds a 'Time Index' ordinal and the ``Time [d]`` simulation-time column.
    """
    if verbose:
        sample_filepath = os.path.join(data_dir, file_name_template.format(0))
        variable_names = extract_plfotran_tecplot_variable_names(sample_filepath)
        logger.debug("Variables found in tecplot files:")
        for i, var in enumerate(variable_names):
            logger.debug("  %2d: %s", i, var)

    all_data = []
    for i in range(n_files):
        filepath = os.path.join(data_dir, file_name_template.format(i))
        df = read_tec_file(filepath)
        df["Time Index"] = i
        df[TIME_COL] = parse_tecplot_time_days(filepath)
        all_data.append(df)

    return pd.concat(all_data, ignore_index=True)


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
