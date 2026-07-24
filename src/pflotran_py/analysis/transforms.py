"""Spatial DataFrame transforms shared by the 2D presentation modules.

Pure compute (no plotting): surface-cell identification and single-point
time-series extraction, pulled out of the former step3/step4 plot files.
"""

import numpy as np

from .columns import time_axis_column


def identify_surface_cells(df):
    """Mark the top cell of each (X, Y) column as a surface cell.

    Adds a boolean 'is_surface' column: True where Z is maximal for that
    (X, Y) at that time step.
    """
    surface_df = df.copy()
    surface_df["is_surface"] = False

    for time_idx in df["Time Index"].unique():
        time_data = df[df["Time Index"] == time_idx]

        max_z_by_xy = time_data.groupby(["X [m]", "Y [m]"])["Z [m]"].max().reset_index()
        max_z_by_xy.columns = ["X [m]", "Y [m]", "max_z"]

        merged = time_data.merge(max_z_by_xy, on=["X [m]", "Y [m]"])
        surface_indices = merged[merged["Z [m]"] == merged["max_z"]].index
        surface_df.loc[surface_indices, "is_surface"] = True

    return surface_df


def find_target_point(df):
    """Return the top-left surface cell coordinates (min X, min Y, max Z)."""
    max_z = df["Z [m]"].max()
    surface = df[df["Z [m]"] == max_z]
    min_x = surface["X [m]"].min()
    min_y = surface["Y [m]"].min()
    target = surface[(surface["X [m]"] == min_x) & (surface["Y [m]"] == min_y)]
    return target.iloc[0][["X [m]", "Y [m]", "Z [m]"]].values


def extract_point_time_series(df, target_x, target_y, target_z, tolerance=1e-6):
    """Return all timesteps for a single grid point, ordered by time."""
    point_data = df[
        (np.abs(df["X [m]"] - target_x) < tolerance)
        & (np.abs(df["Y [m]"] - target_y) < tolerance)
        & (np.abs(df["Z [m]"] - target_z) < tolerance)
    ].copy()
    sort_col = time_axis_column(point_data)
    return point_data.sort_values(sort_col)
