#!/usr/bin/env python3

# PFLOTRAN Flux Time Series Visualization Script
# Plots CO2 and CH4 flux time series for the top-left surface point.


import os
import sys
import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm


########################################################################
# 1) Install dependencies
########################################################################
def install_dependencies():
    packages = ["pandas", "numpy", "bokeh", "tqdm", "scipy"]

    for package in packages:
        try:
            if package == "bokeh":
                import bokeh
            elif package == "scipy":
                import scipy
            else:
                __import__(package)
            print(f"{package} is already installed")
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
            print(f"{package} installed successfully")


########################################################################
# 2) Load in the processed data (CSV) Note, you need step1_extract.py done
########################################################################
def load_data(filename="pflotran_data.pkl"):
    if os.path.exists(filename):
        print(f"Loading data from: {filename}")
        df = pd.read_pickle(filename)
        print(f"Data loaded successfully! Shape: {df.shape}")
        return df
    elif os.path.exists(filename.replace(".pkl", ".csv")):
        csv_filename = filename.replace(".pkl", ".csv")
        df = pd.read_csv(csv_filename)
        print(f"Data loaded successfully! Shape: {df.shape}")
        return df
    else:
        print(f"Data file not found: {filename}")
        print("Please run step1_extract.py first to process the data.")
        return None


########################################################################
# 3) Calculate concentration gradients and fluxes (for Co2 and CH4 conc)
########################################################################
def calculate_concentration_gradients(df):
    from scipy.spatial.distance import cdist

    # Create a copy of the dataframe
    flux_df = df.copy()

    # Species of interest for flux calculations
    co2_col = "CO2(aq) [M]"
    ch4_free_col = "Free CH4(aq) [M]"

    # Initialize gradient columns
    flux_df["CO2_flux_x"] = 0.0
    flux_df["CO2_flux_y"] = 0.0
    flux_df["CO2_flux_z"] = 0.0
    flux_df["CO2_flux_magnitude"] = 0.0
    flux_df["CH4_flux_x"] = 0.0
    flux_df["CH4_flux_y"] = 0.0
    flux_df["CH4_flux_z"] = 0.0
    flux_df["CH4_flux_magnitude"] = 0.0

    # Get unique time indices
    time_indices = flux_df["Time Index"].unique()

    print("Calculating concentration gradients and fluxes...")

    for time_idx in tqdm(time_indices, desc="Processing time steps"):
        # Get data for this time step
        time_data = flux_df[flux_df["Time Index"] == time_idx].copy()

        if len(time_data) < 8:  # Need enough points for gradient calculation
            continue

        # Create coordinate arrays
        coords = time_data[["X [m]", "Y [m]", "Z [m]"]].values

        # Calculate gradients for each point
        for idx, row in time_data.iterrows():
            point = row[["X [m]", "Y [m]", "Z [m]"]].values

            # Find neighboring points (within some distance threshold)
            distances = np.sqrt(np.sum((coords - point) ** 2, axis=1))
            neighbors_idx = np.where((distances > 0) & (distances <= 0.5))[
                0
            ]  # Adjust threshold as needed

            if len(neighbors_idx) >= 3:  # Need at least 3 neighbors for gradient
                neighbor_coords = coords[neighbors_idx]
                neighbor_co2 = time_data.iloc[neighbors_idx][co2_col].values
                neighbor_ch4 = time_data.iloc[neighbors_idx][ch4_free_col].values

                # Simple finite difference approximation
                try:
                    # Calculate gradients using least squares fit
                    A = np.column_stack(
                        [neighbor_coords - point, np.ones(len(neighbor_coords))]
                    )

                    # CO2 gradients
                    co2_gradients, _, _, _ = np.linalg.lstsq(
                        A, neighbor_co2 - row[co2_col], rcond=None
                    )
                    flux_df.loc[idx, "CO2_flux_x"] = (
                        -co2_gradients[0] if len(co2_gradients) > 0 else 0
                    )
                    flux_df.loc[idx, "CO2_flux_y"] = (
                        -co2_gradients[1] if len(co2_gradients) > 1 else 0
                    )
                    flux_df.loc[idx, "CO2_flux_z"] = (
                        -co2_gradients[2] if len(co2_gradients) > 2 else 0
                    )

                    # CH4 gradients
                    ch4_gradients, _, _, _ = np.linalg.lstsq(
                        A, neighbor_ch4 - row[ch4_free_col], rcond=None
                    )
                    flux_df.loc[idx, "CH4_flux_x"] = (
                        -ch4_gradients[0] if len(ch4_gradients) > 0 else 0
                    )
                    flux_df.loc[idx, "CH4_flux_y"] = (
                        -ch4_gradients[1] if len(ch4_gradients) > 1 else 0
                    )
                    flux_df.loc[idx, "CH4_flux_z"] = (
                        -ch4_gradients[2] if len(ch4_gradients) > 2 else 0
                    )

                except np.linalg.LinAlgError:
                    # If linear algebra fails, set fluxes to zero
                    pass

    # Calculate flux magnitudes
    flux_df["CO2_flux_magnitude"] = np.sqrt(
        flux_df["CO2_flux_x"] ** 2
        + flux_df["CO2_flux_y"] ** 2
        + flux_df["CO2_flux_z"] ** 2
    )
    flux_df["CH4_flux_magnitude"] = np.sqrt(
        flux_df["CH4_flux_x"] ** 2
        + flux_df["CH4_flux_y"] ** 2
        + flux_df["CH4_flux_z"] ** 2
    )

    return flux_df


########################################################################
# 4) Find the top-left surface point and extract its time series
########################################################################
# Note that you can really choose any point, this is just an example
# You can just change max_z to min_z for bottom-left, or max_x for top-right, etc.
def find_target_point(df):
    # Find the maximum Z value (surface)
    max_z = df["Z [m]"].max()

    # Get all points at the surface
    surface_points = df[df["Z [m]"] == max_z]

    # Find the top-left point (minimum X and minimum Y)
    min_x = surface_points["X [m]"].min()
    min_y = surface_points["Y [m]"].min()

    # Get the specific point
    target_point = surface_points[
        (surface_points["X [m]"] == min_x) & (surface_points["Y [m]"] == min_y)
    ]

    if len(target_point) > 0:
        x, y, z = target_point.iloc[0][["X [m]", "Y [m]", "Z [m]"]].values
        print(f"Target point found: X={x}, Y={y}, Z={z}")
        return x, y, z
    else:
        print("Could not find target point")
        return None, None, None


########################################################################
# 5) Function that will extract time series for the target point
########################################################################
def extract_point_time_series(flux_df, target_x, target_y, target_z):
    tolerance = 1e-6  # Small tolerance for floating point comparison

    # Filter data for the target point
    point_data = flux_df[
        (np.abs(flux_df["X [m]"] - target_x) < tolerance)
        & (np.abs(flux_df["Y [m]"] - target_y) < tolerance)
        & (np.abs(flux_df["Z [m]"] - target_z) < tolerance)
    ].copy()

    # Sort by time
    point_data = point_data.sort_values("Time Index")

    print(f"Extracted {len(point_data)} time points for target location")

    if len(point_data) > 0:
        print(
            f"Time range: {point_data['Time Index'].min()} to {point_data['Time Index'].max()}"
        )
        print(
            f"CO2 flux range: {point_data['CO2_flux_magnitude'].min():.3e} to {point_data['CO2_flux_magnitude'].max():.3e}"
        )
        print(
            f"CH4 flux range: {point_data['CH4_flux_magnitude'].min():.3e} to {point_data['CH4_flux_magnitude'].max():.3e}"
        )

    return point_data


########################################################################
# 6) Creating the bokeh plot
########################################################################
def create_flux_time_series_plot(point_data, target_x, target_y, target_z):
    from bokeh.plotting import figure, save, output_file
    from bokeh.layouts import column
    from bokeh.models import HoverTool, ColumnDataSource, Legend
    from bokeh.io import curdoc

    if len(point_data) == 0:
        print("No data available for plotting")
        return

    print("Creating flux time series visualization...")

    output_file("flux_time_series.html")

    # Prepare data source
    source = ColumnDataSource(point_data)

    # Create the main plot
    p = figure(
        title=f"CO2 and CH4 Flux Time Series - Point ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})",
        x_axis_label="Time Index",
        y_axis_label="Flux Magnitude [M·m⁻¹]",
        width=1000,
        height=600,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        toolbar_location="above",
    )

    # Plot CO2 flux
    co2_line = p.line(
        "Time Index",
        "CO2_flux_magnitude",
        source=source,
        line_width=3,
        color="blue",
        alpha=0.8,
        legend_label="CO2 Flux Magnitude",
    )

    co2_circles = p.circle(
        "Time Index",
        "CO2_flux_magnitude",
        source=source,
        size=8,
        color="blue",
        alpha=0.8,
    )

    # Plot CH4 flux
    ch4_line = p.line(
        "Time Index",
        "CH4_flux_magnitude",
        source=source,
        line_width=3,
        color="green",
        alpha=0.8,
        legend_label="CH4 Flux Magnitude",
    )

    ch4_circles = p.circle(
        "Time Index",
        "CH4_flux_magnitude",
        source=source,
        size=8,
        color="green",
        alpha=0.8,
    )

    # Add hover tools
    hover_co2 = HoverTool(
        tooltips=[
            ("Time Index", "@{Time Index}"),
            ("CO2 Flux Magnitude [M·m⁻¹]", "@CO2_flux_magnitude{0.000e+0}"),
            ("CO2 Flux X [M·m⁻¹]", "@CO2_flux_x{0.000e+0}"),
            ("CO2 Flux Y [M·m⁻¹]", "@CO2_flux_y{0.000e+0}"),
            ("CO2 Flux Z [M·m⁻¹]", "@CO2_flux_z{0.000e+0}"),
            ("CO2 Concentration [M]", "@{CO2(aq) [M]}{0.000e+0}"),
        ],
        renderers=[co2_line, co2_circles],
    )

    hover_ch4 = HoverTool(
        tooltips=[
            ("Time Index", "@{Time Index}"),
            ("CH4 Flux Magnitude [M·m⁻¹]", "@CH4_flux_magnitude{0.000e+0}"),
            ("CH4 Flux X [M·m⁻¹]", "@CH4_flux_x{0.000e+0}"),
            ("CH4 Flux Y [M·m⁻¹]", "@CH4_flux_y{0.000e+0}"),
            ("CH4 Flux Z [M·m⁻¹]", "@CH4_flux_z{0.000e+0}"),
            ("CH4 Concentration [M]", "@{Free CH4(aq) [M]}{0.000e+0}"),
        ],
        renderers=[ch4_line, ch4_circles],
    )

    p.add_tools(hover_co2, hover_ch4)

    # Customize legend
    p.legend.location = "top_left"
    p.legend.click_policy = "hide"
    p.legend.label_text_font_size = "12pt"

    # Use scientific notation for y-axis if needed
    if (
        point_data["CO2_flux_magnitude"].max() < 1e-3
        or point_data["CH4_flux_magnitude"].max() < 1e-3
    ):
        p.yaxis.formatter.use_scientific = True

    # Create a secondary plot for individual flux components
    p2 = figure(
        title=f"Flux Components - Point ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})",
        x_axis_label="Time Index",
        y_axis_label="Flux Component [M·m⁻¹]",
        width=1000,
        height=400,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        toolbar_location="above",
    )

    # Plot flux components
    p2.line(
        "Time Index",
        "CO2_flux_x",
        source=source,
        line_width=2,
        color="red",
        alpha=0.7,
        legend_label="CO2 Flux X",
    )
    p2.line(
        "Time Index",
        "CO2_flux_y",
        source=source,
        line_width=2,
        color="orange",
        alpha=0.7,
        legend_label="CO2 Flux Y",
    )
    p2.line(
        "Time Index",
        "CO2_flux_z",
        source=source,
        line_width=2,
        color="purple",
        alpha=0.7,
        legend_label="CO2 Flux Z",
    )

    p2.line(
        "Time Index",
        "CH4_flux_x",
        source=source,
        line_width=2,
        color="darkgreen",
        alpha=0.7,
        legend_label="CH4 Flux X",
        line_dash="dashed",
    )
    p2.line(
        "Time Index",
        "CH4_flux_y",
        source=source,
        line_width=2,
        color="lightgreen",
        alpha=0.7,
        legend_label="CH4 Flux Y",
        line_dash="dashed",
    )
    p2.line(
        "Time Index",
        "CH4_flux_z",
        source=source,
        line_width=2,
        color="olive",
        alpha=0.7,
        legend_label="CH4 Flux Z",
        line_dash="dashed",
    )

    p2.legend.location = "top_left"
    p2.legend.click_policy = "hide"
    p2.legend.label_text_font_size = "10pt"

    if (
        point_data["CO2_flux_x"].abs().max() < 1e-3
        or point_data["CH4_flux_x"].abs().max() < 1e-3
    ):
        p2.yaxis.formatter.use_scientific = True

    # Create final layout
    layout = column(p, p2)
    save(layout)

    print("Flux time series visualization saved as: flux_time_series.html")


########################################################################
# 0) Calling main
########################################################################
def main():
    """Main execution function."""
    print("Starting PFLOTRAN Flux Time Series Visualization...")

    # Install dependencies
    install_dependencies()

    # Load data
    df = load_data()

    # Check if we have the required gas species columns
    required_columns = ["CO2(aq) [M]", "Free CH4(aq) [M]"]
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        print(f"Missing required columns: {missing_columns}")
        print("Available columns:")
        for col in df.columns:
            if any(gas in col.lower() for gas in ["co2", "ch4", "methane"]):
                print(f"   - {col}")
        return

    print(
        f"Dataset contains {len(df)} data points across {df['Time Index'].nunique()} time steps"
    )

    # Find the target point (top-left at surface)
    target_x, target_y, target_z = find_target_point(df)
    if target_x is None:
        print("Cannot proceed without target point. Exiting.")
        return

    # Calculate concentration gradients and fluxes
    print("Calculating flux fields...")
    flux_df = calculate_concentration_gradients(df)

    # Extract time series for the target point
    print("Extracting time series for target point...")
    point_data = extract_point_time_series(flux_df, target_x, target_y, target_z)

    if len(point_data) == 0:
        print("No time series data found for target point. Exiting.")
        return

    # Create visualization
    create_flux_time_series_plot(point_data, target_x, target_y, target_z)

    print("\nFlux time series visualization complete!")
    print("\nGenerated files:")
    print("   flux_time_series.html - Interactive flux time series plot")
    print("\nTo view the results:")
    print("   Open the HTML file in your web browser")
    print("   Or use: open flux_time_series.html")


if __name__ == "__main__":
    main()
