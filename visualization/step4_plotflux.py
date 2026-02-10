#!/usr/bin/env python3
"""
PFLOTRAN Gradient & Flux Time Series

Plots concentration gradient and diffusive flux time series for a single
grid point (default: top-left surface cell).

Units are linked to shared_utils — tooltips always match computed columns.
"""

import os
import numpy as np
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import column
from bokeh.models import HoverTool, ColumnDataSource

from shared_utils import (
    load_data,
    calculate_gradients,
    convert_to_flux,
    gradient_col,
    flux_col,
    gradient_tooltips,
    flux_tooltips,
    GRADIENT_UNITS,
    FLUX_UNITS,
)

DEFAULT_SPECIES_MAP = {
    "CO2": "CO2(aq) [M]",
    "CH4": "Free CH4(aq) [M]",
}

DEFAULT_TEMPERATURE_C = 8.0


def find_target_point(df):
    """Find the top-left surface cell (min X, min Y, max Z)."""
    max_z = df["Z [m]"].max()
    surface = df[df["Z [m]"] == max_z]
    min_x = surface["X [m]"].min()
    min_y = surface["Y [m]"].min()
    target = surface[(surface["X [m]"] == min_x) & (surface["Y [m]"] == min_y)]
    return target.iloc[0][["X [m]", "Y [m]", "Z [m]"]].values


def extract_point_time_series(df, target_x, target_y, target_z, tolerance=1e-6):
    """Extract all timesteps for a single grid point."""
    point_data = df[
        (np.abs(df["X [m]"] - target_x) < tolerance)
        & (np.abs(df["Y [m]"] - target_y) < tolerance)
        & (np.abs(df["Z [m]"] - target_z) < tolerance)
    ].copy()
    return point_data.sort_values("Time Index")


def _col(species, component, use_flux):
    """Return the correct column name for gradient or flux."""
    return (
        flux_col(species, component) if use_flux else gradient_col(species, component)
    )


def _units(use_flux):
    return FLUX_UNITS if use_flux else GRADIENT_UNITS


def _make_hover(species_list, conc_cols, renderers, use_flux):
    """Build hover tooltip from linked shared_utils functions."""
    tips = [("Time Index", "@{Time Index}")]
    for species in species_list:
        getter = flux_tooltips if use_flux else gradient_tooltips
        tips.extend(getter(species))
        tips.append(
            (f"{species} Concentration [M]", f"@{{{conc_cols[species]}}}{{0.000e+0}}")
        )
    return HoverTool(tooltips=tips, renderers=renderers)


def create_time_series_plot(
    point_data,
    target_coords,
    species_map,
    output_dir=".",
    output_filename="flux_time_series.html",
    use_flux=False,
):
    """Create magnitude + component time series plots."""
    target_x, target_y, target_z = target_coords
    species_list = list(species_map.keys())
    units = _units(use_flux)
    quantity = "Flux" if use_flux else "Gradient"

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    output_file(output_path)

    source = ColumnDataSource(point_data)

    # Magnitude plot
    p = figure(
        title=f"{quantity} Magnitude — ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})",
        x_axis_label="Time Index",
        y_axis_label=f"{quantity} [{units}]",
        width=1000,
        height=500,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    colors = {"CO2": "blue", "CH4": "green", "SO4": "orange", "O2": "red"}
    renderers = []
    for species in species_list:
        mag = _col(species, "magnitude", use_flux)
        color = colors.get(species, "gray")
        line = p.line(
            "Time Index",
            mag,
            source=source,
            line_width=3,
            color=color,
            alpha=0.8,
            legend_label=f"{species} |{quantity[0]}|",
        )
        circles = p.scatter(
            "Time Index", mag, source=source, size=6, color=color, alpha=0.8
        )
        renderers.extend([line, circles])

    hover = _make_hover(species_list, species_map, renderers, use_flux)
    p.add_tools(hover)
    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    # Check if scientific notation needed
    max_vals = [point_data[_col(s, "magnitude", use_flux)].max() for s in species_list]
    if any(v < 1e-3 for v in max_vals):
        p.yaxis.formatter.use_scientific = True

    # Component plot
    p2 = figure(
        title=f"{quantity} Components — ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})",
        x_axis_label="Time Index",
        y_axis_label=f"{quantity} component [{units}]",
        width=1000,
        height=350,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    component_colors = {
        "CO2": {"x": "red", "y": "orange", "z": "purple"},
        "CH4": {"x": "darkgreen", "y": "lightgreen", "z": "olive"},
    }
    for species in species_list:
        sc = component_colors.get(species, {"x": "gray", "y": "silver", "z": "black"})
        for comp in ["x", "y", "z"]:
            col_name = _col(species, comp, use_flux)
            dash = "dashed" if species == "CH4" else "solid"
            p2.line(
                "Time Index",
                col_name,
                source=source,
                line_width=2,
                color=sc[comp],
                alpha=0.7,
                legend_label=f"{species} {comp}",
                line_dash=dash,
            )

    p2.legend.location = "top_left"
    p2.legend.click_policy = "hide"
    p2.legend.label_text_font_size = "10pt"

    max_component_vals = [
        point_data[_col(s, "x", use_flux)].abs().max() for s in species_list
    ]
    if any(v < 1e-3 for v in max_component_vals):
        p2.yaxis.formatter.use_scientific = True

    layout = column(p, p2)
    save(layout)
    print(f"{quantity} time series saved: {output_path}")


def main(
    species_map=None,
    output_dir=".",
    output_filename="gradient_time_series.html",
    compute_flux=True,
    temperature_c=DEFAULT_TEMPERATURE_C,
):
    species_map = species_map or DEFAULT_SPECIES_MAP
    print("Starting PFLOTRAN Time Series Visualization...")

    df = load_data()
    target_x, target_y, target_z = find_target_point(df)
    print(f"Target point: ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})")

    # Compute gradients
    print("Computing concentration gradients...")
    df = calculate_gradients(df, species_map)

    if compute_flux:
        print(f"Converting to true flux (T={temperature_c}°C)...")
        df = convert_to_flux(df, list(species_map.keys()), temperature_c=temperature_c)

    point_data = extract_point_time_series(df, target_x, target_y, target_z)

    # Always produce gradient plot
    create_time_series_plot(
        point_data,
        (target_x, target_y, target_z),
        species_map,
        output_dir=output_dir,
        output_filename=output_filename,
        use_flux=False,
    )

    # If flux computed, also produce flux plot
    if compute_flux:
        flux_filename = output_filename.replace(".html", "_flux.html")
        create_time_series_plot(
            point_data,
            (target_x, target_y, target_z),
            species_map,
            output_dir=output_dir,
            output_filename=flux_filename,
            use_flux=True,
        )

    print("Step 4 complete.")


if __name__ == "__main__":
    main()
