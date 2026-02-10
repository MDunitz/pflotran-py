#!/usr/bin/env python3
"""
PFLOTRAN Surface Gradient & Flux Visualization

Computes concentration gradients (∇C) and optionally true diffusive flux
(J = -D·∇C) for CO₂ and CH₄, then renders 2D surface maps in Bokeh.

Units are pulled from shared_utils so tooltips always match computed values.
"""

import os
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import column, row
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.transform import linear_cmap
from bokeh.palettes import Viridis256

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

# Default species mapping: short name → DataFrame column
DEFAULT_SPECIES_MAP = {
    "CO2": "CO2(aq) [M]",
    "CH4": "Free CH4(aq) [M]",
}

# Default simulation temperature for Stokes-Einstein correction
DEFAULT_TEMPERATURE_C = 8.0


def create_hover_tool(species, conc_col, renderers, use_flux=False):
    """Create Bokeh hover tool with linked units from shared_utils.

    Parameters
    ----------
    species : str
        Short species name (e.g. 'CO2')
    conc_col : str
        DataFrame column name for concentration
    renderers : list
        Bokeh renderers to attach the hover to
    use_flux : bool
        If True, show flux tooltips; otherwise show gradient tooltips
    """
    tips = gradient_tooltips(species) if not use_flux else flux_tooltips(species)
    return HoverTool(
        tooltips=[
            ("X [m]", "@{X [m]}{0.00}"),
            ("Y [m]", "@{Y [m]}{0.00}"),
            ("Z [m]", "@{Z [m]}{0.00}"),
            (f"{species} Concentration [M]", f"@{{{conc_col}}}"),
        ]
        + tips,
        renderers=renderers,
    )


def identify_surface_cells(df):
    """Mark cells at the top of each (X,Y) column as surface cells."""
    surface_df = df.copy()
    surface_df["is_surface"] = False

    for time_idx in df["Time Index"].unique():
        time_mask = df["Time Index"] == time_idx
        time_data = df[time_mask]

        max_z_by_xy = time_data.groupby(["X [m]", "Y [m]"])["Z [m]"].max().reset_index()
        max_z_by_xy.columns = ["X [m]", "Y [m]", "max_z"]

        merged = time_data.merge(max_z_by_xy, on=["X [m]", "Y [m]"])
        surface_indices = merged[merged["Z [m]"] == merged["max_z"]].index
        surface_df.loc[surface_indices, "is_surface"] = True

    return surface_df


def _magnitude_col(species, use_flux):
    """Return the magnitude column name for gradient or flux."""
    return (
        flux_col(species, "magnitude")
        if use_flux
        else gradient_col(species, "magnitude")
    )


def _units_label(use_flux):
    return FLUX_UNITS if use_flux else GRADIENT_UNITS


def create_surface_visualization(
    df,
    species_map,
    n_timesteps=5,
    output_dir=".",
    output_filename="co2_ch4_surface_visualization.html",
    use_flux=False,
    verbose=False,
):
    """Create 2D surface maps of gradient or flux magnitude.

    Parameters
    ----------
    df : DataFrame
        Must contain gradient (and optionally flux) columns from shared_utils.
    species_map : dict
        Short name → concentration column name
    n_timesteps : int
        Number of time steps to plot
    use_flux : bool
        If True, color by flux magnitude; otherwise by gradient magnitude
    """
    time_indices = sorted(df["Time Index"].unique())
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    output_file(output_path)

    species_list = list(species_map.keys())
    units = _units_label(use_flux)
    quantity = "Flux |J|" if use_flux else "Gradient |∇C|"

    plots = []
    for time_idx in time_indices[:n_timesteps]:
        time_data = df[df["Time Index"] == time_idx]
        surface_data = (
            time_data[time_data["is_surface"]]
            if "is_surface" in time_data.columns
            else time_data
        )
        if len(surface_data) == 0:
            surface_data = time_data

        # Skip if no signal at all
        magnitudes = [
            surface_data[_magnitude_col(s, use_flux)].max() for s in species_list
        ]
        if all(m <= 1e-25 for m in magnitudes):
            if verbose:
                print(f"Skipping time step {time_idx} — no signal above 1e-25")
            continue

        row_plots = []
        for species in species_list:
            mag_col = _magnitude_col(species, use_flux)
            vmin = surface_data[mag_col].min()
            vmax = surface_data[mag_col].max()

            p = figure(
                title=f"{species} {quantity} [{units}] — t={time_idx}",
                x_axis_label="X [m]",
                y_axis_label="Y [m]",
                width=400,
                height=400,
                tools="pan,wheel_zoom,box_zoom,reset,save",
            )

            if vmax > vmin and vmax > 1e-25:
                source = ColumnDataSource(surface_data)
                circles = p.scatter(
                    "X [m]",
                    "Y [m]",
                    size=15,
                    color=linear_cmap(mag_col, Viridis256, vmin, vmax),
                    source=source,
                    alpha=0.7,
                )
                hover = create_hover_tool(
                    species, species_map[species], [circles], use_flux=use_flux
                )
                p.add_tools(hover)

            row_plots.append(p)

        plots.append(row(*row_plots))

    layout = column(*plots)
    save(layout)
    print(f"Surface visualization saved: {output_path}")


def main(
    species_map=None,
    n_timesteps=5,
    output_dir=".",
    output_filename="co2_ch4_surface_visualization.html",
    compute_flux=True,
    temperature_c=DEFAULT_TEMPERATURE_C,
):
    species_map = species_map or DEFAULT_SPECIES_MAP
    print("Starting PFLOTRAN Surface Visualization...")

    df = load_data()

    # Compute gradients
    print("Computing concentration gradients...")
    df = calculate_gradients(df, species_map)

    # Optionally compute true flux with temperature-corrected D
    if compute_flux:
        print(f"Converting to true flux (Fick's law, T={temperature_c}°C)...")
        df = convert_to_flux(df, list(species_map.keys()), temperature_c=temperature_c)

    # Identify surface cells
    print("Identifying surface cells...")
    df = identify_surface_cells(df)

    # Create both gradient and flux visualizations
    create_surface_visualization(
        df,
        species_map,
        n_timesteps=n_timesteps,
        output_dir=output_dir,
        output_filename=output_filename,
        use_flux=False,
    )

    if compute_flux:
        flux_filename = output_filename.replace(".html", "_flux.html")
        create_surface_visualization(
            df,
            species_map,
            n_timesteps=n_timesteps,
            output_dir=output_dir,
            output_filename=flux_filename,
            use_flux=True,
        )

    print("Step 3 complete.")


if __name__ == "__main__":
    main()
