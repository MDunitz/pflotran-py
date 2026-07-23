"""Bokeh 2D presentation layer for gradient / flux visualizations.

Generalized primitives (panels, hovers, colormaps, save) plus the two
top-level renderers used by the pipeline: surface maps and single-point time
series. Both render gradient or flux depending on ``use_flux``. The tooltip
builders live here (not in the physics module) because they are presentation.
"""

import logging
import os

import numpy as np
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import column, row
from bokeh.models import HoverTool, ColumnDataSource, ColorBar, LinearColorMapper
from bokeh.transform import linear_cmap
from bokeh.palettes import Viridis256

from .columns import (
    GRADIENT_UNITS,
    FLUX_UNITS,
    TIME_COL,
    gradient_col,
    flux_col,
    time_axis_column,
)

logger = logging.getLogger(__name__)

# Signal floor: below this a gradient/flux magnitude is treated as no signal.
SIGNAL_FLOOR = 1e-25

SURFACE_MARKER_SIZE = 15
SPECIES_LINE_COLORS = {"CO2": "blue", "CH4": "green", "SO4": "orange", "O2": "red"}
COMPONENT_COLORS = {
    "CO2": {"x": "red", "y": "orange", "z": "purple"},
    "CH4": {"x": "darkgreen", "y": "lightgreen", "z": "olive"},
}
_DEFAULT_COMPONENT_COLORS = {"x": "gray", "y": "silver", "z": "black"}


# ═════════════════════════════════════════════════════════════════════
# Gradient-vs-flux switch helpers
# ═════════════════════════════════════════════════════════════════════


def _magnitude_col(species, use_flux):
    col = flux_col if use_flux else gradient_col
    return col(species, "magnitude")


def _component_col(species, component, use_flux):
    col = flux_col if use_flux else gradient_col
    return col(species, component)


def _units_label(use_flux):
    return FLUX_UNITS if use_flux else GRADIENT_UNITS


# ═════════════════════════════════════════════════════════════════════
# Tooltips
# ═════════════════════════════════════════════════════════════════════


def gradient_tooltips(species):
    """Bokeh tooltip entries for a species' gradient columns."""
    return [
        (
            f"{species} |∇C| [{GRADIENT_UNITS}]",
            f'@{gradient_col(species, "magnitude")}',
        ),
        (f"{species} ∂C/∂x [{GRADIENT_UNITS}]", f'@{gradient_col(species, "x")}'),
        (f"{species} ∂C/∂y [{GRADIENT_UNITS}]", f'@{gradient_col(species, "y")}'),
        (f"{species} ∂C/∂z [{GRADIENT_UNITS}]", f'@{gradient_col(species, "z")}'),
    ]


def flux_tooltips(species):
    """Bokeh tooltip entries for a species' flux columns."""
    return [
        (f"{species} |J| [{FLUX_UNITS}]", f'@{flux_col(species, "magnitude")}'),
        (f"{species} Jx [{FLUX_UNITS}]", f'@{flux_col(species, "x")}'),
        (f"{species} Jy [{FLUX_UNITS}]", f'@{flux_col(species, "y")}'),
        (f"{species} Jz [{FLUX_UNITS}]", f'@{flux_col(species, "z")}'),
    ]


def _species_tips(species, use_flux):
    return flux_tooltips(species) if use_flux else gradient_tooltips(species)


def surface_hover(species, conc_col, renderers, use_flux):
    """Hover tool for a surface map panel (coords + concentration + q-tips)."""
    return HoverTool(
        tooltips=[
            ("X [m]", "@{X [m]}{0.00}"),
            ("Y [m]", "@{Y [m]}{0.00}"),
            ("Z [m]", "@{Z [m]}{0.00}"),
            (f"{species} Concentration [M]", f"@{{{conc_col}}}"),
        ]
        + _species_tips(species, use_flux),
        renderers=renderers,
    )


def timeseries_hover(species_list, conc_cols, renderers, use_flux, time_col):
    """Hover tool for the time-series panels."""
    tips = [(time_col, f"@{{{time_col}}}")]
    for species in species_list:
        tips.extend(_species_tips(species, use_flux))
        tips.append(
            (f"{species} Concentration [M]", f"@{{{conc_cols[species]}}}{{0.000e+0}}")
        )
    return HoverTool(tooltips=tips, renderers=renderers)


# ═════════════════════════════════════════════════════════════════════
# Save
# ═════════════════════════════════════════════════════════════════════


def save_layout_html(layout, output_dir, output_filename):
    """Write a Bokeh layout to output_dir/output_filename, returning the path."""
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    output_file(output_path)
    save(layout)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Surface maps
# ═════════════════════════════════════════════════════════════════════


def _surface_panel(surface_data, species, conc_col, use_flux, time_idx):
    """One 2D surface-map panel for a single species/timestep, or None."""
    quantity = "Flux |J|" if use_flux else "Gradient |∇C|"
    units = _units_label(use_flux)
    mag_col = _magnitude_col(species, use_flux)
    vmin = surface_data[mag_col].min()
    vmax = surface_data[mag_col].max()

    if TIME_COL in surface_data.columns:
        t = surface_data[time_axis_column(surface_data)].iloc[0]
        title = f"{species} {quantity} [{units}] — t={t:g} d"
    else:
        title = f"{species} {quantity} [{units}] — t={time_idx}"

    p = figure(
        title=title,
        x_axis_label="X [m]",
        y_axis_label="Y [m]",
        width=400,
        height=400,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    if vmax > vmin and vmax > SIGNAL_FLOOR:
        source = ColumnDataSource(surface_data)
        circles = p.scatter(
            "X [m]",
            "Y [m]",
            size=SURFACE_MARKER_SIZE,
            color=linear_cmap(mag_col, Viridis256, vmin, vmax),
            source=source,
            alpha=0.7,
        )
        p.add_tools(surface_hover(species, conc_col, [circles], use_flux))

    return p


def render_surface_maps(
    df,
    species_map,
    n_timesteps,
    output_dir,
    output_filename,
    use_flux,
):
    """Render 2D surface maps of gradient or flux magnitude to an HTML file.

    df must contain the gradient (and, for use_flux, the flux) columns and an
    'is_surface' flag from transforms.identify_surface_cells().
    """
    species_list = list(species_map.keys())
    time_indices = sorted(df["Time Index"].unique())

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

        magnitudes = [
            surface_data[_magnitude_col(s, use_flux)].max() for s in species_list
        ]
        if all(m <= SIGNAL_FLOOR for m in magnitudes):
            logger.debug("Skipping time step %s — no signal above floor", time_idx)
            continue

        row_plots = [
            _surface_panel(surface_data, s, species_map[s], use_flux, time_idx)
            for s in species_list
        ]
        plots.append(row(*row_plots))

    output_path = save_layout_html(column(*plots), output_dir, output_filename)
    logger.info("Surface visualization saved: %s", output_path)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Time series at a single point
# ═════════════════════════════════════════════════════════════════════


def _timeseries_magnitude_panel(
    source, point_data, species_list, species_map, use_flux, time_col, target
):
    target_x, target_y, target_z = target
    quantity = "Flux" if use_flux else "Gradient"
    units = _units_label(use_flux)

    p = figure(
        title=f"{quantity} Magnitude — ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})",
        x_axis_label=time_col,
        y_axis_label=f"{quantity} [{units}]",
        width=1000,
        height=500,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    renderers = []
    for species in species_list:
        mag = _magnitude_col(species, use_flux)
        color = SPECIES_LINE_COLORS.get(species, "gray")
        line = p.line(
            time_col,
            mag,
            source=source,
            line_width=3,
            color=color,
            alpha=0.8,
            legend_label=f"{species} |{quantity[0]}|",
        )
        circles = p.scatter(
            time_col, mag, source=source, size=6, color=color, alpha=0.8
        )
        renderers.extend([line, circles])

    p.add_tools(
        timeseries_hover(species_list, species_map, renderers, use_flux, time_col)
    )
    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    max_vals = [point_data[_magnitude_col(s, use_flux)].max() for s in species_list]
    if any(v < 1e-3 for v in max_vals):
        p.yaxis.formatter.use_scientific = True
    return p


def _timeseries_components_panel(
    source, point_data, species_list, use_flux, time_col, target
):
    target_x, target_y, target_z = target
    quantity = "Flux" if use_flux else "Gradient"
    units = _units_label(use_flux)

    p = figure(
        title=f"{quantity} Components — ({target_x:.3f}, {target_y:.3f}, {target_z:.3f})",
        x_axis_label=time_col,
        y_axis_label=f"{quantity} component [{units}]",
        width=1000,
        height=350,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    for species in species_list:
        sc = COMPONENT_COLORS.get(species, _DEFAULT_COMPONENT_COLORS)
        dash = "dashed" if species == "CH4" else "solid"
        for comp in ["x", "y", "z"]:
            p.line(
                time_col,
                _component_col(species, comp, use_flux),
                source=source,
                line_width=2,
                color=sc[comp],
                alpha=0.7,
                legend_label=f"{species} {comp}",
                line_dash=dash,
            )

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"
    p.legend.label_text_font_size = "10pt"

    max_component_vals = [
        point_data[_component_col(s, "x", use_flux)].abs().max() for s in species_list
    ]
    if any(v < 1e-3 for v in max_component_vals):
        p.yaxis.formatter.use_scientific = True
    return p


def render_time_series(
    point_data,
    target_coords,
    species_map,
    output_dir,
    output_filename,
    use_flux,
):
    """Render magnitude + component time-series panels to an HTML file."""
    species_list = list(species_map.keys())
    time_col = time_axis_column(point_data)
    source = ColumnDataSource(point_data)

    magnitude_panel = _timeseries_magnitude_panel(
        source, point_data, species_list, species_map, use_flux, time_col, target_coords
    )
    components_panel = _timeseries_components_panel(
        source, point_data, species_list, use_flux, time_col, target_coords
    )

    output_path = save_layout_html(
        column(magnitude_panel, components_panel), output_dir, output_filename
    )
    quantity = "Flux" if use_flux else "Gradient"
    logger.info("%s time series saved: %s", quantity, output_path)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Domain-mean time series (mean over all cells per snapshot)
# ═════════════════════════════════════════════════════════════════════


def render_mean_timeseries(
    flux_df, series_map, y_label, title, output_dir, filename, sci=False
):
    """Render domain-mean-vs-time lines (one per series) to an HTML file.

    series_map maps a legend label (species short name) to the DataFrame
    column whose per-snapshot spatial mean is plotted.
    """
    time_col = time_axis_column(flux_df)
    grouped = flux_df.groupby("Time Index")
    times = grouped[time_col].first().tolist()

    p = figure(
        title=title,
        x_axis_label=time_col,
        y_axis_label=y_label,
        width=800,
        height=500,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    for label, col in series_map.items():
        means = grouped[col].mean().tolist()
        color = SPECIES_LINE_COLORS.get(label, "gray")
        p.line(times, means, line_width=2, color=color, alpha=0.8, legend_label=label)
        p.scatter(times, means, size=6, color=color, alpha=0.8)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"
    if sci:
        p.yaxis.formatter.use_scientific = True

    return save_layout_html(p, output_dir, filename)


# ═════════════════════════════════════════════════════════════════════
# Arbitrary-z concentration slice (XY map at one depth + timestep)
# ═════════════════════════════════════════════════════════════════════


def render_z_slice(
    df, conc_col, output_dir, filename, z=None, time_idx=None, color_label=None
):
    """Render an XY concentration map at one Z layer and timestep to HTML.

    Defaults to the final timestep and the middle Z layer. Distinct from the
    surface maps (top layer, gradient/flux); this is a raw-concentration slice
    at an arbitrary depth.
    """
    color_label = color_label or conc_col
    time_col = time_axis_column(df)

    if time_idx is None:
        time_idx = sorted(df["Time Index"].unique())[-1]
    snap = df[df["Time Index"] == time_idx]

    if z is None:
        z_levels = sorted(snap["Z [m]"].unique())
        z = z_levels[len(z_levels) // 2]
    layer = snap[np.isclose(snap["Z [m]"], z)]

    t_label = layer[time_col].iloc[0] if time_col in layer.columns else time_idx
    vmin = float(layer[conc_col].min())
    vmax = float(layer[conc_col].max())
    v_high = vmax if vmax > vmin else vmin + SIGNAL_FLOOR
    mapper = LinearColorMapper(palette=Viridis256, low=vmin, high=v_high)

    p = figure(
        title=f"{color_label} at z={z:g} m, t={t_label:g} d",
        x_axis_label="X [m]",
        y_axis_label="Y [m]",
        width=560,
        height=520,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    source = ColumnDataSource(layer)
    renderer = p.scatter(
        "X [m]",
        "Y [m]",
        size=18,
        marker="square",
        color={"field": conc_col, "transform": mapper},
        source=source,
        alpha=0.9,
    )
    p.add_layout(ColorBar(color_mapper=mapper, title=color_label), "right")
    p.add_tools(
        HoverTool(
            tooltips=[
                ("X [m]", "@{X [m]}{0.00}"),
                ("Y [m]", "@{Y [m]}{0.00}"),
                (color_label, f"@{{{conc_col}}}"),
            ],
            renderers=[renderer],
        )
    )
    return save_layout_html(p, output_dir, filename)
