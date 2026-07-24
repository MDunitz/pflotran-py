"""Plotly helpers for PFLOTRAN 3D concentration visualizations."""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..analysis.columns import TIME_COL, time_axis_column

SPATIAL_COLUMNS = ["X [m]", "Y [m]", "Z [m]"]
DEFAULT_COLORSCALE = "Viridis"
NEAR_ZERO_FLOOR = 1e-15


def color_limits(values):
    """Return (vmin, vmax) with a floor for near-zero concentration ranges."""
    var_min = float(values.min())
    var_max = float(values.max())
    if var_max <= NEAR_ZERO_FLOOR:
        return 0.0, NEAR_ZERO_FLOOR
    return var_min, var_max


def colorbar_tickformat(var_max):
    """Pick a compact tick format for very small concentrations."""
    return ".2e" if var_max < 1e-6 else ".3f"


def variable_display_name(variable):
    """Short label for colorbar titles (drops 'Free'/'Total' prefix)."""
    return variable.split(" ")[1] if " " in variable else variable


def filter_available_variables(df, variables_to_plot):
    """Return only variables present in the dataframe."""
    return [var for var in variables_to_plot if var in df.columns]


def subplot_grid(n_plots, n_cols=2):
    """Return (n_rows, n_cols) for a subplot grid."""
    n_rows = (n_plots + n_cols - 1) // n_cols
    return n_rows, n_cols


def colorbar_position(row, col, n_cols):
    """Place colorbars on the right edge of each subplot column."""
    colorbar_x = 1.02 if col == n_cols else 0.48
    colorbar_y = 0.75 if row == 1 else 0.25
    return colorbar_x, colorbar_y


def format_time_label(time_days):
    """Compact label for slider / frame titles (days or years)."""
    t = float(time_days)
    if abs(t) >= 365.25 and abs(t % 365.25) < 1e-6:
        years = t / 365.25
        return f"{years:g} y"
    if abs(t - round(t)) < 1e-9:
        return f"{round(t):g} d"
    return f"{t:.4g} d"


def animation_play_pause_buttons(frame_duration=100, x=0.1, y=0.02):
    """Standard play/pause controls for Plotly frame animations."""
    return [
        {
            "type": "buttons",
            "showactive": False,
            "x": x,
            "y": y,
            "buttons": [
                {
                    "label": "Play",
                    "method": "animate",
                    "args": [None, {"frame": {"duration": frame_duration}}],
                },
                {
                    "label": "Pause",
                    "method": "animate",
                    "args": [
                        [None],
                        {"frame": {"duration": 0}, "mode": "immediate"},
                    ],
                },
            ],
        }
    ]


def time_slider(frame_keys, labels, prefix="Time: "):
    """Slider control stepping through simulation snapshots."""
    return [
        {
            "steps": [
                {
                    "args": [
                        [str(key)],
                        {"frame": {"duration": 0}, "mode": "immediate"},
                    ],
                    "label": label,
                    "method": "animate",
                }
                for key, label in zip(frame_keys, labels)
            ],
            "currentvalue": {"prefix": prefix},
            "len": 0.8,
            "x": 0.1,
            "y": 0.02,
        }
    ]


def snapshot_times(df):
    """Return ordered (Time Index, display time) pairs for animation frames."""
    time_col = time_axis_column(df)
    cols = ["Time Index"]
    if time_col != "Time Index":
        cols.append(time_col)
    rows = df[cols].drop_duplicates(subset=["Time Index"]).sort_values("Time Index")
    return list(zip(rows["Time Index"].tolist(), rows[time_col].tolist()))


def make_3d_scatter_trace(
    data,
    variable,
    var_min=None,
    var_max=None,
    colorbar_x=1.02,
    colorbar_y=0.5,
    name=None,
):
    """Build one 3D scatter trace for a single variable snapshot."""
    if var_min is None or var_max is None:
        var_min, var_max = color_limits(data[variable])

    return go.Scatter3d(
        x=data["X [m]"],
        y=data["Y [m]"],
        z=data["Z [m]"],
        mode="markers",
        marker=dict(
            size=3,
            color=data[variable],
            colorscale=DEFAULT_COLORSCALE,
            cmin=var_min,
            cmax=var_max,
            colorbar=dict(
                title=dict(
                    text=variable_display_name(variable),
                    font=dict(size=10),
                ),
                x=colorbar_x,
                y=colorbar_y,
                len=0.3,
                thickness=15,
                tickfont=dict(size=8),
                tickformat=colorbar_tickformat(var_max),
            ),
            opacity=0.7,
        ),
        name=name or variable,
    )


def create_single_variable_plot(df, variable):
    """Create a single-variable 3D plot with time animation."""
    time_col = time_axis_column(df)
    plot_cols = SPATIAL_COLUMNS + [variable, time_col]
    if time_col != "Time Index" and "Time Index" in df.columns:
        plot_cols.append("Time Index")
    plot_data = df[plot_cols].copy()
    var_min, var_max = color_limits(plot_data[variable])

    fig = px.scatter_3d(
        plot_data,
        x="X [m]",
        y="Y [m]",
        z="Z [m]",
        color=variable,
        animation_frame=time_col,
        color_continuous_scale=DEFAULT_COLORSCALE,
        range_color=[var_min, var_max],
        labels={variable: variable, time_col: time_col},
        title=f"3D Visualization of {variable} Over Time",
        opacity=0.7,
    )

    fig.update_layout(
        scene=dict(xaxis_title="X [m]", yaxis_title="Y [m]", zaxis_title="Z [m]"),
        coloraxis_colorbar=dict(tickformat=colorbar_tickformat(var_max)),
        updatemenus=animation_play_pause_buttons(frame_duration=100),
    )
    return fig


def create_multi_variable_plot(df, variables_to_plot, verbose=False):
    """Create a multi-variable 3D plot with time animation."""
    available_vars = filter_available_variables(df, variables_to_plot)
    meta_cols = SPATIAL_COLUMNS + ["Time Index"]
    if TIME_COL in df.columns:
        meta_cols.append(TIME_COL)
    if verbose:
        if not available_vars:
            print("None of the requested variables found in data!")
            print(f"Requested: {variables_to_plot}")
            print(
                "Available: " f"{[col for col in df.columns if col not in meta_cols]}"
            )
            return None
        print(f"Creating multi-variable plot for: {available_vars}")

    if not available_vars:
        return None

    n_plots = len(available_vars)
    n_cols = 2
    n_rows, n_cols = subplot_grid(n_plots, n_cols=n_cols)

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        specs=[[{"type": "scatter3d"} for _ in range(n_cols)] for _ in range(n_rows)],
        subplot_titles=available_vars,
        vertical_spacing=0.15,
        horizontal_spacing=0.05,
    )

    snapshots = snapshot_times(df)
    var_limits = {var: color_limits(df[var]) for var in available_vars}
    frame_keys = []
    frame_labels = []

    frames = []
    for time_idx, time_value in snapshots:
        frame_key = str(time_value)
        frame_label = format_time_label(time_value)
        frame_keys.append(frame_key)
        frame_labels.append(frame_label)
        frame_data = []
        time_data = df[df["Time Index"] == time_idx]
        for i, var in enumerate(available_vars):
            row = (i // n_cols) + 1
            col = (i % n_cols) + 1
            var_min, var_max = var_limits[var]
            colorbar_x, colorbar_y = colorbar_position(row, col, n_cols)
            frame_data.append(
                make_3d_scatter_trace(
                    time_data,
                    var,
                    var_min,
                    var_max,
                    colorbar_x=colorbar_x,
                    colorbar_y=colorbar_y,
                    name=f"{var} ({frame_label})",
                )
            )
        frames.append(go.Frame(data=frame_data, name=frame_key))

    initial_data = df[df["Time Index"] == snapshots[0][0]]
    for i, var in enumerate(available_vars):
        row = (i // n_cols) + 1
        col = (i % n_cols) + 1
        var_min, var_max = var_limits[var]
        colorbar_x, colorbar_y = colorbar_position(row, col, n_cols)
        fig.add_trace(
            make_3d_scatter_trace(
                initial_data,
                var,
                var_min,
                var_max,
                colorbar_x=colorbar_x,
                colorbar_y=colorbar_y,
                name=var,
            ),
            row=row,
            col=col,
        )

    slider_prefix = "Time: " if time_axis_column(df) == TIME_COL else "Time Index: "
    fig.update_layout(
        title="Multi-Variable 3D Concentration Visualization Over Time",
        height=900,
        width=1200,
        showlegend=False,
        margin=dict(l=50, r=150, t=80, b=100),
        updatemenus=animation_play_pause_buttons(frame_duration=150),
        sliders=time_slider(frame_keys, frame_labels, prefix=slider_prefix),
    )

    for _ in range(1, n_plots + 1):
        fig.update_scenes(
            xaxis_title="X [m]",
            yaxis_title="Y [m]",
            zaxis_title="Z [m]",
            selector=dict(type="scene"),
        )

    fig.frames = frames
    return fig
