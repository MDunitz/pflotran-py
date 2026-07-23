"""Column naming, unit strings, and time-axis labeling.

Single source of truth linking computed columns to their display units, so
gradient/flux column names and Bokeh/Plotly tooltips cannot drift apart.
"""

# ═════════════════════════════════════════════════════════════════════
# Unit strings — displayed in tooltips and axis labels
# ═════════════════════════════════════════════════════════════════════

GRADIENT_UNITS = "M/m"  # mol/L per meter
FLUX_UNITS = "mol/(m²·s)"  # molar flux

# ═════════════════════════════════════════════════════════════════════
# Time-axis labeling
# ═════════════════════════════════════════════════════════════════════

# Canonical simulation-time column written by the extractors (Tecplot + HDF5).
# Values are always converted to days so annual / multi-year runs label axes
# correctly instead of using snapshot ordinals (0, 1, 2, ...).
TIME_COL = "Time [d]"
TIME_UNIT_TO_DAYS = {
    "s": 1.0 / 86400.0,
    "sec": 1.0 / 86400.0,
    "m": 1.0 / 1440.0,
    "min": 1.0 / 1440.0,
    "h": 1.0 / 24.0,
    "hr": 1.0 / 24.0,
    "d": 1.0,
    "day": 1.0,
    "y": 365.25,
    "yr": 365.25,
    "year": 365.25,
}


def time_to_days(value, unit):
    """Convert a PFLOTRAN time value + unit token to days."""
    try:
        scale = TIME_UNIT_TO_DAYS[unit.lower()]
    except KeyError as exc:
        raise ValueError(f"Unsupported PFLOTRAN time unit: {unit!r}") from exc
    return float(value) * scale


def time_axis_column(df):
    """Prefer real simulation time; fall back to snapshot ordinal."""
    return TIME_COL if TIME_COL in df.columns else "Time Index"


# ═════════════════════════════════════════════════════════════════════
# Column naming — links computation to display
# ═════════════════════════════════════════════════════════════════════


def gradient_col(species, component):
    """Column name for a gradient component.

    Examples:
        gradient_col('CO2', 'x')         -> 'CO2_grad_x'
        gradient_col('CO2', 'magnitude') -> 'CO2_grad_magnitude'
    """
    return f"{species}_grad_{component}"


def flux_col(species, component):
    """Column name for a flux component.

    Examples:
        flux_col('CO2', 'x')         -> 'CO2_flux_x'
        flux_col('CO2', 'magnitude') -> 'CO2_flux_magnitude'
    """
    return f"{species}_flux_{component}"
