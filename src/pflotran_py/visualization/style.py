"""Presentation constants for Bokeh visualizations.

Colors and marker sizes live here (plotting-only). Numeric-zero guards are
also presentation helpers: they are not physical detection thresholds and must
not be treated as species- or unit-specific cutoffs (see issue #28).
"""

SURFACE_MARKER_SIZE = 15
SPECIES_LINE_COLORS = {"CO2": "blue", "CH4": "green", "SO4": "orange", "O2": "red"}
COMPONENT_COLORS = {
    "CO2": {"x": "red", "y": "orange", "z": "purple"},
    "CH4": {"x": "darkgreen", "y": "lightgreen", "z": "olive"},
}
DEFAULT_COMPONENT_COLORS = {"x": "gray", "y": "silver", "z": "black"}

# Absolute floor catches true zeros / underflow. Relative tolerance makes
# empty-panel and colormap checks independent of quantity units (M/m vs
# mol/(m²·s)), which differ by many orders of magnitude across species.
NUMERIC_ABS_FLOOR = 0.0
NUMERIC_REL_TOLERANCE = 1e-12


def has_resolvable_span(
    vmin,
    vmax,
    *,
    abs_floor=NUMERIC_ABS_FLOOR,
    rel_tol=NUMERIC_REL_TOLERANCE,
):
    """True if [vmin, vmax] is wide enough to map to a colormap."""
    vmin = float(vmin)
    vmax = float(vmax)
    span = vmax - vmin
    scale = max(abs(vmin), abs(vmax), abs_floor, 1.0)
    return span > max(abs_floor, rel_tol * scale)


def has_resolvable_magnitude(
    magnitude,
    *,
    abs_floor=NUMERIC_ABS_FLOOR,
    rel_tol=NUMERIC_REL_TOLERANCE,
):
    """True if a magnitude is above the numeric-zero guard.

    The relative term uses a unit scale of 1.0 so the same guard applies to
    gradients and fluxes without a per-quantity constant. This is a floating-
    point hygiene check, not a physical detection limit.
    """
    return abs(float(magnitude)) > max(abs_floor, rel_tol)


def colormap_high(
    vmin,
    vmax,
    *,
    abs_floor=NUMERIC_ABS_FLOOR,
    rel_tol=NUMERIC_REL_TOLERANCE,
):
    """Return a ColorMapper high that is strictly above ``vmin`` when flat."""
    vmin = float(vmin)
    vmax = float(vmax)
    if has_resolvable_span(vmin, vmax, abs_floor=abs_floor, rel_tol=rel_tol):
        return vmax
    scale = max(abs(vmin), abs(vmax), 1.0)
    return vmin + max(abs_floor, rel_tol * scale)
