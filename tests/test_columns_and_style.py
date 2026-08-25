"""Unit tests for analysis column helpers and viz style guards."""

import pytest

from pflotran_py.analysis.columns import TIME_UNIT_TO_DAYS, time_to_days
from pflotran_py.visualization.style import (
    has_resolvable_magnitude,
    has_resolvable_span,
    colormap_high,
)


def test_time_unit_table_matches_astropy_subday_units():
    assert TIME_UNIT_TO_DAYS["s"] == pytest.approx(1.0 / 86400.0)
    assert TIME_UNIT_TO_DAYS["h"] == pytest.approx(1.0 / 24.0)
    assert TIME_UNIT_TO_DAYS["d"] == 1.0


def test_time_unit_year_matches_pflotran_365_day_year():
    # PFLOTRAN units.F90: y/yr/year -> 365 * 24 * 3600 seconds (not Julian 365.25).
    assert TIME_UNIT_TO_DAYS["y"] == 365.0
    assert time_to_days(1.0, "y") == pytest.approx(365.0)
    assert time_to_days(2.0, "yr") == pytest.approx(730.0)


def test_has_resolvable_span_rejects_flat_fields():
    assert not has_resolvable_span(1.0, 1.0)
    assert has_resolvable_span(0.0, 1.0)


def test_has_resolvable_magnitude_is_relative_not_species_specific():
    # Same guard for flux-scale and gradient-scale values.
    assert has_resolvable_magnitude(3.8e-11)
    assert not has_resolvable_magnitude(1e-25)
    assert not has_resolvable_magnitude(0.0)


def test_colormap_high_avoids_zero_width():
    assert colormap_high(1.0, 2.0) == 2.0
    assert colormap_high(5.0, 5.0) > 5.0
