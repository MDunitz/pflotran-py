"""Upstream a_w scaling for fermentation / hydrolysis on closed-batch decks."""

import pytest

from pflotran_py.comparison.decks import (
    _apply_aw_upstream_inhibition,
    one_minus_aw_factor,
)
from pflotran_py.generator.pflotran_generator import DEFAULT_RATE_CONSTANTS


def test_one_minus_aw_factor_matches_sandbox_shape():
    assert one_minus_aw_factor(1.0, 0.85) == 1.0
    assert one_minus_aw_factor(0.85, 0.85) == 0.0
    assert one_minus_aw_factor(0.80, 0.85) == 0.0
    assert one_minus_aw_factor(0.925, 0.85) == pytest.approx(0.5)


def test_upstream_inhibition_scales_fermentation_and_hydrolysis():
    kwargs = _apply_aw_upstream_inhibition(
        {
            "fixed_water_activity": 0.91,
            "aw_threshold_fermentation": 0.85,
            "aw_threshold_hydrolysis": 0.85,
            "cellulose_hydrolysis": {},
        }
    )
    expected = one_minus_aw_factor(0.91, 0.85)
    assert kwargs["rate_constants"]["fermentation"] == pytest.approx(
        DEFAULT_RATE_CONSTANTS["fermentation"] * expected
    )
    hydro_rate = float(
        kwargs["cellulose_hydrolysis"]["rate_constant"].lower().replace("d", "e")
    )
    assert hydro_rate == pytest.approx(2.0e-7 * expected)


def test_control_aw_leaves_upstream_rates_unchanged():
    kwargs = _apply_aw_upstream_inhibition(
        {
            "fixed_water_activity": 1.0,
            "cellulose_hydrolysis": {},
        }
    )
    assert kwargs["rate_constants"]["fermentation"] == DEFAULT_RATE_CONSTANTS[
        "fermentation"
    ]
    hydro_rate = float(
        kwargs["cellulose_hydrolysis"]["rate_constant"].lower().replace("d", "e")
    )
    assert hydro_rate == pytest.approx(2.0e-7)
