"""Tests for the bottle starting-carbon inventory helper."""

import pytest

from pflotran_py.comparison.carbon_inventory import (
    CARBONS_PER_DOM1,
    MEASURED_INCUBATION_STARTING_CARBON_MOLES,
    cellulose_volume_fraction_for_starting_carbon,
    default_comparison_starting_carbon_moles,
    parse_constraint_value,
    starting_carbon_moles,
)
from pflotran_py.generator.bottle_generator import LIQUID_VOLUME_L
from pflotran_py.generator.pflotran_generator import DEFAULT_CELLULOSE_HYDROLYSIS


def test_parse_constraint_value_handles_fortran_exponents():
    assert parse_constraint_value("5.00 T") == pytest.approx(5.0)
    assert parse_constraint_value("1.00d-03 T") == pytest.approx(1e-3)
    assert parse_constraint_value("2.d-4 T") == pytest.approx(2e-4)


def test_dissolved_only_inventory_is_glucose_times_liquid_volume():
    moles = starting_carbon_moles(cellulose_hydrolysis=False)
    assert moles == pytest.approx(5.0 * LIQUID_VOLUME_L * CARBONS_PER_DOM1)


def test_default_cellulose_pool_matches_measured_incubation_starting_c():
    """Deck inventory normalisation: model C ≈ recipe-derived biomass C."""
    expected_vf = cellulose_volume_fraction_for_starting_carbon(
        MEASURED_INCUBATION_STARTING_CARBON_MOLES
    )
    assert DEFAULT_CELLULOSE_HYDROLYSIS["volume_fraction"] == pytest.approx(
        expected_vf, rel=1e-5
    )

    moles = starting_carbon_moles(cellulose_hydrolysis=True)
    assert moles == pytest.approx(MEASURED_INCUBATION_STARTING_CARBON_MOLES, rel=1e-5)
    assert default_comparison_starting_carbon_moles() == pytest.approx(moles)


def test_volume_fraction_helper_round_trips():
    vf = cellulose_volume_fraction_for_starting_carbon(0.0565)
    moles = starting_carbon_moles(
        cellulose_hydrolysis={"volume_fraction": vf},
    )
    assert moles == pytest.approx(0.0565, rel=1e-6)
