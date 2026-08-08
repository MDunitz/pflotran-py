"""Tests for the model-to-instrument unit conversion.

Three things are worth protecting. The Henry constants must keep their
literature values, because they are duplicated from the saltyBiomass repository
and a silent divergence between the two copies would make the model and the
measurement disagree for an invisible reason. The physics must keep its
direction: colder water holds more gas, salt holds less. And the bottle
geometry must keep matching the deck, or the comparison is between two
differently-sized bottles.
"""

import pytest
from astropy import units as u

from pflotran_py.comparison.headspace import (
    DEFAULT_BOTTLE,
    GASES,
    BottleGeometry,
    aqueous_concentration_to_dissolved_moles,
    aqueous_concentration_to_headspace_moles,
    gas_water_partition_coefficient,
    headspace_fraction,
    headspace_moles_to_total_moles,
    henry_solubility,
    henry_solubility_in_brine,
    model_headspace_moles,
    salting_out_factor,
    setschenow_salts_from_composition,
    validity_warnings,
)

# ─────────────────────────────────────────────────────────────────────
# Constants pinned against the literature
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "gas_name,expected_mol_per_l_atm",
    [("CH4", 1.4e-3), ("CO2", 3.3e-2), ("N2O", 2.4e-2)],
)
def test_henry_constants_match_the_published_values(gas_name, expected_mol_per_l_atm):
    """Sander (2023), 298.15 K, converted to the conventional mol/(L*atm).

    These are duplicated from saltyBiomass's brine_gas.py. This test is what
    stops the two copies drifting apart unnoticed.
    """
    solubility = henry_solubility(GASES[gas_name], 25 * u.deg_C)
    conventional = (solubility * (101325 * u.Pa)).to(u.mol / u.L).value
    assert conventional == pytest.approx(expected_mol_per_l_atm, rel=0.02)


def test_methane_is_far_less_soluble_than_carbon_dioxide():
    """A factor of roughly 24. If this ever inverted, every headspace fraction
    in the comparison would be wrong in a way that still looked plausible."""
    ch4 = henry_solubility(GASES["CH4"], 25 * u.deg_C)
    co2 = henry_solubility(GASES["CO2"], 25 * u.deg_C)
    assert float(co2 / ch4) == pytest.approx(23.6, rel=0.1)


# ─────────────────────────────────────────────────────────────────────
# Direction of the physics
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("gas_name", ["CH4", "CO2", "N2O"])
def test_cold_water_holds_more_gas(gas_name):
    cold = henry_solubility(GASES[gas_name], 8 * u.deg_C)
    warm = henry_solubility(GASES[gas_name], 25 * u.deg_C)
    assert cold > warm


@pytest.mark.parametrize("gas_name", ["CH4", "CO2", "N2O"])
def test_salt_reduces_solubility(gas_name):
    """Salting out, not salting in. A factor above 1 would put gas into
    solution as salt is added, which is backwards."""
    assert salting_out_factor(GASES[gas_name], nacl_molarity=3.0) < 1.0
    assert salting_out_factor(GASES[gas_name], mgcl2_molarity=2.0) < 1.0


def test_fresh_water_has_no_salting_out():
    assert salting_out_factor(GASES["CH4"]) == pytest.approx(1.0)


def test_salting_out_is_monotone_in_salt():
    factors = [
        salting_out_factor(GASES["CH4"], nacl_molarity=c) for c in (0, 1, 2, 3, 4)
    ]
    assert factors == sorted(factors, reverse=True)


def test_brine_solubility_is_below_fresh_water_solubility():
    fresh = henry_solubility_in_brine(GASES["CO2"])
    briny = henry_solubility_in_brine(GASES["CO2"], nacl_molarity=3.59)
    assert briny < fresh


# ─────────────────────────────────────────────────────────────────────
# Partition
# ─────────────────────────────────────────────────────────────────────


def test_partition_coefficient_matches_the_textbook_value_for_carbon_dioxide():
    """The dimensionless Henry constant for carbon dioxide at 25 C is about
    1.2 -- gas-phase concentration slightly above aqueous."""
    assert gas_water_partition_coefficient(GASES["CO2"], 25 * u.deg_C) == pytest.approx(
        1.22, rel=0.05
    )


def test_almost_all_methane_ends_up_in_the_headspace():
    """About 99 percent, so a headspace measurement is nearly the whole story
    for methane."""
    assert headspace_fraction(GASES["CH4"]) == pytest.approx(0.99, abs=0.01)


def test_a_fifth_of_carbon_dioxide_stays_dissolved():
    """Roughly 80 percent in the headspace, so a fifth of what was produced is
    invisible to the detector."""
    assert headspace_fraction(GASES["CO2"]) == pytest.approx(0.80, abs=0.02)


def test_salt_pushes_gas_into_the_headspace():
    fresh = headspace_fraction(GASES["CO2"])
    briny = headspace_fraction(GASES["CO2"], nacl_molarity=3.59)
    assert briny > fresh


def test_headspace_fraction_stays_between_zero_and_one():
    for gas_name in GASES:
        for nacl in (0.0, 1.0, 3.6):
            for mgcl2 in (0.0, 1.0, 1.94):
                fraction = headspace_fraction(
                    GASES[gas_name], nacl_molarity=nacl, mgcl2_molarity=mgcl2
                )
                assert 0.0 < fraction < 1.0


def test_a_bottle_that_is_all_liquid_puts_nothing_in_the_headspace():
    no_headspace = BottleGeometry(liquid_volume=0.125 * u.L, headspace_volume=0.0 * u.L)
    assert headspace_fraction(GASES["CH4"], bottle=no_headspace) == 0.0


# ─────────────────────────────────────────────────────────────────────
# Mass conservation
# ─────────────────────────────────────────────────────────────────────


def test_headspace_and_dissolved_moles_sum_to_the_total():
    """The partition must not create or destroy gas."""
    concentration = 1e-3 * u.mol / u.L
    in_gas = aqueous_concentration_to_headspace_moles(concentration, GASES["CH4"])
    in_liquid = aqueous_concentration_to_dissolved_moles(concentration)
    fraction = headspace_fraction(GASES["CH4"])
    assert float(in_gas / (in_gas + in_liquid)) == pytest.approx(fraction, rel=1e-9)


def test_total_moles_round_trips_from_headspace_moles():
    concentration = 5e-4 * u.mol / u.L
    in_gas = aqueous_concentration_to_headspace_moles(concentration, GASES["CO2"])
    in_liquid = aqueous_concentration_to_dissolved_moles(concentration)
    recovered = headspace_moles_to_total_moles(in_gas, GASES["CO2"])
    assert float(recovered / (in_gas + in_liquid)) == pytest.approx(1.0, rel=1e-9)


def test_a_bare_number_is_read_as_molar():
    with_units = aqueous_concentration_to_headspace_moles(
        1e-3 * u.mol / u.L, GASES["CH4"]
    )
    without_units = aqueous_concentration_to_headspace_moles(1e-3, GASES["CH4"])
    assert with_units.to_value(u.mol) == pytest.approx(without_units.to_value(u.mol))


# ─────────────────────────────────────────────────────────────────────
# The ebullition proxy
# ─────────────────────────────────────────────────────────────────────


def test_degassed_methane_is_added_to_the_headspace():
    """The reaction network has no methane gas phase; it moves dissolved
    methane into Tracer2 to represent bubbles. In a sealed vial that methane
    can only be in the headspace, so it must be counted."""
    dissolved_only = model_headspace_moles(1e-4, GASES["CH4"])
    with_bubbles = model_headspace_moles(
        1e-4, GASES["CH4"], degassed_concentration=1e-4
    )
    assert with_bubbles > dissolved_only


def test_degassed_methane_is_not_passed_through_henry_twice():
    """It has already left solution, so it enters as moles directly. Applying
    the partition to it again would invent methane that never existed."""
    degassed = 1e-4 * u.mol / u.L
    total = model_headspace_moles(0.0, GASES["CH4"], degassed_concentration=degassed)
    expected = (degassed * DEFAULT_BOTTLE.liquid_volume).to(u.mol)
    assert total.to_value(u.mol) == pytest.approx(expected.to_value(u.mol), rel=1e-9)


def test_omitting_the_proxy_leaves_the_result_unchanged():
    """Gases with no ebullition proxy must be unaffected by its existence."""
    without = model_headspace_moles(1e-4, GASES["CO2"]).to_value(u.mol)
    explicit_none = model_headspace_moles(
        1e-4, GASES["CO2"], degassed_concentration=None
    ).to_value(u.mol)
    assert without == pytest.approx(explicit_none)


# ─────────────────────────────────────────────────────────────────────
# Geometry agreement with the deck
# ─────────────────────────────────────────────────────────────────────


def test_bottle_geometry_matches_the_generated_deck():
    """If these drift apart, the model and the measurement describe bottles of
    different sizes and the comparison is meaningless."""
    from pflotran_py.generator import bottle_generator

    assert DEFAULT_BOTTLE.total_volume.to_value(u.L) == pytest.approx(
        bottle_generator.VIAL_VOLUME_L
    )
    assert DEFAULT_BOTTLE.headspace_volume.to_value(u.L) == pytest.approx(
        bottle_generator.HEADSPACE_VOLUME_L
    )
    assert DEFAULT_BOTTLE.liquid_volume.to_value(u.L) == pytest.approx(
        bottle_generator.LIQUID_VOLUME_L
    )


# ─────────────────────────────────────────────────────────────────────
# Composition mapping and validity
# ─────────────────────────────────────────────────────────────────────


def test_composition_maps_sodium_and_magnesium_onto_the_two_salts():
    nacl, mgcl2 = setschenow_salts_from_composition({"Na+": 3.59, "Mg++": 0.0})
    assert nacl == pytest.approx(3.59)
    assert mgcl2 == 0.0


def test_missing_ions_map_to_zero():
    assert setschenow_salts_from_composition({}) == (0.0, 0.0)


def test_the_strongest_magnesium_batch_is_flagged_as_extrapolated():
    """Mg_H sits at 1.94 mol/L magnesium, ionic strength 5.8. Both the ionic
    strength and the divalent coefficient are past where they were fitted, and
    a figure should be able to say so."""
    warnings = validity_warnings(nacl_molarity=0.0, mgcl2_molarity=1.94)
    assert warnings


def test_fresh_water_raises_no_warnings():
    assert validity_warnings() == []


def test_a_moderate_sodium_chloride_brine_raises_no_warnings():
    assert validity_warnings(nacl_molarity=3.59) == []
