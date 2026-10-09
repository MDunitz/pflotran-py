"""Tests for deriving ion concentrations from measured brine recipes.

These run offline against fixture tables rather than the live sheets, so they
test the chemistry and the consistency checks rather than the network.
"""

import pandas as pd
import pytest

from pflotran_py.comparison.brines import (
    ALL_IONS,
    brine_ion_molarities,
    check_against_water_activity,
    check_recipe_consistency,
    incubation_dilution_factor,
    ionic_strength,
    salt_loading_g_per_ml,
    to_pflotran_constraints,
)


def _brine(**overrides):
    """A Brines-tab row with everything absent unless named."""
    row = {
        "Brine ID": 1.0,
        "Brine Name": "test",
        "Brine Density (g/mL)": 1.0,
        "Final Volume (mL)": 1000.0,
        "NaCl (g)": 0.0,
        "MgCl2*6H2O (g)": 0.0,
        "Artificial Sea Salt (g)": 0.0,
        "NaSO4 (g)": 0.0,
    }
    row.update(overrides)
    return pd.Series(row)


# ─────────────────────────────────────────────────────────────────────
# Salt dissolution
# ─────────────────────────────────────────────────────────────────────


def test_nacl_gives_equal_sodium_and_chloride():
    molarities = brine_ion_molarities(_brine(**{"NaCl (g)": 58.44}))
    assert molarities["Na+"] == pytest.approx(1.0, rel=1e-3)
    assert molarities["Cl-"] == pytest.approx(1.0, rel=1e-3)
    assert molarities["Mg++"] == 0.0


def test_magnesium_chloride_gives_two_chloride_per_magnesium():
    """MgCl2 dissociates into one Mg2+ and two Cl-."""
    molarities = brine_ion_molarities(_brine(**{"MgCl2*6H2O (g)": 203.30}))
    assert molarities["Mg++"] == pytest.approx(1.0, rel=1e-3)
    assert molarities["Cl-"] == pytest.approx(2.0, rel=1e-3)


def test_hexahydrate_molar_mass_is_used():
    """Using the anhydrous mass (95.2) instead of the hexahydrate (203.3) would
    overstate magnesium by 2.1x, which is the kind of error that silently
    doubles the modelled salt stress."""
    molarities = brine_ion_molarities(_brine(**{"MgCl2*6H2O (g)": 95.21}))
    assert molarities["Mg++"] == pytest.approx(0.468, abs=0.01)


def test_sodium_sulfate_gives_two_sodium_per_sulfate():
    molarities = brine_ion_molarities(_brine(**{"NaSO4 (g)": 142.04}))
    assert molarities["Na+"] == pytest.approx(2.0, rel=1e-3)
    assert molarities["SO4--"] == pytest.approx(1.0, rel=1e-3)


def test_sea_salt_decomposes_into_the_major_seawater_ions():
    molarities = brine_ion_molarities(_brine(**{"Artificial Sea Salt (g)": 35.0}))
    # Roughly one litre of seawater's worth of salt: chloride should dominate,
    # then sodium, and every major ion should be present.
    assert molarities["Cl-"] > molarities["Na+"] > molarities["SO4--"]
    assert molarities["Mg++"] > 0
    assert molarities["Ca++"] > 0
    assert molarities["K+"] > 0
    # Standard seawater chlorinity is about 0.55 M at 35 g/L.
    assert molarities["Cl-"] == pytest.approx(0.543, abs=0.02)


def test_salts_are_additive():
    both = brine_ion_molarities(_brine(**{"NaCl (g)": 58.44, "MgCl2*6H2O (g)": 203.30}))
    assert both["Cl-"] == pytest.approx(3.0, rel=1e-3)
    assert both["Na+"] == pytest.approx(1.0, rel=1e-3)


def test_water_control_has_no_salt():
    """Brine C has no recorded make-up volume, which is correct: it is water."""
    molarities = brine_ion_molarities(_brine(**{"Final Volume (mL)": float("nan")}))
    assert all(value == 0.0 for value in molarities.values())
    assert set(molarities) == set(ALL_IONS)


# ─────────────────────────────────────────────────────────────────────
# Dilution
# ─────────────────────────────────────────────────────────────────────


def test_dilution_matches_the_laboratory_convention():
    """135 g brine into 45 g sludge is the Exp004 recipe, and the sheet's own
    incubation molarity is exactly 0.75 of the brine molarity."""
    assert incubation_dilution_factor(135.0, 45.0) == pytest.approx(0.75)


def test_dilution_reproduces_the_sheet_value_for_the_strongest_brine():
    """Na_H: 140 g NaCl in 500 mL is 4.791 M, and the sheet records 3.593 M in
    the incubation. Getting this wrong scales every modelled concentration."""
    brine = brine_ion_molarities(
        _brine(**{"NaCl (g)": 140.0, "Final Volume (mL)": 500.0})
    )
    assert brine["Na+"] == pytest.approx(4.791, abs=0.01)
    incubation = brine["Na+"] * incubation_dilution_factor(135.0, 45.0)
    assert incubation == pytest.approx(3.593, abs=0.01)


def test_dilution_of_nothing_is_not_a_division_by_zero():
    assert incubation_dilution_factor(0.0, 0.0) == 0.0


# ─────────────────────────────────────────────────────────────────────
# Ionic strength
# ─────────────────────────────────────────────────────────────────────


def test_ionic_strength_of_a_one_one_salt_equals_its_molarity():
    assert ionic_strength({"Na+": 1.0, "Cl-": 1.0}) == pytest.approx(1.0)


def test_divalent_ions_count_four_times():
    """I = 0.5 * sum(c * z^2), so 1 M MgCl2 gives 0.5*(1*4 + 2*1) = 3."""
    assert ionic_strength({"Mg++": 1.0, "Cl-": 2.0}) == pytest.approx(3.0)


# ─────────────────────────────────────────────────────────────────────
# PFLOTRAN constraint rendering
# ─────────────────────────────────────────────────────────────────────


def test_one_ion_carries_the_charge_balance_code():
    row = {ion: 0.0 for ion in ALL_IONS}
    row.update({"Na+": 1.0, "Cl-": 1.0})
    constraints = to_pflotran_constraints(row)
    assert constraints["Cl-"].endswith(" Z")
    assert constraints["Na+"].endswith(" T")


def test_absent_ions_are_omitted_not_zeroed():
    """An explicit zero is worse than an omission: PFLOTRAN dislikes zero
    concentrations, and omitting lets the generator's trace default apply."""
    row = {ion: 0.0 for ion in ALL_IONS}
    row.update({"Na+": 1.0, "Cl-": 1.0})
    constraints = to_pflotran_constraints(row)
    assert "Mg++" not in constraints
    assert "K+" not in constraints


# ─────────────────────────────────────────────────────────────────────
# Consistency checks
# ─────────────────────────────────────────────────────────────────────


def test_consistent_recipes_raise_nothing():
    """Denser brine, more salt. Nothing to report."""
    brines = pd.DataFrame(
        [
            _brine(
                **{
                    "Brine Name": "Na_L",
                    "NaCl (g)": 20.0,
                    "Brine Density (g/mL)": 1.05,
                }
            ),
            _brine(
                **{
                    "Brine Name": "Na_H",
                    "NaCl (g)": 140.0,
                    "Brine Density (g/mL)": 1.18,
                }
            ),
        ]
    )
    assert check_recipe_consistency(brines) == []


def test_density_check_catches_an_inverted_recipe():
    """The lighter brine is recorded as carrying more salt, which cannot be."""
    brines = pd.DataFrame(
        [
            _brine(
                **{
                    "Brine Name": "SW_H",
                    "Artificial Sea Salt (g)": 300.0,
                    "Brine Density (g/mL)": 1.1744,
                }
            ),
            _brine(
                **{
                    "Brine Name": "SW_M",
                    "Artificial Sea Salt (g)": 200.0,
                    "Final Volume (mL)": 200.0,
                    "Brine Density (g/mL)": 1.1097,
                }
            ),
        ]
    )
    problems = check_recipe_consistency(brines)
    assert len(problems) == 1
    assert {problems[0]["brine_a"], problems[0]["brine_b"]} == {"SW_H", "SW_M"}


def test_density_check_tolerates_near_identical_twins():
    """A brine and its sulfate-spiked twin differ by a fraction of a percent;
    flagging that pair would bury the real problems in noise."""
    brines = pd.DataFrame(
        [
            _brine(
                **{
                    "Brine Name": "SW_H",
                    "Artificial Sea Salt (g)": 300.0,
                    "Brine Density (g/mL)": 1.1744,
                }
            ),
            _brine(
                **{
                    "Brine Name": "SWSu_H",
                    "Artificial Sea Salt (g)": 300.0,
                    "NaSO4 (g)": 4.0,
                    "Brine Density (g/mL)": 1.1647,
                }
            ),
        ]
    )
    assert check_recipe_consistency(brines) == []


def test_water_activity_check_catches_an_inverted_recipe():
    """Water activity is measured independently, so it outranks the recipe."""
    table = pd.DataFrame(
        [
            {
                "Brine Name": "SW_H",
                "Measured Water Activity": 0.8664,
                "Ionic Strength": 4.46,
            },
            {
                "Brine Name": "SW_M",
                "Measured Water Activity": 0.9196,
                "Ionic Strength": 14.83,
            },
        ]
    )
    problems = check_against_water_activity(table)
    assert len(problems) == 1
    assert problems[0]["overstated"] == "SW_M"


def test_water_activity_check_passes_a_consistent_series():
    table = pd.DataFrame(
        [
            {
                "Brine Name": "Na_L",
                "Measured Water Activity": 0.962,
                "Ionic Strength": 1.20,
            },
            {
                "Brine Name": "Na_M",
                "Measured Water Activity": 0.903,
                "Ionic Strength": 2.40,
            },
            {
                "Brine Name": "Na_H",
                "Measured Water Activity": 0.845,
                "Ionic Strength": 3.59,
            },
        ]
    )
    assert check_against_water_activity(table) == []


def test_salt_loading_handles_a_missing_volume():
    import math

    assert math.isnan(salt_loading_g_per_ml(_brine(**{"Final Volume (mL)": None})))
