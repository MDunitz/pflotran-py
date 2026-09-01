"""Tests for the documented corrections to recorded brine data.

Two things are worth protecting here. The correction itself must stay scoped to
the four brines it argues for and must not drift onto anything else. And the
mass-balance guard must keep catching the class of error that produced it,
since that is the check that would have caught it on day one.
"""

import logging

import pandas as pd
import pytest

from pflotran_py.comparison.corrections import (
    SOLUBILITY_G_PER_G_WATER,
    check_mass_balance,
    correct_sea_salt_makeup_volume,
)

AFFECTED = ["SW_M", "SWSu_M", "SW_L", "SWSu_L"]
UNAFFECTED = ["SW_H", "SWSu_H", "Na_H", "Na_M", "Na_L", "Mg_H", "Mg_M", "Mg_L", "Su_C"]


def _sheet():
    """The Brines tab as recorded, reduced to the columns under test."""
    return pd.DataFrame(
        [
            {"Brine Name": "Su_C", "Final Volume (mL)": 200.0},
            {"Brine Name": "SW_H", "Final Volume (mL)": 1000.0},
            {"Brine Name": "SWSu_H", "Final Volume (mL)": 400.0},
            {"Brine Name": "SW_M", "Final Volume (mL)": 200.0},
            {"Brine Name": "SWSu_M", "Final Volume (mL)": 200.0},
            {"Brine Name": "SW_L", "Final Volume (mL)": 200.0},
            {"Brine Name": "SWSu_L", "Final Volume (mL)": 200.0},
            {"Brine Name": "Na_H", "Final Volume (mL)": 500.0},
            {"Brine Name": "Na_M", "Final Volume (mL)": 200.0},
            {"Brine Name": "Na_L", "Final Volume (mL)": 200.0},
            {"Brine Name": "Mg_H", "Final Volume (mL)": 400.0},
            {"Brine Name": "Mg_M", "Final Volume (mL)": 200.0},
            {"Brine Name": "Mg_L", "Final Volume (mL)": 200.0},
        ]
    )


# ─────────────────────────────────────────────────────────────────────
# Scope of the correction
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("brine_name", AFFECTED)
def test_affected_brines_are_corrected_to_a_litre(brine_name):
    corrected = correct_sea_salt_makeup_volume(_sheet())
    volume = corrected.loc[
        corrected["Brine Name"] == brine_name, "Final Volume (mL)"
    ].iloc[0]
    assert volume == 1000.0


@pytest.mark.parametrize("brine_name", UNAFFECTED)
def test_every_other_brine_is_left_alone(brine_name):
    """Including Na_M, Na_L, Mg_M and Mg_L, which also record 200 mL. The
    correction is scoped by name, not by matching the erroneous value, so a
    brine that legitimately holds 200 mL is never touched."""
    original = _sheet()
    corrected = correct_sea_salt_makeup_volume(original)
    before = original.loc[
        original["Brine Name"] == brine_name, "Final Volume (mL)"
    ].iloc[0]
    after = corrected.loc[
        corrected["Brine Name"] == brine_name, "Final Volume (mL)"
    ].iloc[0]
    assert before == after


def test_correction_does_not_mutate_its_input():
    original = _sheet()
    correct_sea_salt_makeup_volume(original)
    assert (
        original.loc[original["Brine Name"] == "SW_M", "Final Volume (mL)"].iloc[0]
        == 200.0
    )


def test_correction_is_idempotent():
    once = correct_sea_salt_makeup_volume(_sheet())
    twice = correct_sea_salt_makeup_volume(once)
    pd.testing.assert_frame_equal(once, twice)


def test_correction_is_a_no_op_on_a_sheet_fixed_at_source():
    """When the sheet itself is corrected, this should stop firing rather than
    keep announcing a change it did not make."""
    fixed = _sheet()
    fixed.loc[fixed["Brine Name"].isin(AFFECTED), "Final Volume (mL)"] = 1000.0
    with_logging = correct_sea_salt_makeup_volume(fixed)
    pd.testing.assert_frame_equal(with_logging, fixed)


def test_correction_announces_itself(caplog):
    """A silent rewrite of laboratory data would be worse than no correction."""
    with caplog.at_level(logging.INFO):
        correct_sea_salt_makeup_volume(_sheet())
    assert sum("make-up volume" in record.message for record in caplog.records) == len(
        AFFECTED
    )


def test_correction_tolerates_a_missing_brine():
    partial = pd.DataFrame([{"Brine Name": "SW_H", "Final Volume (mL)": 1000.0}])
    assert len(correct_sea_salt_makeup_volume(partial)) == 1


# ─────────────────────────────────────────────────────────────────────
# Mass-balance guard
# ─────────────────────────────────────────────────────────────────────


def _brine_row(name, salt_column, salt_g, volume_ml, density):
    row = {
        "Brine Name": name,
        "Final Volume (mL)": volume_ml,
        "Brine Density (g/mL)": density,
        "NaCl (g)": 0.0,
        "MgCl2*6H2O (g)": 0.0,
        "Artificial Sea Salt (g)": 0.0,
        "NaSO4 (g)": 0.0,
    }
    row[salt_column] = salt_g
    return row


def test_the_original_error_is_caught():
    """SW_M as recorded: 200 g of sea salt in 200 mL at 1.1097 g/mL leaves
    21.9 g of water, a ratio of 9.1 against a limit of 0.46."""
    brines = pd.DataFrame(
        [_brine_row("SW_M", "Artificial Sea Salt (g)", 200.0, 200.0, 1.1097)]
    )
    problems = check_mass_balance(brines)
    assert len(problems) == 1
    assert problems[0]["brine"] == "SW_M"
    assert problems[0]["salt_per_g_water"] == pytest.approx(9.12, abs=0.05)


def test_the_corrected_recipe_passes():
    brines = pd.DataFrame(
        [_brine_row("SW_M", "Artificial Sea Salt (g)", 200.0, 1000.0, 1.1097)]
    )
    assert check_mass_balance(brines) == []


def test_the_undisputed_high_brine_passes():
    """SW_H sits at 0.343 g salt per g water, just under the 0.46 limit, which
    is what a brine named 'high' should look like."""
    brines = pd.DataFrame(
        [_brine_row("SW_H", "Artificial Sea Salt (g)", 300.0, 1000.0, 1.1744)]
    )
    assert check_mass_balance(brines) == []


def test_magnesium_chloride_is_not_judged_by_the_sea_salt_limit():
    """Mg_H sits at 0.79 g salt per g water. That is above the sea-salt limit
    of 0.46 but well under magnesium chloride hexahydrate's 1.17, and a single
    shared threshold would raise a false alarm on a perfectly good brine."""
    brines = pd.DataFrame([_brine_row("Mg_H", "MgCl2*6H2O (g)", 210.0, 400.0, 1.1901)])
    assert check_mass_balance(brines) == []
    assert (
        SOLUBILITY_G_PER_G_WATER["MgCl2*6H2O (g)"]
        > SOLUBILITY_G_PER_G_WATER["Artificial Sea Salt (g)"]
    )


def test_a_negative_water_content_is_impossible_not_a_crash():
    """Density times volume below the salt mass means the recipe implies less
    than no water. It must be reported, not divided by."""
    brines = pd.DataFrame([_brine_row("absurd", "NaCl (g)", 500.0, 100.0, 1.05)])
    problems = check_mass_balance(brines)
    assert len(problems) == 1


def test_rows_without_a_density_are_skipped():
    brines = pd.DataFrame(
        [_brine_row("no_density", "NaCl (g)", 50.0, 200.0, float("nan"))]
    )
    assert check_mass_balance(brines) == []
