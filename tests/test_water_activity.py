"""Tests for pflotran_py.geochem.water_activity (Pitzer osmotic-coefficient a_w)."""

import astropy.units as u
import pytest

from pflotran_py.geochem.water_activity import water_activity

MK = u.mol / u.kg


def aw(**ions):
    return water_activity({ion: molality * MK for ion, molality in ions.items()})


def test_nacl_halite_saturation():
    # Accepted a_w at NaCl saturation (6.14 mol/kg) is 0.753.
    assert aw(Na=6.14, Cl=6.14) == pytest.approx(0.753, abs=0.002)


def test_dilute_approaches_unity():
    assert aw(Na=0.001, Cl=0.001) == pytest.approx(1.0, abs=1e-4)


def test_seawater():
    sw = dict(Na=0.4861, K=0.0106, Mg=0.0547, Ca=0.0107, Cl=0.5657, SO4=0.0293)
    assert aw(**sw) == pytest.approx(0.981, abs=0.003)


def test_mgcl2_below_nacl_saturation():
    # A 3 mol/kg MgCl2 brine is already below NaCl's saturation a_w (0.753):
    # low a_w is reachable with MgCl2 well before NaCl can precipitate.
    assert aw(Mg=3.0, Cl=6.0) < 0.753


def test_magnesium_lowers_aw_at_matched_ion_molality():
    # At identical total ion molality (sum m_i = 4), the divalent MgCl2 brine
    # depresses a_w more than NaCl -- the ion-specific (Hofmeister) signal that
    # a salinity-only axis would miss.
    nacl = aw(Na=2.0, Cl=2.0)
    mgcl2 = aw(Mg=4.0 / 3.0, Cl=8.0 / 3.0)
    assert mgcl2 < nacl


def test_returns_dimensionless_float():
    value = aw(Na=1.0, Cl=1.0)
    assert isinstance(value, float)
    assert 0.0 < value <= 1.0


def test_pitzer_from_batch_matches_recipe_ions():
    from pflotran_py.geochem.water_activity import pitzer_water_activity_from_batch

    control = {
        "Na+": 0.0,
        "Cl-": 0.0,
        "Mg++": 0.0,
        "SO4--": 0.0,
        "Ca++": 0.0,
        "K+": 0.0,
    }
    assert pitzer_water_activity_from_batch(control) == pytest.approx(1.0)

    # MgCl2-like high salt: Pitzer must sit below the ideal Raoult estimate and
    # closer to the meter reading (~0.82) than 1 - 0.017*sum(c) (~0.90).
    mg_h = {
        "Na+": 0.0,
        "Cl-": 3.87,
        "Mg++": 1.94,
        "SO4--": 0.0,
        "Ca++": 0.0,
        "K+": 0.0,
    }
    pitzer = pitzer_water_activity_from_batch(mg_h)
    ideal = 1.0 - 0.017 * (3.87 + 1.94)
    assert pitzer < ideal
    assert pitzer == pytest.approx(0.82, abs=0.05)


def test_deck_emits_fixed_water_activity(tmp_path):
    from pflotran_py.comparison.decks import generate_deck_for_batch

    batch = {
        "Experiment": "Exp004",
        "Batch ID": 5,
        "Brine Name": "Mg_H",
        "Measured Water Activity": 0.824,
        "Ionic Strength": 5.81,
        "Na+": 0.0,
        "Cl-": 3.873585833743236,
        "Mg++": 1.936792916871618,
        "SO4--": 0.0,
        "Ca++": 0.0,
        "K+": 0.0,
    }
    path = generate_deck_for_batch(
        batch,
        output_dir=str(tmp_path),
        final_time_days=1,
        cellulose_hydrolysis={},
    )
    text = open(path).read()
    assert text.count("FIXED_WATER_ACTIVITY") == 3
