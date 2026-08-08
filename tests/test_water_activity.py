"""Tests for geochem.water_activity (Pitzer osmotic-coefficient a_w)."""

import astropy.units as u
import pytest

from geochem.water_activity import water_activity

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
