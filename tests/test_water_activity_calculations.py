"""Unit tests for the pure calculation utilities in geochem.water_activity.

The Pitzer osmotic-coefficient terms are now computed inside PHREEQC
(pitzer.dat), so the remaining pure-Python pieces are the molal ionic-strength
utility and the a_w <-> osmotic-coefficient identity.
"""

import math

import astropy.units as u
import pytest

from pflotran_py.geochem.constants import M_WATER
from pflotran_py.geochem.water_activity import (
    ionic_strength,
    osmotic_coefficient,
    water_activity,
)

MK = u.mol / u.kg


def test_ionic_strength_1_1_salt():
    # I = 1/2 (m_Na*1 + m_Cl*1) = m for a 1:1 salt.
    assert ionic_strength({"Na": 2.0, "Cl": 2.0}) == pytest.approx(2.0)


def test_ionic_strength_2_1_salt():
    # MgCl2 at m: I = 1/2 (m*4 + 2m*1) = 3m.
    assert ionic_strength({"Mg": 1.0, "Cl": 2.0}) == pytest.approx(3.0)


def test_osmotic_coefficient_inverts_water_activity():
    # phi = -ln(a_w)/(M_w sum m) must reproduce a_w = exp(-M_w phi sum m).
    molalities = {"Na": 3.0 * MK, "Cl": 3.0 * MK}
    a_w = water_activity(molalities)
    phi = osmotic_coefficient(molalities)
    sum_m = 6.0
    m_w = M_WATER.to_value(u.kg / u.mol)
    assert math.exp(-m_w * phi * sum_m) == pytest.approx(a_w, rel=1e-9)


def test_osmotic_coefficient_concentrated_nacl_above_one():
    # NaCl phi rises well above 1 by ~6 mol/kg (strong positive deviation).
    phi = osmotic_coefficient({"Na": 6.14 * MK, "Cl": 6.14 * MK})
    assert phi > 1.1
