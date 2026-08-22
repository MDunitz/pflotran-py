"""Unit tests for the extracted Pitzer osmotic-coefficient calculation pieces.

These cover the individual functions that ``osmotic_coefficient`` composes, so a
regression can be localised to a single term rather than only observed at the
a_w level.
"""

import math

import astropy.units as u
import pytest

from pflotran_py.geochem.constants import A_PHI, B_PITZER, PITZER
from pflotran_py.geochem.water_activity import (
    b_phi,
    debye_huckel_osmotic_term,
    ionic_strength,
    osmotic_coefficient,
    pair_interaction_sum,
    total_charge_molality,
)

MK = u.mol / u.kg


def test_ionic_strength_1_1_salt():
    # I = 1/2 (m_Na*1 + m_Cl*1) = m for a 1:1 salt.
    assert ionic_strength({"Na": 2.0, "Cl": 2.0}) == pytest.approx(2.0)


def test_ionic_strength_2_1_salt():
    # MgCl2 at m: I = 1/2 (m*4 + 2m*1) = 3m.
    assert ionic_strength({"Mg": 1.0, "Cl": 2.0}) == pytest.approx(3.0)


def test_total_charge_molality():
    # Z = sum m_i |z_i| = 1*2 + 2*1 = 4 for MgCl2 at m=1.
    assert total_charge_molality({"Mg": 1.0, "Cl": 2.0}) == pytest.approx(4.0)


def test_debye_huckel_term_matches_closed_form():
    strength = 3.0
    expected = -A_PHI * strength**1.5 / (1.0 + B_PITZER * math.sqrt(strength))
    assert debye_huckel_osmotic_term(strength) == pytest.approx(expected)


def test_debye_huckel_term_is_negative():
    assert debye_huckel_osmotic_term(1.0) < 0.0


def test_b_phi_at_zero_ionic_strength_sums_betas():
    # At sqrt(I)=0 every exponential is 1, so B^phi = b0 + b1 (+ b2 for 2:2).
    p = PITZER[("Na", "Cl")]
    assert b_phi(("Na", "Cl"), 0.0) == pytest.approx(p["b0"] + p["b1"])
    q = PITZER[("Mg", "SO4")]
    assert b_phi(("Mg", "SO4"), 0.0) == pytest.approx(q["b0"] + q["b1"] + q["b2"])


def test_b_phi_decays_with_ionic_strength():
    # The beta^1 exponential decays, so B^phi drops toward b0 as I grows.
    p = PITZER[("Na", "Cl")]
    assert b_phi(("Na", "Cl"), 3.0) < b_phi(("Na", "Cl"), 0.0)
    assert b_phi(("Na", "Cl"), 3.0) > p["b0"]


def test_pair_interaction_sum_only_counts_present_pairs():
    # A pure-NaCl brine has a single (Na, Cl) pair contribution; adding a cation
    # with no partner anion present contributes nothing on its own.
    molal = {"Na": 1.0, "Cl": 1.0}
    strength = ionic_strength(molal)
    z_sum = total_charge_molality(molal)
    got = pair_interaction_sum(molal, z_sum, math.sqrt(strength))
    cphi = PITZER[("Na", "Cl")]["cphi"]
    c_ca = cphi / 2.0
    expected = 1.0 * 1.0 * (b_phi(("Na", "Cl"), math.sqrt(strength)) + z_sum * c_ca)
    assert got == pytest.approx(expected)


def test_osmotic_coefficient_dilute_near_one_from_below():
    # The Debye-Huckel term makes phi dip just below 1 at low ionic strength
    # (phi has a minimum before rising), so it is near 1 but not above it.
    phi = osmotic_coefficient({"Na": 0.001 * MK, "Cl": 0.001 * MK})
    assert phi == pytest.approx(1.0, abs=2e-2)
    assert phi < 1.0
