"""Water activity of multi-component brines via the Pitzer ion-interaction model.

Computes the solvent water activity ``a_w`` of an aqueous electrolyte mixture
from ion molalities, using the Pitzer osmotic coefficient. Binary cation-anion
interaction terms only; higher-order cation-cation / anion-anion mixing
(``theta``, ``psi``) is omitted -- a first-order approximation adequate for
placing seawater-derived brines on the a_w scale (accurate to ~1-2% for
Na/K/Mg/Ca/Cl/SO4 brines up to their validated concentration ranges).

Validated against NaCl: a_w = 0.7525 at halite saturation (6.14 mol/kg),
matching the accepted 0.753.

Constants and Pitzer parameters live in ``constants.py``; the molarity ->
molality conversion lives in ``conversions.py``. The osmotic-coefficient sum is
decomposed into small pure functions (ionic strength, total charge molality,
B^phi, Debye-Huckel term, pair-interaction sum) so each piece is unit-testable
in isolation.

Temperature dependence of A_phi and the beta/C parameters is not implemented
(parameters are 25 degC); this is the main known limitation.
"""

import math

import astropy.units as u

from .constants import (
    A_PHI,
    ANIONS,
    B_PITZER,
    BATCH_ION_TO_PITZER,
    CATIONS,
    ION_CHARGE,
    M_WATER,
    PITZER,
)
from .conversions import molarities_to_molalities

_MOLAL = u.mol / u.kg


def _as_molal_floats(molalities):
    """Strip a {ion: molality Quantity} mapping to plain floats in mol/kg."""
    return {ion: q.to_value(_MOLAL) for ion, q in molalities.items()}


def ionic_strength(molal):
    """Molal ionic strength of a brine.

    Lewis-Randall (1921):  I = 1/2 * sum_i m_i * z_i**2

    ``molal`` maps ion label -> molality [mol/kg] (plain float). Returns I
    [mol/kg].
    """
    return 0.5 * sum(m * ION_CHARGE[ion] ** 2 for ion, m in molal.items())


def total_charge_molality(molal):
    """Pitzer Z term:  Z = sum_i m_i * |z_i|   [mol/kg]."""
    return sum(m * abs(ION_CHARGE[ion]) for ion, m in molal.items())


def b_phi(pair, sqrt_i):
    """Pitzer B^phi_ca ionic-strength function for one cation-anion pair.

    B^phi_ca = beta^0 + beta^1 * exp(-alpha_1 * sqrt(I))
                       + beta^2 * exp(-alpha_2 * sqrt(I))

    The beta^2 term is present only for 2:2 pairs (b2 != 0). ``sqrt_i`` is
    sqrt(ionic strength).
    """
    p = PITZER[pair]
    val = p["b0"] + p["b1"] * math.exp(-p["a1"] * sqrt_i)
    if p["b2"]:
        val += p["b2"] * math.exp(-p["a2"] * sqrt_i)
    return val


def debye_huckel_osmotic_term(ionic_strength_value):
    """Debye-Huckel contribution to the Pitzer osmotic coefficient.

    -A_phi * I**1.5 / (1 + b * sqrt(I))

    with A_phi the Debye-Huckel osmotic parameter and b the Pitzer
    closest-approach parameter. ``ionic_strength_value`` is I [mol/kg].
    """
    sqrt_i = math.sqrt(ionic_strength_value)
    return -A_PHI * ionic_strength_value**1.5 / (1.0 + B_PITZER * sqrt_i)


def pair_interaction_sum(molal, z_sum, sqrt_i):
    """Cation-anion pair contribution to the Pitzer osmotic coefficient.

    sum_c sum_a  m_c * m_a * (B^phi_ca + Z * C_ca)

    with C_ca = C^phi_ca / (2 * sqrt(|z_c * z_a|)). ``z_sum`` is the Z term and
    ``sqrt_i`` is sqrt(ionic strength).
    """
    total = 0.0
    for cation in CATIONS:
        if cation not in molal:
            continue
        for anion in ANIONS:
            if anion not in molal:
                continue
            cphi = PITZER[(cation, anion)]["cphi"]
            c_ca = cphi / (2.0 * math.sqrt(abs(ION_CHARGE[cation] * ION_CHARGE[anion])))
            total += (
                molal[cation]
                * molal[anion]
                * (b_phi((cation, anion), sqrt_i) + z_sum * c_ca)
            )
    return total


def osmotic_coefficient(molalities):
    """Pitzer osmotic coefficient phi of an electrolyte mixture.

    Pitzer osmotic coefficient (Pitzer 1973, J. Phys. Chem. 77:268), binary
    cation-anion terms only:

        phi = 1 + (2 / sum_i m_i) * [ D-H term + pair-interaction sum ]

    ``molalities`` maps ion label -> molality Quantity [mol/kg].
    """
    molal = _as_molal_floats(molalities)
    sum_m = sum(molal.values())
    strength = ionic_strength(molal)
    z_sum = total_charge_molality(molal)
    sqrt_i = math.sqrt(strength)

    debye_huckel = debye_huckel_osmotic_term(strength)
    pair_sum = pair_interaction_sum(molal, z_sum, sqrt_i)

    return 1.0 + (2.0 / sum_m) * (debye_huckel + pair_sum)


def water_activity(molalities):
    """Solvent water activity a_w of a brine from its ion molalities.

    Water activity from osmotic coefficient:

        ln(a_w) = -M_w * phi * sum_i m_i

    where M_w is the molar mass of water [kg/mol], phi the Pitzer osmotic
    coefficient, and sum_i m_i the total molality of all dissolved ions.

    ``molalities`` maps ion label -> molality Quantity [mol/kg]. Returns a
    dimensionless float in (0, 1].
    """
    phi = osmotic_coefficient(molalities)
    sum_m = sum(q.to_value(_MOLAL) for q in molalities.values())
    m_w = M_WATER.to_value(u.kg / u.mol)
    return math.exp(-m_w * phi * sum_m)


def pitzer_water_activity_from_molarities(molarities):
    """Pitzer a_w from a mapping of PFLOTRAN ion names to molarity [mol/L].

    Returns 1.0 when no salt ions are present.
    """
    molalities = molarities_to_molalities(molarities)
    if not molalities:
        return 1.0
    pitzer_molalities = {
        BATCH_ION_TO_PITZER[ion]: quantity
        for ion, quantity in molalities.items()
        if ion in BATCH_ION_TO_PITZER
    }
    if not pitzer_molalities:
        return 1.0
    return float(water_activity(pitzer_molalities))


def pitzer_water_activity_from_batch(batch_row):
    """Pitzer a_w for one incubation batch composition row."""
    molarities = {
        ion: float(batch_row.get(ion, 0.0) or 0.0) for ion in BATCH_ION_TO_PITZER
    }
    return pitzer_water_activity_from_molarities(molarities)
