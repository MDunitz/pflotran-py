"""Water activity of multi-component brines via the Pitzer ion-interaction model.

Computes the solvent water activity ``a_w`` of an aqueous electrolyte mixture
from ion molalities, using the Pitzer osmotic coefficient. Binary cation-anion
interaction terms only; higher-order cation-cation / anion-anion mixing
(``theta``, ``psi``) is omitted -- a first-order approximation adequate for
placing seawater-derived brines on the a_w scale (accurate to ~1-2% for
Na/K/Mg/Ca/Cl/SO4 brines up to their validated concentration ranges).

Validated against NaCl: a_w = 0.7525 at halite saturation (6.14 mol/kg),
matching the accepted 0.753.

Parameters are 25 degC values (Pitzer, 1991, "Activity Coefficients in
Electrolyte Solutions", CRC). Temperature dependence of A_phi and the
beta/C parameters is not yet implemented (brines are commonly 30-40 degC);
this is the main known limitation.
"""

import math

import astropy.units as u

# Debye-Huckel osmotic-coefficient parameter, 25 degC. Units (kg/mol)**0.5.
A_PHI = 0.3915
# Pitzer universal closest-approach parameter b. Units (kg/mol)**0.5.
B_PITZER = 1.2
# Molar mass of water (constant tracked with astropy units).
M_WATER = 18.0153 * u.g / u.mol

ION_CHARGE = {"Na": 1, "K": 1, "Mg": 2, "Ca": 2, "Cl": -1, "SO4": -2}
CATIONS = ("Na", "K", "Mg", "Ca")
ANIONS = ("Cl", "SO4")

# Pitzer binary parameters per (cation, anion) pair, 25 degC (Pitzer 1991).
#   b0, b1, b2 : Pitzer beta^0, beta^1, beta^2   [kg/mol]
#   cphi       : Pitzer C^phi                    [kg^2/mol^2]
#   a1, a2     : alpha_1, alpha_2                [(kg/mol)^0.5]
# 2:2 pairs (Mg-SO4, Ca-SO4) carry a nonzero b2 with a1=1.4, a2=12.0; all
# other charge types use a single beta^1 term (a1=2.0, b2=0, a2 unused).
PITZER = {
    ("Na", "Cl"): dict(b0=0.07650, b1=0.2664, b2=0.0, cphi=0.00127, a1=2.0, a2=0.0),
    ("K", "Cl"): dict(b0=0.04835, b1=0.2122, b2=0.0, cphi=-0.00084, a1=2.0, a2=0.0),
    ("Mg", "Cl"): dict(b0=0.35235, b1=1.6815, b2=0.0, cphi=0.00519, a1=2.0, a2=0.0),
    ("Ca", "Cl"): dict(b0=0.31590, b1=1.6140, b2=0.0, cphi=-0.00034, a1=2.0, a2=0.0),
    ("Na", "SO4"): dict(b0=0.01958, b1=1.1130, b2=0.0, cphi=0.00497, a1=2.0, a2=0.0),
    ("K", "SO4"): dict(b0=0.04995, b1=0.7793, b2=0.0, cphi=0.0, a1=2.0, a2=0.0),
    ("Mg", "SO4"): dict(b0=0.22100, b1=3.3430, b2=-37.23, cphi=0.0250, a1=1.4, a2=12.0),
    ("Ca", "SO4"): dict(b0=0.20000, b1=3.1973, b2=-54.24, cphi=0.0, a1=1.4, a2=12.0),
}


def _b_phi(pair, sqrt_i):
    """Pitzer B^phi_ca ionic-strength function for one cation-anion pair.

    B^phi_ca = beta^0 + beta^1 * exp(-alpha_1 * sqrt(I))
                       + beta^2 * exp(-alpha_2 * sqrt(I))

    The beta^2 term is present only for 2:2 pairs (b2 != 0).
    """
    p = PITZER[pair]
    val = p["b0"] + p["b1"] * math.exp(-p["a1"] * sqrt_i)
    if p["b2"]:
        val += p["b2"] * math.exp(-p["a2"] * sqrt_i)
    return val


def osmotic_coefficient(molalities):
    """Pitzer osmotic coefficient phi of an electrolyte mixture.

    Pitzer osmotic coefficient (Pitzer 1973, J. Phys. Chem. 77:268), binary
    cation-anion terms only:

        phi = 1 + (2 / sum_i m_i) * [ -A_phi * I**1.5 / (1 + b*sqrt(I))
                  + sum_c sum_a  m_c * m_a * (B^phi_ca + Z * C_ca) ]

    with
        C_ca = C^phi_ca / (2 * sqrt(|z_c * z_a|))
        I    = 0.5 * sum_i m_i * z_i**2        (ionic strength)
        Z    = sum_i m_i * |z_i|
        b    = 1.2

    ``molalities`` maps ion label -> molality Quantity in mol/kg (water).
    """
    m = {ion: q.to_value(u.mol / u.kg) for ion, q in molalities.items()}
    sum_m = sum(m.values())
    ionic_strength = 0.5 * sum(m[i] * ION_CHARGE[i] ** 2 for i in m)
    z_sum = sum(m[i] * abs(ION_CHARGE[i]) for i in m)
    sqrt_i = math.sqrt(ionic_strength)

    debye_huckel = -A_PHI * ionic_strength**1.5 / (1.0 + B_PITZER * sqrt_i)

    pair_sum = 0.0
    for cation in CATIONS:
        if cation not in m:
            continue
        for anion in ANIONS:
            if anion not in m:
                continue
            cphi = PITZER[(cation, anion)]["cphi"]
            c_ca = cphi / (2.0 * math.sqrt(abs(ION_CHARGE[cation] * ION_CHARGE[anion])))
            pair_sum += (
                m[cation] * m[anion] * (_b_phi((cation, anion), sqrt_i) + z_sum * c_ca)
            )

    return 1.0 + (2.0 / sum_m) * (debye_huckel + pair_sum)


def water_activity(molalities):
    """Solvent water activity a_w of a brine from its ion molalities.

    Water activity from osmotic coefficient:

        ln(a_w) = -M_w * phi * sum_i m_i

    where M_w is the molar mass of water [kg/mol], phi the Pitzer osmotic
    coefficient, and sum_i m_i the total molality of all dissolved ions.

    ``molalities`` maps ion label -> molality Quantity in mol/kg (water).
    Returns a dimensionless float in (0, 1].
    """
    phi = osmotic_coefficient(molalities)
    sum_m = sum(q.to_value(u.mol / u.kg) for q in molalities.values())
    m_w = M_WATER.to_value(u.kg / u.mol)
    return math.exp(-m_w * phi * sum_m)
