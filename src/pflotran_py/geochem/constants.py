"""Physical constants and ion/Pitzer parameter tables for brine water activity.

Data only -- no calculations. Values are 25 degC (Pitzer, 1991, "Activity
Coefficients in Electrolyte Solutions", CRC). Constants carrying a physical
dimension are astropy Quantities; the Pitzer interaction parameters are kept as
plain indices since they are only ever combined through the dimensionless
osmotic-coefficient algebra, never multiplied through a dimensional calculation.
"""

import astropy.units as u

# Debye-Huckel osmotic-coefficient parameter A_phi, 25 degC. Units (kg/mol)**0.5.
A_PHI = 0.3915
# Pitzer universal closest-approach parameter b. Units (kg/mol)**0.5.
B_PITZER = 1.2
# Molar mass of water.
M_WATER = 18.0153 * u.g / u.mol

# Ionic charge number (signed).
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

# PFLOTRAN / batch-table ion names -> Pitzer labels.
BATCH_ION_TO_PITZER = {
    "Na+": "Na",
    "K+": "K",
    "Mg++": "Mg",
    "Ca++": "Ca",
    "Cl-": "Cl",
    "SO4--": "SO4",
}

# Approximate ionic molar masses for density / molality conversion.
ION_MOLAR_MASS = {
    "Na+": 22.9898 * u.g / u.mol,
    "K+": 39.0983 * u.g / u.mol,
    "Mg++": 24.305 * u.g / u.mol,
    "Ca++": 40.078 * u.g / u.mol,
    "Cl-": 35.453 * u.g / u.mol,
    "SO4--": 96.06 * u.g / u.mol,
}
