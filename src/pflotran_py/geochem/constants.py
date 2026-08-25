"""Physical constants and ion tables for brine water activity.

Data only -- no calculations. Water activity is computed by PHREEQC with the
``pitzer.dat`` database (Harvie-Weare-Pitzer, Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-
CO2-H2O at 25 degC), so the ion-interaction parameters live in that database.
What is defined here is the ion metadata needed to build a PHREEQC input and to
convert between concentration scales.
"""

import astropy.units as u

# Molar mass of water, for converting between a_w and osmotic coefficient.
M_WATER = 18.0153 * u.g / u.mol

# Signed ionic charge number, for the molal ionic-strength utility.
ION_CHARGE = {"Na": 1, "K": 1, "Mg": 2, "Ca": 2, "Cl": -1, "SO4": -2}

# PFLOTRAN / batch-table ion names -> internal Pitzer labels.
BATCH_ION_TO_PITZER = {
    "Na+": "Na",
    "K+": "K",
    "Mg++": "Mg",
    "Ca++": "Ca",
    "Cl-": "Cl",
    "SO4--": "SO4",
}

# Internal ion label -> PHREEQC master-species element name. Sulfate enters a
# PHREEQC SOLUTION block as the redox element S(6).
PITZER_TO_PHREEQC_ELEMENT = {
    "Na": "Na",
    "K": "K",
    "Mg": "Mg",
    "Ca": "Ca",
    "Cl": "Cl",
    "SO4": "S(6)",
}

# Ionic molar masses, for the density-based molarity -> molality conversion.
ION_MOLAR_MASS = {
    "Na+": 22.9898 * u.g / u.mol,
    "K+": 39.0983 * u.g / u.mol,
    "Mg++": 24.305 * u.g / u.mol,
    "Ca++": 40.078 * u.g / u.mol,
    "Cl-": 35.453 * u.g / u.mol,
    "SO4--": 96.06 * u.g / u.mol,
}

# ═════════════════════════════════════════════════════════════════════
# Hand-rolled binary NaCl Pitzer parameters (legacy a_w <-> molality invert)
# ═════════════════════════════════════════════════════════════════════
#
# Only ``generator.bottle_generator`` uses these today, for the optional
# NaCl-only water-activity inversion that builds thought-experiment decks.
# The primary multi-salt a_w path is PHREEQC ``pitzer.dat`` via
# ``water_activity.py``; do not extend this table for new science work.
#
# Refs: Pitzer & Mayorga (1973) J. Phys. Chem. 77(19), 2300-2308, Table I.

NACL_PITZER_A_PHI = 0.3915  # Debye-Huckel osmotic coefficient, 25 C [kg^0.5/mol^0.5]
NACL_PITZER_B = 1.2  # universal Pitzer constant [kg^0.5/mol^0.5]
NACL_PITZER_ALPHA = 2.0  # universal for 1:1 electrolytes [kg^0.5/mol^0.5]
NACL_BETA0 = 0.0765
NACL_BETA1 = 0.2664
NACL_CPHI = 0.00127
NACL_NU = 2  # ions per formula unit: Na+ and Cl-
# Bare float twin of M_WATER for the bottle inversion (no astropy required there).
NACL_M_WATER_G_PER_MOL = 18.0153
