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
