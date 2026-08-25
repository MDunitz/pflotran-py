"""Physical constants and ion tables for brine chemistry / water activity.

Data only -- no calculations. Water activity is computed by PHREEQC with the
``pitzer.dat`` database (Harvie-Weare-Pitzer, Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-
CO2-H2O at 25 degC), so the ion-interaction parameters live in that database.
What is defined here is the ion and salt metadata needed to:

* build a PHREEQC input and convert between concentration scales
* dissolve weighed salts into PFLOTRAN ion molarities (``comparison.brines``)

Sheet locations and column names stay in ``comparison.brines`` (data plumbing).
"""

import astropy.units as u

# Molar mass of water, for converting between a_w and osmotic coefficient.
M_WATER = 18.0153 * u.g / u.mol

# Signed ionic charge number, for the molal ionic-strength utility.
ION_CHARGE = {"Na": 1, "K": 1, "Mg": 2, "Ca": 2, "Cl": -1, "SO4": -2}

# PFLOTRAN / batch-table ion names -> signed charge (ionic strength on molarity).
BATCH_ION_CHARGE = {
    "Na+": 1,
    "Cl-": -1,
    "Mg++": 2,
    "SO4--": -2,
    "Ca++": 2,
    "K+": 1,
}

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

# Bare-float twin for recipe arithmetic that does not use astropy.
ION_MOLAR_MASS_G_PER_MOL = {
    ion: float(mass.to_value(u.g / u.mol)) for ion, mass in ION_MOLAR_MASS.items()
}

# Ions the incubation / PFLOTRAN batch table tracks.
ALL_IONS = ["Na+", "Cl-", "Mg++", "SO4--", "Ca++", "K+"]

# ═════════════════════════════════════════════════════════════════════
# Salt recipes -> ions (comparison.brines)
# ═════════════════════════════════════════════════════════════════════

MOLAR_MASS_G_PER_MOL = {
    "NaCl": 58.44,
    "MgCl2*6H2O": 203.30,  # hexahydrate; anhydrous would be wrong by ~2.1x
    "Na2SO4": 142.04,
}

# Ions each salt releases per formula unit.
SALT_DISSOCIATION = {
    "NaCl": {"Na+": 1, "Cl-": 1},
    "MgCl2*6H2O": {"Mg++": 1, "Cl-": 2},
    "Na2SO4": {"Na+": 2, "SO4--": 1},
}

# Mass fraction of each major ion in seawater salt (Millero 2013, Chemical
# Oceanography 4th ed.). These sum to 0.9927; the balance is bicarbonate,
# bromide and strontium, which the reaction network does not track.
SEA_SALT_ION_MASS_FRACTION = {
    "Cl-": 0.5503,
    "Na+": 0.3059,
    "SO4--": 0.0768,
    "Mg++": 0.0368,
    "Ca++": 0.0118,
    "K+": 0.0111,
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
NACL_M_WATER_G_PER_MOL = float(M_WATER.to_value(u.g / u.mol))
