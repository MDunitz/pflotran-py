"""Corrections to the recorded brine data, each with its evidence.

Every correction here is explicit and self-scoped: one named function per
finding, applying to named rows only, with the reasoning written out. There is
deliberately no generic "apply this table of fixes" mechanism, because a
correction to laboratory data is a claim about what happened at the bench and
should have to argue for itself in prose.

The pattern follows the ``data_corrections`` package in the saltyBiomass
repository, which handles the equivalent problem for pressure records.
"""

import logging

logger = logging.getLogger(__name__)

# Recorded and inferred sea-salt masses in grams for the 200 mL recipes.
_SEA_SALT_MASS_CORRECTIONS = {
    "SW_M": (200.0, 40.0),
    "SWSu_M": (200.0, 40.0),
    "SW_L": (100.0, 20.0),
    "SWSu_L": (100.0, 20.0),
}


def correct_sea_salt_mass(brines):
    """Apply the inferred sea-salt masses while retaining volume and sulfate.

    Density and measured water activity support a fivefold reduction in
    sea-salt concentration. The SWSu_M stock and top-up recorded in issue
    #42 support the 200 mL interpretation. SWSu_L retains the sheet's
    nominal sulfate mass because the actual addition is unresolved.
    """
    corrected = brines.copy()
    salt_column = "Artificial Sea Salt (g)"
    if salt_column not in corrected:
        return corrected

    for name, (recorded_mass, inferred_mass) in _SEA_SALT_MASS_CORRECTIONS.items():
        match = (
            (corrected["Brine Name"] == name)
            & (corrected["Final Volume (mL)"] == 200.0)
            & (corrected[salt_column] == recorded_mass)
        )
        if match.any():
            corrected.loc[match, salt_column] = inferred_mass
            logger.info(
                "Corrected %s sea-salt mass: %.0f g -> %.0f g at 200 mL",
                name,
                recorded_mass,
                inferred_mass,
            )

    return corrected


# ═════════════════════════════════════════════════════════════════════
# Guard against the next one
# ═════════════════════════════════════════════════════════════════════

# Grams of salt that will dissolve in one gram of water at 20 C. Per family,
# because a single threshold is wrong: magnesium chloride hexahydrate is far
# more soluble than a sea-salt mix, and applying the sea-salt limit to it
# produces a false alarm on a perfectly good brine.
SOLUBILITY_G_PER_G_WATER = {
    "Artificial Sea Salt (g)": 0.46,  # ~78% NaCl mix
    "NaCl (g)": 0.36,
    # Convert anhydrous MgCl2 saturation (0.546 g/g total water) to
    # hydrate mass per gram of added water, accounting for crystal water.
    "MgCl2*6H2O (g)": (0.546 * 203.30 / 95.21) / (1 - 0.546 * 108.09 / 95.21),
}

_SALT_COLUMNS_FOR_MASS = [
    "NaCl (g)",
    "MgCl2*6H2O (g)",
    "Artificial Sea Salt (g)",
    "NaSO4 (g)",
]


def check_mass_balance(brines):
    """Flag recipes that would need more salt dissolved than water can hold.

    The check that would have caught the sea-salt volume error immediately, and
    the one most worth keeping: it needs no calibration, no reference brine and
    no fitted parameter. Given a measured density and a recorded volume, the
    water content follows by arithmetic, and a salt-to-water ratio above
    saturation means one of the two recorded numbers is wrong.

    Returns
    -------
    list of dict
        One entry per impossible recipe; empty when every recipe is physically
        realisable.
    """
    import pandas as pd

    problems = []
    for _, row in brines.iterrows():
        family = next(
            (
                column
                for column in SOLUBILITY_G_PER_G_WATER
                if not pd.isna(row.get(column)) and float(row.get(column) or 0) > 0
            ),
            None,
        )
        if family is None:
            continue

        volume_ml = row.get("Final Volume (mL)")
        density = row.get("Brine Density (g/mL)")
        if pd.isna(volume_ml) or pd.isna(density):
            continue

        salt_g = sum(
            float(row.get(column) or 0)
            for column in _SALT_COLUMNS_FOR_MASS
            if not pd.isna(row.get(column))
        )
        water_g = density * volume_ml - salt_g
        if water_g <= 0:
            ratio = float("inf")
        else:
            ratio = salt_g / water_g

        limit = SOLUBILITY_G_PER_G_WATER[family]
        if ratio > limit:
            problems.append(
                {
                    "brine": row.get("Brine Name"),
                    "family": family.replace(" (g)", ""),
                    "salt_per_g_water": ratio,
                    "solubility_limit": limit,
                    "message": (
                        f"{row.get('Brine Name')}: {salt_g:.1f} g of "
                        f"{family.replace(' (g)', '')} in "
                        f"{volume_ml:.0f} mL at density {density:.4f} g/mL leaves "
                        f"only {water_g:.1f} g of water, a ratio of {ratio:.2f} g "
                        f"salt per g water against a solubility limit of {limit}. "
                        f"This solution cannot be made."
                    ),
                }
            )
    return problems
