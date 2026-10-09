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

# Brines whose recorded make-up volume is contradicted by three independent
# measurements. See correct_sea_salt_makeup_volume for the evidence.
_SEA_SALT_VOLUME_CORRECTION = {
    "SW_M": 1000.0,
    "SWSu_M": 1000.0,
    "SW_L": 1000.0,
    "SWSu_L": 1000.0,
}

_RECORDED_SEA_SALT_VOLUME_ML = 200.0


def correct_sea_salt_makeup_volume(brines):
    """Correct the make-up volume of four sea-salt brines from 200 to 1000 mL.

    The Brines tab records SW_M, SWSu_M, SW_L and SWSu_L as made up to 200 mL.
    Three independent lines of evidence say the volume was 1000 mL and that the
    recorded 200 is a transcription error.

    **1. The recorded recipe is physically impossible.** Mass balance fixes the
    water content of a solution without any fitted parameter: a solution of
    density rho occupying volume V has total mass rho*V, so the water in it is
    rho*V minus the salt. For SW_M that is 1.1097 g/mL times 200 mL minus 200 g
    of salt, leaving 21.9 g of water to dissolve 200 g of sea salt -- a ratio of
    9.1 g of salt per gram of water. A sea-salt mix is roughly 78 percent sodium
    chloride, which alone saturates at 0.36 g per gram of water, so the mix
    saturates near 0.46. The recorded SW_L and SWSu_L recipes are likewise
    impossible at about 0.90. The two undisputed sea-salt brines, SW_H and
    SWSu_H, sit at 0.343 and 0.353, just under saturation, which is what a brine
    named "high" should look like.

    **2. The measured density implies a volume near 1000 mL.** Over this range
    density is linear in concentration, ``rho = rho_water + c*(1 - rho_water*v)``
    where ``v`` is the salt's apparent specific volume. Calibrating ``v`` on the
    two undisputed sea-salt brines gives 0.433 mL/g, and inverting for the four
    disputed ones yields make-up volumes of 1018, 1003, 943 and 905 mL against a
    recorded 200.

    That method was validated on the six brines whose volumes are not in
    question. Calibrated on the most concentrated member of each family, it
    reproduces the others' measured densities to within 0.62 percent for sodium
    chloride and 0.28 percent for magnesium chloride. A method accurate to under
    one percent on brines of known volume, disagreeing by a factor of five on
    these four, is not being defeated by measurement error.

    **3. Independently measured water activity agrees, and the recorded
    recipe contradicts it.** Water activity is measured on a different
    instrument from both the balance and the densitometer. Ranking the six
    sea-salt brines by measured water activity against the density-implied
    concentration gives a Spearman correlation of -0.986: more salt, less
    available water, as it must be. Ranking against the *recorded*
    concentration gives +0.406 -- the wrong sign, which would mean adding salt
    raises water activity.

    The series also contains three matched pairs, each a plain brine and its
    sulfate-spiked twin. Under the corrected volumes the members of each pair
    agree to within 4 percent in density-implied concentration and to within
    0.001 in measured water activity. Two independent instruments agreeing
    pair-by-pair across three pairs is not something a wrong volume produces by
    chance.

    Together the corrected volumes complete the concentration series the
    experiment was evidently designed around: 0.3, 0.2 and 0.1 g/mL, from
    300, 200 and 100 g of salt made up to a common litre.

    What this evidence does *not* settle
    ------------------------------------
    Density and mass balance both constrain the *concentration*, not the split
    between mass and volume that produced it. A recipe of 200 g in 1000 mL and
    one of 40 g in 200 mL are the same solution and cannot be told apart this
    way. The volume is corrected rather than the mass for two reasons: 300, 200
    and 100 g made up to a common volume is the obvious reading of the design,
    and the note recorded against SWSu_M describes adding 133 mL of a sulfate
    sub-solution, which cannot fit inside a 200 mL final volume alongside 200 g
    of salt.

    For the model this distinction does not matter at all -- only the resulting
    concentration enters a PFLOTRAN deck. It matters for the bench record,
    which is why it is stated rather than glossed.

    Parameters
    ----------
    brines : pandas.DataFrame
        The Brines tab as read from the sheet.

    Returns
    -------
    pandas.DataFrame
        A copy with ``Final Volume (mL)`` corrected for the four affected
        brines. Rows are matched by ``Brine Name`` and only rows still carrying
        the erroneous 200 mL are touched, so re-running is a no-op and a
        corrected sheet is left alone.
    """
    corrected = brines.copy()

    for brine_name, volume_ml in _SEA_SALT_VOLUME_CORRECTION.items():
        match = corrected["Brine Name"] == brine_name
        if not match.any():
            continue

        already_right = corrected.loc[match, "Final Volume (mL)"] != (
            _RECORDED_SEA_SALT_VOLUME_ML
        )
        if already_right.all():
            logger.debug(
                "Brine %s already carries a corrected volume; leaving it alone.",
                brine_name,
            )
            continue

        corrected.loc[match, "Final Volume (mL)"] = volume_ml
        logger.info(
            "Corrected %s make-up volume: %.0f mL -> %.0f mL "
            "(recorded volume is contradicted by mass balance, measured "
            "density and measured water activity)",
            brine_name,
            _RECORDED_SEA_SALT_VOLUME_ML,
            volume_ml,
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
    "MgCl2*6H2O (g)": 1.17,  # 54.6 g anhydrous MgCl2 per 100 g water
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
