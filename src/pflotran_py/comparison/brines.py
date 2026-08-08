"""Measured brine compositions for the sealed-bottle incubations.

A PFLOTRAN deck needs the concentration of each individual ion. The laboratory
records something different: the *mass of each salt* weighed out into a brine,
and the volume that brine was made up to. This module converts one into the
other, so that a deck can be built from what was actually in the bottle rather
than from a guess.

Why this matters more than it might sound. The deck generator can invert a
target water activity into a NaCl concentration, but three of the lowest water
activities in the older exported data would need 7 to 12 mol/kg NaCl, and NaCl
saturates near 6.1 mol/kg. Those batches were never NaCl. Reading the real
recipe shows what they were: magnesium chloride and artificial sea salt, which
reach much lower water activities at attainable concentrations. Building decks
from an inverted water activity alone would have modelled a solution that
cannot exist.

Data source
-----------
Two tabs of the batch-details workbook, published as CSV, no credentials
needed:

* **Brines** -- one row per brine, giving the mass of each salt and the volume
  the brine was made up to.
* **Batches** (one tab per experiment) -- one row per batch, linking a batch to
  its brine by ``Brine ID`` and recording the *measured* water activity, which
  is the quantity the model is ultimately being asked to reproduce.

Dilution into the incubation
----------------------------
A bottle is not filled with brine. It receives a mass of brine plus a mass of
sludge, and the sludge carries its own water, so every ion is diluted. The
laboratory's own convention, inferred from the ``Molarity Salt in Incubations``
column on the Brines tab and reproduced here, is the simple mass ratio::

    dilution = brine_g / (brine_g + sludge_g)

For the Exp004 batches that is 135 / (135 + 45) = 0.75, and applying it to the
brine molarity reproduces the sheet's own incubation molarity exactly for all
six brines that carry one. Using the sheet's convention rather than a more
elaborate volumetric one keeps the model consistent with how the experiment is
recorded, which matters more here than the third decimal place.

Sea salt
--------
Artificial sea salt is a mixture, so unlike NaCl and MgCl2 it has no single
molar mass. It is decomposed into ions using the standard mass fractions of the
major ions in seawater salt (Millero 2013). Those six ions account for about
99.3 percent of the mass; the remainder is bicarbonate, bromide and strontium,
which this model does not track. This decomposition is the least certain step
in the module and is the first thing to check if a sea-salt batch behaves oddly.
"""

import os

import pandas as pd

# ═════════════════════════════════════════════════════════════════════
# Sheet locations
# ═════════════════════════════════════════════════════════════════════
#
# Mirrors saltyBiomass experiments/analysis/incubations/constants.py. Kept as
# literals rather than imported because this repository does not depend on that
# one; if the sheet is ever re-keyed, both copies must move together.

BATCH_WORKBOOK_ID = "12zRu4bYAjCBPrT-hq38L6ESGVg6FcMWl9kv9kGb911I"
BRINES_TAB_GID = "742598617"
BATCH_TAB_GIDS = {
    "Exp003": "1377824489",
    "Exp004": "559062879",
}

# ═════════════════════════════════════════════════════════════════════
# Chemistry
# ═════════════════════════════════════════════════════════════════════

MOLAR_MASS_G_PER_MOL = {
    "NaCl": 58.44,
    "MgCl2*6H2O": 203.30,  # hexahydrate; the anhydrous value would be wrong by 2.1x
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

ION_MOLAR_MASS_G_PER_MOL = {
    "Cl-": 35.45,
    "Na+": 22.99,
    "SO4--": 96.06,
    "Mg++": 24.31,
    "Ca++": 40.08,
    "K+": 39.10,
}

# Column on the Brines tab -> the salt it records. The sheet writes the sulfate
# salt as "NaSO4"; the compound weighed out is sodium sulfate, Na2SO4, which is
# what the molar mass and dissociation above assume.
SALT_COLUMNS = {
    "NaCl (g)": "NaCl",
    "MgCl2*6H2O (g)": "MgCl2*6H2O",
    "NaSO4 (g)": "Na2SO4",
}
SEA_SALT_COLUMN = "Artificial Sea Salt (g)"

ALL_IONS = ["Na+", "Cl-", "Mg++", "SO4--", "Ca++", "K+"]


# ═════════════════════════════════════════════════════════════════════
# Fetching
# ═════════════════════════════════════════════════════════════════════


def _sheet_url(workbook_id, gid):
    return (
        f"https://docs.google.com/spreadsheets/d/{workbook_id}/export"
        f"?gid={gid}&format=csv"
    )


def fetch_brines():
    """Read the Brines tab: one row per brine with its salt masses."""
    brines = pd.read_csv(_sheet_url(BATCH_WORKBOOK_ID, BRINES_TAB_GID))
    return brines.dropna(subset=["Brine ID"])


def fetch_batches(experiment_id):
    """Read one experiment's Batches tab: one row per batch.

    Parameters
    ----------
    experiment_id : str
        Key of :data:`BATCH_TAB_GIDS`, e.g. ``"Exp004"``.
    """
    batches = pd.read_csv(_sheet_url(BATCH_WORKBOOK_ID, BATCH_TAB_GIDS[experiment_id]))
    batches = batches.dropna(subset=["Batch ID"])
    batches["Experiment"] = experiment_id
    return batches


# ═════════════════════════════════════════════════════════════════════
# Composition
# ═════════════════════════════════════════════════════════════════════


def brine_ion_molarities(brine_row):
    """Ion concentrations in a brine, before dilution into the incubation.

    Parameters
    ----------
    brine_row : pandas.Series
        One row of the Brines tab.

    Returns
    -------
    dict
        Ion name to concentration [mol/L]. Ions absent from the brine are
        present with a value of zero, so callers can index without checking.

    Notes
    -----
    Returns all-zero concentrations when the brine has no recorded make-up
    volume. That is the case for the water control (brine ``C``), which is
    correct: it has no salt.
    """
    molarities = {ion: 0.0 for ion in ALL_IONS}

    volume_ml = brine_row.get("Final Volume (mL)")
    if pd.isna(volume_ml) or volume_ml == 0:
        return molarities
    volume_l = volume_ml / 1000.0

    for column, salt in SALT_COLUMNS.items():
        mass_g = brine_row.get(column, 0.0)
        if pd.isna(mass_g) or mass_g == 0:
            continue
        moles = mass_g / MOLAR_MASS_G_PER_MOL[salt]
        for ion, stoichiometry in SALT_DISSOCIATION[salt].items():
            molarities[ion] += stoichiometry * moles / volume_l

    sea_salt_g = brine_row.get(SEA_SALT_COLUMN, 0.0)
    if not pd.isna(sea_salt_g) and sea_salt_g != 0:
        for ion, mass_fraction in SEA_SALT_ION_MASS_FRACTION.items():
            ion_moles = sea_salt_g * mass_fraction / ION_MOLAR_MASS_G_PER_MOL[ion]
            molarities[ion] += ion_moles / volume_l

    return molarities


def incubation_dilution_factor(brine_g, sludge_g):
    """Fraction of the incubation liquid contributed by the brine.

    Uses the laboratory's own mass-ratio convention, which reproduces the
    ``Molarity Salt in Incubations`` column on the Brines tab exactly. See the
    module docstring for why that convention is preferred over a volumetric one.
    """
    total = brine_g + sludge_g
    if total == 0:
        return 0.0
    return brine_g / total


def build_batch_table(experiment_ids=None):
    """One row per batch: measured water activity plus ion concentrations.

    This is the table the deck generator should be driven from. Each row is a
    real bottle condition, with the salt content that was actually weighed out
    and the water activity that was actually measured.

    Parameters
    ----------
    experiment_ids : list of str, optional
        Defaults to every experiment in :data:`BATCH_TAB_GIDS`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``Experiment``, ``Batch ID``, ``Brine ID``, ``Brine Name``,
        ``Measured Water Activity``, ``Measured pH``, ``Dilution``, one column
        per ion in :data:`ALL_IONS` (mol/L in the incubation), and
        ``Ionic Strength`` (mol/L).
    """
    experiment_ids = experiment_ids or list(BATCH_TAB_GIDS)
    brines = fetch_brines().set_index("Brine ID")

    rows = []
    for experiment_id in experiment_ids:
        for _, batch in fetch_batches(experiment_id).iterrows():
            brine_id = batch["Brine ID"]
            if brine_id not in brines.index:
                continue

            brine_molarities = brine_ion_molarities(brines.loc[brine_id])
            dilution = incubation_dilution_factor(
                float(batch.get("Brine (g)", 0) or 0),
                float(batch.get("Sludge (g)", 0) or 0),
            )

            record = {
                "Experiment": experiment_id,
                "Batch ID": int(batch["Batch ID"]),
                "Brine ID": int(brine_id),
                "Brine Name": batch.get("Brine Name"),
                "Measured Water Activity": batch.get("Measured Water Activity"),
                "Measured pH": batch.get("Measured pH"),
                "Dilution": dilution,
            }
            for ion in ALL_IONS:
                record[ion] = brine_molarities[ion] * dilution
            record["Ionic Strength"] = ionic_strength(
                {ion: record[ion] for ion in ALL_IONS}
            )
            rows.append(record)

    return pd.DataFrame(rows)


_ION_CHARGE = {"Na+": 1, "Cl-": -1, "Mg++": 2, "SO4--": -2, "Ca++": 2, "K+": 1}


# ═════════════════════════════════════════════════════════════════════
# Consistency checking
# ═════════════════════════════════════════════════════════════════════


def salt_loading_g_per_ml(brine_row):
    """Total mass of salt per unit volume of brine [g/mL].

    A single number summarising how concentrated a brine is, regardless of
    which salt it was made from. Used only for the consistency check below.
    """
    volume_ml = brine_row.get("Final Volume (mL)")
    if pd.isna(volume_ml) or volume_ml == 0:
        return float("nan")
    total_g = sum(
        (brine_row.get(column, 0.0) or 0.0)
        for column in list(SALT_COLUMNS) + [SEA_SALT_COLUMN]
        if not pd.isna(brine_row.get(column, 0.0))
    )
    return total_g / volume_ml


# Two brines are only compared when they differ by more than this much, so that
# nominally-identical pairs (a brine and its sulfate-spiked twin) are not
# flagged over differences that are really measurement noise.
_LOADING_TOLERANCE_FRACTION = 0.05  # 5 percent of the larger loading
_DENSITY_TOLERANCE_G_PER_ML = 0.005


def check_recipe_consistency(brines=None):
    """Cross-check each recipe against the brine's own measured density.

    Dissolving salt in water makes it denser, always and monotonically. So
    within a family of brines made from the same salt, ranking them by the salt
    loading computed from the recipe must give the same order as ranking them by
    the density that was measured in the laboratory. Where those two orders
    disagree, one of the two recorded numbers is wrong.

    This check exists because the recipe is the input the model is built from.
    A wrong make-up volume does not fail loudly -- it produces a deck at
    entirely the wrong salinity, which then disagrees with the measurements for
    a reason that has nothing to do with the science being tested.

    Returns
    -------
    list of dict
        One entry per detected inconsistency, naming the two brines whose
        density order and recipe order disagree. Empty when every family is
        self-consistent.

    Notes
    -----
    Deliberately reports rather than repairs. Which of the two numbers is wrong
    is a question about what happened at the bench, and belongs to whoever ran
    it.
    """
    brines = fetch_brines() if brines is None else brines

    families = {
        "sea salt": SEA_SALT_COLUMN,
        "NaCl": "NaCl (g)",
        "MgCl2*6H2O": "MgCl2*6H2O (g)",
    }

    problems = []
    for family_name, column in families.items():
        members = []
        for _, row in brines.iterrows():
            mass = row.get(column, 0.0)
            density = row.get("Brine Density (g/mL)")
            if pd.isna(mass) or mass == 0 or pd.isna(density):
                continue
            loading = salt_loading_g_per_ml(row)
            if pd.isna(loading):
                continue
            members.append((row["Brine Name"], loading, float(density)))

        for i, (name_a, loading_a, density_a) in enumerate(members):
            for name_b, loading_b, density_b in members[i + 1 :]:
                loading_gap = abs(loading_a - loading_b)
                density_gap = abs(density_a - density_b)
                if (
                    loading_gap
                    < _LOADING_TOLERANCE_FRACTION * max(loading_a, loading_b)
                    or density_gap < _DENSITY_TOLERANCE_G_PER_ML
                ):
                    continue

                denser = density_a > density_b
                more_concentrated = loading_a > loading_b
                if denser != more_concentrated:
                    problems.append(
                        {
                            "family": family_name,
                            "brine_a": name_a,
                            "brine_b": name_b,
                            "loading_a": loading_a,
                            "loading_b": loading_b,
                            "density_a": density_a,
                            "density_b": density_b,
                            "message": (
                                f"{family_name}: recipe says {name_a} "
                                f"({loading_a:.3f} g/mL) is "
                                f"{'more' if more_concentrated else 'less'} "
                                f"concentrated than {name_b} ({loading_b:.3f} g/mL), "
                                f"but the measured density says the opposite "
                                f"({density_a:.4f} vs {density_b:.4f} g/mL). "
                                f"One of the two recorded values is wrong."
                            ),
                        }
                    )
    return problems


def check_against_water_activity(table, tolerance=0.01):
    """Cross-check derived concentrations against measured water activity.

    This is the stronger of the two checks, because water activity is measured
    on a separate instrument from anything used to make the brine, so it is a
    genuinely independent witness. Adding salt lowers water activity, without
    exception. Within a family of brines made from the same salt, ranking by
    ionic strength must therefore give the reverse of ranking by measured water
    activity.

    Where it does not, believe the water activity. It is a direct measurement
    of the quantity that actually drives the inhibition being modelled, whereas
    the ionic strength is derived from a recipe through several arithmetic steps
    any one of which can carry a transcription error.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`build_batch_table`.
    tolerance : float
        Water activities closer together than this are treated as equal, since
        the batches concerned were intended to match.

    Returns
    -------
    list of dict
        One entry per contradiction, empty when consistent.
    """
    problems = []

    def family_of(brine_name):
        name = str(brine_name)
        if name.startswith("SW"):
            return "sea salt"
        if name.startswith("Na_"):
            return "NaCl"
        if name.startswith("Mg_"):
            return "MgCl2"
        return None

    working = table.copy()
    working["Family"] = working["Brine Name"].map(family_of)

    for family, group in working.dropna(subset=["Family"]).groupby("Family"):
        rows = [
            row
            for _, row in group.dropna(subset=["Measured Water Activity"]).iterrows()
        ]
        for i, row_a in enumerate(rows):
            for row_b in rows[i + 1 :]:
                aw_a = float(row_a["Measured Water Activity"])
                aw_b = float(row_b["Measured Water Activity"])
                strength_a = float(row_a["Ionic Strength"])
                strength_b = float(row_b["Ionic Strength"])

                if abs(aw_a - aw_b) < tolerance:
                    continue
                if abs(strength_a - strength_b) < 1e-6:
                    continue

                saltier = strength_a > strength_b
                drier = aw_a < aw_b
                if saltier != drier:
                    name_a = row_a["Brine Name"]
                    name_b = row_b["Brine Name"]
                    # Name the brine the recipe overstates, so the message points
                    # at the row to go and check.
                    if saltier:
                        overstated, other = name_a, name_b
                        overstated_i, other_i = strength_a, strength_b
                        overstated_aw, other_aw = aw_a, aw_b
                    else:
                        overstated, other = name_b, name_a
                        overstated_i, other_i = strength_b, strength_a
                        overstated_aw, other_aw = aw_b, aw_a
                    problems.append(
                        {
                            "family": family,
                            "brine_a": name_a,
                            "brine_b": name_b,
                            "overstated": overstated,
                            "message": (
                                f"{family}: the recipe makes {overstated} saltier "
                                f"than {other} (ionic strength {overstated_i:.2f} vs "
                                f"{other_i:.2f} mol/L), but {overstated} has the "
                                f"higher measured water activity "
                                f"({overstated_aw:.4f} vs {other_aw:.4f}), so it is "
                                f"really the weaker brine. The recipe for "
                                f"{overstated} overstates its concentration."
                            ),
                        }
                    )
    return problems


def ionic_strength(molarities):
    """Ionic strength, I = 0.5 * sum(c_i * z_i^2), in mol/L.

    Reported alongside each batch because it is the single number that best
    orders these conditions by how far outside seawater-calibrated territory
    they sit.
    """
    return 0.5 * sum(
        concentration * _ION_CHARGE[ion] ** 2
        for ion, concentration in molarities.items()
    )


def to_pflotran_constraints(batch_row, charge_balance_on="Cl-"):
    """Convert a batch row into PFLOTRAN constraint strings.

    Suitable for the ``brine`` argument of
    :class:`~pflotran_py.generator.bottle_generator.BottleGenerator`.

    One ion carries the ``Z`` charge-balance code rather than a fixed value, so
    PFLOTRAN closes the charge balance itself. Chloride is the default choice
    because it is the most abundant anion in every one of these brines and
    because the deck's inherited default already balances on it.

    Ions at zero concentration are omitted, so the generator's inherited trace
    defaults apply instead of an explicit zero, which PFLOTRAN dislikes.
    """
    constraints = {}
    for ion in ALL_IONS:
        concentration = batch_row[ion]
        if concentration <= 0:
            continue
        code = "Z" if ion == charge_balance_on else "T"
        constraints[ion] = f"{concentration:.4e} {code}"
    return constraints


# ═════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════

DEFAULT_OUTPUT = os.path.join("data", "incubation_batch_composition.csv")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Pull measured brine recipes and derive per-batch ion concentrations."
    )
    parser.add_argument(
        "--output", default=DEFAULT_OUTPUT, help="CSV path to write the batch table to."
    )
    parser.add_argument(
        "--experiments",
        nargs="+",
        default=list(BATCH_TAB_GIDS),
        help="Experiment IDs to include.",
    )
    args = parser.parse_args()

    table = build_batch_table(args.experiments)

    density_problems = check_recipe_consistency()
    activity_problems = check_against_water_activity(table)

    if density_problems or activity_problems:
        print("RECIPE CONSISTENCY WARNINGS")
        print("A deck built from an affected brine will sit at the wrong salinity,")
        print("and will then disagree with the measurements for a reason that has")
        print("nothing to do with the science. Resolve at the bench record first.")
        print()
        for problem in density_problems:
            print(f"  [density]  {problem['message']}")
        for problem in activity_problems:
            print(f"  [activity] {problem['message']}")
        print()

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    table.to_csv(args.output, index=False)

    pd.set_option("display.width", 200)
    display_columns = [
        "Experiment",
        "Batch ID",
        "Brine Name",
        "Measured Water Activity",
        "Na+",
        "Cl-",
        "Mg++",
        "SO4--",
        "Ionic Strength",
    ]
    print(table[display_columns].to_string(index=False, float_format="%.4f"))
    print()
    print(f"Wrote {len(table)} batches to {args.output}")


if __name__ == "__main__":
    main()
