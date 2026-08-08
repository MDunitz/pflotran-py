"""Build one closed-batch deck per measured incubation batch.

The deck generator can build a bottle at any composition, and the brine module
knows what composition each real batch actually had. This joins the two: one
deck per batch, built from the salt that was weighed out rather than from a
target water activity inverted through an idealised salt.

Why not drive the decks from water activity
-------------------------------------------
It is tempting to take each batch's measured water activity, invert it to a
salt concentration, and build a deck from that. The generator can do this, and
for a quick sweep it is reasonable. It is the wrong thing to do here, for two
reasons.

The inversion assumes a single salt. Half of these batches are magnesium
chloride or artificial sea salt, and at a given water activity those need quite
different ion concentrations from sodium chloride -- at the three matched water
activities in this experiment, the ionic strength needed runs Na < sea salt <
Mg, spanning roughly a factor of 1.6. A deck built from an inverted water
activity would put the right water activity on the wrong solution.

More importantly, water activity is what the model is being asked to *predict*.
PFLOTRAN computes it from the solution composition at each timestep, and the
inhibition sandboxes act on the value it computes. Feeding in a composition
reverse-engineered from the measured water activity would hand the model part
of its own answer. Building from the weighed salt keeps the measured water
activity as an independent check: after a run, the water activity PFLOTRAN
computed can be compared against the one the meter read, and that agreement or
disagreement is informative rather than circular.
"""

import logging
import os

import pandas as pd

from ..generator.bottle_generator import BOTTLE_FINAL_TIME_DAYS, BottleGenerator
from .brines import to_pflotran_constraints

logger = logging.getLogger(__name__)


def deck_filename(batch_row):
    """A filename that identifies a batch without needing the table to hand."""
    return (
        f"{batch_row['Experiment']}_B{int(batch_row['Batch ID']):02d}"
        f"_{batch_row['Brine Name']}.in"
    )


def generate_deck_for_batch(
    batch_row,
    output_dir=".",
    final_time_days=BOTTLE_FINAL_TIME_DAYS,
    **generator_kwargs,
):
    """Build one closed-batch deck for one measured batch.

    Parameters
    ----------
    batch_row : mapping
        A row of the batch composition table from
        :func:`pflotran_py.comparison.brines.build_batch_table`.
    output_dir : str
        Directory to write into; created if absent.
    final_time_days : int
        Simulated duration. The default spans the whole measured record.

    Returns
    -------
    str
        Path to the written deck.
    """
    constraints = to_pflotran_constraints(batch_row)

    water_activity = batch_row.get("Measured Water Activity")
    label = (
        f"{batch_row['Experiment']} batch {int(batch_row['Batch ID'])}, brine "
        f"{batch_row['Brine Name']}, measured water activity "
        f"{water_activity:.4f}, ionic strength "
        f"{batch_row['Ionic Strength']:.2f} mol/L"
    )

    generator = BottleGenerator(
        brine=constraints,
        label=label,
        final_time_days=final_time_days,
        **generator_kwargs,
    )

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, deck_filename(batch_row))
    generator.generate(path)
    return path


def generate_decks_from_batch_table(
    table,
    output_dir="decks",
    final_time_days=BOTTLE_FINAL_TIME_DAYS,
    **generator_kwargs,
):
    """Build one deck per row of the batch composition table.

    Batches whose brine carries no salt at all -- the water controls -- still
    get a deck. They are the most informative single condition in the set,
    because they show what the reaction network does with no salt inhibition of
    any kind, which is the ceiling every salted condition is measured against.

    Returns
    -------
    pandas.DataFrame
        The input table with a ``Deck`` column naming each written file.
    """
    written = []
    for _, batch in table.iterrows():
        path = generate_deck_for_batch(
            batch,
            output_dir=output_dir,
            final_time_days=final_time_days,
            **generator_kwargs,
        )
        written.append(path)

    result = table.copy()
    result["Deck"] = written
    return result


def main():
    import argparse

    from .brines import build_batch_table

    parser = argparse.ArgumentParser(
        description="Generate one closed-batch PFLOTRAN deck per measured incubation batch."
    )
    parser.add_argument(
        "--composition",
        default=os.path.join("data", "incubation_batch_composition.csv"),
        help="Batch composition CSV. Rebuilt from the sheets if absent.",
    )
    parser.add_argument("--output-dir", default="decks")
    parser.add_argument(
        "--final-time-days",
        type=int,
        default=BOTTLE_FINAL_TIME_DAYS,
        help="Simulated duration in days.",
    )
    parser.add_argument(
        "--aw-threshold",
        type=float,
        default=0.5,
        help=(
            "Water activity below which the sandboxes inhibit methanogenesis. "
            "Note that every measured batch sits at 0.824 or above, so at the "
            "default of 0.5 the sandboxes never engage."
        ),
    )
    parser.add_argument(
        "--experiments",
        nargs="+",
        default=None,
        help="Restrict to these experiment IDs.",
    )
    parser.add_argument(
        "--cellulose-hydrolysis",
        action="store_true",
        help=(
            "Hold the substrate carbon in a solid pool that hydrolyses into "
            "solution, instead of as 5 mol/L of dissolved glucose. Stops the "
            "organic pool from setting the water activity."
        ),
    )
    parser.add_argument(
        "--salinity-threshold",
        type=float,
        default=None,
        help=(
            "Chloride concentration [mol/L] at which to centre a sigmoidal "
            "inhibition on the network's methanogenesis reactions. Omit to "
            "leave the network's inhibition as it is."
        ),
    )
    parser.add_argument(
        "--salinity-interval",
        type=float,
        default=0.5,
        help="Width of that transition, in decades.",
    )
    args = parser.parse_args()

    extra = {}
    if args.cellulose_hydrolysis:
        extra["cellulose_hydrolysis"] = {}
    if args.salinity_threshold is not None:
        extra["salinity_inhibition"] = {
            "species": "Cl-",
            "threshold": args.salinity_threshold,
            "interval": args.salinity_interval,
        }

    if os.path.exists(args.composition):
        table = pd.read_csv(args.composition)
    else:
        logger.info("No composition file at %s; rebuilding it.", args.composition)
        table = build_batch_table()

    if args.experiments:
        table = table[table["Experiment"].isin(args.experiments)]

    result = generate_decks_from_batch_table(
        table,
        output_dir=args.output_dir,
        final_time_days=args.final_time_days,
        aw_threshold=args.aw_threshold,
        **extra,
    )

    print()
    print(f"Generated {len(result)} decks in {args.output_dir}/")
    print(f"Each runs to {args.final_time_days} days.")
    print()
    for _, row in result.iterrows():
        print(
            f"  {os.path.basename(row['Deck']):32s} a_w={row['Measured Water Activity']:.4f}  "
            f"I={row['Ionic Strength']:5.2f} mol/L"
        )

    lowest = result["Measured Water Activity"].min()
    if lowest > args.aw_threshold:
        print()
        print(
            f"Note: the driest batch sits at water activity {lowest:.3f}, above the "
            f"sandbox threshold of {args.aw_threshold}. No batch in this set will "
            f"trigger water-activity inhibition; any modelled salt effect comes from "
            f"the chloride Monod term instead."
        )


if __name__ == "__main__":
    main()
