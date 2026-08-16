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

More importantly, the *composition* is what the model is being asked to carry.
Building from the weighed salt keeps the measured water activity as an
independent check against the Pitzer a_w the sandboxes use for inhibition
(see ``FIXED_WATER_ACTIVITY``), and against PFLOTRAN's own ideal-Raoult
estimate if that is left to run.
"""

import logging
import os

import pandas as pd

from ..geochem.water_activity import pitzer_water_activity_from_batch
from ..generator.bottle_generator import BOTTLE_FINAL_TIME_DAYS, BottleGenerator
from ..generator.pflotran_generator import (
    DEFAULT_CELLULOSE_HYDROLYSIS,
    DEFAULT_RATE_CONSTANTS,
)
from .brines import to_pflotran_constraints

logger = logging.getLogger(__name__)


def one_minus_aw_factor(water_activity, a_crit):
    """Same continuous factor the AWINHIBIT sandboxes use: max(0,(a_w-a_crit)/(1-a_crit))."""
    a_w = float(water_activity)
    crit = float(a_crit)
    if crit >= 1.0:
        return 1.0
    if a_w <= crit:
        return 0.0
    return (a_w - crit) / (1.0 - crit)


def _fortran_float(value):
    """Parse a PFLOTRAN-style float string such as ``2.d-7``."""
    return float(str(value).strip().lower().replace("d", "e"))


def _to_fortran_rate(value):
    """Emit a PFLOTRAN-friendly scientific literal."""
    return f"{float(value):.1e}".replace("e", "d")


def _apply_aw_upstream_inhibition(kwargs):
    """Scale fermentation and cellulose hydrolysis by ONE_MINUS_AW on fixed a_w.

    Closed-batch decks pass a constant FIXED_WATER_ACTIVITY, so multiplying the
    upstream rate constants by the same factor the methanogenesis sandboxes use
    is equivalent to inhibiting those steps with a_w -- without a new Fortran
    sandbox. Controls (a_w = 1) are unchanged.
    """
    fixed = kwargs.get("fixed_water_activity")
    if fixed is None:
        return kwargs

    ferment_crit = kwargs.pop("aw_threshold_fermentation", 0.90)
    hydro_crit = kwargs.pop("aw_threshold_hydrolysis", 0.85)
    f_ferment = one_minus_aw_factor(fixed, ferment_crit)
    f_hydro = one_minus_aw_factor(fixed, hydro_crit)

    rates = dict(kwargs.get("rate_constants") or {})
    base_ferment = rates.get(
        "fermentation", DEFAULT_RATE_CONSTANTS["fermentation"]
    )
    rates["fermentation"] = base_ferment * f_ferment
    kwargs["rate_constants"] = rates

    if kwargs.get("cellulose_hydrolysis") is not None:
        hydro = {
            **DEFAULT_CELLULOSE_HYDROLYSIS,
            **dict(kwargs["cellulose_hydrolysis"]),
        }
        base_hydro = _fortran_float(hydro["rate_constant"])
        hydro["rate_constant"] = _to_fortran_rate(base_hydro * f_hydro)
        kwargs["cellulose_hydrolysis"] = hydro

    return kwargs


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
    kwargs = dict(generator_kwargs)
    source = kwargs.pop("water_activity_source", "pitzer")
    if "fixed_water_activity" not in kwargs:
        if source == "pitzer":
            kwargs["fixed_water_activity"] = pitzer_water_activity_from_batch(
                batch_row
            )
        elif source == "measured":
            kwargs["fixed_water_activity"] = float(water_activity)
        elif source == "pflotran":
            kwargs["fixed_water_activity"] = None
        else:
            raise ValueError(
                f"Unknown water_activity_source {source!r}; "
                "expected 'pitzer', 'measured', or 'pflotran'"
            )

    aw_upstream = kwargs.pop("aw_upstream_inhibition", False)
    if not aw_upstream:
        kwargs.pop("aw_threshold_fermentation", None)
        kwargs.pop("aw_threshold_hydrolysis", None)
    else:
        kwargs = _apply_aw_upstream_inhibition(kwargs)

    pitzer_aw = kwargs.get("fixed_water_activity")
    if pitzer_aw is not None:
        label = (
            f"{batch_row['Experiment']} batch {int(batch_row['Batch ID'])}, brine "
            f"{batch_row['Brine Name']}, measured water activity "
            f"{water_activity:.4f}, sandbox a_w {float(pitzer_aw):.4f}, "
            f"ionic strength {batch_row['Ionic Strength']:.2f} mol/L"
        )
    else:
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
        **kwargs,
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
        default=0.91,
        help=(
            "Critical water activity for hydrogenotrophic methanogenesis "
            "(and the fallback if pathway-specific flags are omitted). With "
            "ONE_MINUS_AW this is a_crit where that pathway's rate hits zero."
        ),
    )
    parser.add_argument(
        "--aw-threshold-acetate",
        type=float,
        default=0.92,
        help=(
            "a_crit for acetoclastic methanogenesis. Default 0.92 (most "
            "salt-sensitive of the three); literature ordering, not a fit."
        ),
    )
    parser.add_argument(
        "--aw-threshold-methyl",
        type=float,
        default=0.91,
        help=(
            "a_crit for methylotrophic methanogenesis. Default 0.91 "
            "(with hydrogenotrophic; acetoclastic is higher)."
        ),
    )
    parser.add_argument(
        "--aw-inhibition-type",
        choices=("ONE_MINUS_AW", "SMOOTHSTEP", "THRESHOLD"),
        default="ONE_MINUS_AW",
        help=(
            "Shape of the a_w rate factor. ONE_MINUS_AW is continuous in "
            "(1-a_w); SMOOTHSTEP is a log10 logistic (interval 0.20)."
        ),
    )
    parser.add_argument(
        "--aw-upstream-inhibition",
        action="store_true",
        help=(
            "Also scale fermentation and cellulose hydrolysis by ONE_MINUS_AW "
            "on each batch's fixed a_w (default a_crit 0.85). Osmotic stress "
            "hits the whole cascade, not only methanogens."
        ),
    )
    parser.add_argument(
        "--aw-threshold-fermentation",
        type=float,
        default=0.90,
        help="a_crit for upstream fermentation rate scaling (with --aw-upstream-inhibition).",
    )
    parser.add_argument(
        "--aw-threshold-hydrolysis",
        type=float,
        default=0.85,
        help="a_crit for cellulose hydrolysis rate scaling (with --aw-upstream-inhibition).",
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
    parser.add_argument(
        "--no-cl-inhibition",
        action="store_true",
        help=(
            "Drop the reaction network's own chloride Monod inhibition. "
            "Use when the AWINHIBIT sandboxes own methanogenesis so salt is "
            "keyed on water activity rather than double-counted via Cl-."
        ),
    )
    parser.add_argument(
        "--keep-network-methanogenesis",
        action="store_true",
        help=(
            "Leave the network's three methanogenesis reactions in place "
            "alongside the sandboxes. Attribution only -- with network-rate "
            "sandboxes this double-produces methane."
        ),
    )
    parser.add_argument(
        "--disable-reactions",
        nargs="+",
        default=None,
        metavar="RATE_KEY",
        help=(
            "Omit these reactions from the deck by rate key, e.g. "
            "sulfate_reduction methane_so4_oxidation. For mechanism "
            "attribution only; a deck built this way is not a prediction."
        ),
    )
    parser.add_argument(
        "--use-measured-aw",
        action="store_true",
        help=(
            "Pass each batch's meter-read water activity to the sandboxes "
            "instead of the Pitzer value from the weighed recipe."
        ),
    )
    parser.add_argument(
        "--no-fixed-aw",
        action="store_true",
        help=(
            "Do not pass FIXED_WATER_ACTIVITY; sandboxes use PFLOTRAN's "
            "ideal Raoult a_w. Attribution only."
        ),
    )
    args = parser.parse_args()

    extra = {}
    if args.cellulose_hydrolysis:
        extra["cellulose_hydrolysis"] = {}
        # Bottle controls drift to pH ~7.9 while acetate accumulates. The
        # network's H+_below Monod (Ki = 2.88e-7, half at pH ~6.5) then
        # throttles acetoclastic methanogenesis to ~0.04, which is why control
        # CH4 floors near 1e-4 mol with a full acetate pool. Half-inhibition
        # near pH 7.5 (Ki = 3.16e-8) matches the upper edge of the usual
        # acetoclast optimum without fitting methane.
        extra["thresholds"] = {"h_plus_inhibition_3": 3.16e-8}
        # Acetate half-saturation in the inherited network is 40 mM, far above
        # literature acetoclastic Ks (typically ~0.2-5 mM). After raising
        # hydrolysis, controls bank ~80 mM acetate. Use 2 mM (within the
        # Methanosaeta/Methanosarcina range) so Monod stays near saturation
        # while the pool is drawn down; not a methane-endpoint fit. Keep the
        # network acetoclastic rate constant (1.5e-8); a 3x rate lift overshot
        # control CH4.
        extra["half_saturation"] = {"acetate": 2.0e-3}
    if args.salinity_threshold is not None:
        extra["salinity_inhibition"] = {
            "species": "Cl-",
            "threshold": args.salinity_threshold,
            "interval": args.salinity_interval,
        }
        # Cl- smoothstep attaches to network methanogenesis reactions.
        extra["aw_sandbox_replaces_network_methanogenesis"] = False
    if args.no_cl_inhibition:
        extra["enable_cl_inhibition"] = False
    if args.keep_network_methanogenesis:
        extra["aw_sandbox_replaces_network_methanogenesis"] = False
    if args.disable_reactions:
        extra["disabled_rate_keys"] = set(args.disable_reactions)
    if args.no_fixed_aw:
        extra["water_activity_source"] = "pflotran"
    elif args.use_measured_aw:
        extra["water_activity_source"] = "measured"
    if args.aw_upstream_inhibition:
        extra["aw_upstream_inhibition"] = True
        extra["aw_threshold_fermentation"] = args.aw_threshold_fermentation
        extra["aw_threshold_hydrolysis"] = args.aw_threshold_hydrolysis

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
        aw_threshold_acetate=args.aw_threshold_acetate,
        aw_threshold_methyl=args.aw_threshold_methyl,
        aw_inhibition_type=args.aw_inhibition_type,
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
    if args.aw_inhibition_type == "ONE_MINUS_AW":
        print()
        print(
            f"Sandbox ONE_MINUS_AW a_crit by pathway: "
            f"H2={args.aw_threshold}, methyl={args.aw_threshold_methyl}, "
            f"acetate={args.aw_threshold_acetate} "
            f"(driest measured = {lowest:.3f})."
        )
    elif lowest > args.aw_threshold:
        print()
        print(
            f"Note: the driest batch sits at water activity {lowest:.3f}, above the "
            f"sandbox threshold of {args.aw_threshold}. No batch in this set will "
            f"trigger water-activity inhibition."
        )
    else:
        print()
        print(
            f"Sandbox a_w threshold {args.aw_threshold}: engages for batches at "
            f"or below that water activity (driest measured = {lowest:.3f})."
        )


if __name__ == "__main__":
    main()
