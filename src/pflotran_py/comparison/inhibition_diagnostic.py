"""Attribute the modelled salt suppression to the terms that produce it.

The comparison figures show modelled methane collapsing more abruptly than the
measurements do, and every salted bottle landing on much the same floor. There
are four candidate causes, and at present they are all switched on at once, so
the figures cannot say which is responsible:

1. **A chloride Monod inhibition baked into the reaction network.** Threshold
   0.2 mol/L, ``INHIBIT_ABOVE_THRESHOLD``, attached to the three methanogenesis
   pathways *and to fermentation*. Every salted batch here carries at least
   1.2 mol/L chloride, so the factor ``K / (K + C)`` runs from about a seventh
   in the weakest brine to a thirtieth in the strongest -- and because it also
   throttles fermentation, that suppression compounds through the carbon chain.
   ``REFERENCES.md`` records the threshold as empirical, with no citation.
2. **A fitted chloride smoothstep.** Added by ``salinity_inhibition``, on the
   same three methanogenesis reactions, deliberately *not* on fermentation. It
   was fitted on Exp004, whose brines carry no sulfate, so it had to absorb the
   entire salt effect of a sulfate-free solution.
3. **Sulfate competition.** Sulfate reducers outcompete methanogens for acetate
   and hydrogen, and sulfate-dependent anaerobic methane oxidation destroys
   methane after it is made. Both are real chemistry and both act only on the
   Exp003 sea-salt brines -- which is how a term fitted on sulfate-free brines
   ends up over-suppressing sea salt when carried across.
4. **Carbon supply.** The solid carbon pool was rescaled to the recipe-derived
   starting inventory (~0.0565 mol C). If hydrolysis is now the binding
   constraint, yield would be flat across brines for reasons having nothing to
   do with salt, which would confound any conclusion about inhibition.

This module runs the five decks that separate those, reads the final headspace
methane out of each, and reports the fold change each term is responsible for.
It fits nothing and introduces no parameter: every variant is the existing deck
with a term switched off.

Variants ``no_salt_terms`` and ``ceiling`` are diagnostics rather than physical
models. Sulfate reduction is something these incubations certainly do, and a
deck without it is not a claim about the world.

Usage::

    python -m pflotran_py.comparison.inhibition_diagnostic          # all stages
    python -m pflotran_py.comparison.inhibition_diagnostic --stage decks
    python -m pflotran_py.comparison.inhibition_diagnostic --stage summarise
"""

import dataclasses
import logging
import os
import shutil

import numpy as np
import pandas as pd

from ..generator.bottle_generator import BOTTLE_FINAL_TIME_DAYS
from .carbon_inventory import default_comparison_starting_carbon_moles
from .decks import deck_filename, generate_decks_from_batch_table
from .run_decks import run_deck

logger = logging.getLogger(__name__)

# The chloride smoothstep as currently fitted, so the variants that keep it
# reproduce the committed figures rather than some other parameterisation.
FITTED_SALINITY_THRESHOLD = 0.75
FITTED_SALINITY_INTERVAL = 1.0

# Sulfate pathways that act on methanogenesis, directly or via its substrate.
# Sulfate reduction diverts acetate away from acetoclastic methanogenesis;
# methane_so4_oxidation is anaerobic methane oxidation, which removes methane
# already produced.
SULFATE_RATE_KEYS = frozenset({"sulfate_reduction", "methane_so4_oxidation"})


@dataclasses.dataclass(frozen=True)
class Variant:
    """One deck configuration, and the question turning it off answers."""

    key: str
    label: str
    question: str
    cl_monod: bool
    cl_smoothstep: bool
    sulfate: bool

    def generator_kwargs(
        self,
        *,
        salinity_threshold=FITTED_SALINITY_THRESHOLD,
        salinity_interval=FITTED_SALINITY_INTERVAL,
    ):
        """Keyword arguments that build this variant's decks."""
        kwargs = {
            "enable_cl_inhibition": self.cl_monod,
            "cellulose_hydrolysis": {},
            # Attribution variants that keep a Cl- smoothstep need the network
            # methanogenesis reactions present to attach it to.
            "aw_sandbox_replaces_network_methanogenesis": False,
        }
        if self.cl_smoothstep:
            kwargs["salinity_inhibition"] = {
                "species": "Cl-",
                "threshold": salinity_threshold,
                "interval": salinity_interval,
            }
        if not self.sulfate:
            kwargs["disabled_rate_keys"] = set(SULFATE_RATE_KEYS)
        return kwargs


VARIANTS = (
    Variant(
        key="baseline",
        label="Baseline (everything on)",
        question="Reproduces the committed comparison figures.",
        cl_monod=True,
        cl_smoothstep=True,
        sulfate=True,
    ),
    Variant(
        key="no_cl_monod",
        label="Legacy Cl- Monod off",
        question="How much suppression comes from the uncited 0.2 M Monod term?",
        cl_monod=False,
        cl_smoothstep=True,
        sulfate=True,
    ),
    Variant(
        key="no_smoothstep",
        label="Fitted Cl- smoothstep off",
        question="How much comes from the term fitted on Exp004?",
        cl_monod=True,
        cl_smoothstep=False,
        sulfate=True,
    ),
    Variant(
        key="no_salt_terms",
        label="Both Cl- terms off",
        question="What does sulfate competition alone do to the sea-salt bottles?",
        cl_monod=False,
        cl_smoothstep=False,
        sulfate=True,
    ),
    Variant(
        key="ceiling",
        label="No Cl- terms, no sulfate pathways",
        question="Is carbon supply now the binding constraint?",
        cl_monod=False,
        cl_smoothstep=False,
        sulfate=False,
    ),
)

VARIANTS_BY_KEY = {variant.key: variant for variant in VARIANTS}

# Fold changes worth naming, each isolating one term by comparing the variant
# that lacks it against the variant that has it and is otherwise identical.
ATTRIBUTIONS = (
    ("cl_monod_fold", "no_cl_monod", "baseline", "legacy Cl- Monod"),
    ("cl_smoothstep_fold", "no_smoothstep", "baseline", "fitted Cl- smoothstep"),
    ("sulfate_fold", "ceiling", "no_salt_terms", "sulfate competition + AOM"),
)


# ═════════════════════════════════════════════════════════════════════
# Stage 1 -- decks
# ═════════════════════════════════════════════════════════════════════


def generate_variant_decks(
    table,
    variant,
    deck_root,
    *,
    final_time_days=BOTTLE_FINAL_TIME_DAYS,
    aw_threshold=0.5,
    **kwargs,
):
    """Build one deck per batch for a single variant."""
    output_dir = os.path.join(deck_root, variant.key)
    return generate_decks_from_batch_table(
        table,
        output_dir=output_dir,
        final_time_days=final_time_days,
        aw_threshold=aw_threshold,
        **variant.generator_kwargs(**kwargs),
    )


def generate_all_decks(table, deck_root, *, variants=VARIANTS, **kwargs):
    """Build every variant's decks. Returns ``{variant key: deck count}``."""
    counts = {}
    for variant in variants:
        result = generate_variant_decks(table, variant, deck_root, **kwargs)
        counts[variant.key] = len(result)
        logger.info("%s: %d decks", variant.key, len(result))
    return counts


# ═════════════════════════════════════════════════════════════════════
# Stage 2 -- runs
# ═════════════════════════════════════════════════════════════════════


def run_all_variants(
    deck_root,
    run_root,
    repo_root,
    *,
    variants=VARIANTS,
    timeout_seconds=3600,
    clean=False,
):
    """Run every variant's decks, keeping each variant in its own subtree.

    Sequential, for the reason ``run_decks.run_all`` gives: each run is
    chemistry-heavy, and running several at once mostly buys memory pressure.
    A failure is recorded and the set continues, so one non-convergent deck
    does not cost the whole diagnostic.
    """
    if clean and os.path.isdir(run_root):
        shutil.rmtree(run_root)

    results = []
    for variant in variants:
        variant_decks = os.path.join(deck_root, variant.key)
        if not os.path.isdir(variant_decks):
            logger.warning("No decks for %s at %s; skipping", variant.key, variant_decks)
            continue

        decks = sorted(
            os.path.join(variant_decks, name)
            for name in os.listdir(variant_decks)
            if name.endswith(".in")
        )
        for deck in decks:
            result = run_deck(
                deck,
                os.path.join(run_root, variant.key),
                repo_root,
                timeout_seconds=timeout_seconds,
            )
            result["variant"] = variant.key
            results.append(result)
            status = "ok" if result["returncode"] == 0 else "FAILED"
            logger.info("  %-12s %-32s %s", variant.key, result["name"], status)
    return results


# ═════════════════════════════════════════════════════════════════════
# Stage 3 -- read the runs
# ═════════════════════════════════════════════════════════════════════


def final_methane_moles(run_dir, batch_row):
    """Headspace methane at the end of one run, in moles, or None."""
    from .figures import load_model_run, model_methane_headspace

    run_frame = load_model_run(run_dir)
    if run_frame is None or run_frame.empty:
        return None
    _, moles = model_methane_headspace(run_frame, batch_row)
    if len(moles) == 0:
        return None
    return float(moles[-1])


def collect(table, run_root, *, variants=VARIANTS, starting_carbon_moles=None):
    """Final modelled methane for every (variant, batch) pair.

    Returns
    -------
    pandas.DataFrame
        One row per variant per batch, with ``ch4_moles`` and
        ``ch4_per_starting_c``. Missing runs are dropped with a warning rather
        than raising, so a partially complete set can still be inspected.
    """
    if starting_carbon_moles is None:
        starting_carbon_moles = default_comparison_starting_carbon_moles()

    rows = []
    for variant in variants:
        for _, batch in table.iterrows():
            name = os.path.splitext(deck_filename(batch))[0]
            run_dir = os.path.join(run_root, variant.key, name)
            if not os.path.isdir(run_dir):
                logger.warning("Missing run: %s", run_dir)
                continue
            moles = final_methane_moles(run_dir, batch)
            if moles is None:
                logger.warning("No readable output: %s", run_dir)
                continue
            rows.append(
                {
                    "variant": variant.key,
                    "variant_label": variant.label,
                    "batch_name": name,
                    "experiment": batch["Experiment"],
                    "batch_id": int(batch["Batch ID"]),
                    "brine": batch["Brine Name"],
                    "water_activity": batch["Measured Water Activity"],
                    "ionic_strength": batch["Ionic Strength"],
                    "chloride": batch.get("Cl-", np.nan),
                    "sulfate": batch.get("SO4--", np.nan),
                    "ch4_moles": moles,
                    "ch4_per_starting_c": moles / starting_carbon_moles,
                }
            )
    return pd.DataFrame(rows)


def attribute(collected):
    """Fold change each term is responsible for, per batch.

    Inhibition factors multiply, so a fold change is the natural unit: a value
    of 30 for ``cl_monod_fold`` means switching that term off multiplies
    modelled methane by thirty. Values are reported per batch rather than
    averaged, because the whole question is whether a term's effect has the
    right *shape* across water activity, not whether it has the right size
    somewhere.
    """
    if collected.empty:
        return pd.DataFrame()

    wide = collected.pivot_table(
        index=[
            "batch_name",
            "experiment",
            "batch_id",
            "brine",
            "water_activity",
            "ionic_strength",
        ],
        columns="variant",
        values="ch4_moles",
    ).reset_index()

    for column, without, with_, _ in ATTRIBUTIONS:
        if without in wide.columns and with_ in wide.columns:
            wide[column] = wide[without] / wide[with_].replace(0.0, np.nan)

    return wide.sort_values(["experiment", "batch_id"])


def measured_suppression(table, ecsv_glob, molecule="CH4_FID"):
    """Measured final methane per batch, and its ratio to the experiment control.

    The control is the salt-free bottle in the same experiment, which is the
    ceiling every salted condition is measured against. Returned separately
    from the model so that a measured/modelled comparison of *suppression*
    (rather than of absolute moles) is possible even where the absolute
    inventory is uncertain.
    """
    from .figures import load_measured

    measured = load_measured(ecsv_glob, molecule)
    rows = []
    for _, batch in table.iterrows():
        points = measured[
            (measured["Experiment"] == batch["Experiment"])
            & (measured["Batch ID"] == batch["Batch ID"])
        ]
        if points.empty:
            continue
        latest = points.sort_values("Days since start")
        final_day = latest["Days since start"].iloc[-1]
        at_end = latest[latest["Days since start"] == final_day]
        rows.append(
            {
                "experiment": batch["Experiment"],
                "batch_id": int(batch["Batch ID"]),
                "brine": batch["Brine Name"],
                "water_activity": batch["Measured Water Activity"],
                "measured_ch4_moles": float(at_end["Cumulative Moles"].mean()),
                "measured_day": float(final_day),
            }
        )

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame

    # The salt-free control is named "C". "Su_C" is a separate substrate
    # control and is not the salt-free reference, so the match is exact.
    controls = (
        frame[frame["brine"] == "C"]
        .groupby("experiment")["measured_ch4_moles"]
        .max()
    )
    frame["measured_suppression"] = frame.apply(
        lambda row: (
            controls.get(row["experiment"], np.nan) / row["measured_ch4_moles"]
            if row["measured_ch4_moles"] > 0
            else np.nan
        ),
        axis=1,
    )
    return frame


# ═════════════════════════════════════════════════════════════════════
# Stage 4 -- figure
# ═════════════════════════════════════════════════════════════════════


def plot_decomposition(
    collected, output_path, measured=None, starting_carbon_moles=None
):
    """Two panels: yield against water activity per variant, and attribution.

    The left panel is the diagnostic's point. If one variant's curve has the
    measured *shape* -- a gradual decline rather than a cliff onto a floor --
    then the terms absent from that variant are what produced the cliff. The
    right panel says how large each term's effect is per batch, which is what
    decides whether a term is doing real work or is redundant.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if starting_carbon_moles is None:
        starting_carbon_moles = default_comparison_starting_carbon_moles()

    figure, (left, right) = plt.subplots(1, 2, figsize=(15, 6.5))

    # Two variants can land on top of each other -- that is itself a result,
    # meaning the term between them does nothing -- so the styles are varied
    # and the lines drawn semi-transparent to keep a hidden line visible.
    colours = plt.get_cmap("viridis")(np.linspace(0.05, 0.8, len(VARIANTS)))
    styles = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
    markers = ["o", "s", "^", "D", "v"]
    for index, (colour, variant) in enumerate(zip(colours, VARIANTS)):
        subset = collected[collected["variant"] == variant.key]
        if subset.empty:
            continue
        subset = subset.sort_values("water_activity")
        left.plot(
            subset["water_activity"],
            subset["ch4_per_starting_c"],
            marker=markers[index % len(markers)],
            markersize=5,
            linewidth=1.8,
            linestyle=styles[index % len(styles)],
            color=colour,
            alpha=0.8,
            label=variant.label,
        )

    if measured is not None and not measured.empty:
        points = measured.dropna(subset=["measured_ch4_moles"])
        left.scatter(
            points["water_activity"],
            points["measured_ch4_moles"] / starting_carbon_moles,
            marker="x",
            s=70,
            color="crimson",
            zorder=5,
            label="measured",
        )

    left.set_yscale("log")
    left.invert_xaxis()
    left.set_xlabel("Measured water activity")
    left.set_ylabel("Final CH$_4$ [mol per mol starting C]")
    left.set_title("Yield against water activity, one line per variant")
    left.grid(alpha=0.3, which="both")
    # The x axis is inverted, so the upper right corner is low water activity
    # with high yield -- a region nothing occupies.
    left.legend(fontsize=8, loc="upper right")

    attribution = attribute(collected)
    fold_columns = [
        (column, description)
        for column, _, _, description in ATTRIBUTIONS
        if column in attribution.columns
    ]
    if fold_columns and not attribution.empty:
        labels = [
            f"{row.experiment[-3:]}B{row.batch_id:02d} {row.brine}"
            for row in attribution.itertuples()
        ]
        positions = np.arange(len(labels))
        width = 0.8 / max(len(fold_columns), 1)
        for index, (column, description) in enumerate(fold_columns):
            right.bar(
                positions + index * width,
                attribution[column].to_numpy(),
                width=width,
                label=description,
            )
        right.set_xticks(positions + width * (len(fold_columns) - 1) / 2)
        right.set_xticklabels(labels, rotation=90, fontsize=7)
        right.axhline(1.0, color="black", linewidth=0.8)
        right.set_yscale("log")
        right.set_ylabel("Fold change in CH$_4$ when the term is removed")
        right.set_title("What each term contributes, per batch")
        right.grid(axis="y", alpha=0.3, which="both")
        right.legend(fontsize=8)

    figure.suptitle(
        "Where the modelled salt suppression comes from -- "
        f"starting C = {starting_carbon_moles:.4f} mol per bottle",
        fontsize=11,
    )
    figure.text(
        0.5,
        0.005,
        "Diagnostic only. The two right-hand variants omit sulfate reduction and "
        "anaerobic methane oxidation, which these incubations certainly perform; "
        "they bound the carbon-supply ceiling and are not predictions.",
        ha="center",
        fontsize=7.5,
        style="italic",
    )
    figure.tight_layout(rect=(0, 0.03, 1, 0.96))
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Command line
# ═════════════════════════════════════════════════════════════════════


def summarise(
    table,
    run_root,
    output_dir,
    *,
    ecsv_glob=None,
    starting_carbon_moles=None,
):
    """Read the runs, write the tables and the figure, return the paths."""
    if starting_carbon_moles is None:
        starting_carbon_moles = default_comparison_starting_carbon_moles()

    collected = collect(table, run_root, starting_carbon_moles=starting_carbon_moles)
    if collected.empty:
        raise SystemExit(
            f"No readable runs under {run_root}/. Run the runs stage first."
        )

    measured = None
    if ecsv_glob:
        try:
            measured = measured_suppression(table, ecsv_glob)
        except FileNotFoundError:
            logger.warning("No measured tables matched %s; plotting model only", ecsv_glob)

    os.makedirs(output_dir, exist_ok=True)
    written = {}

    written["per_variant"] = os.path.join(output_dir, "per_variant.csv")
    collected.to_csv(written["per_variant"], index=False)

    attribution = attribute(collected)
    written["attribution"] = os.path.join(output_dir, "attribution.csv")
    attribution.to_csv(written["attribution"], index=False)

    if measured is not None and not measured.empty:
        written["measured"] = os.path.join(output_dir, "measured.csv")
        measured.to_csv(written["measured"], index=False)

    written["figure"] = plot_decomposition(
        collected,
        os.path.join(output_dir, "inhibition_decomposition.png"),
        measured=measured,
        starting_carbon_moles=starting_carbon_moles,
    )
    return written, collected, attribution


def _print_attribution(attribution):
    """Human-readable summary of the fold changes, for the terminal."""
    fold_columns = [
        (column, description)
        for column, _, _, description in ATTRIBUTIONS
        if column in attribution.columns
    ]
    if not fold_columns:
        print("Not enough variants completed to attribute anything.")
        return

    print()
    print("Fold change in final methane when a term is removed")
    print("(a value of 30 means removing that term multiplies methane by thirty)")
    print()
    header = f"  {'batch':22s} {'a_w':>6s}"
    for _, description in fold_columns:
        header += f" {description[:22]:>24s}"
    print(header)
    for row in attribution.itertuples():
        line = f"  {row.batch_name[:22]:22s} {row.water_activity:6.3f}"
        for column, _ in fold_columns:
            value = getattr(row, column, np.nan)
            line += f" {value:24.1f}" if np.isfinite(value) else f" {'--':>24s}"
        print(line)

    # Controls carry no salt, so a salt term cannot act on them; including
    # them in the summary would dilute the effect toward one.
    salted = attribution[attribution["water_activity"] < 0.99]
    print()
    print("Across the salted batches only:")
    for column, description in fold_columns:
        values = salted[column].dropna()
        if values.empty:
            continue
        print(
            f"  {description:34s} median {values.median():8.1f}x   "
            f"range {values.min():.1f}x to {values.max():.1f}x"
        )


def main():
    import argparse

    from .brines import build_batch_table

    parser = argparse.ArgumentParser(
        description=(
            "Decompose the modelled salt suppression into the terms that cause "
            "it, by running the same decks with one term switched off at a time."
        )
    )
    parser.add_argument(
        "--composition",
        default=os.path.join("data", "incubation_batch_composition.csv"),
        help="Batch composition CSV. Rebuilt from the sheets if absent.",
    )
    parser.add_argument("--deck-root", default=os.path.join("decks", "diagnostic"))
    parser.add_argument("--run-root", default=os.path.join("runs", "diagnostic"))
    parser.add_argument(
        "--output-dir", default=os.path.join("output", "comparison", "diagnostic")
    )
    parser.add_argument("--repo-root", default=os.getcwd())
    parser.add_argument(
        "--stage",
        choices=("decks", "runs", "summarise", "all"),
        default="all",
        help="Run one stage only. The runs stage needs a working container.",
    )
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=sorted(VARIANTS_BY_KEY),
        default=None,
        help="Restrict to these variants.",
    )
    parser.add_argument("--experiments", nargs="+", default=None)
    parser.add_argument(
        "--final-time-days", type=int, default=BOTTLE_FINAL_TIME_DAYS
    )
    parser.add_argument("--aw-threshold", type=float, default=0.5)
    parser.add_argument(
        "--salinity-threshold", type=float, default=FITTED_SALINITY_THRESHOLD
    )
    parser.add_argument(
        "--salinity-interval", type=float, default=FITTED_SALINITY_INTERVAL
    )
    parser.add_argument(
        "--ecsv-glob",
        default=os.path.expanduser(
            "~/Documents/GitHub/saltyBiomass/data/transformed/*Exp00[34]*.ecsv"
        ),
        help="Measured pipeline tables, overlaid on the figure.",
    )
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if os.path.exists(args.composition):
        table = pd.read_csv(args.composition)
    else:
        logger.info("No composition file at %s; rebuilding it.", args.composition)
        table = build_batch_table()
    if args.experiments:
        table = table[table["Experiment"].isin(args.experiments)]

    variants = (
        tuple(VARIANTS_BY_KEY[key] for key in args.variants)
        if args.variants
        else VARIANTS
    )

    if args.stage in ("decks", "all"):
        print(f"==> Generating decks for {len(variants)} variants")
        counts = generate_all_decks(
            table,
            args.deck_root,
            variants=variants,
            final_time_days=args.final_time_days,
            aw_threshold=args.aw_threshold,
            salinity_threshold=args.salinity_threshold,
            salinity_interval=args.salinity_interval,
        )
        total = sum(counts.values())
        print(f"    {total} decks under {args.deck_root}/")

    if args.stage in ("runs", "all"):
        print("==> Running decks (this is the slow part)")
        results = run_all_variants(
            args.deck_root,
            args.run_root,
            args.repo_root,
            variants=variants,
            timeout_seconds=args.timeout,
            clean=args.clean,
        )
        failed = [r for r in results if r["returncode"] != 0]
        print(f"    {len(results) - len(failed)} of {len(results)} runs completed.")
        for result in failed:
            print(f"    FAILED {result['variant']}/{result['name']}")

    if args.stage in ("summarise", "all"):
        print("==> Reading runs and attributing")
        written, collected, attribution = summarise(
            table,
            args.run_root,
            args.output_dir,
            ecsv_glob=args.ecsv_glob,
        )
        _print_attribution(attribution)
        print()
        for label, path in written.items():
            print(f"  wrote {path}")


if __name__ == "__main__":
    main()
