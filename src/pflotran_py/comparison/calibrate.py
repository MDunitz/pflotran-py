"""Fit the salinity inhibition on one experiment and test it on another.

RETIRED FROM THE DEFAULT COMPARISON
-----------------------------------
This module fits a chloride smoothstep (``salinity_inhibition``) that the live
Exp003/Exp004 comparison no longer uses. Salt stress there is owned by the
AWINHIBIT sandboxes (meter a_w × ``AW_CRIT_*``). Keep this file for historical
replay and attribution; do not treat its grids or winners as current science
defaults. See ``inhibition_diagnostic`` and the README "What the comparison
currently shows" section.

Two parameters of the salinity inhibition term -- the chloride concentration at
which it centres, and the width of the transition -- have no independent source.
They were chosen by comparing modelled methane against measured methane, which
makes every agreement statistic computed on the same data circular.

This module separates the two steps so that the second one means something. The
grid of candidate parameters and the objective are both fixed below, before any
run happens. The winner is selected mechanically on the fitting experiment
alone. Only then is it applied to the held-out experiment, and the result is
whatever it is.

What this is and is not
-----------------------
It is a pre-registered protocol: the objective, the grid and the split are
written down in advance and the selection is arithmetic rather than judgement.

It is not a blind prediction. The held-out experiment was examined during the
work that led here, so its behaviour is known to the person running this. That
knowledge cannot be unlearned, and it can leak into choices as ordinary as which
grid to sweep. Read the held-out result as a check on whether two parameters
generalise across salt systems, which it can genuinely answer, rather than as a
forecast of unseen data, which it cannot.

The split is a real one in the way that matters most. The fitting experiment
uses sodium chloride and magnesium chloride brines; the held-out experiment uses
artificial sea salt, a different composition with sulfate and a different
divalent balance. Parameters fitted on one and applied to the other are being
asked to transfer across salt chemistry, not merely across replicates.
"""

import contextlib
import io
import logging
import os

import numpy as np
import pandas as pd

from .brines import build_batch_table
from .decks import generate_deck_for_batch
from .figures import load_measured, load_model_run, model_methane_headspace
from .run_decks import run_deck
from .scoring import score

logger = logging.getLogger(__name__)

# ═════════════════════════════════════════════════════════════════════
# Fixed in advance
# ═════════════════════════════════════════════════════════════════════

# Chloride concentrations [mol/L] at which the inhibition may centre. The range
# brackets the measured brines, whose incubation chloride runs from zero to
# about 3.9 mol/L.
THRESHOLD_GRID = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)

# Width of the transition, in decades. Narrower is a sharper cliff.
INTERVAL_GRID = (0.3, 0.5, 0.75, 1.0)

FIT_EXPERIMENT = "Exp004"
HELD_OUT_EXPERIMENT = "Exp003"


# ═════════════════════════════════════════════════════════════════════
# Running one configuration
# ═════════════════════════════════════════════════════════════════════


def evaluate(
    batches,
    measured,
    threshold,
    interval,
    work_root,
    repo_root,
    tag,
):
    """Run every batch at one parameter pair and return its score.

    Returns
    -------
    tuple
        ``(score, per_batch)`` where ``per_batch`` lists one record per
        condition.
    """
    deck_dir = os.path.join(work_root, f"decks_{tag}")
    run_root = os.path.join(work_root, f"runs_{tag}")

    modelled, observed, records = [], [], []
    for _, batch in batches.iterrows():
        name = (
            f"{batch['Experiment']}_B{int(batch['Batch ID']):02d}_{batch['Brine Name']}"
        )
        # The generator narrates each deck it writes, which is helpful once and
        # unreadable across a grid of twenty-four configurations.
        with contextlib.redirect_stdout(io.StringIO()):
            deck = generate_deck_for_batch(
                batch,
                output_dir=deck_dir,
                cellulose_hydrolysis={},
                aw_sandbox_replaces_network_methanogenesis=False,
                salinity_inhibition={
                    "species": "Cl-",
                    "threshold": threshold,
                    "interval": interval,
                },
            )
        target = os.path.join(deck_dir, f"{name}.in")
        if deck != target:
            os.replace(deck, target)

        result = run_deck(target, run_root, repo_root)
        if result["returncode"] != 0:
            logger.warning(
                "%s failed at threshold=%s interval=%s", name, threshold, interval
            )
            continue

        frame = load_model_run(result["workdir"])
        if frame is None:
            continue
        _, series = model_methane_headspace(frame, batch)

        rows = measured[
            (measured["Experiment"] == batch["Experiment"])
            & (measured["Batch ID"] == batch["Batch ID"])
        ]
        if not len(rows):
            continue
        final = (
            rows.sort_values("Days since start")
            .groupby("Replicate ID")["Cumulative Moles"]
            .last()
        )
        observed_value = float(np.median(final))

        modelled.append(series[-1])
        observed.append(observed_value)
        records.append(
            {
                "batch": name,
                "water activity": batch["Measured Water Activity"],
                "modelled": series[-1],
                "measured": observed_value,
                "log10 ratio": np.log10(series[-1] / observed_value),
            }
        )

    return score(modelled, observed), records


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "RETIRED from the default comparison. Fit the Cl- smoothstep on one "
            "experiment, then apply it to a held-out one (historical replay)."
        )
    )
    parser.add_argument("--work-root", default="calibration")
    parser.add_argument("--repo-root", default=os.getcwd())
    parser.add_argument(
        "--composition",
        default=os.path.join("data", "incubation_batch_composition.csv"),
    )
    parser.add_argument(
        "--ecsv-glob",
        default=os.path.expanduser(
            "~/Documents/GitHub/saltyBiomass/data/transformed/*Exp00[34]*.ecsv"
        ),
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    print(
        "NOTE: calibrate is RETIRED from the default comparison "
        "(Cl- smoothstep fit; live path uses AWINHIBIT). Historical replay only.\n"
    )

    table = (
        pd.read_csv(args.composition)
        if os.path.exists(args.composition)
        else build_batch_table()
    )
    measured = load_measured(args.ecsv_glob, "CH4_FID")

    fit_batches = table[table["Experiment"] == FIT_EXPERIMENT]
    held_out_batches = table[table["Experiment"] == HELD_OUT_EXPERIMENT]

    print("=" * 76)
    print(f"STEP 1  fit on {FIT_EXPERIMENT} only ({len(fit_batches)} batches)")
    print("=" * 76)
    print("Objective: root mean squared log10(modelled/measured). Lower is better.")
    print(f"Grid: {len(THRESHOLD_GRID)} thresholds x {len(INTERVAL_GRID)} intervals")
    print()

    results = []
    for threshold in THRESHOLD_GRID:
        row = []
        for interval in INTERVAL_GRID:
            value, _ = evaluate(
                fit_batches,
                measured,
                threshold,
                interval,
                args.work_root,
                args.repo_root,
                tag=f"fit_t{threshold}_i{interval}",
            )
            results.append((value, threshold, interval))
            row.append(value)
        print(
            f"  threshold {threshold:4.2f} mol/L   "
            + "   ".join(
                f"interval {i:.2f}: {v:5.2f}" for i, v in zip(INTERVAL_GRID, row)
            )
        )

    best_score, best_threshold, best_interval = min(results)
    print()
    print(
        f"Selected: threshold {best_threshold} mol/L, interval {best_interval} decades"
    )
    print(f"  score on {FIT_EXPERIMENT}: {best_score:.2f}")

    print()
    print("=" * 76)
    print(
        f"STEP 2  apply unchanged to {HELD_OUT_EXPERIMENT} ({len(held_out_batches)} batches)"
    )
    print("=" * 76)
    print("No parameter is adjusted below this line.")
    print()

    held_out_score, held_out_records = evaluate(
        held_out_batches,
        measured,
        best_threshold,
        best_interval,
        args.work_root,
        args.repo_root,
        tag="heldout",
    )

    frame = pd.DataFrame(held_out_records)
    pd.set_option("display.width", 200)
    print(frame.to_string(index=False, float_format=lambda v: f"{v:.3e}"))
    print()
    print(f"score on {FIT_EXPERIMENT} (fitted)  : {best_score:.2f}")
    print(f"score on {HELD_OUT_EXPERIMENT} (held out): {held_out_score:.2f}")
    print()
    ratios = frame["log10 ratio"].abs()
    print(
        f"held-out: median miss a factor of {10**ratios.median():.1f}, "
        f"{100 * (ratios <= 1).mean():.0f}% within a factor of ten"
    )


if __name__ == "__main__":
    main()
