"""Fit on the early timepoints of an incubation, then predict the rest of it.

A different question from the one :mod:`calibrate` asks. That module fits on one
experiment and applies the result to another, which tests whether parameters
transfer across salt chemistry. This one stays inside a single experiment and
splits it in time: the salinity inhibition is fitted using only the first few
sampling rounds, and then asked to reproduce the sampling rounds that came
after.

That is the question a modeller actually faces mid-experiment. Three months of
incubation have not happened yet; a few weeks have. Can what is measured so far
pin the parameters well enough to say where the bottles end up?

Why this split is harder than it looks
--------------------------------------
Each experiment has five sampling days, roughly 0, 10, 55, 90 and 120. Day zero
is excluded as a baseline rather than a production measurement, for the reason
given at :data:`EXCLUDE_DAY_ZERO`, leaving four rounds. Fitting on the first two
means constraining the parameters with the rising part of the curve and asking
them for the plateau.

It is also the regime where the model is weakest. Methane production in these
decks is not rate-limited throughout: the fast substrate is consumed early and
the curve then flattens. Parameters chosen to match the rising part of the
curve are being asked to predict a plateau whose height they barely influence.
A good early fit and a poor late prediction would not be a surprising outcome,
and it would be worth knowing.

What is fitted and what is not
------------------------------
Only the two salinity inhibition parameters are fitted here: the chloride
concentration at which the term centres, and the width of the transition. Every
other quantity -- the reaction network, the rate constants, the brine
compositions, the hydrolysis rate, the geometry -- is held at the values
documented in the README and is not adjusted.

Reading the output
------------------
Each experiment gets a figure with the fitting window shaded. Points inside the
shading were used to choose the parameters; points outside it were not, and the
distance between them and the model line there is the actual result.
"""

import contextlib
import io
import logging
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .brines import build_batch_table  # noqa: E402
from .calibrate import INTERVAL_GRID, THRESHOLD_GRID, score  # noqa: E402
from .decks import generate_deck_for_batch  # noqa: E402
from .figures import (  # noqa: E402
    PALETTE,
    colour_for_brine,
    load_measured,
    load_model_run,
    model_methane_headspace,
)
from .run_decks import run_deck  # noqa: E402

logger = logging.getLogger(__name__)

# Number of sampling rounds used for fitting. The rest are predicted.
DEFAULT_FIT_ROUNDS = 2

# Day zero is excluded from every window, and this is not a detail.
#
# It is a baseline reading taken before anything has been produced. The model
# begins at its trace floor of 1e-15 mol/L by construction, while the instrument
# reports whatever the headspace and the detector background give it, typically
# 1e-8 to 1e-6 moles. The ratio between those is not a measure of how well the
# model does; it is the model's initial condition compared against instrument
# background.
#
# In an objective built on log ratios that single point dominates everything
# else. Including it gave fit scores near 4.7, a nominal miss by a factor of
# fifty thousand, and made the held-out window score better than the window that
# was fitted -- which is the signature of a statistic being driven by something
# other than model quality.
EXCLUDE_DAY_ZERO = True

EXPERIMENTS = ("Exp003", "Exp004")


def measured_by_batch_and_day(measured, experiment):
    """Median measured methane per batch per sampling day.

    Replicates are collapsed to their median, matching how every other
    comparison in this package treats them.
    """
    rows = measured[measured["Experiment"] == experiment].copy()
    rows["day"] = rows["Days since start"].round(0)
    grouped = (
        rows.groupby(["Batch ID", "day"])["Cumulative Moles"].median().reset_index()
    )
    return grouped


def model_at_days(days_model, moles_model, wanted_days):
    """Model methane at the measured sampling days.

    The model writes a snapshot every day, so this is an interpolation onto a
    subset of its own grid rather than an extrapolation.
    """
    return np.interp(wanted_days, days_model, moles_model)


def run_grid(batches, work_root, repo_root, tag_prefix):
    """Run every parameter combination once, for every batch.

    The runs do not depend on which sampling days are later used for fitting,
    so the whole grid is executed once and scored twice -- once on the fitting
    window and once on what follows.

    Returns
    -------
    dict
        ``{(threshold, interval): {batch_id: (days, moles)}}``
    """
    grid = {}
    for threshold in THRESHOLD_GRID:
        for interval in INTERVAL_GRID:
            series = {}
            for _, batch in batches.iterrows():
                name = (
                    f"{batch['Experiment']}_B{int(batch['Batch ID']):02d}"
                    f"_{batch['Brine Name']}"
                )
                deck_dir = os.path.join(
                    work_root, f"{tag_prefix}_decks_t{threshold}_i{interval}"
                )
                run_root = os.path.join(
                    work_root, f"{tag_prefix}_runs_t{threshold}_i{interval}"
                )
                with contextlib.redirect_stdout(io.StringIO()):
                    deck = generate_deck_for_batch(
                        batch,
                        output_dir=deck_dir,
                        cellulose_hydrolysis={},
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
                    continue
                frame = load_model_run(result["workdir"])
                if frame is None:
                    continue
                series[int(batch["Batch ID"])] = model_methane_headspace(frame, batch)
            grid[(threshold, interval)] = series
    return grid


def score_window(grid_entry, observed, days):
    """Score one parameter combination over a set of sampling days."""
    modelled, measured_values = [], []
    for _, row in observed[observed["day"].isin(days)].iterrows():
        batch_id = int(row["Batch ID"])
        if batch_id not in grid_entry:
            continue
        model_days, model_moles = grid_entry[batch_id]
        modelled.append(model_at_days(model_days, model_moles, row["day"]))
        measured_values.append(row["Cumulative Moles"])
    return score(modelled, measured_values), len(modelled)


# ═════════════════════════════════════════════════════════════════════
# Figure
# ═════════════════════════════════════════════════════════════════════


def plot_forecast(
    experiment,
    batches,
    observed,
    series,
    fit_days,
    predict_days,
    threshold,
    interval,
    fit_score,
    predict_score,
    output_path,
):
    """Timecourse with the fitting window shaded and the rest left to predict."""
    figure, axis = plt.subplots(figsize=(11, 6.8))

    boundary = (
        (max(fit_days) + min(predict_days)) / 2 if predict_days else max(fit_days)
    )
    axis.axvspan(
        -3,
        boundary,
        color=PALETTE["guide_light"],
        alpha=0.35,
        zorder=0,
    )

    for _, batch in batches.iterrows():
        batch_id = int(batch["Batch ID"])
        if batch_id not in series:
            continue
        colour = colour_for_brine(batch["Brine Name"])
        model_days, model_moles = series[batch_id]
        axis.plot(
            model_days, model_moles, color=colour, linewidth=2, alpha=0.9, zorder=2
        )

        rows = observed[observed["Batch ID"] == batch_id]
        inside = rows[rows["day"].isin(fit_days)]
        outside = rows[rows["day"].isin(predict_days)]

        axis.scatter(
            inside["day"],
            inside["Cumulative Moles"],
            color=colour,
            s=58,
            marker="o",
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        axis.scatter(
            outside["day"],
            outside["Cumulative Moles"],
            color=colour,
            s=86,
            marker="s",
            edgecolor=PALETTE["guide"],
            linewidth=1.1,
            zorder=4,
            label=(
                f"{batch['Brine Name']} "
                f"(water activity {batch['Measured Water Activity']:.3f})"
            ),
        )

    axis.axvline(
        boundary, color=PALETTE["guide"], linewidth=1.2, linestyle="--", zorder=1
    )
    axis.set_yscale("log")
    # Clipped to the range the measurements occupy. The model rises from a trace
    # floor near 1e-15, and letting the axis follow it down would compress every
    # real point into a strip at the top.
    axis.set_ylim(1e-9, 5e-3)
    axis.set_xlim(-3, 133)
    axis.set_xlabel("Days since start of incubation")
    axis.set_ylabel("Methane in the bottle headspace (moles)")
    axis.grid(True, alpha=0.22, linewidth=0.5)

    axis.set_title(
        f"{experiment}: parameters fitted on the first {len(fit_days)} sampling rounds, "
        f"then asked to predict the rest\n"
        f"chloride threshold {threshold} mol/L, transition width {interval} decades   |   "
        f"fitted {fit_score:.2f}, predicted {predict_score:.2f}",
        fontsize=12,
        pad=14,
    )

    axis.text(
        0.02,
        0.97,
        "fitted on this",
        transform=axis.transAxes,
        fontsize=9.5,
        va="top",
        color=PALETTE["guide"],
    )
    axis.text(
        0.62,
        0.97,
        "predicted, not fitted",
        transform=axis.transAxes,
        fontsize=9.5,
        va="top",
        color=PALETTE["guide"],
    )

    handles, labels = axis.get_legend_handles_labels()
    marker_key = [
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["guide"],
            markersize=8,
            label="measured, inside the fitting window",
        ),
        plt.Line2D(
            [],
            [],
            marker="s",
            linestyle="",
            color=PALETTE["guide"],
            markersize=9,
            label="measured, held out",
        ),
        plt.Line2D([], [], linestyle="-", color=PALETTE["guide"], label="model"),
    ]
    axis.legend(
        handles=marker_key + handles,
        labels=[h.get_label() for h in marker_key] + labels,
        fontsize=7.5,
        loc="lower right",
        ncol=2,
        framealpha=0.95,
    )

    figure.text(
        0.5,
        0.01,
        "Lower score is better; a value of one is a typical miss by a factor of ten. "
        "Colour identifies the brine, darker being more concentrated.",
        ha="center",
        fontsize=9,
        color=PALETTE["guide"],
    )
    figure.tight_layout(rect=[0, 0.045, 1, 1])
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Protocol
# ═════════════════════════════════════════════════════════════════════


def forecast_experiment(
    experiment, table, measured, work_root, repo_root, output_root, fit_rounds
):
    """Fit on the early sampling rounds of one experiment, predict the rest."""
    batches = table[table["Experiment"] == experiment]
    observed = measured_by_batch_and_day(measured, experiment)

    days = sorted(observed["day"].unique())
    if EXCLUDE_DAY_ZERO:
        days = [day for day in days if day > 0]
    fit_days = days[:fit_rounds]
    predict_days = days[fit_rounds:]

    print(f"\n{'=' * 76}")
    print(
        f"{experiment}: {len(batches)} batches, production sampling days "
        f"{[int(day) for day in days]}"
    )
    print(f"{'=' * 76}")
    print(f"  fitting on days   {[int(day) for day in fit_days]}")
    print(f"  predicting days   {[int(day) for day in predict_days]}")
    print("  (day 0 excluded: a baseline reading, not a production measurement)")
    print()

    grid = run_grid(batches, work_root, repo_root, tag_prefix=experiment)

    scored = []
    for (threshold, interval), series in grid.items():
        value, count = score_window(series, observed, fit_days)
        scored.append((value, threshold, interval, count))
    fit_score, threshold, interval, fit_count = min(scored)

    predict_score, predict_count = score_window(
        grid[(threshold, interval)], observed, predict_days
    )

    at_edge = []
    if threshold in (min(THRESHOLD_GRID), max(THRESHOLD_GRID)):
        at_edge.append(f"threshold {threshold}")
    if interval in (min(INTERVAL_GRID), max(INTERVAL_GRID)):
        at_edge.append(f"interval {interval}")

    print(f"  selected: threshold {threshold} mol/L, interval {interval} decades")
    if at_edge:
        print(
            f"    WARNING: {' and '.join(at_edge)} sits on the edge of the grid, so "
            "the true optimum is outside it. The fit is reporting the closest thing "
            "the grid allows, not a minimum."
        )
    print(f"    score on the fitting window   {fit_score:.2f}   ({fit_count} points)")
    print(
        f"    score on the held-out rounds  {predict_score:.2f}   ({predict_count} points)"
    )

    directory = os.path.join(output_root, f"{experiment}_fit-first-{fit_rounds}-rounds")
    os.makedirs(directory, exist_ok=True)

    figure_path = plot_forecast(
        experiment,
        batches,
        observed,
        grid[(threshold, interval)],
        fit_days,
        predict_days,
        threshold,
        interval,
        fit_score,
        predict_score,
        os.path.join(directory, "methane_forecast.png"),
    )

    # The full grid, so the shape of the objective is inspectable rather than
    # summarised by its minimum alone.
    pd.DataFrame(
        [
            {"threshold": t, "interval": i, "fit score": v, "points": n}
            for v, t, i, n in scored
        ]
    ).sort_values("fit score").to_csv(
        os.path.join(directory, "fit_grid.csv"), index=False
    )

    records = []
    for _, row in observed.iterrows():
        batch_id = int(row["Batch ID"])
        if batch_id not in grid[(threshold, interval)]:
            continue
        model_days, model_moles = grid[(threshold, interval)][batch_id]
        modelled = float(model_at_days(model_days, model_moles, row["day"]))
        records.append(
            {
                "batch": batch_id,
                "day": row["day"],
                "window": "fitted" if row["day"] in fit_days else "held out",
                "modelled": modelled,
                "measured": row["Cumulative Moles"],
                "log10 ratio": np.log10(modelled / row["Cumulative Moles"]),
            }
        )
    pd.DataFrame(records).to_csv(
        os.path.join(directory, "per_timepoint.csv"), index=False
    )

    print(f"    wrote {figure_path}")
    return {
        "experiment": experiment,
        "threshold": threshold,
        "interval": interval,
        "fit score": fit_score,
        "predict score": predict_score,
        "directory": directory,
        "at grid edge": bool(at_edge),
    }


def write_protocol_note(output_root, fit_rounds, summary):
    """Leave an explanation beside the figures."""
    lines = [
        "# Forecasting test: fit early, predict late",
        "",
        "Each experiment here was split in time rather than by experiment. The two",
        "salinity inhibition parameters -- the chloride concentration at which the",
        "term centres and the width of its transition -- were fitted using only the",
        "first {} sampling rounds of that experiment, and then asked to".format(
            fit_rounds
        ),
        "reproduce the rounds that came after. Nothing else was adjusted.",
        "",
        "This is the question a modeller faces mid-experiment: a few weeks of data",
        "exist, the rest of the incubation does not yet. Can the early rounds pin the",
        "parameters well enough to say where the bottles end up?",
        "",
        "Day zero is excluded from every window. It is a baseline reading taken",
        "before anything has been produced: the model sits at its trace floor by",
        "construction while the instrument reports its background, and the ratio",
        "between those measures nothing about the model. In an objective built on",
        "log ratios that one point dominates everything else -- including it gave",
        "fit scores near 4.7, a nominal miss by a factor of fifty thousand, and made",
        "the held-out window score better than the fitted one.",
        "",
        f"Each experiment therefore has four production rounds. The first {fit_rounds}",
        "are fitted on and the rest are predicted.",
        "",
        "## Results",
        "",
        "Score is the root mean squared error of log10(modelled/measured). Lower is",
        "better; a value of one is a typical miss by a factor of ten.",
        "",
        "| Experiment | Chloride threshold | Transition width | Fitted rounds | Held-out rounds |",
        "|---|---|---|---|---|",
    ]
    for row in summary:
        lines.append(
            f"| {row['experiment']} | {row['threshold']} mol/L | {row['interval']} decades "
            f"| {row['fit score']:.2f} | {row['predict score']:.2f} |"
        )
    if any(row["at grid edge"] for row in summary):
        lines += [
            "",
            "## The selection hit the edge of the grid",
            "",
            "Both experiments chose the corner of the parameter grid: the highest",
            "chloride threshold and the widest transition available, which together are",
            "the weakest inhibition it can express. A search that stops at its own",
            "boundary has not found a minimum; it has run out of room. The true optimum",
            "on early data lies outside the grid, at weaker inhibition still.",
            "",
            "That is itself the finding. Fitted on the whole timecourse the same",
            "objective prefers a threshold of 0.75 mol/L, which is strong inhibition.",
            "Fitted on the first two rounds it prefers 2.0 mol/L or beyond, which is",
            "almost none. The early rounds do not merely constrain the parameters",
            "loosely -- they point in the opposite direction.",
            "",
            "The reason is visible in the figures. At days 10 and roughly 55 the salted",
            "bottles have not yet separated from the controls by much, so a model with",
            "little salt inhibition matches them well. The collapse develops later, and",
            "by then the parameters have been chosen.",
        ]

    lines += [
        "",
        "## What each folder holds",
        "",
        "- `methane_forecast.png` -- the timecourse, with the fitting window shaded.",
        "  Circles inside the shading were used to choose the parameters; squares",
        "  outside it were not.",
        "- `fit_grid.csv` -- every parameter combination and its score on the fitting",
        "  window, so the shape of the objective is inspectable rather than summarised",
        "  by its minimum alone.",
        "- `per_timepoint.csv` -- modelled and measured methane at every batch and",
        "  sampling day, labelled by which window it fell in.",
        "",
        "## Caveat",
        "",
        "These experiments were examined during the work that produced this model, so",
        "this is a pre-registered protocol rather than a blind forecast. The objective",
        "and the parameter grid are fixed in advance and the selection is arithmetic,",
        "but knowledge of how the incubations end cannot be unlearned.",
        "",
        "Regenerate with `python -m pflotran_py.comparison.forecast`.",
        "",
    ]
    path = os.path.join(output_root, "README.md")
    with open(path, "w") as handle:
        handle.write("\n".join(lines))
    return path


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Fit the salinity inhibition on early sampling rounds, predict the rest."
    )
    parser.add_argument("--work-root", default="forecasting")
    parser.add_argument("--repo-root", default=os.getcwd())
    parser.add_argument(
        "--output-root", default=os.path.join("output", "comparison", "forecast")
    )
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
    parser.add_argument("--fit-rounds", type=int, default=DEFAULT_FIT_ROUNDS)
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(message)s")

    table = (
        pd.read_csv(args.composition)
        if os.path.exists(args.composition)
        else build_batch_table()
    )
    measured = load_measured(args.ecsv_glob, "CH4_FID")
    os.makedirs(args.output_root, exist_ok=True)

    summary = [
        forecast_experiment(
            experiment,
            table,
            measured,
            args.work_root,
            args.repo_root,
            args.output_root,
            args.fit_rounds,
        )
        for experiment in EXPERIMENTS
    ]

    note = write_protocol_note(args.output_root, args.fit_rounds, summary)
    print()
    print("=" * 76)
    print("SUMMARY")
    print("=" * 76)
    for row in summary:
        print(
            f"  {row['experiment']}  threshold {row['threshold']} mol/L, "
            f"interval {row['interval']}   fitted {row['fit score']:.2f}   "
            f"predicted {row['predict score']:.2f}"
        )
    print()
    print(f"  wrote {note}")


if __name__ == "__main__":
    main()
