"""Fit on the early timepoints of an incubation, then predict the rest of it.

A different question from the one :mod:`calibrate` asks. That module fits on one
experiment and applies the result to another, which tests whether parameters
transfer across salt chemistry. This one stays inside a single experiment and
splits it in time: the water-activity inhibition thresholds are fitted using only
the first few sampling rounds, and then asked to reproduce the sampling rounds
that came after.

Decks match the live comparison path: cellulose hydrolysis, AWINHIBIT
``ONE_MINUS_AW`` sandboxes, and the network's chloride Monod inhibition off so
salt is not double-counted. See :mod:`decks` and the README "Running it"
section.

Implementation is split across :mod:`forecast_split`, :mod:`forecast_grid`,
:mod:`forecast_scoring`, and :mod:`scoring` so the live comparison path and the
retired calibrate replay do not share a monolithic module.
"""

import logging
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .brines import build_batch_table  # noqa: E402
from .figures import (  # noqa: E402
    PALETTE,
    colour_for_brine,
    load_measured,
)
from .forecast_grid import (  # noqa: E402
    AW_ACETATE_OFFSET,
    AW_METHYL_OFFSET,
    AW_THRESHOLD_GRID,
    measured_by_batch_and_day,
    model_at_days,
    pathway_aw_thresholds,
    run_grid,
)
from .forecast_anchor import (  # noqa: E402
    DEFAULT_ANCHOR_DAY,
    anchor_table,
    model_days_to_physical,
)
from .forecast_scoring import (  # noqa: E402
    observed_by_batch,
    score_flux_window,
    score_window,
)
from .forecast_split import (  # noqa: E402
    DEFAULT_FIT_ROUNDS,
    DEFAULT_FLUX_SKIP_DAYS,
    DEFAULT_HOLDOUT_EARLY_DAYS,
    EXCLUDE_DAY_ZERO,
    flux_interval_pairs,
    forecast_output_dir,
    partition_flux_pairs,
    split_forecast_days,
    window_label,
)
from .scoring import interval_production_rate  # noqa: E402

logger = logging.getLogger(__name__)

EXPERIMENTS = ("Exp003", "Exp004")

# Re-exported for tests and callers that imported from this module.
__all__ = [
    "AW_THRESHOLD_GRID",
    "DEFAULT_ANCHOR_DAY",
    "DEFAULT_FIT_ROUNDS",
    "DEFAULT_FLUX_SKIP_DAYS",
    "DEFAULT_HOLDOUT_EARLY_DAYS",
    "EXCLUDE_DAY_ZERO",
    "EXPERIMENTS",
    "flux_interval_pairs",
    "forecast_deck_kwargs",
    "forecast_experiment",
    "forecast_output_dir",
    "partition_flux_pairs",
    "pathway_aw_thresholds",
    "split_forecast_days",
    "window_label",
]

from .forecast_grid import forecast_deck_kwargs  # noqa: E402


def plot_forecast(
    experiment,
    batches,
    observed,
    series,
    fit_days,
    predict_days,
    aw_threshold,
    aw_threshold_methyl,
    aw_threshold_acetate,
    fit_score,
    predict_score,
    output_path,
    early_holdout_days=(),
    holdout_early_days=0,
    anchor_day=0,
):
    """Timecourse with the fitting window shaded and the rest left to predict."""
    figure, axis = plt.subplots(figsize=(11, 6.8))

    if anchor_day > 0:
        axis.axvline(
            anchor_day,
            color=PALETTE["guide"],
            linewidth=1.0,
            linestyle=":",
            zorder=1,
        )

    if early_holdout_days:
        axis.axvspan(
            -3,
            max(early_holdout_days),
            color=PALETTE["guide_light"],
            alpha=0.15,
            zorder=0,
        )

    boundary = (
        (max(fit_days) + min(predict_days)) / 2 if predict_days else max(fit_days)
    )
    axis.axvspan(
        max(early_holdout_days) if early_holdout_days else -3,
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
        if anchor_day > 0:
            model_days = model_days_to_physical(model_days, anchor_day)
        axis.plot(
            model_days, model_moles, color=colour, linewidth=2, alpha=0.9, zorder=2
        )

        rows = observed[observed["Batch ID"] == batch_id]
        inside = rows[rows["day"].isin(fit_days)]
        outside = rows[rows["day"].isin(predict_days)]
        early = rows[rows["day"].isin(early_holdout_days)]

        if len(early):
            axis.scatter(
                early["day"],
                early["Cumulative Moles"],
                color=colour,
                s=46,
                marker="x",
                linewidth=1.4,
                zorder=2,
            )
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
    axis.set_ylim(1e-9, 5e-3)
    axis.set_xlim(-3, 133)
    axis.set_xlabel("Days since start of incubation")
    axis.set_ylabel("Methane in the bottle headspace (moles)")
    axis.grid(True, alpha=0.22, linewidth=0.5)

    holdout_note = (
        f", first {int(holdout_early_days)} days withheld from fitting"
        if holdout_early_days > 0
        else ""
    )
    anchor_note = (
        f", warm-started at day {int(anchor_day)} from measured CH4 and CO2"
        if anchor_day > 0
        else ""
    )
    axis.set_title(
        f"{experiment}: parameters fitted on the first {len(fit_days)} eligible "
        f"sampling rounds{holdout_note}{anchor_note}, then asked to predict the rest\n"
        f"ONE_MINUS_AW a_crit: H2={aw_threshold:.2f}, methyl={aw_threshold_methyl:.2f}, "
        f"acetate={aw_threshold_acetate:.2f}   |   "
        f"fitted {fit_score:.2f}, predicted {predict_score:.2f}",
        fontsize=12,
        pad=14,
    )

    axis.text(
        0.02,
        0.97,
        "early holdout" if early_holdout_days else "fitted on this",
        transform=axis.transAxes,
        fontsize=9.5,
        va="top",
        color=PALETTE["guide"],
    )
    if early_holdout_days:
        axis.text(
            0.28,
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
    marker_key = []
    if early_holdout_days:
        marker_key.append(
            plt.Line2D(
                [],
                [],
                marker="x",
                linestyle="",
                color=PALETTE["guide"],
                markersize=8,
                label="measured, early holdout (not fitted)",
            )
        )
    marker_key.extend(
        [
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
    )
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


def forecast_experiment(
    experiment,
    table,
    measured,
    work_root,
    repo_root,
    output_root,
    fit_rounds,
    holdout_early_days=DEFAULT_HOLDOUT_EARLY_DAYS,
    aw_threshold_grid=None,
    score_flux=False,
    flux_skip_days=DEFAULT_FLUX_SKIP_DAYS,
    anchor_day=0,
    measured_co2=None,
):
    """Fit on the early sampling rounds of one experiment, predict the rest."""
    batches = table[table["Experiment"] == experiment]
    observed = measured_by_batch_and_day(measured, experiment)

    fit_days, predict_days, early_holdout_days, production_days = split_forecast_days(
        observed["day"].unique(),
        fit_rounds,
        holdout_early_days=holdout_early_days,
        anchor_day=anchor_day,
    )
    if not fit_days:
        raise ValueError(
            f"{experiment}: no sampling days left to fit after excluding day zero"
            + (f", anchor day {anchor_day}" if anchor_day else "")
            + f", and the first {holdout_early_days} days"
        )
    if not predict_days:
        raise ValueError(
            f"{experiment}: no sampling days left to predict with fit_rounds="
            f"{fit_rounds} and holdout_early_days={holdout_early_days}"
        )

    all_flux_pairs = []
    fit_flux_pairs = []
    predict_flux_pairs = []

    print(f"\n{'=' * 76}")
    print(
        f"{experiment}: {len(batches)} batches, production sampling days "
        f"{[int(day) for day in production_days]}"
    )
    print(f"{'=' * 76}")
    if anchor_day > 0:
        print(
            f"  anchor day        {int(anchor_day)} "
            "(decks warm-started from interpolated measured CH4 and CO2 headspace)"
        )
    if early_holdout_days:
        print(
            f"  early holdout     {[int(day) for day in early_holdout_days]} "
            f"(first {int(holdout_early_days)} days: fast-substrate methane, "
            "not used for fitting)"
        )
    print(f"  fitting on days   {[int(day) for day in fit_days]}")
    print(f"  predicting days   {[int(day) for day in predict_days]}")
    if score_flux:
        all_flux_pairs = flux_interval_pairs(production_days, flux_skip_days)
        fit_flux_pairs, predict_flux_pairs = partition_flux_pairs(
            all_flux_pairs, fit_days, predict_days
        )
        if not fit_flux_pairs:
            raise ValueError(
                f"{experiment}: no production-rate intervals left to fit after "
                f"skipping intervals starting at or below day {flux_skip_days}"
            )
        if not predict_flux_pairs:
            raise ValueError(
                f"{experiment}: no production-rate intervals left to predict"
            )
        print(
            f"  scoring on flux   intervals starting after day {flux_skip_days} "
            "(mol/day, average over each sampling interval)"
        )
        print(
            f"  fit intervals     "
            f"{[(int(t0), int(t1)) for t0, t1 in fit_flux_pairs]}"
        )
        print(
            f"  predict intervals "
            f"{[(int(t0), int(t1)) for t0, t1 in predict_flux_pairs]}"
        )
    print("  (day 0 excluded: a baseline reading, not a production measurement)")
    print()

    if aw_threshold_grid is None:
        aw_threshold_grid = AW_THRESHOLD_GRID

    grid = run_grid(
        batches,
        work_root,
        repo_root,
        tag_prefix=experiment,
        aw_threshold_grid=aw_threshold_grid,
        anchor_day=anchor_day,
        measured_ch4=measured if anchor_day > 0 else None,
        measured_co2=measured_co2,
    )

    scored = []
    for (aw_h2, aw_methyl, aw_acetate), series in grid.items():
        if score_flux:
            value, count = score_flux_window(
                series, observed, fit_flux_pairs, anchor_day=anchor_day
            )
        else:
            value, count = score_window(
                series, observed, fit_days, anchor_day=anchor_day
            )
        scored.append((value, aw_h2, aw_methyl, aw_acetate, count))
    fit_score, aw_h2, aw_methyl, aw_acetate, fit_count = min(scored)

    if score_flux:
        predict_score, predict_count = score_flux_window(
            grid[(aw_h2, aw_methyl, aw_acetate)],
            observed,
            predict_flux_pairs,
            anchor_day=anchor_day,
        )
    else:
        predict_score, predict_count = score_window(
            grid[(aw_h2, aw_methyl, aw_acetate)],
            observed,
            predict_days,
            anchor_day=anchor_day,
        )

    at_edge = []
    if aw_h2 in (min(aw_threshold_grid), max(aw_threshold_grid)):
        at_edge.append(f"hydrogenotrophic a_crit {aw_h2}")

    print(
        "  selected: ONE_MINUS_AW a_crit "
        f"H2={aw_h2:.2f}, methyl={aw_methyl:.2f}, acetate={aw_acetate:.2f}"
    )
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
    if score_flux:
        print(
            "    (scores are on interval-averaged production rates, not cumulative moles)"
        )

    directory = forecast_output_dir(
        output_root,
        experiment,
        fit_rounds,
        holdout_early_days,
        score_flux=score_flux,
        anchor_day=anchor_day,
    )
    os.makedirs(directory, exist_ok=True)

    if anchor_day > 0:
        anchors = anchor_table(measured, measured_co2, batches, anchor_day)
        anchors.drop(columns=["concentrations"]).to_csv(
            os.path.join(directory, "anchor_initial.csv"), index=False
        )

    figure_path = plot_forecast(
        experiment,
        batches,
        observed,
        grid[(aw_h2, aw_methyl, aw_acetate)],
        fit_days,
        predict_days,
        aw_h2,
        aw_methyl,
        aw_acetate,
        fit_score,
        predict_score,
        os.path.join(directory, "methane_forecast.png"),
        early_holdout_days=early_holdout_days,
        holdout_early_days=holdout_early_days,
        anchor_day=anchor_day,
    )

    pd.DataFrame(
        [
            {
                "aw_threshold_h2": h2,
                "aw_threshold_methyl": methyl,
                "aw_threshold_acetate": acetate,
                "fit score": value,
                "points": count,
            }
            for value, h2, methyl, acetate, count in scored
        ]
    ).sort_values("fit score").to_csv(
        os.path.join(directory, "fit_grid.csv"), index=False
    )

    records = []
    for _, row in observed.iterrows():
        if anchor_day > 0 and row["day"] <= anchor_day:
            continue
        batch_id = int(row["Batch ID"])
        if batch_id not in grid[(aw_h2, aw_methyl, aw_acetate)]:
            continue
        model_days, model_moles = grid[(aw_h2, aw_methyl, aw_acetate)][batch_id]
        sim_day = float(row["day"]) - float(anchor_day)
        modelled = float(model_at_days(model_days, model_moles, sim_day))
        measured_value = float(row["Cumulative Moles"])
        records.append(
            {
                "batch": batch_id,
                "day": row["day"],
                "window": window_label(
                    row["day"], fit_days, predict_days, early_holdout_days
                ),
                "modelled": modelled,
                "measured": measured_value,
                "log10 ratio": (
                    np.log10(modelled / measured_value)
                    if modelled > 0 and measured_value > 0
                    else np.nan
                ),
            }
        )
    pd.DataFrame(records).to_csv(
        os.path.join(directory, "per_timepoint.csv"), index=False
    )

    if score_flux:
        interval_records = []
        lookup = observed_by_batch(observed)
        winning = grid[(aw_h2, aw_methyl, aw_acetate)]
        for batch_id, (model_days, model_moles) in winning.items():
            batch_measured = lookup.get(batch_id, {})
            for t0, t1 in all_flux_pairs:
                if t0 not in batch_measured or t1 not in batch_measured:
                    continue
                model_rate = interval_production_rate(
                    model_at_days(
                        model_days, model_moles, float(t0) - float(anchor_day)
                    ),
                    model_at_days(
                        model_days, model_moles, float(t1) - float(anchor_day)
                    ),
                    t0,
                    t1,
                )
                measured_rate = interval_production_rate(
                    batch_measured[t0], batch_measured[t1], t0, t1
                )
                if (t0, t1) in fit_flux_pairs:
                    window = "fitted"
                elif (t0, t1) in predict_flux_pairs:
                    window = "held out"
                else:
                    window = "excluded"
                interval_records.append(
                    {
                        "batch": batch_id,
                        "day_start": t0,
                        "day_end": t1,
                        "window": window,
                        "modelled_rate": model_rate,
                        "measured_rate": measured_rate,
                        "log10 ratio": (
                            np.log10(model_rate / measured_rate)
                            if model_rate > 0 and measured_rate > 0
                            else np.nan
                        ),
                    }
                )
        pd.DataFrame(interval_records).to_csv(
            os.path.join(directory, "per_interval.csv"), index=False
        )

    print(f"    wrote {figure_path}")
    return {
        "experiment": experiment,
        "aw_threshold": aw_h2,
        "aw_threshold_methyl": aw_methyl,
        "aw_threshold_acetate": aw_acetate,
        "fit score": fit_score,
        "predict score": predict_score,
        "directory": directory,
        "at grid edge": bool(at_edge),
        "fit days": fit_days,
        "predict days": predict_days,
        "early holdout days": early_holdout_days,
        "holdout early days": holdout_early_days,
        "score flux": score_flux,
        "flux skip days": flux_skip_days,
        "anchor day": anchor_day,
    }


def write_protocol_note(
    output_root,
    fit_rounds,
    summary,
    holdout_early_days=0,
    score_flux=False,
    flux_skip_days=0,
    anchor_day=0,
):
    """Leave an explanation beside the figures."""
    lines = [
        "# Forecasting test: fit early, predict late",
        "",
        "Each experiment here was split in time rather than by experiment. The",
        "hydrogenotrophic AWINHIBIT ``a_crit`` was swept on a fixed grid; methyl",
        "and acetoclastic thresholds kept the live comparison spacing (+0.05 and",
        "+0.10 in water activity). The winner was chosen using only the first "
        f"{fit_rounds} eligible sampling rounds of that experiment, and then asked "
        "to reproduce the rounds that came after. Nothing else was adjusted.",
        "",
        "Decks match the live comparison: cellulose hydrolysis, ONE_MINUS_AW",
        "sandboxes, and the network chloride Monod inhibition off.",
        "",
    ]
    if score_flux:
        lines += [
            f"The objective used interval-averaged headspace methane production rates "
            f"(mol/day) rather than cumulative moles. Intervals that start at or below "
            f"day {int(flux_skip_days)} were excluded so the trace-to-first-sample "
            "step-up does not dominate the fit.",
            "",
        ]
    if anchor_day > 0:
        lines += [
            f"Decks were warm-started at day {int(anchor_day)} using measured headspace "
            "CH4 and CO2 interpolated onto that day. Substrate and redox pools were "
            "left at their default t=0 values, so this is a gas-anchored sensitivity "
            "run rather than a full checkpoint restart.",
            "",
        ]
    lines += [
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
    ]
    if holdout_early_days > 0:
        lines += [
            f"Sampling days in the first {int(holdout_early_days)} days are also withheld "
            "from fitting. Methane released that early is dominated by fast-substrate "
            "turnover rather than salt inhibition, so those points are shown on the "
            "figures but do not enter the objective.",
            "",
        ]
    lines += [
        f"After those exclusions, the first {fit_rounds} remaining production rounds "
        "are fitted on and the rest are predicted.",
        "",
        "## Results",
        "",
    ]
    if score_flux:
        lines += [
            "Score is the root mean squared error of log10(modelled/measured) on "
            "interval-averaged production rates (mol/day). Lower is better; a value "
            "of one is a typical miss by a factor of ten.",
            "",
        ]
    else:
        lines += [
            "Score is the root mean squared error of log10(modelled/measured). Lower is",
            "better; a value of one is a typical miss by a factor of ten.",
            "",
        ]
    lines += [
        "| Experiment | H2 a_crit | Methyl a_crit | Acetate a_crit | Fit days | Predict days | Fitted | Held-out |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in summary:
        fit_days = ", ".join(str(int(day)) for day in row["fit days"])
        predict_days = ", ".join(str(int(day)) for day in row["predict days"])
        lines.append(
            f"| {row['experiment']} | {row['aw_threshold']:.2f} | "
            f"{row['aw_threshold_methyl']:.2f} | {row['aw_threshold_acetate']:.2f} "
            f"| {fit_days} | {predict_days} | {row['fit score']:.2f} | "
            f"{row['predict score']:.2f} |"
        )
    if any(row["at grid edge"] for row in summary):
        lines += [
            "",
            "## The selection hit the edge of the grid",
            "",
            "At least one experiment chose a hydrogenotrophic ``a_crit`` on the",
            "boundary of the search grid. A search that stops at its own boundary",
            "has not found a minimum; it has run out of room.",
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
    ]
    if score_flux:
        lines += [
            "- `per_interval.csv` -- modelled and measured interval-averaged production",
            "  rates for each batch and sampling interval.",
        ]
    lines += [
        "",
        "## Caveat",
        "",
        "These experiments were examined during the work that produced this model, so",
        "this is a pre-registered protocol rather than a blind forecast. The objective",
        "and the parameter grid are fixed in advance and the selection is arithmetic,",
        "but knowledge of how the incubations end cannot be unlearned.",
        "",
        "Regenerate with "
        "`python -m pflotran_py.comparison.forecast "
        f"--fit-rounds {fit_rounds}"
        + (" --score-flux" if score_flux else "")
        + (
            f" --flux-skip-days {int(flux_skip_days)}"
            if score_flux and flux_skip_days
            else ""
        )
        + (f" --anchor-day {int(anchor_day)}" if anchor_day > 0 else "")
        + (
            f" --holdout-early-days {int(holdout_early_days)}`."
            if holdout_early_days > 0
            else "`."
        ),
        "",
    ]
    path = os.path.join(output_root, "README.md")
    with open(path, "w") as handle:
        handle.write("\n".join(lines))
    return path


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Fit AWINHIBIT water-activity thresholds on early sampling rounds, "
            "then predict the rest of each incubation."
        )
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
    parser.add_argument(
        "--holdout-early-days",
        type=int,
        default=DEFAULT_HOLDOUT_EARLY_DAYS,
        help=(
            "exclude sampling days at or below this age from fitting; early methane "
            "reflects fast-substrate turnover rather than salt inhibition"
        ),
    )
    parser.add_argument(
        "--aw-threshold-grid",
        type=float,
        nargs="+",
        default=None,
        metavar="A_W",
        help=(
            "hydrogenotrophic a_crit values to sweep; methyl and acetate follow "
            f"the live comparison offsets (+{AW_METHYL_OFFSET:.2f} / "
            f"+{AW_ACETATE_OFFSET:.2f}). Default: {AW_THRESHOLD_GRID}"
        ),
    )
    parser.add_argument(
        "--score-flux",
        action="store_true",
        help=(
            "score interval-averaged headspace methane production rates (mol/day) "
            "instead of cumulative moles at each sampling day"
        ),
    )
    parser.add_argument(
        "--flux-skip-days",
        type=int,
        default=DEFAULT_FLUX_SKIP_DAYS,
        help=(
            "when using --score-flux, omit intervals that start at or below this day "
            f"(default {DEFAULT_FLUX_SKIP_DAYS})"
        ),
    )
    parser.add_argument(
        "--anchor-day",
        type=int,
        default=0,
        help=(
            "warm-start decks from measured headspace CH4 and CO2 interpolated to "
            f"this incubation day (0 disables; try {DEFAULT_ANCHOR_DAY})"
        ),
    )
    args = parser.parse_args()

    if args.anchor_day < 0:
        parser.error("--anchor-day must be non-negative")

    aw_threshold_grid = tuple(args.aw_threshold_grid or AW_THRESHOLD_GRID)
    for value in aw_threshold_grid:
        pathway_aw_thresholds(value)

    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    print(
        "Forecast uses live-comparison decks (cellulose hydrolysis, ONE_MINUS_AW "
        "AWINHIBIT, no network Cl- Monod)."
    )
    if args.score_flux:
        print(
            f"Objective: interval-averaged production rates, skipping intervals "
            f"starting at or below day {args.flux_skip_days}."
        )
    if args.anchor_day > 0:
        print(
            f"Warm start: measured CH4 and CO2 headspace at day {args.anchor_day} "
            "written into initial constraints."
        )
    if args.score_flux or args.anchor_day > 0:
        print()
    elif not args.score_flux:
        print()

    table = (
        pd.read_csv(args.composition)
        if os.path.exists(args.composition)
        else build_batch_table()
    )
    measured = load_measured(args.ecsv_glob, "CH4_FID")
    measured_co2 = (
        load_measured(args.ecsv_glob, "CO2") if args.anchor_day > 0 else None
    )
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
            args.holdout_early_days,
            aw_threshold_grid=aw_threshold_grid,
            score_flux=args.score_flux,
            flux_skip_days=args.flux_skip_days,
            anchor_day=args.anchor_day,
            measured_co2=measured_co2,
        )
        for experiment in EXPERIMENTS
    ]

    note = write_protocol_note(
        args.output_root,
        args.fit_rounds,
        summary,
        holdout_early_days=args.holdout_early_days,
        score_flux=args.score_flux,
        flux_skip_days=args.flux_skip_days,
        anchor_day=args.anchor_day,
    )
    print()
    print("=" * 76)
    print("SUMMARY")
    print("=" * 76)
    for row in summary:
        print(
            f"  {row['experiment']}  a_crit H2={row['aw_threshold']:.2f} "
            f"methyl={row['aw_threshold_methyl']:.2f} "
            f"acetate={row['aw_threshold_acetate']:.2f}   "
            f"fitted {row['fit score']:.2f}   "
            f"predicted {row['predict score']:.2f}"
        )
    print()
    print(f"  wrote {note}")


if __name__ == "__main__":
    main()
