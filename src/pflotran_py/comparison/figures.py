"""Figures comparing the closed-batch model against the measured incubations.

Three figures, each answering a different question.

**Methane over time** asks whether the model reproduces the measured trajectory
in each bottle. Measured replicates are drawn as points, the model as a line,
one panel per salt family, on a logarithmic vertical axis because the measured
values span more than four orders of magnitude across conditions.

**Methane against water activity** asks the question the experiment was designed
around: does production fall as salt rises, and does the model fall with it?
This is the figure that shows whether the inhibition in the model matches the
inhibition in the bottles.

**Carbon dioxide** is presented separately and with a caveat rather than as a
like-for-like comparison, because the reaction network as configured cannot
produce a headspace carbon dioxide prediction. See :func:`plot_carbon_dioxide`.

Every colour comes from ``palette.json`` in this directory, whose anchors are
the colourblind-safe values adopted across the measurement repository. Hue
carries the salt condition; measured and modelled values are told apart by mark
rather than colour, so a reader learns one colour key for the whole set.
"""

import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy import units as u  # noqa: E402

from ..analysis.extract import extract_pflotran_data_hdf5  # noqa: E402
from ..analysis.extract import find_hdf5_output  # noqa: E402
from .headspace import (  # noqa: E402
    GASES,
    model_headspace_moles,
    setschenow_salts_from_composition,
)

PALETTE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "palette.json")

with open(PALETTE_PATH) as _handle:
    PALETTE = json.load(_handle)


def colour_for_brine(brine_name):
    """Map a brine to its palette role: hue by salt, intensity by strength."""
    name = str(brine_name)
    if name in ("C", "Su_C"):
        return PALETTE["control"]

    if name.startswith("Na_"):
        family = "nacl"
    elif name.startswith("Mg_"):
        family = "mgcl2"
    elif name.startswith("SW"):
        family = "seasalt"
    else:
        return PALETTE["guide"]

    if name.endswith("_H"):
        level = "high"
    elif name.endswith("_M"):
        level = "mid"
    else:
        level = "low"
    return PALETTE[f"{family}_{level}"]


def salt_family(brine_name):
    """Plain-language name of the salt a brine was made from."""
    name = str(brine_name)
    if name in ("C", "Su_C"):
        return "No added salt"
    if name.startswith("Na_"):
        return "Sodium chloride"
    if name.startswith("Mg_"):
        return "Magnesium chloride"
    if name.startswith("SW"):
        return "Artificial sea salt"
    return "Other"


# ═════════════════════════════════════════════════════════════════════
# Loading
# ═════════════════════════════════════════════════════════════════════


def load_measured(ecsv_glob, molecule):
    """Measured cumulative headspace moles per replicate over time."""
    from astropy.table import Table

    paths = sorted(glob.glob(ecsv_glob))
    if not paths:
        raise FileNotFoundError(f"No pipeline output matched {ecsv_glob}")

    frames = [Table.read(path).to_pandas() for path in paths]
    measured = pd.concat(frames, ignore_index=True)
    measured = measured[measured["Molecule"] == molecule]
    return measured.dropna(subset=["Cumulative Moles", "Days since start"])


def load_model_run(run_dir):
    """One PFLOTRAN run as a time-ordered frame."""
    h5_path = find_hdf5_output(run_dir, prefix="sim")
    if h5_path is None:
        return None
    return extract_pflotran_data_hdf5(h5_path).sort_values("Time [d]")


def model_methane_headspace(run_frame, batch_row):
    """Model methane in the headspace over time, in moles.

    Combines the Henry-partitioned dissolved methane with the moles the
    ebullition proxy has already moved out of solution.
    """
    nacl, mgcl2 = setschenow_salts_from_composition(batch_row)
    moles = [
        model_headspace_moles(
            row["Total CH4(aq) [M]"],
            GASES["CH4"],
            degassed_concentration=row["Total Tracer2 [M]"],
            nacl_molarity=nacl,
            mgcl2_molarity=mgcl2,
        ).to_value(u.mol)
        for _, row in run_frame.iterrows()
    ]
    return run_frame["Time [d]"].to_numpy(), np.array(moles)


def assemble(composition, run_root, ecsv_glob, molecule="CH4_FID"):
    """Pair every batch's model run with its measured replicates."""
    measured = load_measured(ecsv_glob, molecule)

    paired = []
    for _, batch in composition.iterrows():
        name = (
            f"{batch['Experiment']}_B{int(batch['Batch ID']):02d}_{batch['Brine Name']}"
        )
        run_dir = os.path.join(run_root, name)
        if not os.path.isdir(run_dir):
            continue
        run_frame = load_model_run(run_dir)
        if run_frame is None:
            continue

        rows = measured[
            (measured["Experiment"] == batch["Experiment"])
            & (measured["Batch ID"] == batch["Batch ID"])
        ]
        paired.append(
            {"batch": batch, "model": run_frame, "measured": rows, "name": name}
        )
    return paired


# ═════════════════════════════════════════════════════════════════════
# Figure 1 -- methane over time
# ═════════════════════════════════════════════════════════════════════


def plot_methane_timeseries(paired, output_path):
    """Measured and modelled methane in the headspace, over the incubation."""
    families = ["Sodium chloride", "Magnesium chloride", "Artificial sea salt"]
    figure, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)

    for axis, family in zip(axes, families):
        # One control is shown per panel as the no-salt reference. Both
        # experiments have one; the panel's own experiment is the right
        # comparison, since the two ran at different times.
        members = [
            entry
            for entry in paired
            if salt_family(entry["batch"]["Brine Name"]) == family
        ]
        experiments = {entry["batch"]["Experiment"] for entry in members}
        controls = [
            entry
            for entry in paired
            if entry["batch"]["Brine Name"] == "C"
            and entry["batch"]["Experiment"] in experiments
        ]

        for entry in controls + members:
            batch = entry["batch"]
            colour = colour_for_brine(batch["Brine Name"])
            name = (
                f"No salt ({batch['Experiment']})"
                if batch["Brine Name"] == "C"
                else batch["Brine Name"]
            )
            label = f"{name}, water activity {batch['Measured Water Activity']:.3f}"

            days, moles = model_methane_headspace(entry["model"], batch)
            axis.plot(days, moles, color=colour, linewidth=2, alpha=0.9, zorder=2)

            points = entry["measured"]
            if len(points):
                axis.scatter(
                    points["Days since start"],
                    points["Cumulative Moles"],
                    color=colour,
                    s=26,
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=3,
                    label=label,
                )

        axis.set_yscale("log")
        # Clipped to the range the data occupies. The model starts from a
        # numerical floor near 1e-15, and letting the axis chase it would
        # compress every measured point into a sliver at the top.
        axis.set_ylim(1e-9, 5e-2)
        axis.set_xlim(-3, 133)
        axis.set_xlabel("Days since start of incubation")
        axis.set_title(family, fontsize=12)
        axis.grid(True, alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=7.5, loc="lower right", framealpha=0.95)

    axes[0].set_ylabel("Methane in the bottle headspace (moles)")

    figure.suptitle(
        "Modelled and measured methane in sealed incubation bottles",
        fontsize=14,
        y=0.99,
    )
    figure.text(
        0.5,
        0.005,
        "Lines are the reactive-transport model; filled circles are gas-chromatograph "
        "measurements of individual bottles.\n"
        "Colour identifies the brine; darker shades are more concentrated. Vertical axis "
        "is logarithmic.",
        ha="center",
        fontsize=9,
        color=PALETTE["guide"],
    )
    figure.tight_layout(rect=[0, 0.06, 1, 0.96])
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Figure 2 -- methane against water activity
# ═════════════════════════════════════════════════════════════════════


def plot_methane_against_water_activity(paired, output_path):
    """Final methane production against water activity, model and measurement.

    The figure the experiment was designed to produce. Water activity falls to
    the right, so the horizontal axis reads as increasing salt stress.
    """
    figure, axis = plt.subplots(figsize=(9, 6.5))

    # Trends are drawn within a salt family, never across families. Two brines
    # made from different salts can share a water activity while differing in
    # ionic strength by more than half again, so a line joining them would
    # imply a dose-response that was never measured.
    by_family = {}
    for entry in paired:
        batch = entry["batch"]
        water_activity = batch["Measured Water Activity"]
        colour = colour_for_brine(batch["Brine Name"])

        _, moles = model_methane_headspace(entry["model"], batch)
        model_final = moles[-1]

        points = entry["measured"]
        measured_final = None
        if len(points):
            final = (
                points.sort_values("Days since start")
                .groupby("Replicate ID")["Cumulative Moles"]
                .last()
            )
            measured_final = float(np.median(final))
            axis.scatter(
                [water_activity] * len(final),
                final.to_numpy(),
                color=colour,
                s=55,
                edgecolor="white",
                linewidth=0.8,
                zorder=3,
            )

        axis.scatter(
            water_activity,
            model_final,
            color=colour,
            s=170,
            marker="_",
            linewidth=3,
            zorder=4,
        )

        family = salt_family(batch["Brine Name"])
        if family != "No added salt":
            by_family.setdefault(family, []).append(
                (water_activity, model_final, measured_final)
            )

    for family, records in by_family.items():
        records.sort()
        activities = [r[0] for r in records]
        colour = colour_for_brine(
            {"Sodium chloride": "Na_M", "Magnesium chloride": "Mg_M"}.get(
                family, "SW_M"
            )
        )
        axis.plot(
            activities,
            [r[1] for r in records],
            color=colour,
            linewidth=1.5,
            linestyle="--",
            alpha=0.75,
            zorder=1,
        )
        measured_series = [(a, m) for a, _, m in records if m is not None]
        if len(measured_series) > 1:
            axis.plot(
                *zip(*measured_series),
                color=colour,
                linewidth=1.5,
                linestyle=":",
                alpha=0.75,
                zorder=1,
            )

    axis.set_yscale("log")
    axis.invert_xaxis()
    axis.set_xlabel(
        "Water activity of the brine  (falling to the right, so salt stress increases rightwards)"
    )
    axis.set_ylabel("Methane in the headspace at the end of the incubation (moles)")
    axis.set_title(
        "The model does not reproduce the measured collapse of methane production under salt",
        fontsize=13,
        pad=14,
    )
    axis.grid(True, alpha=0.25, linewidth=0.5)

    handles = [
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["guide"],
            markersize=8,
            label="Measured, one point per bottle",
        ),
        plt.Line2D(
            [],
            [],
            marker="_",
            linestyle="",
            color=PALETTE["guide"],
            markersize=13,
            markeredgewidth=3,
            label="Model prediction",
        ),
        plt.Line2D(
            [],
            [],
            linestyle="--",
            color=PALETTE["guide"],
            label="Model, within one salt",
        ),
        plt.Line2D(
            [],
            [],
            linestyle=":",
            color=PALETTE["guide"],
            label="Measured, within one salt",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["nacl_mid"],
            markersize=8,
            label="Sodium chloride",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["mgcl2_mid"],
            markersize=8,
            label="Magnesium chloride",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["seasalt_mid"],
            markersize=8,
            label="Artificial sea salt",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["control"],
            markersize=8,
            label="No added salt",
        ),
    ]
    axis.legend(
        handles=handles, fontsize=8.5, loc="lower left", framealpha=0.95, ncol=2
    )

    figure.text(
        0.5,
        0.005,
        "Darker shades are more concentrated brines. Vertical axis is logarithmic; trends are "
        "drawn only within one salt, never across salts.\n"
        "Across this range the measurements fall by roughly four orders of magnitude while the "
        "model stays nearly flat.",
        ha="center",
        fontsize=9,
        color=PALETTE["guide"],
    )
    figure.tight_layout(rect=[0, 0.07, 1, 1])
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Figure 3 -- carbon dioxide, with its caveat
# ═════════════════════════════════════════════════════════════════════


def plot_carbon_dioxide(paired, ecsv_glob, output_path):
    """Measured carbon dioxide, beside the model's carbon production.

    These two panels are deliberately not overlaid, because they are not the
    same quantity and drawing them on one axis would imply a comparison that
    cannot presently be made.

    The reason is in the reaction network. No reaction in it produces dissolved
    carbon dioxide: every carbon-oxidising step yields bicarbonate instead, and
    dissolved carbon dioxide is additionally held out of carbonate equilibrium
    by the deck's decoupling list. Its concentration therefore never moves from
    its initial value, and converting it through Henry's law -- the conversion
    used for methane -- would predict no carbon dioxide accumulation at all.

    What the model does predict is the rise in bicarbonate, which is where the
    respired carbon actually goes. That is shown on the right. Turning it into a
    headspace prediction needs carbonate speciation at brine ionic strength,
    which is exactly the extrapolation the conversion module declines to make.

    So the left panel is the measurement, the right panel is the model's carbon
    production, and the gap between them is a task, not a result.
    """
    measured = load_measured(ecsv_glob, "CO2")

    figure, (left, right) = plt.subplots(1, 2, figsize=(14, 5.5))

    for entry in paired:
        batch = entry["batch"]
        colour = colour_for_brine(batch["Brine Name"])
        label = f"{batch['Brine Name']} ({batch['Measured Water Activity']:.3f})"

        points = measured[
            (measured["Experiment"] == batch["Experiment"])
            & (measured["Batch ID"] == batch["Batch ID"])
        ]
        if len(points):
            left.scatter(
                points["Days since start"],
                points["Cumulative Moles"],
                color=colour,
                s=22,
                edgecolor="white",
                linewidth=0.4,
                label=label,
            )

        run_frame = entry["model"]
        right.plot(
            run_frame["Time [d]"],
            run_frame["Total HCO3- [M]"],
            color=colour,
            linewidth=2,
            label=label,
        )

    left.set_yscale("log")
    left.set_xlabel("Days since start of incubation")
    left.set_ylabel("Carbon dioxide in the headspace (moles)")
    left.set_title("Measured: carbon dioxide reaching the detector", fontsize=12)
    left.grid(True, alpha=0.25, linewidth=0.5)

    right.set_xlabel("Days since start of simulation")
    right.set_ylabel("Bicarbonate in the liquid (moles per litre)")
    right.set_title("Model: where the respired carbon goes", fontsize=12)
    right.grid(True, alpha=0.25, linewidth=0.5)

    # One shared key rather than the same legend twice, so neither panel has
    # its data covered by a box repeating what the other already said.
    handles, labels = left.get_legend_handles_labels()
    seen, unique = set(), []
    for handle, label in zip(handles, labels):
        if label not in seen:
            seen.add(label)
            unique.append((handle, label))
    figure.legend(
        [h for h, _ in unique],
        [text for _, text in unique],
        fontsize=8,
        ncol=8,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.115),
        framealpha=0.95,
    )

    figure.suptitle(
        "Carbon dioxide: the two panels are not the same quantity and are not overlaid",
        fontsize=13,
        y=0.99,
    )
    figure.text(
        0.5,
        0.005,
        "No reaction in the network produces dissolved carbon dioxide; every carbon-oxidising "
        "step yields bicarbonate, and dissolved\ncarbon dioxide is held out of carbonate "
        "equilibrium by the deck. It never moves from its starting value, so the model cannot "
        "yet\npredict a headspace carbon dioxide concentration. The right panel shows the carbon "
        "the model does produce.",
        ha="center",
        fontsize=9,
        color=PALETTE["guide"],
    )
    figure.tight_layout(rect=[0, 0.22, 1, 0.95])
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Build the model-versus-measurement comparison figures."
    )
    parser.add_argument(
        "--composition",
        default=os.path.join("data", "incubation_batch_composition.csv"),
    )
    parser.add_argument("--run-root", default="runs")
    parser.add_argument(
        "--ecsv-glob",
        default=os.path.expanduser(
            "~/Documents/GitHub/saltyBiomass/data/transformed/*Exp00[34]*.ecsv"
        ),
    )
    parser.add_argument("--output-dir", default=os.path.join("output", "comparison"))
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    composition = pd.read_csv(args.composition)
    paired = assemble(composition, args.run_root, args.ecsv_glob)

    if not paired:
        raise SystemExit("No model runs paired with measurements; nothing to plot.")

    written = [
        plot_methane_timeseries(
            paired, os.path.join(args.output_dir, "methane_over_time.png")
        ),
        plot_methane_against_water_activity(
            paired, os.path.join(args.output_dir, "methane_vs_water_activity.png")
        ),
        plot_carbon_dioxide(
            paired, args.ecsv_glob, os.path.join(args.output_dir, "carbon_dioxide.png")
        ),
    ]

    print(f"Paired {len(paired)} batches.")
    for path in written:
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
