"""Figures comparing the closed-batch model against the measured incubations.

Four absolute-mole figures, each answering a different question, plus companion
figures that divide both sides by the model's starting carbon inventory.

**Methane over time** asks whether the model reproduces the measured trajectory
in each bottle. Measured replicates are drawn as points, the model as a line,
one panel per salt family, on a logarithmic vertical axis because the measured
values span more than four orders of magnitude across conditions.

**Methane against water activity** asks the question the experiment was designed
around: does production fall as salt rises, and does the model fall with it?
This is the figure that shows whether the inhibition in the model matches the
inhibition in the bottles.

**Carbon dioxide** overlays measured and modelled on shared axes, in the same
layout as the methane figure. Decks built without coupled carbonate cannot be
drawn this way and fall back to the measurements alone with the reason stated
on the figure. See :func:`plot_carbon_dioxide`.

**Measured against modelled** puts the two on opposite axes with the line of
equality drawn, one point per batch condition. It is the most direct reading of
how well the model does and in which direction it errs.

**Per starting C companions** re-express the methane and carbon dioxide
timeseries, and the parity figure, as moles of headspace gas per mole of
starting carbon. The denominator matches the incubations' recipe-derived
biomass C (~0.0565 mol; see :mod:`.carbon_inventory` and the README "Starting
carbon" section). Prefer these for yield questions; absolute-mole overlays
remain useful for trajectory shape and salt ranking.

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
from .carbon_inventory import (  # noqa: E402
    default_comparison_starting_carbon_moles,
)
from .headspace import (  # noqa: E402
    GASES,
    aqueous_concentration_to_headspace_moles,
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


def _per_starting_c(values, starting_carbon_moles):
    """Divide headspace moles by the model's starting carbon inventory."""
    return np.asarray(values, dtype=float) / float(starting_carbon_moles)


def _gas_amount_label(gas_name, per_starting_c):
    if per_starting_c:
        return f"{gas_name} in the bottle headspace (mol per mol starting C)"
    return f"{gas_name} in the bottle headspace (moles)"


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


GAS_PHASE_METHANE_COLUMN = "Active_Gas_CH4(g) [mol_m^3 gas]"


def model_methane_headspace(run_frame, batch_row):
    """Model methane in the headspace over time, in moles.

    Reads the gas phase directly when the deck carries one, which is the
    honest route: PFLOTRAN has done the partition itself and we are simply
    reporting it. Falls back to partitioning the total inventory for older
    decks that model methane departure with the ebullition proxy instead.
    """
    nacl, mgcl2 = setschenow_salts_from_composition(batch_row)
    has_gas_phase = GAS_PHASE_METHANE_COLUMN in run_frame.columns

    moles = []
    for _, row in run_frame.iterrows():
        moles.append(
            model_headspace_moles(
                row["Total CH4(aq) [M]"],
                GASES["CH4"],
                degassed_concentration=(
                    None if has_gas_phase else row.get("Total Tracer2 [M]")
                ),
                gas_concentration=(
                    row[GAS_PHASE_METHANE_COLUMN] if has_gas_phase else None
                ),
                nacl_molarity=nacl,
                mgcl2_molarity=mgcl2,
            ).to_value(u.mol)
        )
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


def plot_methane_timeseries(
    paired,
    output_path,
    *,
    per_starting_c=False,
    starting_carbon_moles=None,
):
    """Measured and modelled methane in the headspace, over the incubation."""
    if per_starting_c and starting_carbon_moles is None:
        starting_carbon_moles = default_comparison_starting_carbon_moles()

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
            y_model = (
                _per_starting_c(moles, starting_carbon_moles)
                if per_starting_c
                else moles
            )
            axis.plot(days, y_model, color=colour, linewidth=2, alpha=0.9, zorder=2)

            points = entry["measured"]
            if len(points):
                y_meas = points["Cumulative Moles"]
                if per_starting_c:
                    y_meas = _per_starting_c(y_meas, starting_carbon_moles)
                axis.scatter(
                    points["Days since start"],
                    y_meas,
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
        y_lo, y_hi = (1e-9, 5e-3)
        if per_starting_c:
            y_lo, y_hi = y_lo / starting_carbon_moles, y_hi / starting_carbon_moles
        axis.set_ylim(y_lo, y_hi)
        axis.set_xlim(-3, 133)
        axis.set_xlabel("Days since start of incubation")
        axis.set_title(family, fontsize=12)
        axis.grid(True, alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=7.5, loc="lower right", framealpha=0.95)

    axes[0].set_ylabel(_gas_amount_label("Methane", per_starting_c))

    if per_starting_c:
        title = (
            "Modelled and measured methane per mole of starting carbon "
            "in sealed incubation bottles"
        )
        footnote = (
            "Lines are the reactive-transport model; filled circles are gas-chromatograph "
            "measurements of individual bottles.\n"
            f"Both sides are divided by the model's starting carbon inventory "
            f"({starting_carbon_moles:.3f} mol C; cellulose hydrolysis pool). "
            "Vertical axis is logarithmic."
        )
    else:
        title = "Modelled and measured methane in sealed incubation bottles"
        footnote = (
            "Lines are the reactive-transport model; filled circles are gas-chromatograph "
            "measurements of individual bottles.\n"
            "Colour identifies the brine; darker shades are more concentrated. Vertical axis "
            "is logarithmic."
        )
    figure.suptitle(title, fontsize=14, y=0.99)
    figure.text(
        0.5,
        0.005,
        footnote,
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
        "The model now collapses under salt, but too abruptly and onto a floor",
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
        "The measured decline is gradual across the whole range. The modelled one is steeper "
        "between water activity 0.96 and 0.90, then flattens onto a residual floor near "
        "5e-8 moles.",
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


COUPLED_CARBON_DIOXIDE_COLUMN = "CO2(aq) [M]"


def model_carbon_dioxide_series(run_frame, batch_row):
    """Model headspace carbon dioxide over time, or None if unavailable.

    Returns None when the deck was built without coupled carbonate. In that
    case the model carries dissolved carbon dioxide as an independent primary
    species that no reaction produces, so it never moves from its initial value
    and there is nothing to plot against a measurement.
    """
    if COUPLED_CARBON_DIOXIDE_COLUMN not in run_frame.columns:
        return None, None
    nacl, mgcl2 = setschenow_salts_from_composition(batch_row)
    moles = [
        aqueous_concentration_to_headspace_moles(
            row[COUPLED_CARBON_DIOXIDE_COLUMN],
            GASES["CO2"],
            nacl_molarity=nacl,
            mgcl2_molarity=mgcl2,
        ).to_value(u.mol)
        for _, row in run_frame.iterrows()
    ]
    return run_frame["Time [d]"].to_numpy(), np.array(moles)


def model_carbon_dioxide_headspace(run_frame, batch_row):
    """Model headspace carbon dioxide at the end of the run, or None."""
    _, series = model_carbon_dioxide_series(run_frame, batch_row)
    return None if series is None else series[-1]


def plot_carbon_dioxide(
    paired,
    ecsv_glob,
    output_path,
    *,
    per_starting_c=False,
    starting_carbon_moles=None,
):
    """Measured and modelled carbon dioxide in the headspace, over the incubation.

    Overlaid on shared axes, in the same layout as the methane figure, because
    with coupled carbonate the two are now the same quantity: PFLOTRAN speciates
    internally and reports the neutral dissolved carbon dioxide, which Henry's
    law turns into the headspace moles a gas chromatograph would sample.

    Decks built without coupled carbonate cannot be drawn this way. There the
    model holds dissolved carbon dioxide as an independent primary species that
    no reaction produces -- every carbon-oxidising step in the network yields
    bicarbonate -- so its concentration never moves and an overlay would show a
    flat line that means nothing. Those runs fall back to a single panel of the
    measurements alone, with the reason stated on the figure rather than left
    for the reader to infer from a suspiciously horizontal curve.
    """
    if per_starting_c and starting_carbon_moles is None:
        starting_carbon_moles = default_comparison_starting_carbon_moles()

    measured = load_measured(ecsv_glob, "CO2")
    predictable = [
        entry
        for entry in paired
        if COUPLED_CARBON_DIOXIDE_COLUMN in entry["model"].columns
    ]

    def _y(values):
        if per_starting_c:
            return _per_starting_c(values, starting_carbon_moles)
        return values

    if not predictable:
        figure, axis = plt.subplots(figsize=(8, 5.5))
        for entry in paired:
            batch = entry["batch"]
            points = measured[
                (measured["Experiment"] == batch["Experiment"])
                & (measured["Batch ID"] == batch["Batch ID"])
            ]
            if len(points):
                axis.scatter(
                    points["Days since start"],
                    _y(points["Cumulative Moles"]),
                    color=colour_for_brine(batch["Brine Name"]),
                    s=24,
                    edgecolor="white",
                    linewidth=0.4,
                )
        axis.set_yscale("log")
        axis.set_xlabel("Days since start of incubation")
        axis.set_ylabel(_gas_amount_label("Carbon dioxide", per_starting_c))
        if per_starting_c:
            axis.set_title(
                "Measured carbon dioxide per mole of starting carbon",
                fontsize=12,
            )
            caveat = (
                f"Measurements only — modelled CO2 curves need the HDF5 runs "
                f"(re-run comparison.figures after run_decks).\n"
                f"Values are divided by the model's starting carbon inventory "
                f"({starting_carbon_moles:.3f} mol C; cellulose hydrolysis pool)."
            )
        else:
            axis.set_title(
                "Measured carbon dioxide; the model cannot predict it", fontsize=12
            )
            caveat = (
                "These decks were built without coupled carbonate, so dissolved carbon dioxide "
                "never moves from its\ninitial value and there is no modelled curve to draw. "
                "Rebuild with couple_carbonate to compare."
            )
        axis.grid(True, alpha=0.25, linewidth=0.5)
        figure.text(
            0.5,
            0.01,
            caveat,
            ha="center",
            fontsize=9,
            color=PALETTE["guide"],
        )
        figure.tight_layout(rect=[0, 0.09, 1, 1])
        figure.savefig(output_path, dpi=200)
        plt.close(figure)
        return output_path

    families = ["Sodium chloride", "Magnesium chloride", "Artificial sea salt"]
    figure, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)

    for axis, family in zip(axes, families):
        members = [
            entry
            for entry in predictable
            if salt_family(entry["batch"]["Brine Name"]) == family
        ]
        experiments = {entry["batch"]["Experiment"] for entry in members}
        controls = [
            entry
            for entry in predictable
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

            days, moles = model_carbon_dioxide_series(entry["model"], batch)
            if days is not None:
                axis.plot(
                    days, _y(moles), color=colour, linewidth=2, alpha=0.9, zorder=2
                )

            points = measured[
                (measured["Experiment"] == batch["Experiment"])
                & (measured["Batch ID"] == batch["Batch ID"])
            ]
            if len(points):
                axis.scatter(
                    points["Days since start"],
                    _y(points["Cumulative Moles"]),
                    color=colour,
                    s=26,
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=3,
                    label=label,
                )

        axis.set_yscale("log")
        y_lo, y_hi = (1e-7, 3e-3)
        if per_starting_c:
            y_lo, y_hi = y_lo / starting_carbon_moles, y_hi / starting_carbon_moles
        axis.set_ylim(y_lo, y_hi)
        axis.set_xlim(-3, 133)
        axis.set_xlabel("Days since start of incubation")
        axis.set_title(family, fontsize=12)
        axis.grid(True, alpha=0.25, linewidth=0.5)
        axis.legend(fontsize=7.5, loc="lower right", framealpha=0.95)

    axes[0].set_ylabel(_gas_amount_label("Carbon dioxide", per_starting_c))

    if per_starting_c:
        title = (
            "Modelled and measured carbon dioxide per mole of starting carbon "
            "in sealed incubation bottles"
        )
        footnote = (
            "Lines are the reactive-transport model; filled circles are gas-chromatograph "
            "measurements of individual bottles.\n"
            f"Both sides are divided by the model's starting carbon inventory "
            f"({starting_carbon_moles:.3f} mol C; cellulose hydrolysis pool). "
            "No parameter was fitted against carbon dioxide."
        )
    else:
        title = "Modelled and measured carbon dioxide in sealed incubation bottles"
        footnote = (
            "Lines are the reactive-transport model; filled circles are gas-chromatograph "
            "measurements of individual bottles.\n"
            "No parameter anywhere in this model was fitted against carbon dioxide, so both the "
            "level and the shape here are predictions."
        )
    figure.suptitle(title, fontsize=14, y=0.99)
    figure.text(
        0.5,
        0.005,
        footnote,
        ha="center",
        fontsize=9,
        color=PALETTE["guide"],
    )
    figure.tight_layout(rect=[0, 0.06, 1, 0.96])
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


# ═════════════════════════════════════════════════════════════════════
# Figure 4 -- measured against modelled, directly
# ═════════════════════════════════════════════════════════════════════


def plot_measured_against_modelled(
    paired,
    ecsv_glob,
    output_path,
    *,
    per_starting_c=False,
    starting_carbon_moles=None,
):
    """Measured on one axis, modelled on the other, with the line of equality.

    The most direct reading of the comparison. A point on the diagonal is a
    condition the model gets right; distance from the diagonal is how wrong it
    is, and in which direction. Because both axes are logarithmic and span four
    orders of magnitude, the shaded band marks agreement within a factor of ten
    -- generous, but the honest resolution of a comparison whose inputs carry
    the uncertainties described in the module docstring.
    """
    if per_starting_c and starting_carbon_moles is None:
        starting_carbon_moles = default_comparison_starting_carbon_moles()

    measured_co2 = load_measured(ecsv_glob, "CO2")

    figure, (left, right) = plt.subplots(1, 2, figsize=(13.5, 6.4))

    panels = [
        (left, "Methane", None),
        (right, "Carbon dioxide", measured_co2),
    ]

    for axis, gas_label, co2_frame in panels:
        pairs = []
        for entry in paired:
            batch = entry["batch"]
            colour = colour_for_brine(batch["Brine Name"])

            if gas_label == "Methane":
                _, series = model_methane_headspace(entry["model"], batch)
                modelled = series[-1]
                points = entry["measured"]
            else:
                modelled = model_carbon_dioxide_headspace(entry["model"], batch)
                points = co2_frame[
                    (co2_frame["Experiment"] == batch["Experiment"])
                    & (co2_frame["Batch ID"] == batch["Batch ID"])
                ]
            if modelled is None or not len(points):
                continue

            final = (
                points.sort_values("Days since start")
                .groupby("Replicate ID")["Cumulative Moles"]
                .last()
            )
            measured = float(np.median(final))
            if measured <= 0 or modelled <= 0:
                continue

            if per_starting_c:
                measured = measured / starting_carbon_moles
                modelled = modelled / starting_carbon_moles

            pairs.append((measured, modelled))
            axis.scatter(
                measured,
                modelled,
                color=colour,
                s=95,
                edgecolor="white",
                linewidth=0.9,
                zorder=3,
            )

        if not pairs:
            continue

        values = np.array(pairs)
        low = min(values.min() * 0.3, 1e-8 if not per_starting_c else 1e-8 / starting_carbon_moles)
        high = values.max() * 3
        line = np.array([low, high])

        axis.fill_between(
            line,
            line / 10,
            line * 10,
            color=PALETTE["guide_light"],
            alpha=0.28,
            zorder=0,
        )
        axis.plot(line, line, color=PALETTE["guide"], linewidth=1.4, zorder=1)

        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(low, high)
        axis.set_ylim(low, high)
        axis.set_aspect("equal")
        unit = "mol per mol starting C" if per_starting_c else "moles"
        axis.set_xlabel(f"Measured {gas_label.lower()} in the headspace ({unit})")
        axis.set_ylabel(f"Modelled {gas_label.lower()} in the headspace ({unit})")

        ratios = np.log10(values[:, 1] / values[:, 0])
        within = np.mean(np.abs(ratios) <= 1)
        axis.set_title(
            f"{gas_label}\n{within:.0%} of conditions within a factor of ten",
            fontsize=12,
        )
        axis.grid(True, alpha=0.22, linewidth=0.5)

        axis.text(
            0.04,
            0.94,
            "model over-predicts",
            transform=axis.transAxes,
            fontsize=8,
            color=PALETTE["guide"],
        )
        axis.text(
            0.55,
            0.05,
            "model under-predicts",
            transform=axis.transAxes,
            fontsize=8,
            color=PALETTE["guide"],
        )

    handles = [
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["control"],
            markersize=9,
            label="No added salt",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["nacl_mid"],
            markersize=9,
            label="Sodium chloride",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["mgcl2_mid"],
            markersize=9,
            label="Magnesium chloride",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=PALETTE["seasalt_mid"],
            markersize=9,
            label="Artificial sea salt",
        ),
        plt.Line2D(
            [], [], linestyle="-", color=PALETTE["guide"], label="Exact agreement"
        ),
    ]
    figure.legend(
        handles=handles,
        fontsize=9,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
    )

    if per_starting_c:
        title = (
            "Measured against modelled production per mole of starting carbon, "
            "one point per batch"
        )
    else:
        title = "Measured against modelled production, one point per batch condition"
    figure.suptitle(title, fontsize=13.5, y=0.98)
    figure.tight_layout(rect=[0, 0.09, 1, 0.95])
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
    parser.add_argument(
        "--starting-carbon-moles",
        type=float,
        default=None,
        help=(
            "Denominator for the per-starting-C companions. Defaults to the "
            "cellulose-hydrolysis inventory used by the comparison pipeline."
        ),
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    composition = pd.read_csv(args.composition)
    paired = assemble(composition, args.run_root, args.ecsv_glob)

    if not paired:
        raise SystemExit("No model runs paired with measurements; nothing to plot.")

    starting_c = args.starting_carbon_moles
    if starting_c is None:
        starting_c = default_comparison_starting_carbon_moles()

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
        plot_measured_against_modelled(
            paired,
            args.ecsv_glob,
            os.path.join(args.output_dir, "measured_vs_modelled.png"),
        ),
        plot_methane_timeseries(
            paired,
            os.path.join(args.output_dir, "methane_over_time_per_starting_c.png"),
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
        plot_carbon_dioxide(
            paired,
            args.ecsv_glob,
            os.path.join(args.output_dir, "carbon_dioxide_per_starting_c.png"),
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
        plot_measured_against_modelled(
            paired,
            args.ecsv_glob,
            os.path.join(args.output_dir, "measured_vs_modelled_per_starting_c.png"),
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
    ]

    print(f"Paired {len(paired)} batches.")
    print(f"Starting carbon inventory for companions: {starting_c:.6f} mol C")
    for path in written:
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()