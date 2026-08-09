"""Draw the package's data flow, so the diagram cannot drift from the code.

An architecture picture drawn by hand goes stale the first time somebody adds a
module and forgets the slide. This one is generated, lives next to the code, and
names real files, so a wrong box is a wrong line of Python rather than a
forgotten redraw.

The flow has two inputs that stay separate for a long time and then meet once.
The left column is what the laboratory measured; the right column is what the
model predicts. They are deliberately not joined until the units agree, which
happens in exactly one place, and that place is worth being able to point at.

Run with:  python scripts/draw_architecture.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

# Roles, not decoration. Measurement, model, the join, and the questions asked
# of the joined result each get their own hue.
COLOURS = {
    "measured": "#0072B2",
    "model": "#D55E00",
    "join": "#009E73",
    "question": "#CC79A7",
    "external": "#666666",
    "unused": "#BBBBBB",
}

BOX_WIDTH = 3.5
BOX_HEIGHT = 0.82


def box(axis, x, y, title, subtitle, role, width=BOX_WIDTH, height=BOX_HEIGHT):
    """One step in the flow: what it is, and which file does it."""
    colour = COLOURS[role]
    axis.add_patch(
        FancyBboxPatch(
            (x - width / 2, y - height / 2),
            width,
            height,
            boxstyle="round,pad=0.06",
            linewidth=1.6,
            edgecolor=colour,
            facecolor=colour,
            alpha=0.13,
            zorder=2,
        )
    )
    axis.text(
        x,
        y + 0.14,
        title,
        ha="center",
        va="center",
        fontsize=9.5,
        weight="bold",
        color=colour,
        zorder=3,
    )
    axis.text(
        x,
        y - 0.17,
        subtitle,
        ha="center",
        va="center",
        fontsize=7.6,
        color="#333333",
        zorder=3,
        family="monospace",
    )
    return (x, y, width, height)


def arrow(axis, start, end, colour="#666666", style="-|>", dashed=False):
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=13,
            linewidth=1.4,
            color=colour,
            linestyle="--" if dashed else "-",
            shrinkA=2,
            shrinkB=2,
            zorder=1,
        )
    )


def draw(output_path):
    figure, axis = plt.subplots(figsize=(15, 11.5))
    axis.set_xlim(0, 15)
    axis.set_ylim(0, 11.5)
    axis.axis("off")

    left, right = 3.4, 11.0
    centre = (left + right) / 2

    axis.text(
        left,
        11.15,
        "WHAT THE LABORATORY MEASURED",
        ha="center",
        fontsize=10.5,
        weight="bold",
        color=COLOURS["measured"],
    )
    axis.text(
        right,
        11.15,
        "WHAT THE MODEL PREDICTS",
        ha="center",
        fontsize=10.5,
        weight="bold",
        color=COLOURS["model"],
    )

    # ── Measurement track ────────────────────────────────────────────
    m1 = box(
        axis,
        left,
        10.3,
        "Batch workbook + chromatograms",
        "Google Sheets  |  1786 mzML files",
        "external",
    )
    m2 = box(
        axis,
        left,
        9.1,
        "Measurement pipeline",
        "saltyBiomass: make run-mzml-pipeline",
        "measured",
    )
    m3 = box(
        axis,
        left,
        7.9,
        "Measured gas per bottle",
        "*.ecsv  ->  Cumulative Moles",
        "measured",
    )

    # ── Composition, which feeds the model side ──────────────────────
    c1 = box(
        axis,
        left,
        6.5,
        "Brine recipes -> ion concentrations",
        "comparison/brines.py + corrections.py",
        "measured",
    )

    # ── Model track ──────────────────────────────────────────────────
    p1 = box(
        axis,
        right,
        10.3,
        "Reaction network + sandboxes",
        "generator/pflotran_templates.py, sandbox/*.F90",
        "external",
    )
    p2 = box(
        axis,
        right,
        9.1,
        "Closed-batch deck per condition",
        "generator/bottle_generator.py, comparison/decks.py",
        "model",
    )
    p3 = box(
        axis,
        right,
        7.9,
        "PFLOTRAN in a container",
        "comparison/run_decks.py  ->  Containerfile",
        "model",
    )
    p4 = box(
        axis,
        right,
        6.5,
        "Simulated chemistry over time",
        "sim.h5  ->  analysis/extract_hdf5.py",
        "model",
    )

    # ── The join ─────────────────────────────────────────────────────
    j1 = box(
        axis,
        centre,
        5.0,
        "THE ONE PLACE THE UNITS MEET",
        "comparison/headspace.py   aqueous mol/L -> headspace moles",
        "join",
        width=6.6,
        height=1.25,
    )
    axis.text(
        centre,
        j1[1] - 0.42,
        "Henry's law + Setschenow salting-out. Nothing above this line is "
        "comparable; everything below it is.",
        ha="center",
        fontsize=7.6,
        color=COLOURS["join"],
        style="italic",
        zorder=3,
    )

    # ── Questions asked of the joined result ─────────────────────────
    q_y = 3.2
    q1 = box(
        axis, 2.5, q_y, "Does it agree?", "comparison/figures.py", "question", width=3.0
    )
    q2 = box(
        axis,
        6.2,
        q_y,
        "Does it transfer?",
        "comparison/calibrate.py",
        "question",
        width=3.0,
    )
    q3 = box(
        axis,
        9.9,
        q_y,
        "Can it forecast?",
        "comparison/forecast.py",
        "question",
        width=3.0,
    )
    q4 = box(
        axis,
        13.0,
        q_y,
        "Is the fit honest?",
        "fit_grid.csv, README.md",
        "question",
        width=2.6,
    )

    for x, label in [
        (2.5, "4 figures:\ntimecourse, dose-response,\nCO2 overlay, parity"),
        (6.2, "fit Exp004,\npredict Exp003\ncold"),
        (9.9, "fit early rounds,\npredict late ones"),
        (13.0, "objective surface,\nwhat was tuned"),
    ]:
        axis.text(x, q_y - 0.78, label, ha="center", fontsize=7.2, color="#444444")

    # ── The branch that does not apply here ──────────────────────────
    u1 = box(
        axis,
        right + 0.2,
        1.15,
        "Spatial visualisation (column runs only)",
        "pipeline.py: 3D scatter, surface maps, flux",
        "unused",
        width=5.6,
    )
    axis.text(
        right + 0.2,
        0.52,
        "Needs gradients. A bottle is one cell, so every gradient is zero and\n"
        "these produce nothing for the closed-batch workflow.",
        ha="center",
        fontsize=7.2,
        color=COLOURS["external"],
        style="italic",
    )

    # ── Arrows ───────────────────────────────────────────────────────
    down = [(m1, m2), (m2, m3), (p1, p2), (p2, p3), (p3, p4)]
    for (x1, y1, _, h1), (x2, y2, _, h2) in down:
        arrow(axis, (x1, y1 - h1 / 2), (x2, y2 + h2 / 2))

    arrow(axis, (m3[0], m3[1] - m3[3] / 2), (c1[0], c1[1] + c1[3] / 2))

    # Composition crosses to the model side: this is the dependency people miss.
    arrow(
        axis,
        (c1[0] + c1[2] / 2, c1[1] + 0.15),
        (p2[0] - p2[2] / 2, p2[1] - 0.2),
        colour=COLOURS["measured"],
    )
    axis.text(
        centre,
        7.35,
        "decks are built from the weighed salt,\nnot from a target water activity",
        ha="center",
        fontsize=7.4,
        color=COLOURS["measured"],
        style="italic",
    )

    arrow(
        axis,
        (c1[0], c1[1] - c1[3] / 2),
        (j1[0] - 1.6, j1[1] + j1[3] / 2),
        colour=COLOURS["measured"],
    )
    arrow(
        axis,
        (p4[0], p4[1] - p4[3] / 2),
        (j1[0] + 1.6, j1[1] + j1[3] / 2),
        colour=COLOURS["model"],
    )

    for q in (q1, q2, q3, q4):
        arrow(
            axis,
            (j1[0], j1[1] - j1[3] / 2),
            (q[0], q[1] + q[3] / 2),
            colour=COLOURS["join"],
        )

    arrow(
        axis,
        (p4[0] + p4[2] / 2, p4[1]),
        (u1[0] + 1.6, u1[1] + u1[3] / 2),
        colour=COLOURS["unused"],
        dashed=True,
    )

    figure.suptitle(
        "pflotran-py: how a measured bottle and a simulated one are made comparable",
        fontsize=13.5,
        y=0.985,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.965])
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


def main():
    directory = os.path.join("output", "docs")
    os.makedirs(directory, exist_ok=True)
    path = draw(os.path.join(directory, "architecture.png"))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
