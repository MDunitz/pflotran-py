"""Draw what the package does, generically, so the diagram cannot go stale.

A hand-drawn architecture picture goes out of date the first time somebody adds
a module and forgets the slide. This one is generated, lives beside the code and
names real files, so a wrong box is a wrong line of Python rather than a
forgotten redraw.

The diagram is about the wrapper, not about any one study. PFLOTRAN is driven by
hand-written text decks; running a hundred variants of one means editing a
hundred files and remembering which was which. What this package adds is a loop
you can run programmatically: describe a simulation in Python, generate a deck
per condition, execute them reproducibly in a container, read the output back as
a table, and feed what you learn into the next set of conditions.

Everything on the diagram is reusable. The salt-inhibition study is shown only
as a worked example, with its three experiment-specific modules named so a new
user knows exactly which parts to replace.

Run with:  python scripts/draw_architecture.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

# One hue per role, so the reader learns the key once: what you supply, what the
# package does for you, what is optional, and what is somebody else's software.
COLOURS = {
    "supply": "#0072B2",  # you write this
    "core": "#009E73",  # the package does this
    "optional": "#CC79A7",  # opt in when you need it
    "external": "#666666",  # PFLOTRAN, the container, the database
    "example": "#D55E00",  # one instantiation, replaceable
}


def box(axis, x, y, title, subtitle, role, width=3.6, height=0.92, fontsize=9.5):
    colour = COLOURS[role]
    axis.add_patch(
        FancyBboxPatch(
            (x - width / 2, y - height / 2),
            width,
            height,
            boxstyle="round,pad=0.06",
            linewidth=1.7,
            edgecolor=colour,
            facecolor=colour,
            alpha=0.13,
            zorder=2,
        )
    )
    axis.text(
        x,
        y + 0.17,
        title,
        ha="center",
        va="center",
        fontsize=fontsize,
        weight="bold",
        color=colour,
        zorder=3,
    )
    axis.text(
        x,
        y - 0.19,
        subtitle,
        ha="center",
        va="center",
        fontsize=7.3,
        color="#333333",
        zorder=3,
        family="monospace",
    )
    return (x, y, width, height)


def arrow(axis, start, end, colour="#555555", dashed=False, curve=0.0, width=1.5):
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=width,
            color=colour,
            linestyle="--" if dashed else "-",
            connectionstyle=f"arc3,rad={curve}",
            shrinkA=3,
            shrinkB=3,
            zorder=1,
        )
    )


def draw(output_path):
    figure, axis = plt.subplots(figsize=(17, 11.5))
    axis.set_xlim(0, 17)
    axis.set_ylim(0, 11.5)
    axis.axis("off")

    mid = 7.2  # the main vertical flow
    side = 2.2  # things that feed into it from outside

    # ── 1. What you write ────────────────────────────────────────────
    axis.text(
        mid,
        10.85,
        "1.  DESCRIBE THE SIMULATION IN PYTHON",
        ha="center",
        fontsize=11.5,
        weight="bold",
        color=COLOURS["supply"],
    )
    s1 = box(
        axis,
        2.6,
        9.95,
        "Chemistry",
        "species, reactions,\nrate constants, inhibition",
        "supply",
        width=4.0,
    )
    s2 = box(
        axis,
        mid,
        9.95,
        "Domain and conditions",
        "1D / 2D / 3D column, or closed batch\ngeometry, time, initial state",
        "supply",
        width=4.4,
    )
    s3 = box(
        axis,
        11.8,
        9.95,
        "What varies",
        "one deck per condition:\nsalinity, temperature, anything",
        "supply",
        width=4.0,
    )

    # ── 2. Generate ──────────────────────────────────────────────────
    g = box(
        axis,
        mid,
        8.25,
        "2.  GENERATE DECKS",
        "generator/   ->   N x .in files",
        "core",
        width=5.4,
        height=1.05,
        fontsize=10.5,
    )

    e1 = box(
        axis,
        side,
        8.25,
        "Reactions PFLOTRAN lacks",
        "sandbox/*.F90\nscripts/patch_pflotran_sandboxes.py",
        "optional",
        width=3.5,
    )
    axis.text(
        side,
        7.55,
        "optional, compiled into the binary",
        ha="center",
        fontsize=7.2,
        color=COLOURS["optional"],
        style="italic",
    )

    # ── 3. Run ───────────────────────────────────────────────────────
    r = box(
        axis,
        mid,
        6.5,
        "3.  RUN THEM REPRODUCIBLY",
        "run_decks.py   ->   Containerfile (Docker / Podman)\n"
        "one directory per deck, logs kept, failures reported not raised",
        "core",
        width=6.0,
        height=1.25,
        fontsize=10.5,
    )

    ext = box(
        axis,
        side,
        6.5,
        "PFLOTRAN and its database",
        "hanford.dat, or your own",
        "external",
        width=3.5,
    )

    # ── 4. Extract ───────────────────────────────────────────────────
    x = box(
        axis,
        mid,
        4.85,
        "4.  READ THE OUTPUT AS A TABLE",
        "analysis/extract.py (.tec)  |  extract_hdf5.py (.h5)",
        "core",
        width=6.0,
        height=1.05,
        fontsize=10.5,
    )

    # ── 5. Analyse and visualise ─────────────────────────────────────
    # Placed hard left, where the arrow fan below does not reach.
    axis.text(
        0.35,
        4.05,
        "5.  ANALYSE AND\n     VISUALISE",
        ha="left",
        va="center",
        fontsize=11,
        weight="bold",
        color=COLOURS["core"],
    )

    v1 = box(
        axis,
        2.4,
        2.95,
        "Spatial",
        "gradients, Fick flux\n3D scatter, surface maps",
        "core",
        width=3.4,
    )
    v2 = box(
        axis,
        6.2,
        2.95,
        "Temporal",
        "concentration and flux\nover time at a point",
        "core",
        width=3.4,
    )
    v3 = box(
        axis,
        10.0,
        2.95,
        "Across runs",
        "compare/\none condition against another",
        "core",
        width=3.4,
    )
    v4 = box(
        axis, 13.8, 2.95, "Against measurements", "comparison/", "optional", width=3.4
    )

    axis.text(
        2.4,
        2.28,
        "needs more than one cell",
        ha="center",
        fontsize=7.2,
        color=COLOURS["external"],
        style="italic",
    )
    axis.text(
        6.2,
        2.28,
        "works for any domain",
        ha="center",
        fontsize=7.2,
        color=COLOURS["external"],
        style="italic",
    )

    c1 = box(
        axis,
        12.0,
        1.25,
        "Bringing model and instrument into the same units",
        "headspace.py:  aqueous concentration  ->  what an instrument samples\n"
        "calibrate.py / forecast.py:  fit on one split, hold out another",
        "optional",
        width=6.4,
        height=1.15,
        fontsize=8.6,
    )

    # ── The loop that makes it a wrapper rather than a script ────────
    arrow(
        axis,
        (15.6, 2.95),
        (13.85, 9.95),
        colour=COLOURS["supply"],
        curve=-0.34,
        width=2.0,
    )
    axis.text(
        16.85,
        6.4,
        "sweep,\ncalibrate,\niterate",
        ha="center",
        fontsize=8.8,
        color=COLOURS["supply"],
        weight="bold",
    )

    # ── Arrows ───────────────────────────────────────────────────────
    for supplied in (s1, s2, s3):
        arrow(
            axis,
            (supplied[0], supplied[1] - supplied[3] / 2),
            (g[0] + (supplied[0] - mid) * 0.3, g[1] + g[3] / 2),
        )
    arrow(
        axis,
        (e1[0] + e1[2] / 2, e1[1]),
        (g[0] - g[2] / 2, g[1]),
        colour=COLOURS["optional"],
    )
    arrow(axis, (g[0], g[1] - g[3] / 2), (r[0], r[1] + r[3] / 2))
    arrow(
        axis,
        (ext[0] + ext[2] / 2, ext[1]),
        (r[0] - r[2] / 2, r[1]),
        colour=COLOURS["external"],
    )
    arrow(axis, (r[0], r[1] - r[3] / 2), (x[0], x[1] + x[3] / 2))
    for view in (v1, v2, v3, v4):
        arrow(axis, (x[0], x[1] - x[3] / 2), (view[0], view[1] + view[3] / 2))
    arrow(
        axis,
        (v4[0], v4[1] - v4[3] / 2),
        (c1[0] + 2.0, c1[1] + c1[3] / 2),
        colour=COLOURS["optional"],
    )

    # ── Worked example, clearly marked as the replaceable part ───────
    axis.add_patch(
        FancyBboxPatch(
            (0.3, 0.18),
            8.3,
            1.35,
            boxstyle="round,pad=0.08",
            linewidth=1.5,
            edgecolor=COLOURS["example"],
            facecolor=COLOURS["example"],
            alpha=0.08,
            zorder=2,
        )
    )
    axis.text(
        0.62,
        1.28,
        "Worked example shipped with the package: salt inhibition of methanogenesis",
        fontsize=8.8,
        weight="bold",
        color=COLOURS["example"],
        va="center",
    )
    axis.text(
        0.62,
        0.72,
        "Three modules are specific to that study and are the ones to replace:\n"
        "  comparison/brines.py     reads its sample sheet\n"
        "  comparison/corrections.py  its documented data fixes\n"
        "  comparison/decks.py      turns its sample table into decks",
        fontsize=7.4,
        color="#333333",
        va="center",
        family="monospace",
    )

    figure.suptitle(
        "pflotran-py: describe a simulation in Python, run many of them, get tables and figures back",
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
