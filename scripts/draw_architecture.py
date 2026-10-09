"""Render the package architecture from ``docs/architecture.csv``.

Edit the CSV to update box text, positions, roles, or connections, then run:
    python scripts/draw_architecture.py
"""

import csv
from pathlib import Path
import textwrap

import matplotlib
from matplotlib.path import Path as MplPath

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "docs" / "architecture.csv"
OUTPUT_PATH = ROOT / "output" / "docs" / "architecture.png"

ROLE_COLOURS = {
    "input": "#0072B2",
    "core": "#009E73",
    "optional": "#CC79A7",
    "external": "#666666",
    "example": "#D55E00",
}
ROLE_LABELS = {
    "input": "Supplied by the user",
    "core": "pflotran-py core",
    "optional": "Optional capability",
    "external": "External software or data",
    "example": "Study-specific example",
}


def read_diagram(path=DATA_PATH):
    """Read node and edge records from the architecture CSV."""
    with path.open(newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))

    nodes = {}
    edges = []
    for row in records:
        if row["kind"] == "node":
            row.update(
                x=float(row["x"]),
                y=float(row["y"]),
                width=float(row["width"]),
                height=float(row["height"]),
                fontsize=float(row["fontsize"]),
            )
            if row["role"] not in ROLE_COLOURS:
                raise ValueError(f"Unknown role {row['role']!r} for {row['id']}")
            nodes[row["id"]] = row
        elif row["kind"] == "edge":
            row["curve"] = float(row["curve"] or 0)
            edges.append(row)
        else:
            raise ValueError(f"Unknown architecture record: {row['kind']!r}")

    for edge in edges:
        for key in ("source", "target"):
            if edge[key] not in nodes:
                raise ValueError(f"Edge {edge['id']} references missing {key} {edge[key]!r}")
    return nodes, edges


def anchor(node, side):
    """Return the center of one side of a box."""
    return {
        "top": (node["x"], node["y"] + node["height"] / 2),
        "bottom": (node["x"], node["y"] - node["height"] / 2),
        "left": (node["x"] - node["width"] / 2, node["y"]),
        "right": (node["x"] + node["width"] / 2, node["y"]),
    }[side]


def draw_node(axis, node):
    colour = ROLE_COLOURS[node["role"]]
    width, height = node["width"], node["height"]
    axis.add_patch(
        FancyBboxPatch(
            (node["x"] - width / 2, node["y"] - height / 2),
            width,
            height,
            boxstyle="round,pad=0.06",
            linewidth=1.5,
            edgecolor=colour,
            facecolor=colour,
            alpha=0.10,
            zorder=2,
        )
    )
    axis.text(
        node["x"],
        node["y"] + 0.16,
        node["title"],
        ha="center",
        va="center",
        fontsize=node["fontsize"],
        weight="bold",
        color=colour,
        zorder=3,
    )
    subtitle = textwrap.fill(node["subtitle"], width=max(22, round(width * 11)))
    axis.text(
        node["x"],
        node["y"] - 0.17,
        subtitle,
        ha="center",
        va="center",
        fontsize=7.1 if node["role"] != "example" else 7.0,
        color="#333333",
        family="monospace",
        zorder=3,
    )


def draw_legend(axis):
    axis.text(10.05, 1.04, "COLOR KEY", fontsize=8, weight="bold", color="#333333")
    positions = {
        "input": (10.1, 0.72),
        "core": (12.9, 0.72),
        "optional": (10.1, 0.39),
        "external": (12.9, 0.39),
        "example": (10.1, 0.06),
    }
    for role, (x, y) in positions.items():
        colour = ROLE_COLOURS[role]
        axis.add_patch(
            FancyBboxPatch(
                (x, y - 0.08), 0.24, 0.16,
                boxstyle="round,pad=0.02",
                linewidth=0,
                facecolor=colour,
                alpha=0.35,
                zorder=2,
            )
        )
        axis.text(x + 0.34, y, ROLE_LABELS[role], va="center", fontsize=7.6, color="#333333")


def draw(output_path=OUTPUT_PATH, data_path=DATA_PATH):
    nodes, edges = read_diagram(data_path)
    figure, axis = plt.subplots(figsize=(16.5, 10.0))
    axis.set_xlim(0, 16.9)
    axis.set_ylim(-0.10, 9.6)
    axis.axis("off")

    for edge in edges:
        start = anchor(nodes[edge["source"]], edge["source_side"])
        end = anchor(nodes[edge["target"]], edge["target_side"])
        colour = ROLE_COLOURS.get(edge["edge_role"] or "", "#666666")
        if edge["route"] == "outside_right":
            route_x = 15.95
            path = MplPath(
                [start, (route_x, start[1]), (route_x, end[1]), end],
                [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO],
            )
            axis.add_patch(
                FancyArrowPatch(
                    path=path,
                    arrowstyle="-|>",
                    mutation_scale=13,
                    linewidth=1.6 if edge["edge_role"] else 1.35,
                    color=colour,
                    linestyle="--",
                    shrinkA=4,
                    shrinkB=5,
                    zorder=1,
                )
            )
            continue
        axis.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.35,
                color=colour,
                linestyle="--" if edge["edge_role"] else "-",
                connectionstyle=f"arc3,rad={edge['curve']}",
                shrinkA=4,
                shrinkB=5,
                zorder=1,
            )
        )

    for node in nodes.values():
        draw_node(axis, node)

    axis.text(
        8.0,
        9.35,
        "How pflotran-py connects a simulation description to PFLOTRAN results",
        ha="center",
        fontsize=14,
        weight="bold",
        color="#222222",
    )
    axis.text(
        16.15,
        5.5,
        "revise inputs",
        ha="center",
        va="center",
        fontsize=8,
        weight="bold",
        color=ROLE_COLOURS["input"],
        rotation=90,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
    )
    draw_legend(axis)

    figure.tight_layout(pad=0.4)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return output_path


if __name__ == "__main__":
    print(f"wrote {draw()}")
