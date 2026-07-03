"""Shared helpers for PFLOTRAN docker end-to-end integration tests."""

import os
import re
import sys
import glob
import shutil
import tempfile
import subprocess

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HANFORD_DB = os.path.join(REPO_ROOT, "sandbox", "hanford.dat")
CUSTOM_PFLOTRAN = os.path.join(REPO_ROOT, "build", "pflotran")

SPECIES_MAP = {
    "CO2": "CO2(aq) [M]",
    "CH4": "Free CH4(aq) [M]",
}

# Decks that list CO2(aq) as a primary species emit "Free CO2(aq) [M]" instead.
CO2_COLUMN_CANDIDATES = ["CO2(aq) [M]", "Free CO2(aq) [M]"]
CH4_COLUMN_CANDIDATES = ["Free CH4(aq) [M]", "CH4(aq) [M]"]
TEMPERATURE_C = 8.0

_STOCK_CANDIDATES = [
    os.environ.get("PFLOTRAN_EXE"),
    "/scratch/pflotran/src/pflotran/pflotran",
    "/scratch/pflotran/src/pflotran/bin/pflotran",
    shutil.which("pflotran"),
]

_CUSTOM_CANDIDATES = [
    os.environ.get("PFLOTRAN_CUSTOM_EXE"),
    CUSTOM_PFLOTRAN,
    "/opt/pflotran-py/pflotran",
]


def _is_runnable(exe):
    """Return True if *exe* can actually be executed on this machine."""
    # Linux binaries bind-mounted from Docker cannot run natively on macOS.
    if sys.platform == "darwin":
        with open(exe, "rb") as fh:
            if fh.read(4) == b"\x7fELF":
                return False
    try:
        proc = subprocess.run(
            [exe],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=15,
        )
        return "PFLOTRAN" in proc.stdout or proc.returncode != 127
    except (OSError, subprocess.TimeoutExpired):
        return False


def _first_executable(candidates):
    for candidate in candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            if _is_runnable(candidate):
                return candidate
    return None


def find_pflotran():
    """Return path to a stock PFLOTRAN executable, or None."""
    return _first_executable(_STOCK_CANDIDATES)


def find_custom_pflotran():
    """Return path to a PFLOTRAN build with AWINHIBIT sandboxes, or None."""
    return _first_executable(_CUSTOM_CANDIDATES)


def _rewrite_database(text):
    return re.sub(
        r"^\s*DATABASE\s+.*$",
        f"  DATABASE {HANFORD_DB}",
        text,
        flags=re.MULTILINE,
    )


def _strip_water_activity_output(text):
    """Remove CHEMISTRY OUTPUT block that breaks on PFLOTRAN v6.

    Decks with ACTIVITY_WATER often include::

        OUTPUT
          WATER_ACTIVITY_COEFFICIENT
        /

    immediately after ACTIVITY_COEFFICIENTS. PFLOTRAN v6 rejects this as an
    unknown chemistry species. The AWINHIBIT sandboxes still receive water
    activity internally via ACTIVITY_WATER + ACTIVITY_COEFFICIENTS TIMESTEP.
    """
    # Match the chemistry-level OUTPUT block that contains ONLY
    # WATER_ACTIVITY_COEFFICIENT, regardless of the preceding line. The main
    # OUTPUT block (PH, FREE_ION, ...) has other entries and is left intact.
    return re.sub(
        r"\n[ \t]*OUTPUT\n[ \t]*WATER_ACTIVITY_COEFFICIENT\n[ \t]*/\n",
        "\n",
        text,
    )


def prepare_input(src_in, workdir, custom_sandboxes=False):
    """Copy deck to workdir/sim.in with portable paths and v6 fixes."""
    dst_in = os.path.join(workdir, "sim.in")
    with open(src_in, "r") as f:
        text = f.read()
    text = _rewrite_database(text)
    if custom_sandboxes:
        text = _strip_water_activity_output(text)
    with open(dst_in, "w") as f:
        f.write(text)
    return dst_in


def run_pflotran(exe, workdir, prefix="sim"):
    proc = subprocess.run(
        [exe, "-input_prefix", prefix],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=1800,
    )
    log = proc.stdout
    if proc.returncode != 0:
        raise RuntimeError(
            f"PFLOTRAN exited {proc.returncode}\n--- last log lines ---\n"
            + "\n".join(log.splitlines()[-40:])
        )
    return log


def discover_snapshots(workdir, prefix="sim"):
    pattern = os.path.join(workdir, f"{prefix}-[0-9][0-9][0-9].tec")
    return sorted(glob.glob(pattern))


def render_images(flux_df, output_dir, species_map=SPECIES_MAP):
    import shared_utils  # imported by callers after sys.path setup

    os.makedirs(output_dir, exist_ok=True)
    times = sorted(flux_df["Time Index"].unique())
    images = {}

    fig, ax = plt.subplots(figsize=(8, 5))
    for species, col in species_map.items():
        means = [flux_df.loc[flux_df["Time Index"] == t, col].mean() for t in times]
        ax.plot(times, means, marker="o", label=f"{species} mean")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Concentration [M]")
    ax.set_title("Mean aqueous concentration vs time")
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, "concentration_timeseries.png")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    images["concentration_timeseries"] = path

    fig, ax = plt.subplots(figsize=(8, 5))
    for species in species_map:
        mag_col = shared_utils.flux_col(species, "magnitude")
        means = [flux_df.loc[flux_df["Time Index"] == t, mag_col].mean() for t in times]
        ax.plot(times, means, marker="s", label=f"{species} |J| mean")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel(f"Diffusive flux magnitude [{shared_utils.FLUX_UNITS}]")
    ax.set_title("Mean diffusive flux magnitude vs time")
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = os.path.join(output_dir, "flux_magnitude_timeseries.png")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    images["flux_magnitude_timeseries"] = path

    final_t = times[-1]
    final = flux_df[flux_df["Time Index"] == final_t]
    z_mid = sorted(final["Z [m]"].unique())[len(final["Z [m]"].unique()) // 2]
    layer = final[np.isclose(final["Z [m]"], z_mid)]
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(
        layer["X [m]"],
        layer["Y [m]"],
        c=layer[species_map["CH4"]],
        s=400,
        marker="s",
        cmap="viridis",
    )
    fig.colorbar(sc, ax=ax, label="Free CH4(aq) [M]")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title(f"CH4 at z={z_mid:g} m, t={final_t} d")
    path = os.path.join(output_dir, "ch4_final_zslice.png")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    images["ch4_final_zslice"] = path

    return images


def resolve_species_map(df, base_map=None):
    """Pick tecplot column names present in *df* (handles deck-to-deck variation)."""
    base_map = base_map or SPECIES_MAP
    resolved = {}
    for species, default_col in base_map.items():
        if species == "CO2":
            candidates = CO2_COLUMN_CANDIDATES
        elif species == "CH4":
            candidates = CH4_COLUMN_CANDIDATES
        else:
            candidates = [default_col]
        for col in candidates:
            if col in df.columns:
                resolved[species] = col
                break
        else:
            raise KeyError(
                f"No column found for {species}; tried {candidates}. "
                f"Available: {sorted(df.columns)}"
            )
    return resolved


def run_full_pipeline(
    src_in,
    output_dir,
    exe,
    custom_sandboxes=False,
    species_map=SPECIES_MAP,
    temperature_c=TEMPERATURE_C,
    keep_workdir=False,
):
    import step1_extract
    import shared_utils

    os.makedirs(output_dir, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="pflotran_e2e_")
    try:
        prepare_input(src_in, workdir, custom_sandboxes=custom_sandboxes)
        log = run_pflotran(exe, workdir)

        snapshots = discover_snapshots(workdir)
        if not snapshots:
            raise RuntimeError("PFLOTRAN produced no snapshot .tec files.")

        df = step1_extract.extract_pflotran_data_tec(
            data_dir=workdir,
            file_name_template="sim-{:03d}.tec",
            n_files=len(snapshots),
        )
        species_map = resolve_species_map(df, base_map=species_map)
        df = shared_utils.calculate_gradients(df, species_map)
        df = shared_utils.convert_to_flux(
            df, list(species_map.keys()), temperature_c=temperature_c
        )

        csv_path = os.path.join(output_dir, "e2e_pflotran_data.csv")
        df.to_csv(csv_path, index=False)
        images = render_images(df, output_dir, species_map=species_map)

        return {
            "input": src_in,
            "executable": exe,
            "log": log,
            "n_snapshots": len(snapshots),
            "dataframe": df,
            "species_map": species_map,
            "csv": csv_path,
            "images": images,
        }
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
