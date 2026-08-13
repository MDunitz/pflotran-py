#!/usr/bin/env python3
"""Regenerate per-starting-C companion figures without a local ``runs/`` tree.

Uses the committed forecast ``per_timepoint.csv`` files for modelled methane
at the sampling days, and the saltyBiomass ``.ecsv`` outputs for measured
CH4 / CO2. Carbon dioxide model curves are omitted here -- they need the
HDF5 runs; the CO2 companion shows measurements only until decks are re-run.

Prefer ``python -m pflotran_py.comparison.figures`` once ``runs/`` is populated.
"""

from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd
from astropy import units as u

from pflotran_py.comparison.carbon_inventory import (
    default_comparison_starting_carbon_moles,
)
from pflotran_py.comparison.figures import (
    assemble,
    plot_carbon_dioxide,
    plot_measured_against_modelled,
    plot_methane_timeseries,
)
from pflotran_py.comparison.headspace import DEFAULT_BOTTLE


def _gas_conc_from_headspace_moles(moles):
    """Invert ``gas_phase_moles`` so figure helpers can read a gas-phase column."""
    volume_m3 = DEFAULT_BOTTLE.headspace_volume.to_value(u.m**3)
    return float(moles) / volume_m3


def paired_from_forecast(
    composition_path,
    ecsv_glob,
    forecast_root=os.path.join("output", "comparison", "forecast"),
):
    """Build ``assemble``-like entries from forecast CSVs + measured tables."""
    composition = pd.read_csv(composition_path)
    measured_ch4 = __import__(
        "pflotran_py.comparison.figures", fromlist=["load_measured"]
    ).load_measured(ecsv_glob, "CH4_FID")

    paired = []
    for _, batch in composition.iterrows():
        experiment = batch["Experiment"]
        batch_id = int(batch["Batch ID"])
        forecast_csv = os.path.join(
            forecast_root,
            f"{experiment}_fit-first-2-rounds",
            "per_timepoint.csv",
        )
        if not os.path.isfile(forecast_csv):
            continue
        forecast = pd.read_csv(forecast_csv)
        rows = forecast[forecast["batch"] == batch_id].sort_values("day")
        if rows.empty:
            continue

        gas_conc = [_gas_conc_from_headspace_moles(m) for m in rows["modelled"]]
        model = pd.DataFrame(
            {
                "Time [d]": rows["day"].to_numpy(),
                "Total CH4(aq) [M]": np.full(len(rows), 1e-12),
                "Active_Gas_CH4(g) [mol_m^3 gas]": gas_conc,
                # No CO2(aq): CO2 companion falls back to measurements only.
            }
        )
        points = measured_ch4[
            (measured_ch4["Experiment"] == experiment)
            & (measured_ch4["Batch ID"] == batch_id)
        ]
        paired.append(
            {
                "batch": batch,
                "model": model,
                "measured": points,
                "name": (
                    f"{experiment}_B{batch_id:02d}_{batch['Brine Name']}"
                ),
            }
        )
    return paired


def main():
    composition = os.path.join("data", "incubation_batch_composition.csv")
    ecsv_glob = os.path.expanduser(
        "~/Documents/GitHub/saltyBiomass/data/transformed/*Exp00[34]*.ecsv"
    )
    output_dir = os.path.join("output", "comparison")
    os.makedirs(output_dir, exist_ok=True)
    starting_c = default_comparison_starting_carbon_moles()

    # Prefer real HDF5 runs when present.
    if os.path.isdir("runs") and any(os.scandir("runs")):
        composition_df = pd.read_csv(composition)
        paired = assemble(composition_df, "runs", ecsv_glob)
        source = "runs/"
    else:
        if not glob.glob(ecsv_glob):
            raise SystemExit(f"No measured tables matched {ecsv_glob}")
        paired = paired_from_forecast(composition, ecsv_glob)
        source = "forecast CSVs (CH4 only; CO2 measured-only)"

    if not paired:
        raise SystemExit("No batches available to plot.")

    written = [
        plot_methane_timeseries(
            paired,
            os.path.join(output_dir, "methane_over_time_per_starting_c.png"),
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
        plot_carbon_dioxide(
            paired,
            ecsv_glob,
            os.path.join(output_dir, "carbon_dioxide_per_starting_c.png"),
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
        plot_measured_against_modelled(
            paired,
            ecsv_glob,
            os.path.join(output_dir, "measured_vs_modelled_per_starting_c.png"),
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
    ]
    print(f"Source: {source}")
    print(f"Starting carbon: {starting_c:.6f} mol C")
    print(f"Paired {len(paired)} batches.")
    for path in written:
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
