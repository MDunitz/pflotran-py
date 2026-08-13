"""Smoke test: companion per-starting-C figures render from synthetic pairs."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pflotran_py.comparison.carbon_inventory import (
    default_comparison_starting_carbon_moles,
)
from pflotran_py.comparison.figures import (
    plot_carbon_dioxide,
    plot_measured_against_modelled,
    plot_methane_timeseries,
)


def _synthetic_paired(tmp_path):
    times = np.linspace(0.0, 120.0, 25)
    # Growing headspace CH4 / dissolved CO2 that Henry will turn into moles.
    ch4 = 1e-6 * (1.0 - np.exp(-times / 40.0))
    co2 = 1e-5 * (1.0 - np.exp(-times / 50.0))
    model = pd.DataFrame(
        {
            "Time [d]": times,
            "Total CH4(aq) [M]": ch4,
            "Active_Gas_CH4(g) [mol_m^3 gas]": ch4 * 400.0,
            "CO2(aq) [M]": co2,
        }
    )
    measured = pd.DataFrame(
        {
            "Experiment": ["Exp004"] * 6,
            "Batch ID": [1] * 6,
            "Replicate ID": [1, 1, 1, 2, 2, 2],
            "Days since start": [10.0, 60.0, 120.0, 10.0, 60.0, 120.0],
            "Cumulative Moles": [1e-6, 5e-5, 2e-4, 1.2e-6, 4.5e-5, 1.8e-4],
            "Molecule": ["CH4_FID"] * 6,
        }
    )
    batch = pd.Series(
        {
            "Experiment": "Exp004",
            "Batch ID": 1,
            "Brine Name": "Na_M",
            "Measured Water Activity": 0.95,
            "Na+": 2.0,
            "Cl-": 2.0,
            "Mg++": 0.0,
            "SO4--": 0.0,
            "Ca++": 0.0,
            "K+": 0.0,
        }
    )
    ecsv = tmp_path / "meas.ecsv"
    # Minimal astropy-readable table via pandas CSV the loader won't use for CO2
    # path below; CO2 plot reads through load_measured, so write a tiny ecsv.
    from astropy.table import Table

    co2_rows = measured.copy()
    co2_rows["Molecule"] = "CO2"
    Table.from_pandas(pd.concat([measured, co2_rows], ignore_index=True)).write(
        ecsv, format="ascii.ecsv", overwrite=True
    )
    paired = [
        {
            "batch": batch,
            "model": model,
            "measured": measured,
            "name": "Exp004_B01_Na_M",
        }
    ]
    return paired, str(ecsv)


def test_companion_figures_write(tmp_path):
    paired, ecsv = _synthetic_paired(tmp_path)
    starting_c = default_comparison_starting_carbon_moles()
    paths = [
        plot_methane_timeseries(
            paired,
            tmp_path / "ch4.png",
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
        plot_carbon_dioxide(
            paired,
            ecsv,
            tmp_path / "co2.png",
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
        plot_measured_against_modelled(
            paired,
            ecsv,
            tmp_path / "parity.png",
            per_starting_c=True,
            starting_carbon_moles=starting_c,
        ),
    ]
    for path in paths:
        assert path.exists()
        assert path.stat().st_size > 1_000
