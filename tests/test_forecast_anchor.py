"""Unit tests for forecast warm-start anchoring."""

from __future__ import annotations

import pandas as pd
import pytest

from pflotran_py.comparison.forecast_anchor import (
    DEFAULT_ANCHOR_DAY,
    anchor_final_time_days,
    build_anchor_concentrations,
    interpolate_cumulative_moles,
    measured_headspace_at_anchor,
    model_days_to_physical,
)
from pflotran_py.comparison.forecast_split import forecast_output_dir, split_forecast_days
from pflotran_py.comparison.headspace import (
    GASES,
    aqueous_concentration_to_headspace_moles,
    headspace_moles_to_aqueous_concentration,
)
from pflotran_py.generator.constants import BOTTLE_FINAL_TIME_DAYS


def test_anchor_final_time_days():
    assert anchor_final_time_days(20) == BOTTLE_FINAL_TIME_DAYS - 20


def test_interpolate_cumulative_moles_between_samples():
    batch = pd.DataFrame(
        {
            "day": [10, 10, 56, 56],
            "Cumulative Moles": [1.0, 3.0, 5.0, 9.0],
        }
    )
    assert interpolate_cumulative_moles(batch, 10) == pytest.approx(2.0)
    assert interpolate_cumulative_moles(batch, 20) == pytest.approx(3.0869565)


def test_measured_headspace_at_anchor_filters_experiment():
    measured = pd.DataFrame(
        {
            "Experiment": ["Exp003", "Exp004"],
            "Batch ID": [1, 1],
            "Days since start": [10.0, 10.0],
            "Cumulative Moles": [1.0, 100.0],
        }
    )
    assert measured_headspace_at_anchor(
        measured, 1, 10, experiment="Exp003"
    ) == pytest.approx(1.0)


def test_build_anchor_concentrations_use_gas_equilibrium():
    batch = pd.Series(
        {
            "Experiment": "Exp003",
            "Batch ID": 1,
            "Na+": 2.0,
            "Mg++": 0.0,
        }
    )
    concentrations = build_anchor_concentrations(batch, 1e-6, 2e-6)
    assert concentrations["CH4(aq)"].endswith("T CH4(g)")
    assert concentrations["HCO3-"].endswith("T CO2(g)")


def test_headspace_moles_to_aqueous_round_trip():
    concentration = 1e-4
    in_gas = aqueous_concentration_to_headspace_moles(concentration, GASES["CH4"])
    recovered = headspace_moles_to_aqueous_concentration(in_gas, GASES["CH4"])
    assert recovered.to_value("mol/L") == pytest.approx(concentration, rel=1e-9)


def test_split_forecast_days_with_anchor():
    assert split_forecast_days(
        [0, 10, 56, 91, 119], fit_rounds=2, anchor_day=20
    ) == ([56, 91], [119], [], [56, 91, 119])


def test_forecast_output_dir_anchor_tag(tmp_path):
    path = forecast_output_dir(tmp_path, "Exp003", 2, 0, anchor_day=20)
    assert path.endswith("Exp003_anchor-20d_fit-first-2-rounds")


def test_model_days_to_physical():
    assert model_days_to_physical([0, 36], 20).tolist() == [20.0, 56.0]


def test_default_anchor_day():
    assert DEFAULT_ANCHOR_DAY == 20
