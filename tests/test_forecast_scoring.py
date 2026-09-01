"""Unit tests for shared and forecast scoring helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pflotran_py.comparison.forecast_scoring import (
    score_flux_window,
    score_window,
)
from pflotran_py.comparison.scoring import interval_production_rate, score


def test_score_perfect_match():
    assert score([1.0, 10.0, 100.0], [1.0, 10.0, 100.0]) == 0.0


def test_score_typical_miss_by_factor_ten():
    assert score([10.0], [1.0]) == pytest.approx(1.0)


def test_score_returns_inf_when_no_usable_points():
    assert score([0.0], [1.0]) == np.inf
    assert score([1.0], [0.0]) == np.inf


def test_interval_production_rate():
    assert interval_production_rate(0.0, 10.0, 0, 10) == pytest.approx(1.0)


def test_interval_production_rate_rejects_non_positive_dt():
    with pytest.raises(ValueError, match="non-positive interval length"):
        interval_production_rate(0.0, 1.0, 5, 5)


def test_score_window_uses_only_requested_days():
    observed = pd.DataFrame(
        {
            "Batch ID": [1, 1, 1],
            "day": [10, 56, 91],
            "Cumulative Moles": [1e-6, 1e-5, 1e-4],
        }
    )
    grid_entry = {
        1: (
            np.array([0.0, 10.0, 56.0, 91.0]),
            np.array([1e-9, 1e-6, 1e-5, 1e-4]),
        )
    }
    fit_score, count = score_window(grid_entry, observed, [56, 91])
    assert count == 2
    assert fit_score == pytest.approx(0.0)


def test_score_flux_window_uses_interval_rates():
    observed = pd.DataFrame(
        {
            "Batch ID": [1, 1],
            "day": [56.0, 91.0],
            "Cumulative Moles": [0.0, 35.0],
        }
    )
    grid_entry = {
        1: (
            np.array([0.0, 56.0, 91.0]),
            np.array([0.0, 0.0, 35.0]),
        )
    }
    flux_score, count = score_flux_window(grid_entry, observed, [(56, 91)])
    assert count == 1
    assert flux_score == pytest.approx(0.0)
