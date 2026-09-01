"""Unit tests for forecast day-window splitting."""

from __future__ import annotations

import pytest

from pflotran_py.comparison.decks import one_minus_aw_factor
from pflotran_py.comparison.forecast_grid import (
    AW_THRESHOLD_GRID,
    forecast_deck_kwargs,
    forecast_grid_points,
    pathway_aw_thresholds,
)
from pflotran_py.comparison.forecast_split import (
    flux_interval_pairs,
    forecast_output_dir,
    partition_flux_pairs,
    split_forecast_days,
    window_label,
)
from pflotran_py.generator.constants import AW_INHIBITION_TYPE


@pytest.mark.parametrize(
    ("days", "fit_rounds", "holdout_early_days", "expected"),
    [
        (
            [0, 10, 56, 91, 119],
            2,
            0,
            ([10, 56], [91, 119], [], [10, 56, 91, 119]),
        ),
        (
            [0, 10, 56, 91, 119],
            2,
            20,
            ([56, 91], [119], [10], [10, 56, 91, 119]),
        ),
        (
            [0, 10, 56, 91, 119],
            1,
            20,
            ([56], [91, 119], [10], [10, 56, 91, 119]),
        ),
    ],
)
def test_split_forecast_days(days, fit_rounds, holdout_early_days, expected):
    assert (
        split_forecast_days(days, fit_rounds, holdout_early_days=holdout_early_days)
        == expected
    )


def test_window_label():
    fit_days = [56, 91]
    predict_days = [119]
    early_holdout_days = [10]

    assert window_label(0, fit_days, predict_days, early_holdout_days) == "baseline"
    assert window_label(10, fit_days, predict_days, early_holdout_days) == "early holdout"
    assert window_label(56, fit_days, predict_days, early_holdout_days) == "fitted"
    assert window_label(119, fit_days, predict_days, early_holdout_days) == "held out"


def test_forecast_output_dir(tmp_path):
    assert forecast_output_dir(tmp_path, "Exp003", 2, 0).endswith(
        "Exp003_fit-first-2-rounds"
    )
    assert forecast_output_dir(tmp_path, "Exp003", 2, 20).endswith(
        "Exp003_holdout-20d_fit-first-2-rounds"
    )
    assert forecast_output_dir(tmp_path, "Exp003", 2, 20, score_flux=True).endswith(
        "Exp003_holdout-20d_flux_fit-first-2-rounds"
    )


def test_flux_interval_pairs_skip_early_days():
    production = [10, 56, 91, 119, 257]
    assert flux_interval_pairs(production, flux_skip_days=20) == [
        (56, 91),
        (91, 119),
        (119, 257),
    ]


def test_partition_flux_pairs_for_forecast_windows():
    pairs = flux_interval_pairs([10, 56, 91, 119, 257], flux_skip_days=20)
    fit_days, predict_days, _, _ = split_forecast_days(
        [0, 10, 56, 91, 119, 257], fit_rounds=2, holdout_early_days=20
    )
    fit_pairs, predict_pairs = partition_flux_pairs(pairs, fit_days, predict_days)
    assert fit_pairs == [(56, 91)]
    assert predict_pairs == [(91, 119), (119, 257)]


def test_pathway_aw_thresholds_match_live_defaults():
    h2, methyl, acetate = pathway_aw_thresholds(0.80)
    assert (h2, methyl, acetate) == (0.80, 0.85, 0.90)


def test_pathway_aw_thresholds_clamp_acetate():
    h2, methyl, acetate = pathway_aw_thresholds(0.92, clamp_acetate=True)
    assert (h2, methyl, acetate) == (0.92, 0.97, 1.0)


@pytest.mark.parametrize("hydrogenotrophic", AW_THRESHOLD_GRID)
def test_pathway_aw_thresholds_keep_ordering(hydrogenotrophic):
    h2, methyl, acetate = pathway_aw_thresholds(hydrogenotrophic)
    assert h2 < methyl < acetate <= 1.0


def test_forecast_deck_kwargs_match_live_comparison():
    kwargs = forecast_deck_kwargs(0.80, 0.85, 0.90)
    assert kwargs["enable_cl_inhibition"] is False
    assert kwargs["aw_inhibition_type"] == AW_INHIBITION_TYPE
    assert kwargs["cellulose_hydrolysis"] == {}
    assert kwargs["aw_threshold"] == 0.80
    assert kwargs["aw_threshold_methyl"] == 0.85
    assert kwargs["aw_threshold_acetate"] == 0.90


def test_forecast_deck_kwargs_upstream():
    kwargs = forecast_deck_kwargs(
        0.80,
        0.85,
        0.90,
        aw_upstream_inhibition=True,
        aw_threshold_fermentation=0.90,
        aw_threshold_hydrolysis=0.85,
    )
    assert kwargs["aw_upstream_inhibition"] is True
    assert kwargs["aw_threshold_fermentation"] == 0.90
    assert kwargs["aw_threshold_hydrolysis"] == 0.85


def test_forecast_grid_points_upstream_cartesian():
    points = forecast_grid_points(
        aw_threshold_grid=(0.80, 0.85),
        aw_upstream_inhibition=True,
        aw_fermentation_grid=(0.90, 0.95),
        aw_hydrolysis_grid=(0.85,),
    )
    assert len(points) == 4
    assert points[0] == (0.80, 0.85, 0.90, 0.90, 0.85)


def test_one_minus_aw_factor_matches_sandbox():
    assert one_minus_aw_factor(0.824, 0.90) == 0.0
    assert one_minus_aw_factor(0.95, 0.90) == pytest.approx(0.5)
