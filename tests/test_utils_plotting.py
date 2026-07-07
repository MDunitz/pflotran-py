"""Unit tests for Plotly visualization helpers."""

import os
import sys

import pandas as pd
import plotly.graph_objects as go
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "visualization"))

import utils_plotting as plotting  # noqa: E402

SAMPLE_CSV = os.path.join(REPO_ROOT, "sample_data", "pflotran_data.csv")


@pytest.fixture(scope="module")
def sample_df():
    return pd.read_csv(SAMPLE_CSV)


def test_color_limits_handles_near_zero():
    series = pd.Series([0.0, 0.0, 1e-20])
    vmin, vmax = plotting.color_limits(series)
    assert vmin == 0.0
    assert vmax == plotting.NEAR_ZERO_FLOOR


def test_color_limits_preserves_range():
    series = pd.Series([1.0, 2.0, 3.0])
    assert plotting.color_limits(series) == (1.0, 3.0)


def test_filter_available_variables(sample_df):
    requested = ["CO2(aq) [M]", "Missing Species", "Free CH4(aq) [M]"]
    found = plotting.filter_available_variables(sample_df, requested)
    assert found == ["CO2(aq) [M]", "Free CH4(aq) [M]"]


def test_make_3d_scatter_trace_returns_scatter3d(sample_df):
    snapshot = sample_df[sample_df["Time Index"] == 0]
    trace = plotting.make_3d_scatter_trace(snapshot, "CO2(aq) [M]")

    assert isinstance(trace, go.Scatter3d)
    assert trace.mode == "markers"
    assert len(trace.x) == len(snapshot)


def test_create_single_variable_plot_has_animation_frames(sample_df):
    fig = plotting.create_single_variable_plot(sample_df, "CO2(aq) [M]")

    assert isinstance(fig, go.Figure)
    assert len(fig.frames) == sample_df["Time Index"].nunique()
    assert "3D Visualization of CO2(aq) [M] Over Time" in fig.layout.title.text


def test_create_multi_variable_plot_builds_subplots(sample_df):
    variables = ["CO2(aq) [M]", "Free CH4(aq) [M]"]
    fig = plotting.create_multi_variable_plot(sample_df, variables)

    assert isinstance(fig, go.Figure)
    assert len(fig.data) == len(variables)
    assert len(fig.frames) == sample_df["Time Index"].nunique()


def test_create_multi_variable_plot_returns_none_for_missing_vars(sample_df):
    fig = plotting.create_multi_variable_plot(
        sample_df, ["Not A Real Column"], verbose=False
    )
    assert fig is None


def test_step2_plot_reexports_helpers():
    import step2_plot

    assert (
        step2_plot.create_single_variable_plot is plotting.create_single_variable_plot
    )
    assert step2_plot.create_multi_variable_plot is plotting.create_multi_variable_plot
