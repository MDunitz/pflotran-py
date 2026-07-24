"""Presentation layer: renders analysis output to HTML figures.

Plotting only. Everything upstream of a figure -- extraction, gradients, flux,
spatial transforms -- lives in ``pflotran_py.analysis``; the naming/unit
contract these modules render against is ``analysis.columns``.

    plotly_plotting  -- 3D scatter (Plotly; Bokeh has no 3D)
    bokeh_plotting   -- 2D surface maps + time series (Bokeh)
"""

from . import bokeh_plotting, plotly_plotting  # noqa: F401
