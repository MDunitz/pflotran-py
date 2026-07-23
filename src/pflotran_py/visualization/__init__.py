"""PFLOTRAN reactive-transport post-processing.

Layered so compute and presentation stay separate:
    extract, extract_hdf5   -- Tecplot / HDF5 -> DataFrame
    physics                 -- gradients + Fick diffusive flux (astropy consts)
    transforms              -- spatial DataFrame ops (surface cells, point series)
    columns                 -- column names, unit strings, time labeling
    config, data_io         -- defaults and pickle/CSV persistence
    plotly_plotting         -- 3D scatter (Plotly)
    bokeh_plotting          -- 2D surface maps + time series (Bokeh)
    pipeline                -- extract -> compute -> render orchestration
"""

from . import (  # noqa: F401
    bokeh_plotting,
    columns,
    config,
    data_io,
    extract,
    extract_hdf5,
    physics,
    pipeline,
    plotly_plotting,
    transforms,
)
