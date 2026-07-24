"""Post-PFLOTRAN analysis: deck output -> DataFrame -> derived quantities.

Import-light by design: nothing here pulls in a plotting backend, so gradient
and flux numbers can be computed (for tables, statistics, or downstream
modelling) without installing Bokeh or Plotly.

    extract, extract_hdf5   -- Tecplot / HDF5 -> DataFrame
    gradients               -- concentration gradients + Fick diffusive flux
    transforms              -- spatial DataFrame ops (surface cells, point series)
    columns                 -- column names, unit strings, time labeling
    constants               -- referenced physical constants (astropy units)
    data_io                 -- pickle / CSV persistence
"""

from . import (  # noqa: F401
    columns,
    constants,
    data_io,
    extract,
    extract_hdf5,
    gradients,
    transforms,
)
