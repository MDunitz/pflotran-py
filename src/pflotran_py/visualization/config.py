"""Default configuration for the visualization pipeline.

Single source for values previously duplicated across the step files and the
orchestrator (species map, temperature, output paths).
"""

# Short species name -> DataFrame concentration column.
DEFAULT_SPECIES_MAP = {
    "CO2": "CO2(aq) [M]",
    "CH4": "Free CH4(aq) [M]",
}

# Simulation temperature [°C] for the Stokes-Einstein diffusion correction.
# Must match the temperature used in the PFLOTRAN .in file.
DEFAULT_TEMPERATURE_C = 8.0

# Number of timesteps rendered in the 2D surface maps.
DEFAULT_N_TIMESTEPS_2D = 5

# Output locations (relative to the current working directory).
DEFAULT_OUTPUT_DIR = "experiments/figures/pflotran_flux_imgs"
DEFAULT_SURFACE_FILENAME = "co2_ch4_surface_visualization.html"
DEFAULT_TIMESERIES_FILENAME = "gradient_time_series.html"

# Pickle written by the extractor and read by the plotting stages.
DEFAULT_DATA_PICKLE = "pflotran_data.pkl"

# Variables rendered by the Plotly 3D stage (step 2).
DEFAULT_3D_VARIABLES = [
    "Total CH4(aq) [M]",
    "Total Acetate- [M]",
    "Total SO4-- [M]",
    "Gamma H2O",
]
