"""
PFLOTRAN Visualization Pipeline Orchestrator

Runs steps 1–4 in sequence:
    1. Extract data from PFLOTRAN output (.tec or .hdf5)
    2. 3D Plotly scatter plots with time animation
    3. 2D surface gradient/flux maps (Bokeh)
    4. Time series gradient/flux at a single point (Bokeh)

For multi-condition comparisons (varying water activity),
see ../compare/comparing_aw.ipynb
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import step1_extract
import step2_plot
import step3_flux
import step4_plotflux

# ─────────────────────────────────────────────────────────────────
# Configuration — edit these as needed
# ─────────────────────────────────────────────────────────────────

# Data source
DATA_FORMAT = 'tec'         # 'tec' or 'hdf5'
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sample_data')
FILE_TEMPLATE = 'test29-{:03d}.tec'
N_FILES = 6

# Species to compute gradients/flux for
SPECIES_MAP = {
    'CO2': 'CO2(aq) [M]',
    'CH4': 'Free CH4(aq) [M]',
}

# Simulation temperature for Stokes-Einstein diffusion correction
# Must match the temperature used in the PFLOTRAN .in file
TEMPERATURE_C = 8.0

# Number of timesteps to render in 2D surface maps (step 3)
N_TIMESTEPS_2D = 5

# Output directories
OUTPUT_DIR_STEP3 = 'experiments/figures/pflotran_flux_imgs'
OUTPUT_DIR_STEP4 = 'experiments/figures/pflotran_flux_imgs'
OUTPUT_FILE_STEP3 = 'co2_ch4_surface_visualization.html'
OUTPUT_FILE_STEP4 = 'gradient_time_series.html'


def main():
    print("Starting PFLOTRAN visualization pipeline...\n")

    # Step 1: Extract data
    step1_extract.main(DATA_FORMAT, DATA_DIR, FILE_TEMPLATE, N_FILES)

    # Step 2: 3D Plotly scatter (keeps Plotly for 3D — good for figures)
    step2_plot.main()

    # Step 3: 2D surface gradient/flux maps (Bokeh)
    step3_flux.main(
        species_map=SPECIES_MAP,
        n_timesteps=N_TIMESTEPS_2D,
        output_dir=OUTPUT_DIR_STEP3,
        output_filename=OUTPUT_FILE_STEP3,
        compute_flux=True,
        temperature_c=TEMPERATURE_C,
    )

    # Step 4: Time series at a single point (Bokeh)
    step4_plotflux.main(
        species_map=SPECIES_MAP,
        output_dir=OUTPUT_DIR_STEP4,
        output_filename=OUTPUT_FILE_STEP4,
        compute_flux=True,
        temperature_c=TEMPERATURE_C,
    )

    print("\nPipeline complete — open the HTML files to view results.")


if __name__ == '__main__':
    main()

