"""PFLOTRAN visualization pipeline: extract -> compute -> render.

Replaces the former step1..4 + orchestrator scripts. Each render stage is
independently runnable from the saved pickle; ``run`` chains them. 3D scatter
uses Plotly (Bokeh has no 3D); the 2D surface maps and time series use Bokeh.

Run the whole pipeline from a checkout with:
    python -m pflotran_py.visualization.pipeline
"""

import logging
import os

from . import bokeh_plotting, physics, plotly_plotting, transforms
from .config import (
    DEFAULT_3D_VARIABLES,
    DEFAULT_DATA_PICKLE,
    DEFAULT_N_TIMESTEPS_2D,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_SPECIES_MAP,
    DEFAULT_SURFACE_FILENAME,
    DEFAULT_TEMPERATURE_C,
    DEFAULT_TIMESERIES_FILENAME,
)
from .data_io import load_data, save_data
from .extract import (
    extract_pflotran_data_tec,
    find_hdf5_output,
    extract_pflotran_data_hdf5,
)

logger = logging.getLogger(__name__)


# ═════════════════════════════════════════════════════════════════════
# Stage 1 — extraction
# ═════════════════════════════════════════════════════════════════════


def extract_and_save(
    data_format, data_dir, file_name_template, n_files, pickle=DEFAULT_DATA_PICKLE
):
    """Extract PFLOTRAN output to a DataFrame and persist it as pickle + CSV."""
    if data_format == "tec":
        logger.info("Extracting data from tecplot files...")
        df = extract_pflotran_data_tec(
            data_dir=data_dir, file_name_template=file_name_template, n_files=n_files
        )
    elif data_format == "hdf5":
        logger.info("Extracting data from hdf5 file...")
        h5_path = find_hdf5_output(data_dir)
        if h5_path is None:
            raise FileNotFoundError(f"No PFLOTRAN .h5 output found in {data_dir}")
        df = extract_pflotran_data_hdf5(h5_path)
    else:
        raise ValueError(f"Unknown data_format: {data_format!r}")

    logger.info("Saving data to %s", pickle)
    save_data(df, pickle)
    return df


# ═════════════════════════════════════════════════════════════════════
# Stage 2 — Plotly 3D scatter
# ═════════════════════════════════════════════════════════════════════


def render_3d(variables_to_plot=None, pickle=DEFAULT_DATA_PICKLE, output_dir="."):
    """Render Plotly 3D multi- and single-variable scatter animations."""
    variables_to_plot = variables_to_plot or DEFAULT_3D_VARIABLES
    df = load_data(pickle)
    os.makedirs(output_dir, exist_ok=True)

    fig = plotly_plotting.create_multi_variable_plot(df, variables_to_plot)
    if fig is None:
        logger.warning("No requested variables found in data; skipping 3D export.")
        return
    multi_path = os.path.join(output_dir, "multi_variable_concentration_3d.html")
    fig.write_html(multi_path)
    logger.info("Multi-variable 3D plot saved: %s", multi_path)

    if variables_to_plot[0] in df.columns:
        single_var = variables_to_plot[0]
        single_fig = plotly_plotting.create_single_variable_plot(df, single_var)
        safe = single_var.replace(" ", "_").replace("[", "").replace("]", "")
        single_path = os.path.join(output_dir, f"single_variable_{safe}.html")
        single_fig.write_html(single_path)
        logger.info("Single-variable 3D plot saved: %s", single_path)


# ═════════════════════════════════════════════════════════════════════
# Stage 3 — Bokeh 2D surface maps
# ═════════════════════════════════════════════════════════════════════


def render_surface_stage(
    species_map=None,
    n_timesteps=DEFAULT_N_TIMESTEPS_2D,
    output_dir=DEFAULT_OUTPUT_DIR,
    output_filename=DEFAULT_SURFACE_FILENAME,
    compute_flux=True,
    temperature_c=DEFAULT_TEMPERATURE_C,
    pickle=DEFAULT_DATA_PICKLE,
):
    """Compute gradients (and flux) and render 2D surface maps."""
    species_map = species_map or DEFAULT_SPECIES_MAP
    df = load_data(pickle)

    logger.info("Computing concentration gradients...")
    df = physics.calculate_gradients(df, species_map)
    if compute_flux:
        logger.info("Converting to flux (Fick's law, T=%s°C)...", temperature_c)
        df = physics.convert_to_flux(df, list(species_map), temperature_c=temperature_c)

    df = transforms.identify_surface_cells(df)

    bokeh_plotting.render_surface_maps(
        df, species_map, n_timesteps, output_dir, output_filename, use_flux=False
    )
    if compute_flux:
        bokeh_plotting.render_surface_maps(
            df,
            species_map,
            n_timesteps,
            output_dir,
            output_filename.replace(".html", "_flux.html"),
            use_flux=True,
        )


# ═════════════════════════════════════════════════════════════════════
# Stage 4 — Bokeh 2D time series at a single point
# ═════════════════════════════════════════════════════════════════════


def render_time_series_stage(
    species_map=None,
    output_dir=DEFAULT_OUTPUT_DIR,
    output_filename=DEFAULT_TIMESERIES_FILENAME,
    compute_flux=True,
    temperature_c=DEFAULT_TEMPERATURE_C,
    pickle=DEFAULT_DATA_PICKLE,
):
    """Compute gradients (and flux) and render single-point time series."""
    species_map = species_map or DEFAULT_SPECIES_MAP
    df = load_data(pickle)

    target = transforms.find_target_point(df)
    logger.info("Target point: (%.3f, %.3f, %.3f)", *target)

    logger.info("Computing concentration gradients...")
    df = physics.calculate_gradients(df, species_map)
    if compute_flux:
        logger.info("Converting to flux (T=%s°C)...", temperature_c)
        df = physics.convert_to_flux(df, list(species_map), temperature_c=temperature_c)

    point_data = transforms.extract_point_time_series(df, *target)

    bokeh_plotting.render_time_series(
        point_data, target, species_map, output_dir, output_filename, use_flux=False
    )
    if compute_flux:
        bokeh_plotting.render_time_series(
            point_data,
            target,
            species_map,
            output_dir,
            output_filename.replace(".html", "_flux.html"),
            use_flux=True,
        )


# ═════════════════════════════════════════════════════════════════════
# Orchestrator
# ═════════════════════════════════════════════════════════════════════

_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
)


def run(
    data_format="tec",
    data_dir=None,
    file_name_template="test29-{:03d}.tec",
    n_files=6,
    species_map=None,
    temperature_c=DEFAULT_TEMPERATURE_C,
    n_timesteps_2d=DEFAULT_N_TIMESTEPS_2D,
    output_dir=DEFAULT_OUTPUT_DIR,
):
    """Run the full pipeline: extract -> 3D -> surface maps -> time series."""
    if data_dir is None:
        data_dir = os.path.join(_REPO_ROOT, "sample_data")
    species_map = species_map or DEFAULT_SPECIES_MAP

    logger.info("Starting PFLOTRAN visualization pipeline...")
    extract_and_save(data_format, data_dir, file_name_template, n_files)
    render_3d()
    render_surface_stage(
        species_map=species_map,
        n_timesteps=n_timesteps_2d,
        output_dir=output_dir,
        compute_flux=True,
        temperature_c=temperature_c,
    )
    render_time_series_stage(
        species_map=species_map,
        output_dir=output_dir,
        compute_flux=True,
        temperature_c=temperature_c,
    )
    logger.info("Pipeline complete — open the HTML files to view results.")


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run()


if __name__ == "__main__":
    main()
