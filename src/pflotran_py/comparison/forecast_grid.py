"""Parameter grid execution for the forecast protocol."""

import contextlib
import io
import os
from itertools import product

import numpy as np

from ..generator.constants import (
    AW_CRIT_ACETOCLASTIC,
    AW_CRIT_HYDROGENOTROPHIC,
    AW_CRIT_METHYLOTROPHIC,
    AW_INHIBITION_TYPE,
    BOTTLE_FINAL_TIME_DAYS,
)
from .decks import (
    DEFAULT_AW_THRESHOLD_FERMENTATION,
    DEFAULT_AW_THRESHOLD_HYDROLYSIS,
)
from .decks import generate_deck_for_batch
from .figures import load_model_run, model_methane_headspace
from .forecast_anchor import anchor_final_time_days, anchor_table
from .run_decks import run_deck

# Hydrogenotrophic a_crit values to sweep. Methyl and acetoclastic thresholds
# follow the live comparison offsets.
AW_THRESHOLD_GRID = (0.75, 0.80, 0.85, 0.90)

# Upstream ONE_MINUS_AW thresholds (fermentation / hydrolysis rate scaling).
AW_FERMENTATION_GRID = (0.85, 0.90, 0.95)
AW_HYDROLYSIS_GRID = (0.80, 0.85, 0.90)

AW_METHYL_OFFSET = AW_CRIT_METHYLOTROPHIC - AW_CRIT_HYDROGENOTROPHIC
AW_ACETATE_OFFSET = AW_CRIT_ACETOCLASTIC - AW_CRIT_HYDROGENOTROPHIC

GridPoint = tuple[float, float, float, float, float]


def pathway_aw_thresholds(hydrogenotrophic, *, clamp_acetate=False):
    """Map one hydrogenotrophic a_crit to the three pathway thresholds.

    By default methyl and acetoclastic thresholds keep the live comparison
    offsets (+0.05 and +0.10). With ``clamp_acetate=True``, acetoclastic
    a_crit is capped at 1.0 so hydrogenotrophic values above 0.90 remain
    valid when the grid needs to extend past the fixed-offset ceiling.
    """
    methyl = hydrogenotrophic + AW_METHYL_OFFSET
    acetate = hydrogenotrophic + AW_ACETATE_OFFSET
    if clamp_acetate:
        acetate = min(acetate, 1.0)
    if not hydrogenotrophic < methyl < acetate <= 1.0:
        raise ValueError(
            f"invalid pathway thresholds from hydrogenotrophic={hydrogenotrophic}: "
            f"need H2 < methyl < acetate <= 1, got {hydrogenotrophic}, {methyl}, "
            f"{acetate}"
        )
    return hydrogenotrophic, methyl, acetate


def forecast_grid_points(
    aw_threshold_grid=None,
    aw_upstream_inhibition=False,
    aw_fermentation_grid=None,
    aw_hydrolysis_grid=None,
    clamp_acetate=False,
):
    """Enumerate parameter combinations for one forecast grid search."""
    if aw_threshold_grid is None:
        aw_threshold_grid = AW_THRESHOLD_GRID

    ferment_grid = aw_fermentation_grid or (DEFAULT_AW_THRESHOLD_FERMENTATION,)
    hydro_grid = aw_hydrolysis_grid or (DEFAULT_AW_THRESHOLD_HYDROLYSIS,)

    points = []
    for aw_h2 in aw_threshold_grid:
        aw_h2, aw_methyl, aw_acetate = pathway_aw_thresholds(
            aw_h2, clamp_acetate=clamp_acetate
        )
        if aw_upstream_inhibition:
            for aw_ferment, aw_hydro in product(ferment_grid, hydro_grid):
                points.append((aw_h2, aw_methyl, aw_acetate, aw_ferment, aw_hydro))
        else:
            points.append(
                (
                    aw_h2,
                    aw_methyl,
                    aw_acetate,
                    DEFAULT_AW_THRESHOLD_FERMENTATION,
                    DEFAULT_AW_THRESHOLD_HYDROLYSIS,
                )
            )
    return points


def forecast_deck_kwargs(
    aw_threshold,
    aw_threshold_methyl,
    aw_threshold_acetate,
    aw_upstream_inhibition=False,
    aw_threshold_fermentation=DEFAULT_AW_THRESHOLD_FERMENTATION,
    aw_threshold_hydrolysis=DEFAULT_AW_THRESHOLD_HYDROLYSIS,
):
    """Keyword arguments for :func:`generate_deck_for_batch` on the live path."""
    kwargs = {
        "cellulose_hydrolysis": {},
        "enable_cl_inhibition": False,
        "aw_inhibition_type": AW_INHIBITION_TYPE,
        "aw_threshold": aw_threshold,
        "aw_threshold_methyl": aw_threshold_methyl,
        "aw_threshold_acetate": aw_threshold_acetate,
    }
    if aw_upstream_inhibition:
        kwargs.update(
            {
                "aw_upstream_inhibition": True,
                "aw_threshold_fermentation": aw_threshold_fermentation,
                "aw_threshold_hydrolysis": aw_threshold_hydrolysis,
            }
        )
    return kwargs


def aw_grid_tag(
    aw_threshold,
    aw_threshold_methyl,
    aw_threshold_acetate,
    aw_threshold_fermentation=None,
    aw_threshold_hydrolysis=None,
    aw_upstream_inhibition=False,
):
    """Filesystem-safe label for one parameter combination."""
    tag = (
        f"aw{aw_threshold:.2f}_{aw_threshold_methyl:.2f}_{aw_threshold_acetate:.2f}"
    )
    if aw_upstream_inhibition:
        tag += (
            f"_ferm{aw_threshold_fermentation:.2f}"
            f"_hyd{aw_threshold_hydrolysis:.2f}"
        )
    return tag


def measured_by_batch_and_day(measured, experiment):
    """Median measured methane per batch per sampling day."""
    rows = measured[measured["Experiment"] == experiment].copy()
    rows["day"] = rows["Days since start"].round(0)
    return (
        rows.groupby(["Batch ID", "day"])["Cumulative Moles"].median().reset_index()
    )


def model_at_days(days_model, moles_model, wanted_days):
    """Model methane at the measured sampling days."""
    return np.interp(wanted_days, days_model, moles_model)


def run_grid(
    batches,
    work_root,
    repo_root,
    tag_prefix,
    aw_threshold_grid=None,
    anchor_day=0,
    measured_ch4=None,
    measured_co2=None,
    aw_upstream_inhibition=False,
    aw_fermentation_grid=None,
    aw_hydrolysis_grid=None,
    clamp_acetate=False,
):
    """Run every AWINHIBIT parameter combination once, for every batch."""
    grid_points = forecast_grid_points(
        aw_threshold_grid=aw_threshold_grid,
        aw_upstream_inhibition=aw_upstream_inhibition,
        aw_fermentation_grid=aw_fermentation_grid,
        aw_hydrolysis_grid=aw_hydrolysis_grid,
        clamp_acetate=clamp_acetate,
    )

    anchor_lookup = {}
    final_time_days = BOTTLE_FINAL_TIME_DAYS
    if anchor_day > 0:
        if measured_ch4 is None or measured_co2 is None:
            raise ValueError(
                "measured_ch4 and measured_co2 are required when anchor_day > 0"
            )
        anchors = anchor_table(measured_ch4, measured_co2, batches, anchor_day)
        anchor_lookup = {
            int(row["Batch ID"]): row["concentrations"]
            for _, row in anchors.iterrows()
        }
        if not anchor_lookup:
            raise ValueError(
                f"no batches with positive CH4 and CO2 headspace at anchor day "
                f"{anchor_day}"
            )
        final_time_days = anchor_final_time_days(anchor_day)

    grid = {}
    for aw_h2, aw_methyl, aw_acetate, aw_ferment, aw_hydro in grid_points:
        tag = aw_grid_tag(
            aw_h2,
            aw_methyl,
            aw_acetate,
            aw_threshold_fermentation=aw_ferment,
            aw_threshold_hydrolysis=aw_hydro,
            aw_upstream_inhibition=aw_upstream_inhibition,
        )
        series = {}
        for _, batch in batches.iterrows():
            batch_id = int(batch["Batch ID"])
            if anchor_day > 0 and batch_id not in anchor_lookup:
                continue
            name = (
                f"{batch['Experiment']}_B{batch_id:02d}"
                f"_{batch['Brine Name']}"
            )
            deck_dir = os.path.join(work_root, f"{tag_prefix}_decks_{tag}")
            run_root = os.path.join(work_root, f"{tag_prefix}_runs_{tag}")
            deck_kwargs = forecast_deck_kwargs(
                aw_h2,
                aw_methyl,
                aw_acetate,
                aw_upstream_inhibition=aw_upstream_inhibition,
                aw_threshold_fermentation=aw_ferment,
                aw_threshold_hydrolysis=aw_hydro,
            )
            if anchor_day > 0:
                deck_kwargs = {
                    **deck_kwargs,
                    "concentrations": anchor_lookup[batch_id],
                }
            with contextlib.redirect_stdout(io.StringIO()):
                deck = generate_deck_for_batch(
                    batch,
                    output_dir=deck_dir,
                    final_time_days=final_time_days,
                    **deck_kwargs,
                )
            target = os.path.join(deck_dir, f"{name}.in")
            if deck != target:
                os.replace(deck, target)

            result = run_deck(target, run_root, repo_root)
            if result["returncode"] != 0:
                continue
            frame = load_model_run(result["workdir"])
            if frame is None:
                continue
            series[batch_id] = model_methane_headspace(frame, batch)
        grid[(aw_h2, aw_methyl, aw_acetate, aw_ferment, aw_hydro)] = series
    return grid
