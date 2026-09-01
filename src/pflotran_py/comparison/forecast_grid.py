"""Parameter grid execution for the forecast protocol."""

import contextlib
import io
import os

import numpy as np

from ..generator.constants import (
    AW_CRIT_ACETOCLASTIC,
    AW_CRIT_HYDROGENOTROPHIC,
    AW_CRIT_METHYLOTROPHIC,
    AW_INHIBITION_TYPE,
    BOTTLE_FINAL_TIME_DAYS,
)
from .decks import generate_deck_for_batch
from .figures import load_model_run, model_methane_headspace
from .forecast_anchor import anchor_final_time_days, anchor_table
from .run_decks import run_deck

# Hydrogenotrophic a_crit values to sweep. Methyl and acetoclastic thresholds
# follow the live comparison offsets.
AW_THRESHOLD_GRID = (0.75, 0.80, 0.85, 0.90)

AW_METHYL_OFFSET = AW_CRIT_METHYLOTROPHIC - AW_CRIT_HYDROGENOTROPHIC
AW_ACETATE_OFFSET = AW_CRIT_ACETOCLASTIC - AW_CRIT_HYDROGENOTROPHIC


def pathway_aw_thresholds(hydrogenotrophic):
    """Map one hydrogenotrophic a_crit to the three pathway thresholds."""
    methyl = hydrogenotrophic + AW_METHYL_OFFSET
    acetate = hydrogenotrophic + AW_ACETATE_OFFSET
    if not hydrogenotrophic < methyl < acetate <= 1.0:
        raise ValueError(
            f"invalid pathway thresholds from hydrogenotrophic={hydrogenotrophic}: "
            f"need H2 < methyl < acetate <= 1, got {hydrogenotrophic}, {methyl}, "
            f"{acetate}"
        )
    return hydrogenotrophic, methyl, acetate


def forecast_deck_kwargs(aw_threshold, aw_threshold_methyl, aw_threshold_acetate):
    """Keyword arguments for :func:`generate_deck_for_batch` on the live path."""
    return {
        "cellulose_hydrolysis": {},
        "enable_cl_inhibition": False,
        "aw_inhibition_type": AW_INHIBITION_TYPE,
        "aw_threshold": aw_threshold,
        "aw_threshold_methyl": aw_threshold_methyl,
        "aw_threshold_acetate": aw_threshold_acetate,
    }


def aw_grid_tag(aw_threshold, aw_threshold_methyl, aw_threshold_acetate):
    """Filesystem-safe label for one parameter triple."""
    return (
        f"aw{aw_threshold:.2f}_{aw_threshold_methyl:.2f}_{aw_threshold_acetate:.2f}"
    )


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
):
    """Run every AWINHIBIT parameter combination once, for every batch."""
    if aw_threshold_grid is None:
        aw_threshold_grid = AW_THRESHOLD_GRID

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
    for aw_h2 in aw_threshold_grid:
        aw_h2, aw_methyl, aw_acetate = pathway_aw_thresholds(aw_h2)
        tag = aw_grid_tag(aw_h2, aw_methyl, aw_acetate)
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
            deck_kwargs = forecast_deck_kwargs(aw_h2, aw_methyl, aw_acetate)
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
        grid[(aw_h2, aw_methyl, aw_acetate)] = series
    return grid
