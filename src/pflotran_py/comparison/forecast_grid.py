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
)
from .decks import generate_deck_for_batch
from .figures import load_model_run, model_methane_headspace
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


def run_grid(batches, work_root, repo_root, tag_prefix, aw_threshold_grid=None):
    """Run every AWINHIBIT parameter combination once, for every batch."""
    if aw_threshold_grid is None:
        aw_threshold_grid = AW_THRESHOLD_GRID

    grid = {}
    for aw_h2 in aw_threshold_grid:
        aw_h2, aw_methyl, aw_acetate = pathway_aw_thresholds(aw_h2)
        tag = aw_grid_tag(aw_h2, aw_methyl, aw_acetate)
        series = {}
        for _, batch in batches.iterrows():
            name = (
                f"{batch['Experiment']}_B{int(batch['Batch ID']):02d}"
                f"_{batch['Brine Name']}"
            )
            deck_dir = os.path.join(work_root, f"{tag_prefix}_decks_{tag}")
            run_root = os.path.join(work_root, f"{tag_prefix}_runs_{tag}")
            with contextlib.redirect_stdout(io.StringIO()):
                deck = generate_deck_for_batch(
                    batch,
                    output_dir=deck_dir,
                    **forecast_deck_kwargs(aw_h2, aw_methyl, aw_acetate),
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
            series[int(batch["Batch ID"])] = model_methane_headspace(frame, batch)
        grid[(aw_h2, aw_methyl, aw_acetate)] = series
    return grid
