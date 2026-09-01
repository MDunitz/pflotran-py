"""Warm-start the forecast from measured headspace gas at an anchor day."""

import numpy as np
import pandas as pd

from ..generator.constants import BOTTLE_FINAL_TIME_DAYS
from .headspace import (
    GASES,
    headspace_moles_to_aqueous_concentration,
    setschenow_salts_from_composition,
)

# Simulation time zero is mapped to this incubation age. Gas concentrations at
# the anchor are interpolated from the measured timecourse when no sample falls
# exactly on this day.
DEFAULT_ANCHOR_DAY = 20


def anchor_final_time_days(anchor_day, total_days=BOTTLE_FINAL_TIME_DAYS):
    """Remaining simulated duration after the anchor."""
    if anchor_day < 0:
        raise ValueError(f"anchor_day must be non-negative, got {anchor_day}")
    if anchor_day >= total_days:
        raise ValueError(
            f"anchor_day {anchor_day} must be before final time {total_days}"
        )
    return total_days - anchor_day


def interpolate_cumulative_moles(batch_measured, target_day):
    """Median cumulative headspace moles at one day, linearly interpolated."""
    by_day = (
        batch_measured.groupby("day")["Cumulative Moles"]
        .median()
        .sort_index()
    )
    if by_day.empty:
        return None
    days = by_day.index.to_numpy(dtype=float)
    values = by_day.to_numpy(dtype=float)
    target_day = float(target_day)
    if target_day in days:
        return float(by_day.loc[target_day])
    if target_day < days.min() or target_day > days.max():
        return None
    return float(np.interp(target_day, days, values))


def measured_headspace_at_anchor(measured, batch_id, anchor_day, experiment=None):
    """Interpolate CH4 (or CO2) cumulative headspace moles at the anchor day."""
    batch_rows = measured[measured["Batch ID"] == batch_id].copy()
    if experiment is not None and "Experiment" in batch_rows.columns:
        batch_rows = batch_rows[batch_rows["Experiment"] == experiment]
    if batch_rows.empty:
        return None
    batch_rows["day"] = batch_rows["Days since start"].round(0)
    return interpolate_cumulative_moles(batch_rows, anchor_day)


def build_anchor_concentrations(batch_row, ch4_headspace_mol, co2_headspace_mol):
    """PFLOTRAN initial constraints from measured headspace CH4 and CO2.

    Carbon dioxide is written through ``HCO3- ... G CO2(g)`` so coupled
    carbonate decks stay self-consistent. Bicarbonate not in the headspace is
    still omitted; see :mod:`headspace` for the caveat on CO2 inversion.
    """
    nacl, mgcl2 = setschenow_salts_from_composition(batch_row)
    ch4_aq = headspace_moles_to_aqueous_concentration(
        ch4_headspace_mol,
        GASES["CH4"],
        nacl_molarity=nacl,
        mgcl2_molarity=mgcl2,
    )
    co2_aq = headspace_moles_to_aqueous_concentration(
        co2_headspace_mol,
        GASES["CO2"],
        nacl_molarity=nacl,
        mgcl2_molarity=mgcl2,
    )
    return {
        "CH4(aq)": f"{ch4_aq.to_value('mol/L'):.6e} T CH4(g)",
        "HCO3-": f"{co2_aq.to_value('mol/L'):.6e} T CO2(g)",
    }


def anchor_table(measured_ch4, measured_co2, batches, anchor_day):
    """Per-batch headspace moles at the anchor and deck initial concentrations."""
    records = []
    for _, batch in batches.iterrows():
        batch_id = int(batch["Batch ID"])
        experiment = batch["Experiment"]
        ch4_mol = measured_headspace_at_anchor(
            measured_ch4, batch_id, anchor_day, experiment=experiment
        )
        co2_mol = measured_headspace_at_anchor(
            measured_co2, batch_id, anchor_day, experiment=experiment
        )
        if ch4_mol is None or co2_mol is None or ch4_mol <= 0 or co2_mol <= 0:
            continue
        records.append(
            {
                "Batch ID": batch_id,
                "anchor day": anchor_day,
                "CH4 headspace (mol)": ch4_mol,
                "CO2 headspace (mol)": co2_mol,
                "concentrations": build_anchor_concentrations(batch, ch4_mol, co2_mol),
            }
        )
    return pd.DataFrame(records)


def filter_observed_after_anchor(observed, anchor_day):
    """Keep post-anchor sampling days and add simulation-time coordinates."""
    rows = observed[observed["day"] > anchor_day].copy()
    rows["sim_day"] = rows["day"] - anchor_day
    return rows


def model_days_to_physical(days_model, anchor_day):
    """Shift simulation time back to incubation age for plotting."""
    return np.asarray(days_model, dtype=float) + float(anchor_day)
