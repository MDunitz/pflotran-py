"""Scoring helpers for the forecast protocol."""

from .forecast_grid import model_at_days
from .scoring import interval_production_rate, score


def _model_day(physical_day, anchor_day=0):
    return float(physical_day) - float(anchor_day)


def observed_by_batch(observed):
    """Median measured methane indexed by batch and sampling day."""
    lookup = {}
    for _, row in observed.iterrows():
        lookup.setdefault(int(row["Batch ID"]), {})[float(row["day"])] = float(
            row["Cumulative Moles"]
        )
    return lookup


def score_window(grid_entry, observed, days, anchor_day=0):
    """Score one parameter combination over cumulative moles at sampling days."""
    modelled, measured_values = [], []
    for _, row in observed[observed["day"].isin(days)].iterrows():
        batch_id = int(row["Batch ID"])
        if batch_id not in grid_entry:
            continue
        model_days, model_moles = grid_entry[batch_id]
        modelled.append(
            model_at_days(
                model_days, model_moles, _model_day(row["day"], anchor_day)
            )
        )
        measured_values.append(row["Cumulative Moles"])
    return score(modelled, measured_values), len(modelled)


def score_flux_window(grid_entry, observed, pairs, anchor_day=0):
    """Score one parameter combination over interval-averaged production rates."""
    lookup = observed_by_batch(observed)
    modelled_rates, measured_rates = [], []
    for batch_id, (model_days, model_moles) in grid_entry.items():
        batch_measured = lookup.get(batch_id)
        if not batch_measured:
            continue
        for t0, t1 in pairs:
            if t0 not in batch_measured or t1 not in batch_measured:
                continue
            measured_rates.append(
                interval_production_rate(
                    batch_measured[t0], batch_measured[t1], t0, t1
                )
            )
            modelled_rates.append(
                interval_production_rate(
                    model_at_days(
                        model_days, model_moles, _model_day(t0, anchor_day)
                    ),
                    model_at_days(
                        model_days, model_moles, _model_day(t1, anchor_day)
                    ),
                    t0,
                    t1,
                )
            )
    return score(modelled_rates, measured_rates), len(modelled_rates)
