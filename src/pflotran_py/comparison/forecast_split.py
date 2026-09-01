"""Temporal splits for the mid-incubation forecast protocol."""

import os

# Number of sampling rounds used for fitting. The rest are predicted.
DEFAULT_FIT_ROUNDS = 2

# Days at or below this threshold are excluded from fitting. Early methane release
# in these incubations is dominated by fast-substrate turnover rather than salt
# inhibition, so including it when choosing inhibition parameters can mislead the
# fit. Set to zero to disable.
DEFAULT_HOLDOUT_EARLY_DAYS = 0

# When scoring on production rates, skip intervals that start at or below this day.
DEFAULT_FLUX_SKIP_DAYS = 20

# Day zero is excluded from every window; see forecast module docstring.
EXCLUDE_DAY_ZERO = True


def split_forecast_days(
    days,
    fit_rounds,
    holdout_early_days=0,
    exclude_day_zero=EXCLUDE_DAY_ZERO,
):
    """Partition sampling days into fit, predict, and early-holdout windows."""
    production_days = sorted(set(days))
    if exclude_day_zero:
        production_days = [day for day in production_days if day > 0]

    if holdout_early_days > 0:
        early_holdout_days = [
            day for day in production_days if day <= holdout_early_days
        ]
        eligible_days = [
            day for day in production_days if day > holdout_early_days
        ]
    else:
        early_holdout_days = []
        eligible_days = production_days

    fit_days = eligible_days[:fit_rounds]
    predict_days = eligible_days[fit_rounds:]
    return fit_days, predict_days, early_holdout_days, production_days


def window_label(day, fit_days, predict_days, early_holdout_days):
    """Label a sampling day for tables and figures."""
    if day in fit_days:
        return "fitted"
    if day in predict_days:
        return "held out"
    if day in early_holdout_days:
        return "early holdout"
    if day == 0:
        return "baseline"
    return "excluded"


def forecast_output_dir(
    output_root, experiment, fit_rounds, holdout_early_days, score_flux=False
):
    """Stable folder name for one forecast configuration."""
    flux_tag = "_flux" if score_flux else ""
    if holdout_early_days > 0:
        return os.path.join(
            output_root,
            f"{experiment}_holdout-{int(holdout_early_days)}d{flux_tag}_"
            f"fit-first-{fit_rounds}-rounds",
        )
    if score_flux:
        return os.path.join(
            output_root, f"{experiment}_flux_fit-first-{fit_rounds}-rounds"
        )
    return os.path.join(output_root, f"{experiment}_fit-first-{fit_rounds}-rounds")


def flux_interval_pairs(production_days, flux_skip_days):
    """Average production-rate intervals between consecutive sampling days."""
    days = sorted(set(production_days))
    pairs = []
    for t0, t1 in zip(days, days[1:]):
        if t0 <= flux_skip_days:
            continue
        pairs.append((t0, t1))
    return pairs


def partition_flux_pairs(pairs, fit_days, predict_days):
    """Split rate intervals into those used for fitting and for prediction."""
    fit_end = max(fit_days)
    predict_start = min(predict_days)
    fit_pairs = [(t0, t1) for t0, t1 in pairs if t1 <= fit_end]
    predict_pairs = [
        (t0, t1)
        for t0, t1 in pairs
        if t1 >= predict_start and t0 >= min(fit_days)
    ]
    return fit_pairs, predict_pairs
