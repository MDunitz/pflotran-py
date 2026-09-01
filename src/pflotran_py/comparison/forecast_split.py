"""Temporal splits for the mid-incubation forecast protocol."""

import os
from datetime import datetime

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


def forecast_run_stamp(when=None):
    """Filesystem-safe UTC-local timestamp for one forecast invocation."""
    moment = when or datetime.now()
    return moment.strftime("%Y%m%d_%H%M%S")


def forecast_figure_basename(run_stamp):
    """Timestamped methane figure filename (without directory)."""
    return f"methane_forecast_{run_stamp}.png"


def split_forecast_days(
    days,
    fit_rounds,
    holdout_early_days=0,
    exclude_day_zero=EXCLUDE_DAY_ZERO,
    anchor_day=0,
):
    """Partition sampling days into fit, predict, and early-holdout windows."""
    production_days = sorted(set(days))
    if exclude_day_zero:
        production_days = [day for day in production_days if day > 0]
    if anchor_day > 0:
        production_days = [day for day in production_days if day > anchor_day]

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
    output_root,
    experiment,
    fit_rounds,
    holdout_early_days,
    score_flux=False,
    anchor_day=0,
    run_stamp=None,
    aw_upstream_inhibition=False,
    multi_parameter_fit=False,
):
    """Stable folder name for one forecast configuration."""
    flux_tag = "_flux" if score_flux else ""
    anchor_tag = f"_anchor-{int(anchor_day)}d" if anchor_day > 0 else ""
    upstream_tag = "_upstream" if aw_upstream_inhibition else ""
    multidof_tag = "_multidof" if multi_parameter_fit else ""
    stamp_tag = f"_{run_stamp}" if run_stamp else ""
    suffix = (
        f"fit-first-{fit_rounds}-rounds{upstream_tag}{multidof_tag}{stamp_tag}"
    )
    if holdout_early_days > 0:
        return os.path.join(
            output_root,
            f"{experiment}_holdout-{int(holdout_early_days)}d"
            f"{flux_tag}{anchor_tag}_{suffix}",
        )
    if score_flux or anchor_day > 0 or aw_upstream_inhibition:
        tags = []
        if score_flux:
            tags.append("flux")
        if anchor_day > 0:
            tags.append(f"anchor-{int(anchor_day)}d")
        if aw_upstream_inhibition and not (score_flux or anchor_day > 0):
            tags.append("upstream")
        prefix = f"{experiment}_{'_'.join(tags)}_" if tags else f"{experiment}_"
        return os.path.join(output_root, f"{prefix}{suffix}")
    return os.path.join(output_root, f"{experiment}_{suffix}")


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
