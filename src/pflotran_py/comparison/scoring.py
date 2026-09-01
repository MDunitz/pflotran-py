"""Shared objectives for comparing modelled and measured methane."""

import numpy as np


def score(modelled, measured):
    """Root mean squared error of log10(modelled / measured).

    Logarithmic because the measurements span four orders of magnitude, so an
    absolute error would be decided entirely by the largest bottle. Squared
    rather than absolute so that a single badly-missed condition is not hidden
    by a good median, which matters when the whole question is whether the
    model can follow a collapse.

    Zero is exact agreement; a value of one means a typical miss by a factor of
    ten.
    """
    modelled = np.asarray(modelled, dtype=float)
    measured = np.asarray(measured, dtype=float)
    usable = (modelled > 0) & (measured > 0)
    if not usable.any():
        return np.inf
    return float(np.sqrt(np.mean(np.log10(modelled[usable] / measured[usable]) ** 2)))


def interval_production_rate(m0, m1, t0, t1):
    """Average headspace methane production rate over one sampling interval."""
    dt = float(t1) - float(t0)
    if dt <= 0:
        raise ValueError(f"non-positive interval length: {t0} to {t1}")
    return (float(m1) - float(m0)) / dt
