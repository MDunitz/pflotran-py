"""Tests for the inhibition-attribution diagnostic.

The runs themselves need a container, so what is checked here is everything
around them: that each variant switches off the term it claims to, that the
decks it produces really lack those blocks, and that the attribution arithmetic
recovers a known fold change.
"""

import numpy as np
import pandas as pd
import pytest

from pflotran_py.comparison.inhibition_diagnostic import (
    ATTRIBUTIONS,
    SULFATE_RATE_KEYS,
    VARIANTS,
    VARIANTS_BY_KEY,
    attribute,
)
from pflotran_py.generator.bottle_generator import BottleGenerator


def test_every_variant_has_a_distinct_switch_combination():
    combinations = {
        (v.cl_monod, v.cl_smoothstep, v.sulfate) for v in VARIANTS
    }
    assert len(combinations) == len(VARIANTS)


def test_baseline_keeps_everything_on():
    baseline = VARIANTS_BY_KEY["baseline"]
    kwargs = baseline.generator_kwargs()
    assert kwargs["enable_cl_inhibition"] is True
    assert kwargs["salinity_inhibition"]["species"] == "Cl-"
    assert "disabled_rate_keys" not in kwargs


@pytest.mark.parametrize(
    "key,expect_monod,expect_smoothstep,expect_sulfate",
    [
        ("baseline", True, True, True),
        ("no_cl_monod", False, True, True),
        ("no_smoothstep", True, False, True),
        ("no_salt_terms", False, False, True),
        ("ceiling", False, False, False),
    ],
)
def test_variant_kwargs_match_their_switches(
    key, expect_monod, expect_smoothstep, expect_sulfate
):
    kwargs = VARIANTS_BY_KEY[key].generator_kwargs()
    assert kwargs["enable_cl_inhibition"] is expect_monod
    assert ("salinity_inhibition" in kwargs) is expect_smoothstep
    if expect_sulfate:
        assert "disabled_rate_keys" not in kwargs
    else:
        assert kwargs["disabled_rate_keys"] == set(SULFATE_RATE_KEYS)


def test_every_variant_uses_the_matched_carbon_inventory():
    """A variant that changed the carbon pool would not be comparable."""
    for variant in VARIANTS:
        assert variant.generator_kwargs()["cellulose_hydrolysis"] == {}


def test_attribution_pairs_differ_in_exactly_one_switch():
    """Each fold change must isolate one term, or it is not an attribution."""
    for _, without_key, with_key, _ in ATTRIBUTIONS:
        without = VARIANTS_BY_KEY[without_key]
        with_ = VARIANTS_BY_KEY[with_key]
        differences = sum(
            getattr(without, field) != getattr(with_, field)
            for field in ("cl_monod", "cl_smoothstep", "sulfate")
        )
        assert differences == 1, f"{without_key} vs {with_key} differ in {differences}"


# ─────────────────────────────────────────────────────────────────────
# Decks really carry (or lack) the blocks the variant names
# ─────────────────────────────────────────────────────────────────────


def _deck_text(tmp_path, variant, name):
    generator = BottleGenerator(
        brine={"Na+": "1.0000e+00 T", "Cl-": "1.0000e+00 Z"},
        final_time_days=1,
        **variant.generator_kwargs(),
    )
    path = generator.generate(str(tmp_path / f"{name}.in"))
    with open(path) as handle:
        return handle.read()


def test_cl_monod_present_only_when_enabled(tmp_path):
    with_monod = _deck_text(tmp_path, VARIANTS_BY_KEY["no_smoothstep"], "with")
    without = _deck_text(tmp_path, VARIANTS_BY_KEY["no_salt_terms"], "without")

    # The legacy term is a MONOD inhibition at the 0.2 mol/L threshold.
    assert "2.00e-01" in with_monod
    assert "2.00e-01" not in without


def test_smoothstep_present_only_when_enabled(tmp_path):
    with_step = _deck_text(tmp_path, VARIANTS_BY_KEY["no_cl_monod"], "with")
    without = _deck_text(tmp_path, VARIANTS_BY_KEY["no_salt_terms"], "without")

    # Cl- smoothstep on the network (SMOOTHSTEP_INTERVAL), not the a_w sandbox.
    assert "SMOOTHSTEP_INTERVAL" in with_step
    assert "SMOOTHSTEP_INTERVAL" not in without


def test_ceiling_deck_drops_the_sulfate_pathways(tmp_path):
    ceiling = _deck_text(tmp_path, VARIANTS_BY_KEY["ceiling"], "ceiling")
    kept = _deck_text(tmp_path, VARIANTS_BY_KEY["no_salt_terms"], "kept")

    assert "sulfate reduction" in kept
    assert "methane oxidation (SO4)" in kept
    assert "sulfate reduction" not in ceiling
    assert "methane oxidation (SO4)" not in ceiling
    # Dropping sulfate must not disturb the other methane sinks.
    assert "methane oxidation (NO3)" in ceiling
    assert "methane oxidation (Fe)" in ceiling


def test_unknown_rate_key_is_rejected():
    with pytest.raises(ValueError, match="Unknown rate key"):
        BottleGenerator(disabled_rate_keys={"not_a_reaction"})


# ─────────────────────────────────────────────────────────────────────
# Attribution arithmetic
# ─────────────────────────────────────────────────────────────────────


def _collected(values_by_variant):
    rows = []
    for variant, moles in values_by_variant.items():
        rows.append(
            {
                "variant": variant,
                "variant_label": variant,
                "batch_name": "Exp003_B03_SW_H",
                "experiment": "Exp003",
                "batch_id": 3,
                "brine": "SW_H",
                "water_activity": 0.95,
                "ionic_strength": 1.2,
                "ch4_moles": moles,
                "ch4_per_starting_c": moles / 0.0565,
            }
        )
    return pd.DataFrame(rows)


def test_attribute_recovers_known_fold_changes():
    collected = _collected(
        {
            "baseline": 1.0e-6,
            "no_cl_monod": 3.0e-5,  # 30x
            "no_smoothstep": 5.0e-6,  # 5x
            "no_salt_terms": 1.0e-4,
            "ceiling": 4.0e-4,  # 4x over no_salt_terms
        }
    )
    result = attribute(collected)
    assert len(result) == 1
    row = result.iloc[0]
    assert row["cl_monod_fold"] == pytest.approx(30.0)
    assert row["cl_smoothstep_fold"] == pytest.approx(5.0)
    assert row["sulfate_fold"] == pytest.approx(4.0)


def test_attribute_handles_a_missing_variant():
    """A non-convergent variant should leave a gap, not raise."""
    collected = _collected({"baseline": 1.0e-6, "no_cl_monod": 3.0e-5})
    result = attribute(collected)
    assert result["cl_monod_fold"].iloc[0] == pytest.approx(30.0)
    assert "sulfate_fold" not in result.columns


def test_attribute_of_nothing_is_empty():
    assert attribute(pd.DataFrame()).empty


def test_zero_denominator_becomes_nan_rather_than_infinity():
    collected = _collected({"baseline": 0.0, "no_cl_monod": 3.0e-5})
    result = attribute(collected)
    assert np.isnan(result["cl_monod_fold"].iloc[0])
