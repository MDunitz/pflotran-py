"""Document the ONE_MINUS_AW rate factor used by the AWINHIBIT sandboxes.

The Fortran implements f = max(0, (a_w - a_crit) / (1 - a_crit)). This test
locks the Python-side expectation so a deck-default change cannot silently
reintroduce a cliff without updating the documented shape.
"""

import pytest

from pflotran_py.generator.bottle_generator import BottleGenerator


def one_minus_aw(water_activity: float, a_crit: float) -> float:
    if a_crit >= 1.0:
        return 1.0
    if water_activity <= a_crit:
        return 0.0
    return (water_activity - a_crit) / (1.0 - a_crit)


@pytest.mark.parametrize(
    "a_w,expected",
    [
        (1.0, 1.0),
        (0.90, 0.5),
        (0.80, 0.0),
        (0.75, 0.0),
        (0.8424, pytest.approx(0.212, abs=1e-3)),  # Mg_H Pitzer
        (0.9127, pytest.approx(0.5635, abs=1e-3)),  # Na_M Pitzer
    ],
)
def test_one_minus_aw_factors_at_a_crit_0_80(a_w, expected):
    assert one_minus_aw(a_w, 0.80) == expected


def test_comparison_default_emits_one_minus_aw():
    gen = BottleGenerator(aw_inhibition_type="ONE_MINUS_AW")
    block = gen._build_reaction_sandbox()
    assert "INHIBITION_TYPE ONE_MINUS_AW" in block
    assert "WATER_ACTIVITY_THRESHOLD 0.9100" in block


def test_pathway_specific_aw_thresholds_are_ordered():
    """Acetoclastic is most salt-sensitive; H2 and methyl share a lower a_crit."""
    gen = BottleGenerator(
        aw_threshold=0.91,
        aw_threshold_methyl=0.91,
        aw_threshold_acetate=0.92,
        aw_inhibition_type="ONE_MINUS_AW",
    )
    block = gen._build_reaction_sandbox()
    import re

    found = re.findall(
        r"(AWINHIBIT(?:ACETATE|METHYL)?)\n\s+WATER_ACTIVITY_THRESHOLD ([0-9.]+)",
        block,
    )
    by_name = {name: float(val) for name, val in found}
    assert by_name["AWINHIBIT"] == 0.91
    assert by_name["AWINHIBITMETHYL"] == 0.91
    assert by_name["AWINHIBITACETATE"] == 0.92
    assert by_name["AWINHIBITACETATE"] > by_name["AWINHIBIT"]
