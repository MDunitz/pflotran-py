"""Document the ONE_MINUS_AW rate factor used by the AWINHIBIT sandboxes.

The Fortran implements f = max(0, (a_w - a_crit) / (1 - a_crit)). This test
locks the Python-side expectation so a deck-default change cannot silently
reintroduce a cliff without updating the documented shape.

The Mg_H / Na_M a_w points below are the computed (PHREEQC/pitzer.dat) values;
the sandbox is fed the meter-read a_w by default (--use-computed-aw to
override), so the factor the deck actually applies uses the meter value.
"""

import re

import pytest

from pflotran_py.generator.bottle_generator import BottleGenerator
from pflotran_py.generator.constants import (
    AW_CRIT_ACETOCLASTIC,
    AW_CRIT_HYDROGENOTROPHIC,
    AW_CRIT_METHYLOTROPHIC,
    AW_INHIBITION_TYPE,
)


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
    """Shape check at a_crit = 0.80, independent of the comparison default."""
    assert one_minus_aw(a_w, 0.80) == expected


def test_comparison_default_emits_one_minus_aw():
    gen = BottleGenerator(aw_inhibition_type=AW_INHIBITION_TYPE)
    block = gen._build_reaction_sandbox()
    assert f"INHIBITION_TYPE {AW_INHIBITION_TYPE}" in block
    assert f"WATER_ACTIVITY_THRESHOLD {AW_CRIT_HYDROGENOTROPHIC:.4f}" in block


def test_pathway_specific_aw_thresholds_are_ordered():
    """Acetoclastic is most salt-sensitive, then methyl, then hydrogenotrophic."""
    gen = BottleGenerator(aw_inhibition_type=AW_INHIBITION_TYPE)
    block = gen._build_reaction_sandbox()

    found = re.findall(
        r"(AWINHIBIT(?:ACETATE|METHYL)?)\n\s+WATER_ACTIVITY_THRESHOLD ([0-9.]+)",
        block,
    )
    by_name = {name: float(val) for name, val in found}
    assert by_name["AWINHIBIT"] == AW_CRIT_HYDROGENOTROPHIC
    assert by_name["AWINHIBITMETHYL"] == AW_CRIT_METHYLOTROPHIC
    assert by_name["AWINHIBITACETATE"] == AW_CRIT_ACETOCLASTIC
    assert by_name["AWINHIBITACETATE"] > by_name["AWINHIBITMETHYL"] > by_name["AWINHIBIT"]


def test_aw_crits_are_ranked_and_match_cli_imports():
    """Constants stay ranked; comparison CLI imports the same objects."""
    from pflotran_py.comparison import decks

    assert AW_CRIT_ACETOCLASTIC > AW_CRIT_METHYLOTROPHIC > AW_CRIT_HYDROGENOTROPHIC
    assert decks.AW_CRIT_HYDROGENOTROPHIC is AW_CRIT_HYDROGENOTROPHIC
    assert decks.AW_CRIT_ACETOCLASTIC is AW_CRIT_ACETOCLASTIC
    assert decks.AW_CRIT_METHYLOTROPHIC is AW_CRIT_METHYLOTROPHIC
    assert decks.AW_INHIBITION_TYPE == AW_INHIBITION_TYPE
