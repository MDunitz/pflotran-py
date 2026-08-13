"""Tests for AWINHIBIT sandboxes owning network methanogenesis."""

from pflotran_py.generator.bottle_generator import BottleGenerator, nacl_brine
from pflotran_py.generator.pflotran_generator import (
    DEFAULT_RATE_CONSTANTS,
    METHANOGENESIS_RATE_KEYS,
)


def test_default_bottle_deck_drops_network_methanogenesis(tmp_path):
    path = tmp_path / "bottle.in"
    BottleGenerator(brine=nacl_brine(molality=1.0)).generate(str(path))
    text = path.read_text()
    assert "# hydrogenotrophic methanogenesis" not in text
    assert "# acetoclastic methanogenesis" not in text
    assert "# methylotrophic methanogenesis" not in text
    assert "AWINHIBIT" in text
    assert "HALF_SATURATION_ACETATE" in text


def test_keep_network_methanogenesis_flag(tmp_path):
    path = tmp_path / "both.in"
    BottleGenerator(
        brine=nacl_brine(molality=1.0),
        aw_sandbox_replaces_network_methanogenesis=False,
    ).generate(str(path))
    text = path.read_text()
    assert "# hydrogenotrophic methanogenesis" in text
    assert "# acetoclastic methanogenesis" in text


def test_sandbox_off_keeps_network_methanogenesis(tmp_path):
    path = tmp_path / "no_sandbox.in"
    BottleGenerator(
        brine=nacl_brine(molality=1.0),
        enable_aw_sandbox=False,
    ).generate(str(path))
    text = path.read_text()
    assert "REACTION_SANDBOX" not in text
    assert "# hydrogenotrophic methanogenesis" in text


def test_sandbox_rates_are_network_defaults():
    for key in METHANOGENESIS_RATE_KEYS:
        assert key in DEFAULT_RATE_CONSTANTS
    assert DEFAULT_RATE_CONSTANTS["hydrogenotrophic_methano"] == 7.20e-09
    assert DEFAULT_RATE_CONSTANTS["acetaclastic_methano"] == 1.50e-08
    assert DEFAULT_RATE_CONSTANTS["methylotrophic_methano"] == 9.10e-06
