"""Tests for the closed-batch (sealed bottle) deck generator.

The point of these tests is to protect the one property that makes a bottle
deck a bottle deck: that it is sealed. The open-column deck this generator
diverges from flushes its entire brine out of the domain within a day, which
silently removes the salt stress the experiment exists to study. A regression
that reintroduced a boundary condition would not crash anything -- it would
just quietly produce wrong science. So the closure assertions below are the
load-bearing ones.
"""

import os

import pytest

from pflotran_py.generator.bottle_generator import (
    BOTTLE_GAS_SATURATION,
    HEADSPACE_VOLUME_L,
    VIAL_VOLUME_L,
    BottleGenerator,
    generate_bottle_series,
    nacl_brine,
    nacl_molality_for_water_activity,
    water_activity_from_nacl_molality,
)


@pytest.fixture
def deck(tmp_path):
    """A generated bottle deck at a moderate salinity, as text."""
    path = tmp_path / "bottle.in"
    BottleGenerator(brine=nacl_brine(molality=2.7)).generate(str(path))
    return path.read_text()


# ─────────────────────────────────────────────────────────────────────
# Closure — the reason this generator exists
# ─────────────────────────────────────────────────────────────────────


def _uncommented_lines(deck_text):
    """Deck lines with comments stripped, so prose about boundary conditions
    is not mistaken for a boundary condition."""
    return [
        line for line in deck_text.splitlines() if not line.lstrip().startswith("#")
    ]


@pytest.mark.parametrize(
    "forbidden",
    [
        "BOUNDARY_CONDITION",
        "FLOW_CONDITION atmospheric",
        "TRANSPORT_CONDITION atmospheric",
        "CONSTRAINT atmospheric",
        "top_surface",
        "bottom_surface",
    ],
)
def test_deck_has_no_opening_to_the_outside(deck, forbidden):
    """No construct that would let mass leave the domain."""
    body = "\n".join(_uncommented_lines(deck))
    assert forbidden not in body


def test_deck_has_exactly_one_coupler(deck):
    """An initial condition and nothing else.

    PFLOTRAN treats a face with no boundary condition as zero-flux, so this is
    what makes the domain closed.
    """
    body = "\n".join(_uncommented_lines(deck))
    assert body.count("INITIAL_CONDITION") == 1


def test_deck_keeps_the_reaction_network(deck):
    """Closure must not have cost us the chemistry inherited from the column
    generator."""
    for expected in ("AWINHIBIT", "AWINHIBITACETATE", "AWINHIBITMETHYL"):
        assert expected in deck
    assert "MICROBIAL_REACTION" in deck
    assert "ACTIVITY_WATER" in deck


# ─────────────────────────────────────────────────────────────────────
# Geometry
# ─────────────────────────────────────────────────────────────────────


def test_domain_is_a_single_cell(deck):
    assert "NXYZ 1 1 1" in deck


def test_gas_saturation_matches_the_real_headspace_ratio(deck):
    """The modelled headspace fraction must equal the real vial's."""
    assert BOTTLE_GAS_SATURATION == pytest.approx(HEADSPACE_VOLUME_L / VIAL_VOLUME_L)
    assert f"GAS_SATURATION {BOTTLE_GAS_SATURATION:.3f}" in deck


def test_cell_volume_equals_the_vial_volume(deck):
    """Parse the cell edge back out of the deck and check the volume."""
    lines = deck.splitlines()
    dxyz_index = next(i for i, line in enumerate(lines) if line.strip() == "DXYZ")
    edges = [float(lines[dxyz_index + n].strip().replace("d0", "")) for n in (1, 2, 3)]
    volume_m3 = edges[0] * edges[1] * edges[2]
    assert volume_m3 == pytest.approx(VIAL_VOLUME_L * 1e-3, rel=1e-3)


def test_temperature_matches_the_post_processing_assumption(deck):
    """config.py corrects diffusion to 18 C; the simulation must agree."""
    assert "TEMPERATURE 18.0d0" in deck


def test_run_spans_the_measured_incubation_window(deck):
    """Both measured series run past 42 days; the default must cover them."""
    final_time = next(
        line for line in deck.splitlines() if line.strip().startswith("FINAL_TIME")
    )
    days = int(final_time.split()[1])
    assert days >= 51


def test_database_path_resolves_in_this_clone(deck):
    """The column generator defaults to a path on one developer's machine."""
    database_line = next(
        line for line in deck.splitlines() if line.strip().startswith("DATABASE")
    )
    path = database_line.split(maxsplit=1)[1].strip()
    assert os.path.exists(path), f"hanford.dat not found at {path}"


# ─────────────────────────────────────────────────────────────────────
# Water activity conversion
# ─────────────────────────────────────────────────────────────────────


def test_pure_water_has_unit_activity():
    assert water_activity_from_nacl_molality(0.0) == 1.0


def test_water_activity_falls_with_salt():
    values = [water_activity_from_nacl_molality(m) for m in (0.5, 1.0, 2.0, 4.0, 6.0)]
    assert values == sorted(values, reverse=True)


def test_saturated_nacl_matches_the_literature_value():
    """Saturated NaCl (about 6.1 mol/kg at 25 C) has a water activity of about
    0.753 -- the standard humidity-calibration value. Agreement here is what
    justifies using this correlation to pick deck compositions at all."""
    assert water_activity_from_nacl_molality(6.1) == pytest.approx(0.753, abs=0.01)


@pytest.mark.parametrize("target", [0.996, 0.905, 0.850, 0.773])
def test_inversion_round_trips(target):
    molality = nacl_molality_for_water_activity(target)
    assert water_activity_from_nacl_molality(molality) == pytest.approx(
        target, abs=1e-5
    )


def test_brine_charge_balances_on_chloride():
    """Sodium is set explicitly; chloride carries the Z code so PFLOTRAN closes
    the charge balance itself."""
    brine = nacl_brine(molality=3.0)
    assert brine["Na+"].endswith(" T")
    assert brine["Cl-"].endswith(" Z")


def test_brine_molarity_is_below_molality():
    """Molarity (per litre of solution) must be below molality (per kg of
    water), because dissolving salt adds volume."""
    molality = 5.0
    brine = nacl_brine(molality=molality)
    molarity = float(brine["Na+"].split()[0])
    assert 0 < molarity < molality


def test_brine_requires_exactly_one_of_its_two_inputs():
    with pytest.raises(ValueError):
        nacl_brine()
    with pytest.raises(ValueError):
        nacl_brine(molality=1.0, water_activity=0.9)


# ─────────────────────────────────────────────────────────────────────
# Series generation
# ─────────────────────────────────────────────────────────────────────


def test_series_writes_one_deck_per_water_activity(tmp_path):
    targets = [0.996, 0.905, 0.773]
    results = generate_bottle_series(targets, output_dir=str(tmp_path))

    assert len(results) == len(targets)
    for (aw, molality, path), expected_aw in zip(results, targets):
        assert aw == expected_aw
        assert molality > 0
        assert os.path.exists(path)


def test_series_decks_differ_only_in_salt(tmp_path):
    """Everything except the brine is held fixed, so a comparison across the
    series isolates the salt effect."""
    results = generate_bottle_series([0.996, 0.905], output_dir=str(tmp_path))
    decks = [open(path).read() for _, _, path in results]

    def strip_variable_lines(text):
        return [
            line
            for line in text.splitlines()
            if not line.lstrip().startswith("#")
            and "Na+" not in line
            and "Cl-" not in line
        ]

    assert strip_variable_lines(decks[0]) == strip_variable_lines(decks[1])
