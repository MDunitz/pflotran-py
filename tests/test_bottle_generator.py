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
    generator. Methanogenesis lives in the AWINHIBIT sandboxes."""
    for expected in ("AWINHIBIT", "AWINHIBITACETATE", "AWINHIBITMETHYL"):
        assert expected in deck
    assert "MICROBIAL_REACTION" in deck
    assert "ACTIVITY_WATER" in deck
    assert "HALF_SATURATION_H2" in deck
    assert "INHIBITION_TYPE SMOOTHSTEP" in deck
    # Network methanogenesis is omitted -- the sandboxes own those pathways.
    assert "# hydrogenotrophic methanogenesis" not in deck
    assert "# acetoclastic methanogenesis" not in deck
    assert "# methylotrophic methanogenesis" not in deck


def test_sandbox_rates_match_the_network_defaults(deck):
    """Sandbox RATE_CONSTANT values are the network methanogenesis rates."""
    assert "RATE_CONSTANT 7.20e-09" in deck  # hydrogenotrophic
    assert "RATE_CONSTANT 1.50e-08" in deck  # acetoclastic
    assert "RATE_CONSTANT 9.10e-06" in deck  # methylotrophic
    assert "WATER_ACTIVITY_THRESHOLD 0.9500" in deck


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
    """The pipeline's output runs to 119 days for Exp003 and 122 for Exp004.
    A shorter simulation would stop partway through the measured record, and
    the comparison would silently cover only part of the experiment."""
    final_time = next(
        line for line in deck.splitlines() if line.strip().startswith("FINAL_TIME")
    )
    days = int(final_time.split()[1])
    assert days >= 122


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


# ─────────────────────────────────────────────────────────────────────
# Salinity inhibition on the reaction network
# ─────────────────────────────────────────────────────────────────────


def test_salinity_inhibition_is_off_by_default(deck):
    """Cl- smoothstep on the network is off unless asked for.

    The a_w sandboxes still use INHIBITION_TYPE SMOOTHSTEP; that is a different
    term. The Cl- one is the one that writes SMOOTHSTEP_INTERVAL.
    """
    assert "SMOOTHSTEP_INTERVAL" not in deck
    assert "INHIBITION_TYPE SMOOTHSTEP" in deck


def _salted_deck(tmp_path, **spec):
    path = tmp_path / "salted.in"
    BottleGenerator(
        brine=nacl_brine(molality=2.7),
        # Extra Cl- smoothstep on the network methanogenesis reactions -- keep
        # those reactions rather than replacing them with the a_w sandboxes.
        aw_sandbox_replaces_network_methanogenesis=False,
        salinity_inhibition={
            "species": "Cl-",
            "threshold": 1.0,
            "interval": 0.5,
            **spec,
        },
    ).generate(str(path))
    return path.read_text()


def test_salinity_inhibition_lands_on_the_methanogenesis_reactions(tmp_path):
    """Three reactions produce methane, and the term has to be on all three."""
    deck = _salted_deck(tmp_path)
    # Sandbox SMOOTHSTEP (a_w) plus three Cl- smoothsteps on methanogenesis.
    assert deck.count("TYPE SMOOTHSTEP") >= 3
    assert deck.count("SMOOTHSTEP_INTERVAL 0.50") == 3


def test_salinity_inhibition_is_sigmoidal_not_hyperbolic(tmp_path):
    """Monod inhibition is hyperbolic and cannot express a collapse; across
    the measured brines it delivers at most about twentyfold, against four
    orders of magnitude in the data. SMOOTHSTEP can."""
    deck = _salted_deck(tmp_path)
    assert "TYPE SMOOTHSTEP" in deck
    assert "SMOOTHSTEP_INTERVAL" in deck


def test_salinity_inhibition_carries_its_parameters(tmp_path):
    deck = _salted_deck(tmp_path, threshold=1.5, interval=1.0)
    assert "SMOOTHSTEP_INTERVAL 1.00" in deck
    assert "THRESHOLD_CONCENTRATION 1.50e+00" in deck
    assert "INHIBIT_ABOVE_THRESHOLD" in deck


def test_salinity_inhibition_does_not_touch_the_oxidation_steps(tmp_path):
    """Salt stress in these incubations is understood to act on the
    methanogens. Inhibiting fermentation or the oxidation steps too would
    suppress the whole carbon chain instead."""
    deck = _salted_deck(tmp_path)
    methane_oxidation = deck[deck.index("methane oxidation (O2)") :][:600]
    assert "TYPE SMOOTHSTEP" not in methane_oxidation


# ─────────────────────────────────────────────────────────────────────
# Carbon inventory
# ─────────────────────────────────────────────────────────────────────


def test_carbon_stays_dissolved_by_default(deck):
    """Existing decks are unchanged unless hydrolysis is asked for."""
    assert "Cellulose_min" not in deck
    assert "DOM1                5.00 T" in deck


def _hydrolysis_deck(tmp_path, **spec):
    path = tmp_path / "hydrolysis.in"
    BottleGenerator(brine=nacl_brine(molality=2.7), cellulose_hydrolysis=spec).generate(
        str(path)
    )
    return path.read_text()


def test_an_empty_dict_switches_hydrolysis_on(tmp_path):
    """An empty dict means "on, with the defaults". Testing the option for
    truthiness rather than for None would read it as "off" and silently do
    nothing, which is how this was first written and first broke."""
    deck = _hydrolysis_deck(tmp_path)
    assert "Cellulose_min" in deck


def test_hydrolysis_moves_carbon_out_of_solution(tmp_path):
    """The dissolved pool drops from molar to millimolar, and the bulk of the
    carbon moves into a solid that dissolves into it."""
    deck = _hydrolysis_deck(tmp_path)
    assert "DOM1                5.00 T" not in deck
    assert "DOM1                1.00d-03 T" in deck


def test_hydrolysis_declares_the_mineral_everywhere_it_is_needed(tmp_path):
    """A mineral has to appear in the MINERALS list, in MINERAL_KINETICS and in
    the constraint, or PFLOTRAN either ignores it or refuses the deck."""
    deck = _hydrolysis_deck(tmp_path)
    assert deck.count("Cellulose_min") >= 3
    assert "RATE_CONSTANT  2.d-8 mol/m^2-sec" in deck


def test_hydrolysis_rate_is_overridable(tmp_path):
    deck = _hydrolysis_deck(tmp_path, rate_constant="5.d-9")
    assert "RATE_CONSTANT  5.d-9 mol/m^2-sec" in deck


def test_the_bottle_database_repairs_the_cellulose_record(tmp_path):
    """The record in hanford.dat declares two species but carries two surplus
    fields, and its second species has a zero coefficient that PFLOTRAN drops.
    The deck is then refused with a species-count mismatch. Left unrepaired,
    hydrolysis cannot run at all."""
    from pflotran_py.generator.bottle_generator import (
        bottle_database_path,
        write_bottle_database,
    )

    path = write_bottle_database(destination_path=str(tmp_path / "db.dat"))
    line = next(line for line in open(path) if line.startswith("'Cellulose_min'"))
    fields = line.split()
    species_count = int(fields[2])
    # name, molar volume, count, one pair per species, eight log K, molar mass
    assert len(fields) == 3 + 2 * species_count + 8 + 1
    assert species_count == 1
    assert "'DOM1'" in line
    assert bottle_database_path().endswith("hanford_bottle.dat")
