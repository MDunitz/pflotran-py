"""Stoichiometric mass-balance checks for FlexBioHill / TempBioHill batch runs.

These live under tests/ (not examples/) so they run in CI. Synthetic .pft
fixtures cover the checker without needing a PFLOTRAN binary.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mass_balance_checks import (
    FLEXBIOHILL_INITIAL,
    FLEXBIOHILL_STOICH,
    check_mass_balance,
    check_mass_balance_from_obs,
    check_species_against_extent,
    parse_obs_pft,
)


def _write_obs_pft(path: Path, finals: dict[str, float]) -> Path:
    names = ["Time [y]"] + [f"Free {sp} [mol/L]" for sp in finals]
    # Two rows: t=0 (unused by checker, which takes the last row) and final.
    row0 = [0.0] + [FLEXBIOHILL_INITIAL.get(sp, 0.0) for sp in finals]
    row1 = [1.0] + [finals[sp] for sp in finals]
    lines = [
        "VARIABLES = " + ", ".join(f'"{n}"' for n in names),
        'ZONE T="obs"',
        " ".join(f"{v:.6e}" for v in row0),
        " ".join(f"{v:.6e}" for v in row1),
        "",
    ]
    path.write_text("\n".join(lines))
    return path


def test_check_species_against_extent_pass_and_fail():
    ok = check_species_against_extent(
        species="Baq",
        actual=4.75e-4,
        initial=5.0e-4,
        stoich=-0.25,
        extent=1.0e-4,
        rtol=0.05,
    )
    assert ok["pass"]
    assert ok["expected"] == pytest.approx(4.75e-4)

    bad = check_species_against_extent(
        species="Baq",
        actual=1.0e-3,
        initial=5.0e-4,
        stoich=-0.25,
        extent=1.0e-4,
        rtol=0.05,
    )
    assert not bad["pass"]


def test_check_mass_balance_flexbiohill_stoich():
    # Extent ξ = 1e-4 from Aaq: A goes 1e-3 -> 9e-4
    finals = {
        "Aaq": 9.0e-4,
        "Baq": 5.0e-4 - 0.25e-4,
        "Caq": 1.0e-10 + 0.33e-4,
        "Daq": 1.0e-10 + 1.0e-4,
    }
    report = check_mass_balance(
        finals,
        initial=FLEXBIOHILL_INITIAL,
        stoich=FLEXBIOHILL_STOICH,
        extent_species="Aaq",
    )
    assert report["pass"]
    assert report["extent_from_Aaq"] == pytest.approx(1.0e-4)
    assert set(report["checks"]) == {"Baq", "Caq", "Daq"}


def test_check_mass_balance_arbitrary_species_set():
    initial = {"Edonor": 2.0, "Eacceptor": 1.0, "Product": 0.0}
    stoich = {"Edonor": -1.0, "Eacceptor": -0.5, "Product": 1.5}
    finals = {
        "Edonor": 1.5,  # ξ = 0.5
        "Eacceptor": 0.75,
        "Product": 0.75,
    }
    report = check_mass_balance(
        finals,
        initial=initial,
        stoich=stoich,
        extent_species="Edonor",
    )
    assert report["pass"]
    assert report["extent_from_Edonor"] == pytest.approx(0.5)
    assert set(report["checks"]) == {"Eacceptor", "Product"}


def test_check_mass_balance_from_obs_roundtrip(tmp_path: Path):
    finals = {
        "Aaq": 9.0e-4,
        "Baq": 5.0e-4 - 0.25e-4,
        "Caq": 1.0e-10 + 0.33e-4,
        "Daq": 1.0e-10 + 1.0e-4,
    }
    obs = _write_obs_pft(tmp_path / "batch-obs-0.pft", finals)
    parsed = parse_obs_pft(obs)
    assert "Free Aaq [mol/L]" in parsed

    report = check_mass_balance_from_obs(obs)
    assert report["pass"]
    assert report["obs"] == str(obs)
