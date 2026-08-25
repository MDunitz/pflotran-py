"""Reusable stoichiometric mass-balance helpers for closed batch cells.

Used by pytest (primary) and the thin TempBioHill example CLI.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

# Hammond 2022 FlexBioHill R1 defaults (mol/L aqueous; Xim mol/m^3)
FLEXBIOHILL_INITIAL = {
    "Aaq": 1.0e-3,
    "Baq": 5.0e-4,
    "Caq": 1.0e-10,
    "Daq": 1.0e-10,
}
# Stoich coefficients relative to extent ξ defined from the donor (Aaq).
# Sign convention: reactants negative, products positive on the reaction
# Aaq + 0.25 Baq -> 0.33 Caq + Daq, with ξ = (A0 - A) / 1.0.
FLEXBIOHILL_STOICH = {
    "Aaq": -1.0,
    "Baq": -0.25,
    "Caq": 0.33,
    "Daq": 1.0,
}
FLEXBIOHILL_EXTENT_SPECIES = "Aaq"


def parse_obs_pft(path: Path) -> dict[str, np.ndarray]:
    """Parse a PFLOTRAN observation .pft (Tecplot-like) file."""
    text = path.read_text()
    lines = text.splitlines()
    var_names: list[str] = []
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().upper().startswith("VARIABLES"):
            var_names = re.findall(r'"([^"]+)"', line)
            data_start = i + 1
            break
    if not var_names:
        raise ValueError(f"Could not parse VARIABLES from {path}")

    while data_start < len(lines) and (
        not lines[data_start].strip()
        or lines[data_start].strip().upper().startswith("ZONE")
    ):
        data_start += 1

    rows = []
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) < len(var_names):
            continue
        try:
            rows.append([float(x) for x in parts[: len(var_names)]])
        except ValueError:
            continue
    arr = np.asarray(rows, dtype=float)
    return {name: arr[:, i] for i, name in enumerate(var_names)}


def find_series(data: Mapping[str, np.ndarray], key: str) -> np.ndarray:
    """Match Free Aaq / Aaq etc. flexibly."""
    key_l = key.lower()
    for name, vals in data.items():
        n = name.lower().replace(" ", "")
        if key_l in n or n.endswith(key_l) or n == key_l:
            return vals
    raise KeyError(f"{key} not found in columns: {list(data)}")


def check_species_against_extent(
    *,
    species: str,
    actual: float,
    initial: float,
    stoich: float,
    extent: float,
    rtol: float,
) -> dict:
    """Compare one species final concentration to c0 + stoich * extent."""
    expected = initial + stoich * extent
    scale = max(abs(expected), abs(stoich * extent), 1e-12)
    err = abs(actual - expected) / scale
    return {
        "expected": expected,
        "actual": actual,
        "rel_err": err,
        "pass": err <= rtol,
    }


def check_mass_balance(
    finals: Mapping[str, float],
    *,
    initial: Mapping[str, float],
    stoich: Mapping[str, float],
    extent_species: str,
    check_species: Iterable[str] | None = None,
    rtol: float = 0.05,
) -> dict:
    """Stoichiometric mass balance for an arbitrary set of aqueous species.

    ``extent`` is derived from ``extent_species`` as
    ``(c0 - c_final) / (-stoich[extent_species])`` so reactant stoich is
    negative. Other species are checked via
    ``c_final ≈ c0 + stoich[sp] * extent``.
    """
    if extent_species not in finals:
        raise KeyError(f"extent species {extent_species!r} missing from finals")
    if extent_species not in initial or extent_species not in stoich:
        raise KeyError(f"extent species {extent_species!r} missing from initial/stoich")

    extent_coeff = stoich[extent_species]
    if extent_coeff == 0:
        raise ValueError(f"stoich[{extent_species!r}] must be non-zero")

    c_ext = finals[extent_species]
    c0_ext = initial[extent_species]
    extent = (c0_ext - c_ext) / (-extent_coeff)

    if check_species is None:
        check_species = [sp for sp in stoich if sp != extent_species]

    report: dict = {
        "final": dict(finals),
        f"extent_from_{extent_species}": extent,
        "checks": {},
    }
    ok = True
    for sp in check_species:
        result = check_species_against_extent(
            species=sp,
            actual=finals[sp],
            initial=initial[sp],
            stoich=stoich[sp],
            extent=extent,
            rtol=rtol,
        )
        report["checks"][sp] = result
        ok = ok and result["pass"]
    report["pass"] = ok
    return report


def check_mass_balance_from_obs(
    obs_path: Path,
    *,
    initial: Mapping[str, float] = FLEXBIOHILL_INITIAL,
    stoich: Mapping[str, float] = FLEXBIOHILL_STOICH,
    extent_species: str = FLEXBIOHILL_EXTENT_SPECIES,
    check_species: Iterable[str] | None = None,
    rtol: float = 0.05,
) -> dict:
    """Load finals from a .pft observation file and run mass-balance checks."""
    data = parse_obs_pft(obs_path)
    species = set(initial) | set(stoich)
    if check_species is not None:
        species |= set(check_species)
    species.add(extent_species)
    finals = {sp: float(find_series(data, sp.lower())[-1]) for sp in species}
    report = check_mass_balance(
        finals,
        initial=initial,
        stoich=stoich,
        extent_species=extent_species,
        check_species=check_species,
        rtol=rtol,
    )
    report["obs"] = str(obs_path)
    return report


def parse_gold(path: Path) -> dict[str, float]:
    """Parse PFLOTRAN .regression.gold free-ion finals."""
    out: dict[str, float] = {}
    text = path.read_text()
    current = None
    for line in text.splitlines():
        m = re.match(r"-- CONCENTRATION:\s*(?:Free\s+)?(\S+)\s*--", line)
        if m:
            current = m.group(1)
            continue
        if current:
            nums = re.findall(r"([-+]?\d+\.\d+E[-+]?\d+)", line, flags=re.I)
            if nums:
                out[current] = float(nums[-1])
                current = None
    return out


def compare_to_gold(
    obs_path: Path,
    gold_path: Path,
    *,
    species: Iterable[str] = ("Aaq", "Baq", "Caq", "Daq", "Xim"),
    rtol: float = 1e-4,
) -> dict:
    data = parse_obs_pft(obs_path)
    gold = parse_gold(gold_path)
    report: dict = {"obs": str(obs_path), "gold": str(gold_path), "checks": {}}
    ok = True
    for sp in species:
        if sp not in gold:
            continue
        act = float(find_series(data, sp.lower())[-1])
        exp = gold[sp]
        scale = max(abs(exp), 1e-30)
        err = abs(act - exp) / scale
        report["checks"][sp] = {
            "gold": exp,
            "actual": act,
            "rel_err": err,
            "pass": err <= rtol,
        }
        ok = ok and err <= rtol
    report["pass"] = ok
    return report
