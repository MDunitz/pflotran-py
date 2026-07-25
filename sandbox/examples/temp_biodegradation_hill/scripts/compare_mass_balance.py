#!/usr/bin/env python3
"""Mass-balance and regression checks for TempBioHill / FlexBioHill batch runs.

Stoichiometry (Hammond 2022 R1):
  Aaq + 0.25 Baq -> 0.33 Caq + Daq

For a closed batch cell (no transport), reaction extent ξ from ΔAaq should
satisfy:
  ΔBaq ≈ 0.25 * ΔAaq
  ΔCaq ≈ 0.33 * ΔAaq
  ΔDaq ≈ 1.00 * ΔAaq
(within numerical tolerance; biomass yield/decay does not affect aqueous stoich).

Also compares final free-ion concentrations against the PFLOTRAN gold file
flexible_biodegradation_hill.regression.gold when available.

Usage:
  python compare_mass_balance.py path/to/run-dir
  python compare_mass_balance.py path/to/foo-obs-0.pft --gold ../gold/...
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

# Initial conditions from the paper / decks (mol/L aqueous; Xim mol/m^3)
A0, B0, C0, D0 = 1.0e-3, 5.0e-4, 1.0e-10, 1.0e-10
STOICH = {"Aaq": 1.0, "Baq": 0.25, "Caq": 0.33, "Daq": 1.0}


def parse_obs_pft(path: Path) -> dict[str, np.ndarray]:
    """Parse a PFLOTRAN observation .pft (Tecplot-like) file."""
    text = path.read_text()
    # Collect all quoted names on the VARIABLES line
    lines = text.splitlines()
    var_names: list[str] = []
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().upper().startswith("VARIABLES"):
            var_names = re.findall(r'"([^"]+)"', line)
            data_start = i + 1
            break
    if not var_names:
        # Fallback: whitespace header
        raise ValueError(f"Could not parse VARIABLES from {path}")

    # Skip ZONE header lines
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


def _find_series(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Match Free Aaq / Aaq etc. flexibly."""
    key_l = key.lower()
    for name, vals in data.items():
        n = name.lower().replace(" ", "")
        if key_l in n or n.endswith(key_l) or n == key_l:
            return vals
    raise KeyError(f"{key} not found in columns: {list(data)}")


def check_mass_balance(obs_path: Path, rtol: float = 0.05) -> dict:
    data = parse_obs_pft(obs_path)
    A = _find_series(data, "aaq")[-1]
    B = _find_series(data, "baq")[-1]
    C = _find_series(data, "caq")[-1]
    D = _find_series(data, "daq")[-1]

    dA = A0 - A
    extent = dA / STOICH["Aaq"]
    expected = {
        "Baq": B0 - STOICH["Baq"] * extent,
        "Caq": C0 + STOICH["Caq"] * extent,
        "Daq": D0 + STOICH["Daq"] * extent,
    }
    actual = {"Baq": B, "Caq": C, "Daq": D}
    report = {
        "obs": str(obs_path),
        "final": {"Aaq": A, "Baq": B, "Caq": C, "Daq": D},
        "extent_from_Aaq": extent,
        "checks": {},
    }
    ok = True
    for sp in ("Baq", "Caq", "Daq"):
        exp, act = expected[sp], actual[sp]
        # relative to characteristic scale of extent contribution
        scale = max(abs(exp), abs(STOICH[sp] * extent), 1e-12)
        err = abs(act - exp) / scale
        report["checks"][sp] = {
            "expected": exp,
            "actual": act,
            "rel_err": err,
            "pass": err <= rtol,
        }
        ok = ok and err <= rtol
    report["pass"] = ok
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


def compare_to_gold(obs_path: Path, gold_path: Path, rtol: float = 1e-4) -> dict:
    data = parse_obs_pft(obs_path)
    gold = parse_gold(gold_path)
    report = {"obs": str(obs_path), "gold": str(gold_path), "checks": {}}
    ok = True
    for sp in ("Aaq", "Baq", "Caq", "Daq", "Xim"):
        if sp not in gold:
            continue
        act = _find_series(data, sp.lower())[-1]
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "target",
        type=Path,
        help="Run directory or *-obs-0.pft file",
    )
    ap.add_argument(
        "--gold",
        type=Path,
        default=None,
        help="Path to flexible_biodegradation_hill.regression.gold",
    )
    ap.add_argument("--rtol-stoich", type=float, default=0.05)
    ap.add_argument("--rtol-gold", type=float, default=1e-3)
    args = ap.parse_args()

    if args.target.is_dir():
        obs_files = sorted(args.target.glob("*-obs-0.pft"))
        if not obs_files:
            raise SystemExit(f"No *-obs-0.pft in {args.target}")
    else:
        obs_files = [args.target]

    example_root = Path(__file__).resolve().parents[1]
    gold = args.gold or (
        example_root / "gold" / "flexible_biodegradation_hill.regression.gold"
    )

    all_ok = True
    for obs in obs_files:
        print(f"\n=== {obs.name} ===")
        mb = check_mass_balance(obs, rtol=args.rtol_stoich)
        print(f"Mass balance (stoich): {'PASS' if mb['pass'] else 'FAIL'}")
        print(f"  extent(Aaq) = {mb['extent_from_Aaq']:.6e}")
        for sp, c in mb["checks"].items():
            print(
                f"  {sp}: act={c['actual']:.6e} exp={c['expected']:.6e} "
                f"rel_err={c['rel_err']:.3e} {'OK' if c['pass'] else 'BAD'}"
            )
        all_ok = all_ok and mb["pass"]

        if gold.exists():
            g = compare_to_gold(obs, gold, rtol=args.rtol_gold)
            print(f"Gold comparison: {'PASS' if g['pass'] else 'FAIL'} ({gold.name})")
            for sp, c in g["checks"].items():
                print(
                    f"  {sp}: act={c['actual']:.6e} gold={c['gold']:.6e} "
                    f"rel_err={c['rel_err']:.3e} {'OK' if c['pass'] else 'BAD'}"
                )
            # Only expect gold match for FlexBioHill / Q10=1 decks
            if "flexbiohill" in obs.name.lower() or "q10_1" in obs.name.lower():
                all_ok = all_ok and g["pass"]
        else:
            print(f"(no gold file at {gold})")

    raise SystemExit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
