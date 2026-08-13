#!/usr/bin/env python3
"""CLI wrapper around tests/mass_balance_checks.py for ad-hoc run-dir checks.

Primary coverage lives in tests/test_tempbiohill_mass_balance.py (pytest).
This script remains for inspecting PFLOTRAN output directories manually:

  python compare_mass_balance.py path/to/run-dir
  python compare_mass_balance.py path/to/foo-obs-0.pft --gold ../gold/...
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[4]
_TESTS = _REPO_ROOT / "tests"
sys.path.insert(0, str(_TESTS))

from mass_balance_checks import (  # noqa: E402
    check_mass_balance_from_obs,
    compare_to_gold,
)

# Back-compat aliases expected by older call sites / docs.
A0 = 1.0e-3
B0 = 5.0e-4
C0 = 1.0e-10
D0 = 1.0e-10
STOICH = {"Aaq": 1.0, "Baq": 0.25, "Caq": 0.33, "Daq": 1.0}


def check_mass_balance(obs_path: Path, rtol: float = 0.05) -> dict:
    """FlexBioHill convenience wrapper (same signature as the old script)."""
    return check_mass_balance_from_obs(obs_path, rtol=rtol)


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
        extent_key = next(k for k in mb if k.startswith("extent_from_"))
        print(f"  {extent_key} = {mb[extent_key]:.6e}")
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
            if "flexbiohill" in obs.name.lower() or "q10_1" in obs.name.lower():
                all_ok = all_ok and g["pass"]
        else:
            print(f"(no gold file at {gold})")

    raise SystemExit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
