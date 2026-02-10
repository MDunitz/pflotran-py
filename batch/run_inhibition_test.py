"""
Three-variant inhibition diagnostic for double-counting analysis.

Tests whether Cl⁻ Monod inhibition and a_w reaction sandbox are
redundantly suppressing the same methanogenesis pathways.

Generates three PFLOTRAN input files at a given seawater multiplier:

  Run A  — Cl⁻ inhibition ON,  a_w sandbox OFF
  Run B  — Cl⁻ inhibition OFF, a_w sandbox ON
  Run C  — Cl⁻ inhibition ON,  a_w sandbox ON   (current default)

Physics reasoning
-----------------
Water activity (a_w) already accounts for the thermodynamic effect of
ALL dissolved ions, including Cl⁻. The Cl⁻ Monod inhibition is a
separate, empirical mechanism that only sees [Cl⁻]. If both are active,
methanogenesis is suppressed twice by the same underlying salinity signal.

Diagnostic: compare CH₄(aq) production across the three runs.
If C ≈ A·B / uninhibited  → independent mechanisms, no double-counting.
If C ≪ A·B / uninhibited  → correlated suppression, redundant.

Expected result: Cl⁻ and a_w are correlated → recommend using a_w only
(more physically complete, captures Mg²⁺, SO₄²⁻, etc.).

Affected reactions (overlap on all three methanogenesis pathways):
  - fermentation              : Cl⁻ only (no a_w sandbox)
  - hydrogenotrophic methano  : Cl⁻ + AWINHIBIT
  - acetoclastic methano      : Cl⁻ + AWINHIBITACETATE
  - methylotrophic methano    : Cl⁻ + AWINHIBITMETHYL

Usage:
    python run_inhibition_test.py                     # default 5× seawater
    python run_inhibition_test.py --multiplier 10     # 10× seawater
    python run_inhibition_test.py --multiplier 1      # baseline, no inhibition expected
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "generator"))
from pflotran_generator import PFLOTRANGenerator  # noqa: E402

# Import seawater scaling from reproducible_varconc
# create_modified_files is in the same directory
from create_modified_files import scale_seawater  # noqa: E402

VARIANTS = {
    "A_cl_only": {"enable_cl_inhibition": True, "enable_aw_sandbox": False},
    "B_aw_only": {"enable_cl_inhibition": False, "enable_aw_sandbox": True},
    "C_both": {"enable_cl_inhibition": True, "enable_aw_sandbox": True},
}


def generate_test_files(multiplier=5, output_dir="inhibition_test", **kwargs):
    """Generate three PFLOTRAN input files for the inhibition diagnostic.

    Parameters
    ----------
    multiplier : float
        Seawater salinity multiplier.
    output_dir : str
        Directory for output files.
    **kwargs
        Additional kwargs passed to PFLOTRANGenerator
        (temperature, dimensions, database_path, etc.)
    """
    os.makedirs(output_dir, exist_ok=True)
    concentrations = scale_seawater(multiplier)

    for variant_name, flags in VARIANTS.items():
        gen = PFLOTRANGenerator(
            concentrations=concentrations,
            **flags,
            **kwargs,
        )
        filename = os.path.join(output_dir, f"{variant_name}_{multiplier}x.in")
        gen.generate(filename)
        print()

    print(f"\nGenerated {len(VARIANTS)} files in {output_dir}/")
    print(f"Seawater multiplier: {multiplier}×")
    print(f"Cl⁻ concentration: {float(concentrations['Cl-'].split()[0]):.3f} M")
    print("\nRun all three with PFLOTRAN, then compare CH₄(aq) profiles.")
    print("If C << A·B/uninhibited → double-counting confirmed → use a_w only.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate inhibition diagnostic test files."
    )
    parser.add_argument(
        "--multiplier", type=float, default=5.0, help="Seawater multiplier (default: 5)"
    )
    parser.add_argument("--output-dir", default="inhibition_test")
    parser.add_argument("--temperature", type=float, default=8.0)
    parser.add_argument("--dimensions", default="1d", choices=["1d", "2d", "3d"])
    args = parser.parse_args()

    generate_test_files(
        multiplier=args.multiplier,
        output_dir=args.output_dir,
        temperature=args.temperature,
        dimensions=args.dimensions,
    )
