"""
Generate PFLOTRAN input files at varying seawater concentrations.

Produces one .in file per multiplier (1×–20× seawater) by passing
scaled ion concentrations directly to PFLOTRANGenerator.

No regex surgery — the generator accepts concentrations as kwargs.

Usage:
    python create_modified_files.py                     # default: 1×–20× seawater
    python create_modified_files.py --max-multiplier 5  # 1×–5× only
"""

import os
import sys
import argparse

# Add the master_input_generator directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'generator'))
from pflotran_generator import PFLOTRANGenerator


# Seawater baseline concentrations [mol/L] at 1× salinity
# Ref: Millero (2013) Chemical Oceanography, Table 2.1
SEAWATER_1X = {
    'Cl-':   5.36e-01,
    'Na+':   4.59e-01,
    'Mg++':  5.23e-02,
    'SO4--': 2.76e-02,
    'Ca++':  1.00e-02,
    'K+':    9.72e-03,
}


def scale_seawater(multiplier):
    """Scale seawater ions by a multiplier and format as PFLOTRAN constraint strings.

    Parameters
    ----------
    multiplier : float
        Salinity multiplier (e.g. 1.0 = normal seawater, 5.0 = 5× seawater).

    Returns
    -------
    dict
        Species → PFLOTRAN constraint string, e.g. {'Cl-': '2.68e+00 Z', ...}
        Cl- uses 'Z' (charge balance), all others use 'T' (total concentration).
    """
    scaled = {}
    for ion, baseline in SEAWATER_1X.items():
        conc = baseline * multiplier
        constraint_type = 'Z' if ion == 'Cl-' else 'T'
        scaled[ion] = f'{conc:.3e} {constraint_type}'
    return scaled


def create_modified_files(
    output_dir='modified_pflotran_files',
    min_multiplier=1,
    max_multiplier=20,
    **generator_kwargs,
):
    """Generate PFLOTRAN input files for a range of seawater multipliers.

    Parameters
    ----------
    output_dir : str
        Directory for output .in files.
    min_multiplier : int
        Starting multiplier (inclusive).
    max_multiplier : int
        Ending multiplier (inclusive).
    **generator_kwargs
        Additional kwargs passed to PFLOTRANGenerator (e.g. aw_threshold,
        dimensions, temperature).
    """
    os.makedirs(output_dir, exist_ok=True)

    for multiplier in range(min_multiplier, max_multiplier + 1):
        concentrations = scale_seawater(multiplier)
        generator = PFLOTRANGenerator(
            concentrations=concentrations,
            **generator_kwargs,
        )
        filename = os.path.join(output_dir, f'9_addnitrogen_{multiplier}M.in')
        generator.generate(filename)

    n_files = max_multiplier - min_multiplier + 1
    print(f'\nGenerated {n_files} files in {output_dir}/')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate PFLOTRAN input files at varying salinity.')
    parser.add_argument('--output-dir', default='modified_pflotran_files')
    parser.add_argument('--min-multiplier', type=int, default=1)
    parser.add_argument('--max-multiplier', type=int, default=20)
    parser.add_argument('--aw-threshold', type=float, default=0.5)
    parser.add_argument('--dimensions', default='1d', choices=['1d', '2d', '3d'])
    parser.add_argument('--temperature', type=float, default=8.0)
    args = parser.parse_args()

    create_modified_files(
        output_dir=args.output_dir,
        min_multiplier=args.min_multiplier,
        max_multiplier=args.max_multiplier,
        aw_threshold=args.aw_threshold,
        dimensions=args.dimensions,
        temperature=args.temperature,
    )

