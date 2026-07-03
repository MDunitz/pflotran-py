"""Smoke tests to verify the package imports and constants are sane."""

import sys
import os

# Allow imports from repo root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from exploratory.constants import T, molar_mass_C


def test_exploratory_constants_import():
    # 0 degC in Kelvin
    assert T == 273.15
    # Carbon molar mass in g/mol
    assert 12.0 < molar_mass_C < 12.02
