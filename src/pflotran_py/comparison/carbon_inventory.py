"""Starting carbon inventory for closed-batch bottle decks.

Comparison figures can use absolute headspace moles or moles per mole of
starting carbon. Both sides of a normalised figure use the **same** inventory
denominator.

Measured starting C
-------------------
Exp003 / Exp004 have no lab TOC / volatile-solids assay. The saltyBiomass
pipeline instead derives dry-biomass carbon from the batch recipe::

    dry biomass (g) = incubation mass
        × (sludge fraction × (1 − 0.732)
           + spirulina fraction × (1 − 0.05))
    moles C = dry biomass × 0.5 / 12.011 g mol⁻¹

Across those experiments that value is nearly constant at about
:data:`MEASURED_INCUBATION_STARTING_CARBON_MOLES` (≈ 0.0565 mol C / bottle).

Model inventory
---------------
With ``--cellulose-hydrolysis``, the deck holds carbon in solid
``Cellulose_min``. The default volume fraction is chosen so model starting C
**matches** that measured inventory (inventory normalisation, not a kinetics
fit). DOM1 is glucose (six carbons); ``Cellulose_min`` dissolves one-to-one to
DOM1, so each mole of mineral is six moles of carbon.
"""

from __future__ import annotations

import re
from typing import Mapping

from ..generator.bottle_generator import LIQUID_VOLUME_L, VIAL_VOLUME_L
from ..generator.pflotran_generator import (
    DEFAULT_CELLULOSE_HYDROLYSIS,
    DEFAULT_INITIAL_CONCENTRATIONS,
)

# Glucose / anhydroglucose unit carbon count.
CARBONS_PER_DOM1 = 6

# Molar volume of Cellulose_min in hanford.dat [m^3 / mol].
# 162.14 cm^3/mol = 162.14e-6 m^3/mol.
CELLULOSE_MOLAR_VOLUME_M3_PER_MOL = 162.14e-6

# Median recipe-derived dry-biomass carbon (mol C / bottle) for Exp003 + Exp004,
# recovered from saltyBiomass pipeline percent-of-starting-C columns (same
# physics as Dry Biomass moles C). Batch-to-batch spread is only a few percent.
MEASURED_INCUBATION_STARTING_CARBON_MOLES = 0.0565

_CONSTRAINT_VALUE_RE = re.compile(
    r"^\s*([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[dDeE][-+]?\d+)?)", re.ASCII
)


def parse_constraint_value(constraint: str) -> float:
    """Leading numeric value from a PFLOTRAN constraint string (``5.00 T``)."""
    match = _CONSTRAINT_VALUE_RE.match(str(constraint))
    if not match:
        raise ValueError(f"Could not parse constraint value from {constraint!r}")
    return float(match.group(1).replace("d", "e").replace("D", "e"))


def cellulose_volume_fraction_for_starting_carbon(
    moles_c: float,
    *,
    liquid_volume_l: float = LIQUID_VOLUME_L,
    bulk_volume_l: float = VIAL_VOLUME_L,
    carbons_per_dom1: int = CARBONS_PER_DOM1,
    dissolved_dom1_mol_l: float = 1.0e-3,
) -> float:
    """``Cellulose_min`` bulk volume fraction that yields ``moles_c`` at t=0.

    Subtracts the millimolar dissolved DOM1 contribution, then converts the
    remaining solid carbon through six carbons per mineral mole and the
    database molar volume.
    """
    aqueous_c = dissolved_dom1_mol_l * liquid_volume_l * carbons_per_dom1
    solid_c = float(moles_c) - aqueous_c
    if solid_c <= 0:
        raise ValueError(
            f"Target starting C ({moles_c} mol) is not larger than the "
            f"dissolved DOM1 contribution ({aqueous_c} mol C)."
        )
    moles_mineral = solid_c / carbons_per_dom1
    return moles_mineral * CELLULOSE_MOLAR_VOLUME_M3_PER_MOL / (bulk_volume_l * 1e-3)


def starting_carbon_moles(
    *,
    cellulose_hydrolysis: Mapping[str, object] | bool | None = True,
    liquid_volume_l: float = LIQUID_VOLUME_L,
    bulk_volume_l: float = VIAL_VOLUME_L,
    carbons_per_dom1: int = CARBONS_PER_DOM1,
    dissolved_dom1_mol_l: float | None = None,
) -> float:
    """Moles of carbon atoms present at ``t = 0`` in one modelled bottle.

    Parameters
    ----------
    cellulose_hydrolysis
        ``True`` (default) uses :data:`DEFAULT_CELLULOSE_HYDROLYSIS`, matching
        the comparison pipeline's ``--cellulose-hydrolysis`` decks. ``None`` /
        ``False`` uses the dissolved-only inventory (5 mol/L DOM1). A mapping
        overrides individual solid-pool fields.
    dissolved_dom1_mol_l
        Override for the aqueous DOM1 concentration [mol/L]. When omitted, the
        value is taken from the cellulose spec or the default concentrations.
    """
    if cellulose_hydrolysis is True:
        spec: dict | None = dict(DEFAULT_CELLULOSE_HYDROLYSIS)
    elif cellulose_hydrolysis is False or cellulose_hydrolysis is None:
        spec = None
    else:
        spec = {**DEFAULT_CELLULOSE_HYDROLYSIS, **dict(cellulose_hydrolysis)}

    if dissolved_dom1_mol_l is None:
        if spec is not None:
            dissolved_dom1_mol_l = parse_constraint_value(str(spec["dom1_initial"]))
        else:
            dissolved_dom1_mol_l = parse_constraint_value(
                DEFAULT_INITIAL_CONCENTRATIONS["DOM1"]
            )

    moles_c = dissolved_dom1_mol_l * liquid_volume_l * carbons_per_dom1

    if spec is not None:
        volume_fraction = float(spec["volume_fraction"])
        moles_mineral = (
            volume_fraction * (bulk_volume_l * 1e-3) / CELLULOSE_MOLAR_VOLUME_M3_PER_MOL
        )
        moles_c += moles_mineral * carbons_per_dom1

    return float(moles_c)


def default_comparison_starting_carbon_moles() -> float:
    """Starting C used by the README comparison pipeline (cellulose on)."""
    return starting_carbon_moles(cellulose_hydrolysis=True)
