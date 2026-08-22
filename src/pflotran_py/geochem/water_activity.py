"""Water activity of multi-component brines via PHREEQC's Pitzer model.

Computes the solvent water activity ``a_w`` of an aqueous electrolyte mixture by
handing its composition to PHREEQC with the ``pitzer.dat`` database and reading
back ``ACT("H2O")`` -- the full Harvie-Weare-Pitzer virial model (including the
cation-cation / anion-anion mixing terms) for
Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O at 25 degC.

Scale
-----
PHREEQC is fed ``units mol/kgw`` (molality). The molarity -> molality step is
applied first, in ``conversions.py``, rather than letting PHREEQC do it:
PHREEQC's own ``mol/l`` conversion assumes rho ~= 1 kg/L for ``pitzer.dat`` and
over-concentrates a multi-molar brine (a 1.94 M MgCl2 batch lands at I ~= 7.1
mol/kg instead of ~= 5.8, dropping a_w from ~0.84 to ~0.80). Supplying molality
directly keeps the density estimate explicit and salt-loading aware.

Water activity is set by the salt ions; H+/OH- at any physical brine pH are
negligible against a multi-molar background, so the PHREEQC input carries no pH
(PHREEQC's neutral default is used) and no carbonate system.

Validation (pitzer.dat, 25 degC):
    NaCl 6.14 mol/kg  -> a_w = 0.7524  (accepted 0.753)
    seawater          -> a_w = 0.9813  (accepted ~0.981)
    1.94 m MgCl2 batch (density-corrected molality) -> a_w = 0.842 (meter 0.824)

Temperature is fixed at 25 degC (``pitzer.dat``'s calibration point); per-batch
incubation temperature is not threaded through.
"""

import math

import astropy.units as u

from .constants import (
    BATCH_ION_TO_PITZER,
    ION_CHARGE,
    M_WATER,
    PITZER_TO_PHREEQC_ELEMENT,
)
from .conversions import molarities_to_molalities

_MOLAL = u.mol / u.kg

# pitzer.dat is calibrated at 25 degC; temperature dependence is not modelled.
TEMPERATURE_C = 25.0

_PHREEQC_INPUT_TEMPLATE = """
SOLUTION 1
    units {units}
    temp {temp}
    {ion_lines}
SELECTED_OUTPUT
    -reset false
USER_PUNCH
    -headings aw
    10 PUNCH ACT("H2O")
END
"""

_phreeqc_instance = None


def _phreeqc():
    """Lazily-created shared PhreeqPython instance bound to pitzer.dat."""
    global _phreeqc_instance
    if _phreeqc_instance is None:
        from phreeqpython import PhreeqPython

        _phreeqc_instance = PhreeqPython(database="pitzer.dat")
    return _phreeqc_instance


def _as_molal_floats(molalities):
    """Strip a {ion: molality Quantity} mapping to plain floats in mol/kg."""
    return {ion: q.to_value(_MOLAL) for ion, q in molalities.items()}


def ionic_strength(molal):
    """Molal ionic strength of a brine.

    Lewis-Randall (1921):  I = 1/2 * sum_i m_i * z_i**2

    ``molal`` maps internal ion label -> molality [mol/kg] (plain float).
    Returns I [mol/kg].
    """
    return 0.5 * sum(m * ION_CHARGE[ion] ** 2 for ion, m in molal.items())


def _water_activity_from_phreeqc(amounts, units):
    """Run PHREEQC (pitzer.dat) on a composition and return ACT("H2O").

    ``amounts`` maps internal ion label -> concentration (plain float) on the
    scale named by ``units`` ("mol/kgw" or "mol/l"). Returns 1.0 for a
    salt-free composition (pure water).
    """
    ion_lines = "\n    ".join(
        f"{PITZER_TO_PHREEQC_ELEMENT[ion]} {value}"
        for ion, value in amounts.items()
        if value
    )
    if not ion_lines:
        return 1.0
    deck = _PHREEQC_INPUT_TEMPLATE.format(
        units=units, temp=TEMPERATURE_C, ion_lines=ion_lines
    )
    phreeqc = _phreeqc()
    phreeqc.ip.run_string(deck)
    output = phreeqc.ip.get_selected_output_array()
    return dict(zip(output[0], output[1]))["aw"]


def water_activity(molalities):
    """Solvent water activity a_w of a brine from its ion molalities.

    Delegates to PHREEQC (pitzer.dat) via ``ACT("H2O")``. ``molalities`` maps
    internal ion label (Na, K, Mg, Ca, Cl, SO4) -> molality Quantity [mol/kg].
    Returns a dimensionless float in (0, 1].
    """
    return _water_activity_from_phreeqc(_as_molal_floats(molalities), "mol/kgw")


def osmotic_coefficient(molalities):
    """Molal osmotic coefficient phi of a brine.

    Recovered from the water activity through the exact thermodynamic identity:

        phi = -ln(a_w) / (M_w * sum_i m_i)

    where M_w is the molar mass of water and sum_i m_i the total ion molality.
    The product M_w * sum_i m_i is dimensionless; astropy carries the units and
    resolves the cancellation (including the g <-> kg scale), so it is converted
    to a dimensionless magnitude rather than stripped by hand. ``molalities``
    maps ion label -> molality Quantity [mol/kg].
    """
    a_w = water_activity(molalities)
    total_molality = u.Quantity(list(molalities.values())).sum()
    scale = (M_WATER * total_molality).to_value(u.dimensionless_unscaled)
    return -math.log(a_w) / scale


def pitzer_water_activity_from_molarities(molarities):
    """Pitzer a_w from a mapping of PFLOTRAN ion names to molarity [mol/L].

    Converts to molality (density-aware, ``conversions.py``) before the PHREEQC
    solve. Returns 1.0 when no salt ions are present.
    """
    molalities = molarities_to_molalities(molarities)
    if not molalities:
        return 1.0
    pitzer_molalities = {
        BATCH_ION_TO_PITZER[ion]: quantity
        for ion, quantity in molalities.items()
        if ion in BATCH_ION_TO_PITZER
    }
    if not pitzer_molalities:
        return 1.0
    return float(water_activity(pitzer_molalities))


def pitzer_water_activity_from_batch(batch_row):
    """Pitzer a_w for one incubation batch composition row."""
    molarities = {
        ion: float(batch_row.get(ion, 0.0) or 0.0) for ion in BATCH_ION_TO_PITZER
    }
    return pitzer_water_activity_from_molarities(molarities)
