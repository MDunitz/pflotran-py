"""Concentration-scale conversions for brine compositions.

Molarity (mol/L, per litre of solution) and molality (mol/kg, per kg of water)
differ by the solution density, which for multi-molar brines departs
substantially from 1 kg/L. This module isolates that conversion so the water
activity model can be fed molality on a stated basis.
"""

import astropy.units as u

from .constants import ION_MOLAR_MASS

# Empirical brine-density slope: kg of solution mass added per kg of dissolved
# salt, i.e. rho ~= rho_water + SALT_DENSITY_SLOPE * (salt mass concentration).
# A single salt-identity-blind coefficient (MgCl2 solutions are actually denser
# per gram than NaCl); adequate only for placing these brines on the a_w scale.
SALT_DENSITY_SLOPE = 0.7
WATER_DENSITY = 1000.0 * u.g / u.L


def molarities_to_molalities(molarities):
    """Convert ion molarities [mol/L] to molalities [mol/kg water].

    Estimates solution density from salt mass loading as
        rho = rho_water + SALT_DENSITY_SLOPE * salt_mass_concentration
    then
        m_i = c_i / (rho - salt_mass_concentration)
    where the denominator is the mass of water per litre of solution.

    ``molarities`` maps a PFLOTRAN/batch ion name to a molarity [mol/L]
    (plain float, interpreted as mol/L). Returns {ion: molality Quantity
    [mol/kg]}; ions with zero/absent amount are dropped.
    """
    concentrations = {
        ion: float(molarities.get(ion, 0.0) or 0.0) * u.mol / u.L
        for ion in ION_MOLAR_MASS
    }
    salt_mass_concentration = sum(
        (concentrations[ion] * ION_MOLAR_MASS[ion]).to(u.g / u.L)
        for ion in ION_MOLAR_MASS
    )
    if salt_mass_concentration.value <= 0:
        return {}

    density = WATER_DENSITY + SALT_DENSITY_SLOPE * salt_mass_concentration
    water_mass_concentration = (density - salt_mass_concentration).to(u.kg / u.L)
    if water_mass_concentration.value <= 0:
        raise ValueError(
            f"Implied water mass is non-positive for salt loading "
            f"{salt_mass_concentration}"
        )

    return {
        ion: (concentrations[ion] / water_mass_concentration).to(u.mol / u.kg)
        for ion in ION_MOLAR_MASS
        if concentrations[ion].value > 0
    }
