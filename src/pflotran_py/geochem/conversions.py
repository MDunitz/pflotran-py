"""Concentration-scale conversions for brine compositions.

Molarity (mol/L, per litre of solution) and molality (mol/kg, per kg of water)
differ by the solution density, which for multi-molar brines departs
substantially from 1 kg/L. The conversion therefore requires a measured
solution density; there is no density estimate baked in.
"""

import astropy.units as u

from .constants import ION_MOLAR_MASS


def molarities_to_molalities(molarities, density):
    """Convert ion molarities [mol/L] to molalities [mol/kg water].

        m_i = c_i / (rho - sum_j c_j * M_j)

    where rho is the measured solution density and the denominator is the mass
    of water per litre of solution (total solution mass minus dissolved salt
    mass). ``molarities`` maps a PFLOTRAN/batch ion name to a molarity [mol/L]
    (plain float, interpreted as mol/L); ``density`` is the measured solution
    density [g/mL]. Returns {ion: molality Quantity [mol/kg]}; returns {} for a
    salt-free composition. Raises if the implied water mass is non-positive.
    """
    solution_density = (density * u.g / u.mL).to(u.g / u.L)
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

    water_mass_concentration = (solution_density - salt_mass_concentration).to(
        u.kg / u.L
    )
    if water_mass_concentration.value <= 0:
        raise ValueError(
            f"Implied water mass is non-positive: density {solution_density} "
            f"minus salt {salt_mass_concentration}"
        )

    return {
        ion: (concentrations[ion] / water_mass_concentration).to(u.mol / u.kg)
        for ion in ION_MOLAR_MASS
        if concentrations[ion].value > 0
    }
