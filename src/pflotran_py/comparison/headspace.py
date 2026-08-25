"""Converting between what the model reports and what the instrument measures.

The two sides of this comparison do not speak the same language.

PFLOTRAN reports **aqueous concentrations**: ``Total CH4(aq) [M]`` and
``Total CO2(aq) [M]``, moles per litre of the liquid in the bottle. The gas
chromatograph reports something else entirely: the composition of the
**headspace**, which the measurement pipeline turns into moles of gas above the
liquid. A number from one side cannot be plotted against a number from the other
until one is converted into the other's terms.

This module does that conversion, in the direction model to instrument. A
sealed vial holds a fixed liquid volume and a fixed gas volume, and at
equilibrium a dissolved gas distributes between them according to Henry's law::

    p = c_aq / H_brine                (Henry, rearranged for partial pressure)
    n_gas = p * V_gas / (R * T)       (ideal gas)

Combining these gives a dimensionless partition coefficient, the ratio of gas
concentration to aqueous concentration at equilibrium::

    K = 1 / (H_brine * R * T)

and from it the fraction of a gas that sits in the headspace::

    f_headspace = K * V_gas / (K * V_gas + V_liquid)

For methane at 18 degrees Celsius in fresh water that fraction is about 0.99 --
methane is poorly soluble, so nearly all of it ends up in the headspace. For
carbon dioxide it is about 0.80. Salt pushes both fractions higher, because
dissolved salt reduces gas solubility.

Salting out
-----------
Gases are less soluble in brine than in fresh water, and at the concentrations
here the effect is large rather than marginal: in the strongest magnesium
chloride batch, methane's solubility falls to about a quarter of its
fresh-water value. Ignoring it would misplace the model by more than the effect
being studied. The correction is the Setschenow relation::

    log10(H_water / H_brine) = ks_NaCl * c_NaCl + ks_MgCl2 * c_MgCl2

What this module deliberately does not do
-----------------------------------------
It does not convert *measured* headspace carbon dioxide back into total carbon
produced. That inversion needs carbonate speciation, because at the pH of these
incubations most dissolved inorganic carbon is bicarbonate rather than dissolved
carbon dioxide gas. The carbonate equilibrium constants available for that job
are calibrated for seawater, valid to an ionic strength of roughly 0.7 mol/L,
and these brines run from 1.2 to 5.8 mol/L. Applying them here would be an
extrapolation of up to eight-fold, and the resulting number would look like a
measurement while being something closer to a guess.

Converting in the other direction -- model to headspace -- avoids the problem
entirely, because PFLOTRAN has already done its own carbonate speciation
internally against its own thermodynamic database, and reports the neutral
dissolved carbon dioxide species separately. Henry's law applies to exactly
that species. So the conversion this module performs needs no carbonate
chemistry of its own.

Provenance of the constants
---------------------------
Henry solubilities and their temperature dependence are from Sander (2023),
Compilation of Henry's law constants, Atmos. Chem. Phys. 23, 10901-12440. The
same values appear in ``software_module/diffusion_calculations/brine_gas.py``
in the saltyBiomass repository, which performs the equivalent calculation for
the measurement side. They are duplicated here rather than imported because
this package is installed and tested standalone, including in a container that
has no access to that repository. If either copy is revised, both must move
together; the test suite pins these values so that a silent divergence fails
loudly.

Setschenow coefficients for magnesium chloride are approximate. Data for
divalent salts is sparse, and :func:`validity_warnings` flags conditions where
that approximation is being pushed.
"""

from dataclasses import dataclass

import numpy as np
from astropy import units as u

from ..generator.constants import (
    BOTTLE_TEMPERATURE_C,
    HEADSPACE_VOLUME_L,
    LIQUID_VOLUME_L,
)

# Universal gas constant.
R = 8.314462618 * u.J / (u.mol * u.K)

# Reference temperature for the tabulated Henry constants.
T_STANDARD = 298.15 * u.K


# ═════════════════════════════════════════════════════════════════════
# Gas properties
# ═════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class GasProperties:
    """Henry solubility and salting-out behaviour for one gas.

    Attributes
    ----------
    name : str
        Gas label as used by the measurement pipeline.
    henry_solubility_standard : astropy.units.Quantity
        Henry solubility constant at 298.15 K, in mol per cubic metre per
        pascal. Larger means more soluble.
    vant_hoff_temperature : astropy.units.Quantity
        The quantity ``-d(solH)/R`` in kelvin, giving the temperature
        dependence of the solubility.
    setschenow_nacl : float
        Salting-out coefficient for sodium chloride, litres per mole.
    setschenow_mgcl2 : float
        Salting-out coefficient for magnesium chloride, litres per mole.
        Approximate; divalent salting-out data is sparse.
    """

    name: str
    henry_solubility_standard: u.Quantity
    vant_hoff_temperature: u.Quantity
    setschenow_nacl: float
    setschenow_mgcl2: float


GASES = {
    "CH4": GasProperties(
        name="CH4",
        henry_solubility_standard=1.4e-5 * u.mol / (u.m**3 * u.Pa),
        vant_hoff_temperature=1900 * u.K,
        setschenow_nacl=0.132,
        setschenow_mgcl2=0.30,
    ),
    "CO2": GasProperties(
        name="CO2",
        henry_solubility_standard=3.3e-4 * u.mol / (u.m**3 * u.Pa),
        vant_hoff_temperature=2400 * u.K,
        setschenow_nacl=0.107,
        setschenow_mgcl2=0.24,
    ),
    "N2O": GasProperties(
        name="N2O",
        henry_solubility_standard=2.4e-4 * u.mol / (u.m**3 * u.Pa),
        vant_hoff_temperature=2700 * u.K,
        setschenow_nacl=0.100,
        setschenow_mgcl2=0.22,
    ),
}


# ═════════════════════════════════════════════════════════════════════
# Bottle geometry
# ═════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class BottleGeometry:
    """The liquid and gas volumes of one sealed incubation vial.

    Defaults come from :mod:`pflotran_py.generator.constants` so the
    comparison and the closed-batch deck describe the same bottle. The test
    suite still asserts agreement in case a future edit re-hardcodes volumes.
    """

    liquid_volume: u.Quantity = LIQUID_VOLUME_L * u.L
    headspace_volume: u.Quantity = HEADSPACE_VOLUME_L * u.L

    @property
    def total_volume(self):
        return self.liquid_volume + self.headspace_volume


DEFAULT_BOTTLE = BottleGeometry()

# Incubation temperature. Same source as the closed-batch deck.
DEFAULT_TEMPERATURE = BOTTLE_TEMPERATURE_C * u.deg_C


# ═════════════════════════════════════════════════════════════════════
# Henry's law
# ═════════════════════════════════════════════════════════════════════


def henry_solubility(gas, temperature=DEFAULT_TEMPERATURE):
    """Henry solubility constant in pure water at a given temperature.

    Uses the van't Hoff form::

        H(T) = H(T_std) * exp[ (-dsolH/R) * (1/T - 1/T_std) ]

    Parameters
    ----------
    gas : GasProperties
    temperature : astropy.units.Quantity
        Temperature, in degrees Celsius or kelvin.

    Returns
    -------
    astropy.units.Quantity
        Solubility in mol per cubic metre per pascal.
    """
    temperature_k = temperature.to(u.K, equivalencies=u.temperature())
    exponent = gas.vant_hoff_temperature * (1.0 / temperature_k - 1.0 / T_STANDARD)
    return gas.henry_solubility_standard * np.exp(
        exponent.to_value(u.dimensionless_unscaled)
    )


def salting_out_factor(gas, nacl_molarity=0.0, mgcl2_molarity=0.0):
    """Factor by which brine reduces a gas's solubility, between 0 and 1.

    The Setschenow relation::

        log10(H_water / H_brine) = ks_NaCl * c_NaCl + ks_MgCl2 * c_MgCl2

    Parameters
    ----------
    gas : GasProperties
    nacl_molarity, mgcl2_molarity : float
        Salt concentrations in mol per litre.

    Returns
    -------
    float
        Multiplier on the fresh-water solubility. A value of 0.25 means the gas
        is four times less soluble in this brine than in fresh water.
    """
    log10_reduction = (
        gas.setschenow_nacl * nacl_molarity + gas.setschenow_mgcl2 * mgcl2_molarity
    )
    return 10.0 ** (-log10_reduction)


def henry_solubility_in_brine(
    gas, temperature=DEFAULT_TEMPERATURE, nacl_molarity=0.0, mgcl2_molarity=0.0
):
    """Henry solubility corrected for both temperature and dissolved salt."""
    return henry_solubility(gas, temperature) * salting_out_factor(
        gas, nacl_molarity, mgcl2_molarity
    )


# ═════════════════════════════════════════════════════════════════════
# Partition between liquid and headspace
# ═════════════════════════════════════════════════════════════════════


def gas_water_partition_coefficient(
    gas, temperature=DEFAULT_TEMPERATURE, nacl_molarity=0.0, mgcl2_molarity=0.0
):
    """Ratio of gas-phase to aqueous concentration at equilibrium.

    ``K = 1 / (H_brine * R * T)``, dimensionless. A value of 25 means that at
    equilibrium a cubic metre of headspace holds 25 times the concentration of
    a cubic metre of liquid.
    """
    temperature_k = temperature.to(u.K, equivalencies=u.temperature())
    henry = henry_solubility_in_brine(gas, temperature, nacl_molarity, mgcl2_molarity)
    return (1.0 / (henry * R * temperature_k)).to_value(u.dimensionless_unscaled)


def headspace_fraction(
    gas,
    bottle=DEFAULT_BOTTLE,
    temperature=DEFAULT_TEMPERATURE,
    nacl_molarity=0.0,
    mgcl2_molarity=0.0,
):
    """Fraction of a gas in the bottle that sits in the headspace.

    ``f = K*V_gas / (K*V_gas + V_liquid)``.

    This is the single number that says how much of what the organisms produced
    the instrument can actually see. For methane it is close to one, so the
    headspace measurement is nearly the whole story. For carbon dioxide it is
    lower, so a meaningful share of the produced gas stays dissolved and never
    reaches the detector.

    Returns
    -------
    float
        Between 0 and 1.
    """
    partition = gas_water_partition_coefficient(
        gas, temperature, nacl_molarity, mgcl2_molarity
    )
    gas_capacity = partition * bottle.headspace_volume.to_value(u.L)
    liquid_capacity = bottle.liquid_volume.to_value(u.L)
    return gas_capacity / (gas_capacity + liquid_capacity)


def aqueous_concentration_to_headspace_moles(
    aqueous_concentration,
    gas,
    bottle=DEFAULT_BOTTLE,
    temperature=DEFAULT_TEMPERATURE,
    nacl_molarity=0.0,
    mgcl2_molarity=0.0,
):
    """Headspace moles in equilibrium with a given aqueous concentration.

    This is the conversion the comparison rests on: it takes a concentration
    PFLOTRAN reports and returns the number of moles a gas chromatograph
    sampling the headspace would find.

    ``n_gas = (c_aq / H_brine) * V_gas / (R*T)``, which reduces to
    ``n_gas = c_aq * K * V_gas``.

    Parameters
    ----------
    aqueous_concentration : astropy.units.Quantity or float
        Dissolved concentration in mol/L. A bare number is taken as mol/L.
    gas : GasProperties

    Returns
    -------
    astropy.units.Quantity
        Moles of gas in the headspace.
    """
    if not isinstance(aqueous_concentration, u.Quantity):
        aqueous_concentration = aqueous_concentration * u.mol / u.L

    partition = gas_water_partition_coefficient(
        gas, temperature, nacl_molarity, mgcl2_molarity
    )
    return (
        aqueous_concentration.to(u.mol / u.L) * partition * bottle.headspace_volume
    ).to(u.mol)


def aqueous_concentration_to_dissolved_moles(
    aqueous_concentration, bottle=DEFAULT_BOTTLE
):
    """Moles still dissolved in the liquid, ``c_aq * V_liquid``."""
    if not isinstance(aqueous_concentration, u.Quantity):
        aqueous_concentration = aqueous_concentration * u.mol / u.L
    return (aqueous_concentration.to(u.mol / u.L) * bottle.liquid_volume).to(u.mol)


def gas_phase_moles(gas_concentration, bottle=DEFAULT_BOTTLE):
    """Headspace moles read straight from a model that carries a gas phase.

    The preferred route. When the deck declares the gas as an active gas
    species, PFLOTRAN performs the partition itself using its own equilibrium
    constants, and reports the gas-phase concentration in moles per cubic metre
    of gas. Multiplying by the headspace volume gives the moles a gas
    chromatograph would sample, with no assumption of ours in between.

    Parameters
    ----------
    gas_concentration : astropy.units.Quantity or float
        Gas-phase concentration, mol per cubic metre of gas. A bare number is
        read in those units.
    """
    if not isinstance(gas_concentration, u.Quantity):
        gas_concentration = gas_concentration * u.mol / u.m**3
    return (gas_concentration.to(u.mol / u.m**3) * bottle.headspace_volume).to(u.mol)


def partition_total_moles(
    total_moles,
    gas,
    bottle=DEFAULT_BOTTLE,
    temperature=DEFAULT_TEMPERATURE,
    nacl_molarity=0.0,
    mgcl2_molarity=0.0,
):
    """Split a known total quantity of gas between headspace and liquid.

    For decks that do **not** carry a gas phase, where the model reports the
    whole inventory of a gas as dissolved because it has nowhere else to put it.

    This function takes that total and divides it, rather than adding a
    headspace on top of it. The distinction is not cosmetic. Multiplying a
    dissolved concentration by the partition coefficient and the headspace
    volume answers the question "if this liquid were equilibrated against a
    headspace, what would the headspace hold?" -- and the answer is larger than
    the methane the model ever made, by a factor of roughly
    ``(K * V_gas + V_liquid) / V_liquid``, about a hundredfold for methane in
    this vial. Doing that and calling the result a prediction inflates the
    model by that factor.

    Parameters
    ----------
    total_moles : astropy.units.Quantity or float
        All of the gas in the bottle, however the model distributes it.

    Returns
    -------
    astropy.units.Quantity
        Moles in the headspace.
    """
    if not isinstance(total_moles, u.Quantity):
        total_moles = total_moles * u.mol
    fraction = headspace_fraction(
        gas,
        bottle=bottle,
        temperature=temperature,
        nacl_molarity=nacl_molarity,
        mgcl2_molarity=mgcl2_molarity,
    )
    return (total_moles * fraction).to(u.mol)


def model_headspace_moles(
    aqueous_concentration,
    gas,
    degassed_concentration=None,
    gas_concentration=None,
    bottle=DEFAULT_BOTTLE,
    temperature=DEFAULT_TEMPERATURE,
    nacl_molarity=0.0,
    mgcl2_molarity=0.0,
):
    """Headspace moles predicted by the model, for one gas at one time.

    Takes whichever route the deck supports.

    When ``gas_concentration`` is given, the deck carries a real gas phase and
    the answer is read from it directly.

    Otherwise the deck has no gas phase, and every mole of the gas is somewhere
    in the liquid: dissolved, plus -- for methane in decks that use the
    ebullition proxy -- whatever the proxy has moved into ``Tracer2``. Those are
    summed into a total and then *partitioned*, so the result can never exceed
    what the model actually produced.

    Parameters
    ----------
    aqueous_concentration : astropy.units.Quantity or float
        The model's dissolved concentration, mol/L.
    degassed_concentration : astropy.units.Quantity or float, optional
        The model's ``Tracer2`` concentration, mol/L, for ebullition-proxy
        decks. Omit where there is no proxy.
    gas_concentration : astropy.units.Quantity or float, optional
        Gas-phase concentration, mol per cubic metre of gas. When present, the
        other two arguments are ignored.
    """
    if gas_concentration is not None:
        return gas_phase_moles(gas_concentration, bottle=bottle)

    if not isinstance(aqueous_concentration, u.Quantity):
        aqueous_concentration = aqueous_concentration * u.mol / u.L
    total = (aqueous_concentration.to(u.mol / u.L) * bottle.liquid_volume).to(u.mol)

    if degassed_concentration is not None:
        if not isinstance(degassed_concentration, u.Quantity):
            degassed_concentration = degassed_concentration * u.mol / u.L
        total = total + (
            degassed_concentration.to(u.mol / u.L) * bottle.liquid_volume
        ).to(u.mol)

    return partition_total_moles(
        total,
        gas,
        bottle=bottle,
        temperature=temperature,
        nacl_molarity=nacl_molarity,
        mgcl2_molarity=mgcl2_molarity,
    )


def headspace_moles_to_total_moles(
    headspace_moles,
    gas,
    bottle=DEFAULT_BOTTLE,
    temperature=DEFAULT_TEMPERATURE,
    nacl_molarity=0.0,
    mgcl2_molarity=0.0,
):
    """Total moles in the bottle implied by a headspace measurement.

    Divides the measured headspace moles by the headspace fraction, which is
    the share of the total the instrument can see.

    Use with care, and read what it assumes. It supposes the liquid and the
    headspace have reached Henry equilibrium, which a static unstirred vial
    may not have done between samplings. For methane the correction is small,
    a few percent, because almost all of it is in the headspace anyway. For
    carbon dioxide the *Henry* part of the correction is larger, and it still
    understates the total, because dissolved inorganic carbon that has
    speciated to bicarbonate is not counted here at all -- see the module
    docstring for why that part is deliberately left alone.
    """
    fraction = headspace_fraction(
        gas,
        bottle=bottle,
        temperature=temperature,
        nacl_molarity=nacl_molarity,
        mgcl2_molarity=mgcl2_molarity,
    )
    return headspace_moles / fraction


# ═════════════════════════════════════════════════════════════════════
# Mapping a measured brine onto the salting-out model
# ═════════════════════════════════════════════════════════════════════


def setschenow_salts_from_composition(batch_row):
    """Reduce a measured ion composition to the two salts Setschenow covers.

    Setschenow coefficients are tabulated per salt, and reliable values exist
    here only for sodium chloride and magnesium chloride. The measured brines
    also contain sulfate, calcium and potassium, chiefly in the sea-salt
    batches.

    The mapping takes sodium molarity as the sodium chloride term and magnesium
    molarity as the magnesium chloride term, which captures the dominant
    salting-out in every batch. The remaining ions are roughly eight percent of
    the ionic strength in the sea-salt batches and absent from the rest, so
    omitting them understates the salting-out slightly. Understating it is the
    conservative direction: it places less gas in the headspace, so the model is
    not flattered.

    Parameters
    ----------
    batch_row : mapping
        A row of the batch composition table, carrying ``Na+`` and ``Mg++``.

    Returns
    -------
    tuple of float
        ``(nacl_molarity, mgcl2_molarity)`` in mol/L.
    """
    return float(batch_row.get("Na+", 0.0)), float(batch_row.get("Mg++", 0.0))


def validity_warnings(nacl_molarity=0.0, mgcl2_molarity=0.0):
    """Conditions under which this conversion is being extrapolated.

    None of these invalidate a result on their own. They mark where the
    uncertainty stops being small, so that a figure can say so rather than
    presenting every point with equal confidence.

    Returns
    -------
    list of str
        Human-readable warnings; empty when the conversion is on solid ground.
    """
    warnings = []
    ionic_strength = nacl_molarity + 3.0 * mgcl2_molarity

    if ionic_strength > 5.0:
        warnings.append(
            f"Ionic strength {ionic_strength:.1f} mol/L is beyond where the "
            f"Setschenow coefficients were fitted; the salting-out correction "
            f"is an extrapolation."
        )
    if mgcl2_molarity > 2.0:
        warnings.append(
            f"Magnesium chloride at {mgcl2_molarity:.1f} mol/L relies on a "
            f"divalent Setschenow coefficient that is itself approximate."
        )
    return warnings
