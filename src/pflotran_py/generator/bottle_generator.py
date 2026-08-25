"""PFLOTRAN input decks for sealed-bottle (closed batch) incubations.

The generator in ``pflotran_generator.py`` builds decks for an open sediment
column: a one-cubic-metre domain with an atmospheric boundary at the top, so
water and solutes move through the domain over the course of the simulation.
That is the right model for coastal sediment, and it is the wrong model for the
laboratory incubations this project measures.

The incubations are sealed 125 mL glass vials holding roughly 25 mL of brine
plus biomass and roughly 100 mL of headspace. Nothing enters and nothing leaves
between sampling. Running an open-column deck against those measurements gives
a badly misleading answer: in the committed five-day test run, chloride falls
from 6.2 M to 7.6e-14 M because the atmospheric boundary flushes the entire
brine out of the domain, which removes exactly the salt stress the experiment
is designed to impose.

This module therefore builds a *closed* deck. It changes four things relative to
the sediment-column generator and deliberately changes nothing else:

1.  **The domain is sealed.** All boundary conditions are removed. PFLOTRAN
    treats a face with no boundary condition as zero-flux, so a domain with only
    an initial condition is a closed batch reactor. This is the change that stops
    the salt washing out.
2.  **The domain is the size of a vial.** A single cell of 1.25e-4 cubic metres
    (a 50 mm cube) rather than one cubic metre, with a gas saturation chosen to
    reproduce the real headspace-to-liquid ratio.
3.  **The run is long enough to compare against.** A hundred and thirty days by
    default, which spans the measured record (119 days for Exp003, 122 for
    Exp004), rather than the 31-day column default.
4.  **The temperature is the incubation temperature.** 18 degrees Celsius, which
    is what the post-processing package in ``config.py`` already assumes, rather
    than the 8 degrees Celsius used for coastal sediment.

Two further changes concern the chemistry rather than the geometry, and are
described where they are implemented: ``couple_carbonate`` lets dissolved carbon
dioxide equilibrate with the carbonate system instead of being frozen as an
independent primary species, and ``methane_gas_phase`` gives the bottle a real
methane headspace in place of the ebullition proxy. Both default on here and
off in the column generator, so existing column results are unaffected.

Everything that constitutes the chemistry -- the microbial reaction network, the
three water-activity inhibition sandboxes, the species and mineral lists, the
mineral kinetics -- is inherited unchanged from ``PFLOTRANGenerator``. This
module subclasses rather than copies so that a change to the reaction network
applies to column and bottle decks alike.

A consequence worth stating plainly: a single-cell domain has no spatial
extent, so every concentration gradient in it is exactly zero, and the
gradient and diffusive-flux stages of the post-processing pipeline produce
nothing meaningful for these runs. That is the correct behaviour for a bottle.
What a bottle gives you is a concentration time series, and that is what should
be compared against the measured headspace data.

Usage::

    from pflotran_py.generator.bottle_generator import BottleGenerator

    gen = BottleGenerator(brine=nacl_brine(molality=3.0))
    gen.generate("bottle_aw0.88.in")

or, to build one deck per measured water activity::

    generate_bottle_series([0.996, 0.905, 0.773, 0.696, 0.609, 0.454],
                           output_dir="decks/")

Sources:
    [1] Pitzer & Mayorga (1973) J. Phys. Chem. 77(19), 2300-2308 -- osmotic
        coefficients for 1:1 electrolytes, used to convert a target water
        activity into a NaCl molality.
    [2] Millero (2013) Chemical Oceanography, 4th ed. -- seawater ion ratios.
    See REFERENCES.md for the reaction-network parameter sourcing, which this
    module inherits unchanged.
"""

import math
import os

from .constants import (
    AW_CRIT_HYDROGENOTROPHIC,
    BOTTLE_FINAL_TIME_DAYS,
    BOTTLE_GAS_PRESSURE_PA,
    BOTTLE_GAS_SATURATION,
    BOTTLE_POROSITY,
    BOTTLE_TEMPERATURE_C,
    HEADSPACE_VOLUME_L,
    LIQUID_VOLUME_L,
    VIAL_VOLUME_L,
)
from .pflotran_generator import PFLOTRANGenerator

# Bottle geometry / run defaults live in ``constants`` and are re-exported
# here so ``from ...bottle_generator import VIAL_VOLUME_L`` keeps working.

# PFLOTRAN works in metres. A cube of this edge length has the vial's volume.
# 0.125 L = 1.25e-4 m^3, and (1.25e-4)^(1/3) = 0.05 m exactly.
_VIAL_EDGE_M = (VIAL_VOLUME_L * 1e-3) ** (1.0 / 3.0)

# ═════════════════════════════════════════════════════════════════════
# Water activity <-> NaCl molality
# ═════════════════════════════════════════════════════════════════════
#
# The AWINHIBIT sandboxes do not take water activity as an input. PFLOTRAN
# computes water activity from the solution composition at each timestep (the
# deck switches this on with ACTIVITY_WATER). So to build a deck that sits at a
# particular water activity, we have to work backwards from the measured water
# activity to a salt concentration that produces it.
#
# The conversion below is a starting estimate, not the answer. It assumes a pure
# NaCl solution at 25 C, whereas the real batches contain NaCl, MgCl2 and sea
# salt at 18 C. After a run completes, the water activity PFLOTRAN actually
# computed should be read back out of the output and compared against the
# measured value -- that round trip, not this function, is what establishes
# whether a deck sits where it was meant to sit.

_PITZER_A_PHI = 0.3915  # Debye-Huckel osmotic coefficient, 25 C [kg^0.5/mol^0.5]
_PITZER_B = 1.2  # universal Pitzer constant [kg^0.5/mol^0.5]
_PITZER_ALPHA = 2.0  # universal for 1:1 electrolytes [kg^0.5/mol^0.5]
_NACL_BETA0 = 0.0765  # Pitzer & Mayorga (1973), Table I
_NACL_BETA1 = 0.2664
_NACL_CPHI = 0.00127

_M_WATER_G_PER_MOL = 18.0153
_NACL_NU = 2  # ions per formula unit: Na+ and Cl-


def nacl_osmotic_coefficient(molality):
    """Pitzer osmotic coefficient for aqueous NaCl at 25 C.

    Pitzer & Mayorga (1973), equation for a 1:1 electrolyte:

        phi - 1 = -A_phi * sqrt(I) / (1 + b*sqrt(I))
                  + m * (beta0 + beta1 * exp(-alpha * sqrt(I)))
                  + m^2 * C_phi

    For a 1:1 electrolyte the ionic strength I equals the molality m.

    Parameters
    ----------
    molality : float
        NaCl molality [mol per kg water].

    Returns
    -------
    float
        Osmotic coefficient, dimensionless.
    """
    ionic_strength = molality
    sqrt_i = math.sqrt(ionic_strength)

    debye_huckel = -_PITZER_A_PHI * sqrt_i / (1.0 + _PITZER_B * sqrt_i)
    second_virial = molality * (
        _NACL_BETA0 + _NACL_BETA1 * math.exp(-_PITZER_ALPHA * sqrt_i)
    )
    third_virial = molality**2 * _NACL_CPHI

    return 1.0 + debye_huckel + second_virial + third_virial


def water_activity_from_nacl_molality(molality):
    """Water activity of an aqueous NaCl solution at 25 C.

    From the definition of the osmotic coefficient:

        ln(a_w) = -nu * m * phi * M_water / 1000

    where nu = 2 for NaCl, m is molality [mol/kg], phi is the osmotic
    coefficient and M_water is 18.0153 g/mol.

    Parameters
    ----------
    molality : float
        NaCl molality [mol per kg water]. Zero returns unit activity.

    Returns
    -------
    float
        Water activity, dimensionless, between 0 and 1.
    """
    if molality <= 0:
        return 1.0
    phi = nacl_osmotic_coefficient(molality)
    ln_aw = -_NACL_NU * molality * phi * _M_WATER_G_PER_MOL / 1000.0
    return math.exp(ln_aw)


def nacl_molality_for_water_activity(target_aw, tolerance=1e-6, max_iterations=200):
    """Invert :func:`water_activity_from_nacl_molality` by bisection.

    Water activity falls monotonically with molality, so bisection is
    sufficient and needs no derivative. The upper bracket is 12 mol/kg, well
    past NaCl saturation (about 6.1 mol/kg at 25 C), so that decks below the
    solubility limit are reachable and ones above it are visibly so.

    Parameters
    ----------
    target_aw : float
        Water activity to hit, between 0 and 1.
    tolerance : float
        Absolute convergence tolerance on water activity.
    max_iterations : int
        Bisection iteration cap.

    Returns
    -------
    float
        NaCl molality [mol per kg water] whose water activity is ``target_aw``.

    Notes
    -----
    A returned molality above roughly 6.1 mol/kg exceeds NaCl solubility at
    25 C. Such a solution cannot be made from NaCl alone -- the real batch at
    that water activity used a different salt or a mixture -- and the deck it
    generates should be treated as a thought experiment until its composition
    is replaced with the measured one.
    """
    low, high = 0.0, 12.0
    for _ in range(max_iterations):
        mid = 0.5 * (low + high)
        aw = water_activity_from_nacl_molality(mid)
        if abs(aw - target_aw) < tolerance:
            return mid
        if aw > target_aw:
            low = mid  # not salty enough yet
        else:
            high = mid
    return 0.5 * (low + high)


def nacl_brine(molality=None, water_activity=None):
    """Build a PFLOTRAN concentration override for a NaCl brine.

    Give either a molality or a target water activity, not both.

    Molality (mol per kg water) is converted to molarity (mol per litre
    solution) using the density of aqueous NaCl, approximated as
    ``rho = 1000 + 40 * m`` kg/m^3, which tracks tabulated NaCl densities to
    better than one percent up to saturation. That approximation is adequate
    here because the deck's chloride is charge-balanced rather than fixed, so
    an error in the sodium molarity is partly absorbed.

    Parameters
    ----------
    molality : float, optional
        NaCl molality [mol per kg water].
    water_activity : float, optional
        Target water activity; converted to molality first.

    Returns
    -------
    dict
        Mapping of species name to PFLOTRAN constraint string, suitable for
        the ``concentrations`` or ``brine`` argument of :class:`BottleGenerator`.
        Sodium is set explicitly; chloride carries the ``Z`` charge-balance
        code so PFLOTRAN closes the charge balance itself.
    """
    if (molality is None) == (water_activity is None):
        raise ValueError("Pass exactly one of molality or water_activity")

    if molality is None:
        molality = nacl_molality_for_water_activity(water_activity)

    # Take one kilogram of water as the basis. Adding m moles of NaCl at
    # 58.44 g/mol gives a solution of known mass; dividing by its density gives
    # its volume; the salt's molarity is then m moles in that volume.
    solution_density_kg_per_l = (1000.0 + 40.0 * molality) / 1000.0
    solution_mass_kg = 1.0 + molality * 58.44 / 1000.0
    solution_volume_l = solution_mass_kg / solution_density_kg_per_l
    molarity = molality / solution_volume_l

    return {
        "Na+": f"{molarity:.4e} T",
        "Cl-": f"{molarity:.4e} Z",
    }


# ═════════════════════════════════════════════════════════════════════
# Generator
# ═════════════════════════════════════════════════════════════════════


def _repo_root():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(package_dir, "..", "..", ".."))


def _default_database_path():
    """Absolute path to the repository's own thermodynamic database.

    The sediment-column generator defaults to a developer-specific absolute
    path that exists on one machine. Resolving against this file's location
    makes a generated deck runnable from any clone.
    """
    return os.path.join(_repo_root(), "sandbox", "hanford.dat")


# ═════════════════════════════════════════════════════════════════════
# Methane gas phase
# ═════════════════════════════════════════════════════════════════════
#
# A sealed vial has a real headspace, and methane moves into it continuously by
# Henry's law. The inherited reaction network models that departure differently:
# with an ebullition proxy, a kinetic reaction CH4(aq) -> Tracer2 that fires only
# once dissolved methane passes 2.5e-3 M, standing in for bubbles nucleating and
# rising away. That is a reasonable picture of submerged sediment, where methane
# must become buoyant enough to escape. It is the wrong picture of a bottle,
# where there is no threshold to cross: the headspace is already there, and
# methane partitions into it from the first molecule produced.
#
# Replacing the proxy with a genuine gas phase requires a gas-phase methane
# species, and the database does not quite provide one. hanford.dat defines
# CH4(g) against Methane(aq), an organic species built on an ethane basis that
# this reaction network does not carry. The deck's methane is CH4(aq), a redox
# species built on bicarbonate and oxygen. The two are the same molecule -- both
# are recorded with molar mass 16.0428 -- but they are written against different
# basis species, so the tabulated gas reaction cannot be used as it stands.
#
# What can be reused is the equilibrium constant. The tabulated log K series for
# CH4(g) is the methane Henry constant: at 25 degrees Celsius it gives
# 1.41e-3 mol/(L*atm), against Sander (2023)'s 1.4e-3, and at 0 degrees the two
# agree to two percent. Volatility is a property of the molecule, not of the
# basis chosen to describe its dissolved form, so the same constants apply when
# the aqueous partner is written as CH4(aq).
#
# The patched database therefore rewrites one line: CH4(g)'s aqueous partner,
# leaving its constants untouched. It is written to a separate file rather than
# edited in place, so the shared database is never modified and the change is
# visible as a diff between two files.

METHANE_GAS_SPECIES = "CH4(g)"
_HANFORD_METHANE_GAS_LINE = "'CH4(g)' 0.0000 1 1.0000 'Methane(aq)'"
_BOTTLE_METHANE_GAS_LINE = "'CH4(g)' 0.0000 1 1.0000 'CH4(aq)'"

# A mineral record in this database reads
#
#     'name' molar_volume n_species [coefficient 'species'] x n  logK x 8  molar_mass
#
# The Cellulose_min record does not. It declares two species but supplies two
# more fields than that implies, and the second species carries a coefficient of
# zero, which PFLOTRAN drops. The parser then finds one species where the header
# promised two and stops with "Number of reaction species does not match
# original: 1 2".
#
# The repair declares the one species the record actually describes and removes
# the two surplus fields. The chemistry is untouched: Cellulose_min still
# dissolves to one DOM1 with a log K of zero at every tabulated temperature.
# This record is unused by the existing decks, which is presumably why the
# malformation has gone unnoticed.
_HANFORD_CELLULOSE_LINE = (
    "'Cellulose_min' 162.14 2 1.0000 'DOM1' 0.0000 'H+' "
    "0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0 0.0 162.1400"
)
_BOTTLE_CELLULOSE_LINE = (
    "'Cellulose_min' 162.14 1 1.0000 'DOM1' "
    "0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 162.1400"
)


def bottle_database_path():
    """Path to the database used by closed-batch decks."""
    return os.path.join(_repo_root(), "sandbox", "hanford_bottle.dat")


def write_bottle_database(source_path=None, destination_path=None):
    """Write a database whose gas-phase methane pairs with ``CH4(aq)``.

    Copies the thermodynamic database and rewrites the single ``CH4(g)`` line so
    that its aqueous partner is the species this reaction network actually
    carries. Equilibrium constants are copied unchanged; see the commentary
    above for why that is sound.

    Returns
    -------
    str
        Path to the written database.

    Raises
    ------
    ValueError
        If the expected ``CH4(g)`` line is absent, rather than writing a
        database that silently lacks a methane gas phase.
    """
    source_path = source_path or _default_database_path()
    destination_path = destination_path or bottle_database_path()

    with open(source_path) as handle:
        text = handle.read()

    if _HANFORD_METHANE_GAS_LINE not in text and _BOTTLE_METHANE_GAS_LINE not in text:
        raise ValueError(
            f"No CH4(g) entry pairing with Methane(aq) found in {source_path}; "
            "the database format may have changed. Refusing to write a bottle "
            "database with no methane gas phase."
        )

    text = text.replace(_HANFORD_METHANE_GAS_LINE, _BOTTLE_METHANE_GAS_LINE)
    text = text.replace(_HANFORD_CELLULOSE_LINE, _BOTTLE_CELLULOSE_LINE)
    with open(destination_path, "w") as handle:
        handle.write(text)
    return destination_path


class BottleGenerator(PFLOTRANGenerator):
    """Generate a closed-batch PFLOTRAN deck for one sealed incubation vial.

    Inherits the full reaction network, sandbox configuration and species lists
    from :class:`PFLOTRANGenerator`. Overrides only the physical setup: domain
    size, closure, run length and temperature.

    Parameters
    ----------
    brine : dict, optional
        Species-to-constraint-string overrides describing the salt content of
        this bottle, e.g. the output of :func:`nacl_brine`. Merged on top of the
        inherited default initial concentrations, so it need only name the ions
        that differ.
    gas_saturation : float
        Fraction of pore volume that is headspace. Defaults to the real vial
        ratio, 100 mL of 125 mL.
    gas_pressure_pa : float
        Initial headspace pressure [Pa].
    porosity : float
        Porosity of the bottle contents. Near 1 because a vial is fluid, not
        packed sediment.
    label : str, optional
        Free-text description written into the deck header, so a directory of
        generated decks is readable without opening each one.

    Other parameters are inherited; see :class:`PFLOTRANGenerator`. The
    defaults for ``temperature``, ``final_time_days`` and ``database_path``
    differ, as described in the module docstring.
    """

    def __init__(
        self,
        brine=None,
        gas_saturation=BOTTLE_GAS_SATURATION,
        gas_pressure_pa=BOTTLE_GAS_PRESSURE_PA,
        porosity=BOTTLE_POROSITY,
        label=None,
        temperature=BOTTLE_TEMPERATURE_C,
        final_time_days=BOTTLE_FINAL_TIME_DAYS,
        database_path=None,
        concentrations=None,
        couple_carbonate=True,
        methane_gas_phase=True,
        **kwargs,
    ):
        merged_concentrations = dict(concentrations or {})
        if brine:
            merged_concentrations.update(brine)

        self.methane_gas_phase = methane_gas_phase

        if database_path is None:
            database_path = (
                write_bottle_database()
                if methane_gas_phase
                else _default_database_path()
            )

        super().__init__(
            concentrations=merged_concentrations,
            temperature=temperature,
            final_time_days=final_time_days,
            database_path=database_path,
            couple_carbonate=couple_carbonate,
            **kwargs,
        )

        self.gas_saturation = gas_saturation
        self.gas_pressure_pa = gas_pressure_pa
        self.porosity = porosity
        self.label = label

        # A bottle has no spatial extent worth resolving, so the inherited grid
        # preset is replaced wholesale with a single cell the size of the vial.
        edge = _VIAL_EDGE_M
        self.grid_cells = "1 1 1"
        self.cell_size_x = f"{edge:.4f}d0"
        self.cell_size_y = f"{edge:.4f}d0"
        self.cell_size_z = f"{edge:.4f}d0"
        self.domain = (edge, edge, edge)

    # ─────────────────────────────────────────────────────────────────
    # Overrides
    # ─────────────────────────────────────────────────────────────────

    def _include_reaction(self, rxn):
        """Drop the ebullition proxy when a real gas phase is present.

        With ``CH4(g)`` declared, PFLOTRAN moves methane between the liquid and
        the headspace itself, using the equilibrium constant from the database.
        Keeping the proxy alongside it would remove the same methane twice --
        once into the gas phase and once into ``Tracer2`` -- and would impose a
        2.5e-3 M threshold on a process that in a sealed vial has none.
        """
        if self.methane_gas_phase and rxn.get("rate_key") == "ebullition":
            return False
        return super()._include_reaction(rxn)

    def _build_mineral_kinetics_and_sorption(self):
        """Add gas-phase methane to the active gas species.

        Active rather than passive: a passive gas species is held at a fixed
        partial pressure, which is what an open system in contact with the
        atmosphere looks like. In a sealed bottle the methane partial pressure
        is free to rise as methane accumulates, which is an active species.
        """
        block = super()._build_mineral_kinetics_and_sorption()
        if not self.methane_gas_phase:
            return block
        return block.replace(
            "    ACTIVE_GAS_SPECIES\n      GAS_TRANSPORT_IS_UNVETTED\n      CO2(g)\n",
            "    ACTIVE_GAS_SPECIES\n      GAS_TRANSPORT_IS_UNVETTED\n      CO2(g)\n"
            f"      {METHANE_GAS_SPECIES}\n",
        )

    def _build_chemistry_output(self):
        """Report gas-phase methane, so the headspace can be read directly.

        Without this the simulation partitions methane into the headspace but
        never writes how much went there, and the comparison would have to infer
        it from the dissolved concentration -- which is the very inference the
        gas phase was added to avoid.
        """
        block = super()._build_chemistry_output()
        if not self.methane_gas_phase:
            return block
        return block.replace(
            "  OUTPUT\n    PH\n",
            f"  OUTPUT\n    PH\n    GAS_CONCENTRATION\n    {METHANE_GAS_SPECIES}\n",
        )

    def _build_constraints(self):
        """Initial constraint only -- a sealed bottle has no boundary.

        The inherited version also emits a ``CONSTRAINT atmospheric`` block for
        the open column's top boundary. A closed deck has no boundary condition
        to attach it to, so emitting it would leave a dangling block that
        suggests an exchange with the atmosphere that does not happen.
        """
        full = super()._build_constraints()
        atmospheric_start = full.find("\nCONSTRAINT atmospheric")
        if atmospheric_start == -1:
            return full
        return full[:atmospheric_start].rstrip()

    def _build_material_properties(self):
        """Material property block for a fluid-filled vial.

        Written out here rather than inherited from the ``MATERIAL_PROPERTIES``
        template because the template describes packed sediment: 0.97 porosity
        and a permeability tensor tuned for a column. Permeability is set very
        low, which costs nothing in a sealed domain with no pressure gradient to
        drive flow, and guards against a spurious circulation if one is ever
        introduced.
        """
        return f"""\
#=========================== material properties ==============================
MATERIAL_PROPERTY bottle
  ID 1
  POROSITY {self.porosity:.3f}d0
  SOIL_COMPRESSIBILITY 1.d-07
  SOIL_REFERENCE_PRESSURE {self.gas_pressure_pa:.5e}
  ROCK_DENSITY 2650.0d0
  SPECIFIC_HEAT 830.0d0
  THERMAL_CONDUCTIVITY_DRY 0.12037926674717922d0
  THERMAL_CONDUCTIVITY_WET 1.6082691464310437d0
  CHARACTERISTIC_CURVES sf01
  PERMEABILITY
    PERM_X 1.d-18
    PERM_Y 1.d-18
    PERM_Z 1.d-18
  /
/

CHARACTERISTIC_CURVES sf01
  SATURATION_FUNCTION VAN_GENUCHTEN
    LIQUID_RESIDUAL_SATURATION 0.d0
    ALPHA 1.d-4
    M 0.5d0
    MAX_CAPILLARY_PRESSURE 1.d6
  /
  PERMEABILITY_FUNCTION MUALEM
    PHASE LIQUID
    LIQUID_RESIDUAL_SATURATION 0.d0
    M 0.5d0
  /
  PERMEABILITY_FUNCTION MUALEM_VG_GAS
    PHASE GAS
    LIQUID_RESIDUAL_SATURATION 0.d0
    GAS_RESIDUAL_SATURATION 1.d-5
    M 0.5d0
  /
END"""

    def _build_regions_and_conditions(self):
        """Regions and conditions for a sealed domain.

        This is the override that closes the bottle. Relative to the inherited
        open-column version it drops:

        * ``FLOW_CONDITION atmospheric`` and ``TRANSPORT_CONDITION atmospheric``
        * ``BOUNDARY_CONDITION top_atm``
        * the ``top_surface`` and ``bottom_surface`` regions those used

        leaving an initial condition applied over the whole domain and nothing
        else. PFLOTRAN treats any face without a boundary condition as
        zero-flux, so the result is a closed batch reactor: no water enters, no
        water leaves, and the brine stays at the concentration it was set to.
        """
        lx, ly, lz = self.domain
        return f"""\
#=========================== regions ==========================================
# A sealed bottle needs exactly one region. There are deliberately no surface
# regions here: a region only matters if something is coupled to it, and
# coupling anything to a face would open the bottle.
REGION all
  COORDINATES
    0.d0 0.d0 0.d0
    {lx:.4f}d0 {ly:.4f}d0 {lz:.4f}d0
  /
END

REGION center_obs
  COORDINATE {lx / 2:.4f}d0 {ly / 2:.4f}d0 {lz / 2:.4f}d0
END

#=========================== observation points ===============================
# With a single cell the observation point is the whole bottle, so the
# observation file is a complete record of the simulation as a time series.
OBSERVATION
  REGION center_obs
END

#=========================== flow and transport conditions ====================
# One condition, used as the initial state. No atmospheric condition is defined
# because there is no boundary for it to act on.
FLOW_CONDITION initial
  TYPE
    GAS_PRESSURE dirichlet
    GAS_SATURATION dirichlet
    TEMPERATURE dirichlet
  /
  GAS_PRESSURE {self.gas_pressure_pa:.5e}
  GAS_SATURATION {self.gas_saturation:.3f}
  TEMPERATURE {self.temperature:.1f}d0
/

TRANSPORT_CONDITION initial
  TYPE dirichlet
  CONSTRAINT_LIST
    0.d0 initial
  /
END

#=========================== condition couplers ===============================
# The initial condition is the only coupler. Adding a BOUNDARY_CONDITION here
# would reopen the domain and let the brine drain out, which is the failure
# mode this deck exists to avoid.
INITIAL_CONDITION
  TRANSPORT_CONDITION initial
  FLOW_CONDITION initial
  REGION all
END

#=========================== stratigraphy couplers ============================
STRATA
  REGION all
  MATERIAL bottle
END

END_SUBSURFACE"""

    def generate(self, filename="bottle_simulation.in"):
        """Write the deck, then substitute the bottle material properties.

        The inherited ``generate`` assembles sections from templates, one of
        which is the sediment ``MATERIAL_PROPERTIES`` block. Rather than
        reimplement the whole assembly to swap one section, the file is written
        by the parent and the block replaced afterwards. The marker comments
        that delimit the block are stable parts of the template.
        """
        path = super().generate(filename)

        with open(path) as handle:
            content = handle.read()

        material_start = content.find(
            "#=========================== material properties"
        )
        material_end = content.find("#=========================== output options")
        if material_start != -1 and material_end != -1:
            content = (
                content[:material_start]
                + self._build_material_properties()
                + "\n\n"
                + content[material_end:]
            )

        header_note = (
            "# Closed-batch (sealed bottle) deck -- no boundary conditions.\n"
            f"# Vial {VIAL_VOLUME_L * 1000:.0f} mL, headspace "
            f"{HEADSPACE_VOLUME_L * 1000:.0f} mL, liquid "
            f"{LIQUID_VOLUME_L * 1000:.0f} mL.\n"
        )
        if self.label:
            header_note += f"# Condition: {self.label}\n"
        content = header_note + content

        with open(path, "w") as handle:
            handle.write(content)

        print(f"  Closed batch: sealed, {self.final_time_days} d")
        print(
            f"  Headspace:    {self.gas_saturation:.0%} of {VIAL_VOLUME_L * 1000:.0f} mL"
        )
        if self.label:
            print(f"  Condition:    {self.label}")
        return path


# ═════════════════════════════════════════════════════════════════════
# Series generation
# ═════════════════════════════════════════════════════════════════════


def generate_bottle_series(
    water_activities,
    output_dir=".",
    filename_template="bottle_aw{aw:.3f}.in",
    **generator_kwargs,
):
    """Generate one closed-batch deck per target water activity.

    Each deck differs only in its NaCl content; everything else is held fixed,
    so a comparison across the series isolates the salt effect.

    Parameters
    ----------
    water_activities : iterable of float
        Target water activities, one deck each. These should be the water
        activities actually measured for the batches being compared against,
        not a convenient round-numbered sweep.
    output_dir : str
        Directory to write into; created if absent.
    filename_template : str
        Format string taking ``aw``.
    **generator_kwargs
        Passed through to :class:`BottleGenerator`.

    Returns
    -------
    list of tuple
        ``(water_activity, molality, path)`` per deck, so the caller can record
        what composition each deck was built from.
    """
    os.makedirs(output_dir, exist_ok=True)
    results = []

    for aw in water_activities:
        molality = nacl_molality_for_water_activity(aw)
        brine = nacl_brine(molality=molality)
        label = f"target a_w = {aw:.3f}, NaCl = {molality:.3f} mol/kg"

        generator = BottleGenerator(brine=brine, label=label, **generator_kwargs)
        path = os.path.join(output_dir, filename_template.format(aw=aw))
        generator.generate(path)
        results.append((aw, molality, path))

    return results


def main():
    """Generate a deck series across the measured water activities.

    The default list is the six water activities present in the measured
    methanogen incubation series.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate closed-batch PFLOTRAN decks for sealed bottle incubations."
    )
    parser.add_argument(
        "--water-activities",
        type=float,
        nargs="+",
        default=[0.996, 0.905, 0.773, 0.696, 0.609, 0.454],
        help="Target water activities, one deck each.",
    )
    parser.add_argument(
        "--output-dir", default="decks", help="Directory to write decks into."
    )
    parser.add_argument(
        "--final-time-days",
        type=int,
        default=BOTTLE_FINAL_TIME_DAYS,
        help="Simulated duration in days.",
    )
    parser.add_argument(
        "--aw-threshold",
        type=float,
        default=AW_CRIT_HYDROGENOTROPHIC,
        help=(
            "Water activity below which the sandboxes inhibit methanogenesis. "
            f"Default {AW_CRIT_HYDROGENOTROPHIC} from generator.constants."
        ),
    )
    args = parser.parse_args()

    results = generate_bottle_series(
        args.water_activities,
        output_dir=args.output_dir,
        final_time_days=args.final_time_days,
        aw_threshold=args.aw_threshold,
    )

    print()
    print(f"Generated {len(results)} closed-batch decks in {args.output_dir}/")
    for aw, molality, path in results:
        saturated = " (above NaCl solubility)" if molality > 6.1 else ""
        print(f"  a_w {aw:.3f}  NaCl {molality:6.3f} mol/kg{saturated}  {path}")


if __name__ == "__main__":
    main()
