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
3.  **The run is long enough to compare against.** Sixty days by default, which
    spans both measured incubation series (42 and 51 days), rather than the
    31-day column default.
4.  **The temperature is the incubation temperature.** 18 degrees Celsius, which
    is what the post-processing package in ``config.py`` already assumes, rather
    than the 8 degrees Celsius used for coastal sediment.

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

from .pflotran_generator import PFLOTRANGenerator

# ═════════════════════════════════════════════════════════════════════
# Vial geometry
# ═════════════════════════════════════════════════════════════════════
#
# These match the constants the measurement pipeline uses to turn headspace
# concentrations into moles (saltyBiomass, experiments/analysis/incubations/
# constants.py). If those change, these must change with them, or the model and
# the measurement will be describing differently-sized bottles.

VIAL_VOLUME_L = 0.125  # total internal volume of the vial [L]
HEADSPACE_VOLUME_L = 0.100  # gas volume above the liquid [L]
LIQUID_VOLUME_L = VIAL_VOLUME_L - HEADSPACE_VOLUME_L  # brine + biomass [L]

# PFLOTRAN works in metres. A cube of this edge length has the vial's volume.
# 0.125 L = 1.25e-4 m^3, and (1.25e-4)^(1/3) = 0.05 m exactly.
_VIAL_EDGE_M = (VIAL_VOLUME_L * 1e-3) ** (1.0 / 3.0)

# Porosity of a bottle is not the porosity of packed sediment. The vial is
# essentially all fluid, so porosity approaches 1. It is held just below 1
# because a porosity of exactly 1 leaves PFLOTRAN with no solid phase for the
# mineral reactions to attach to.
BOTTLE_POROSITY = 0.99

# Fraction of pore space occupied by gas. With porosity ~1 this is just the
# headspace fraction of the vial.
BOTTLE_GAS_SATURATION = HEADSPACE_VOLUME_L / VIAL_VOLUME_L  # 0.8

# Bottles are sealed at roughly local atmospheric pressure. The measurement
# pipeline uses 0.969 atm (Pasadena, ~260 m elevation); the difference from
# 1 atm is well inside the uncertainty on everything else here, so 1 atm is
# used and the value is exposed as a parameter for anyone who needs it exact.
BOTTLE_GAS_PRESSURE_PA = 1.01325e5

# Incubation temperature [deg C]. The post-processing package already assumes
# 18 C for its Stokes-Einstein diffusion correction (see config.py); this makes
# the simulation agree with it.
BOTTLE_TEMPERATURE_C = 18.0

# Long enough to span the measured incubations with margin. The pipeline's own
# output runs to 119 days for Exp003 and 122 for Exp004, so a simulation must
# reach at least 122 days for the comparison to cover the whole measured record
# rather than stopping partway through it.
#
# An earlier value of 60 days was set from the older exported files, which end
# at 42 and 51 days. Those exports turned out to be a stale snapshot; the live
# pipeline output runs twice as long.
BOTTLE_FINAL_TIME_DAYS = 130

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


def _default_database_path():
    """Absolute path to the repository's own thermodynamic database.

    The sediment-column generator defaults to a developer-specific absolute
    path that exists on one machine. Resolving against this file's location
    makes a generated deck runnable from any clone.
    """
    package_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(package_dir, "..", "..", ".."))
    return os.path.join(repo_root, "sandbox", "hanford.dat")


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
        **kwargs,
    ):
        merged_concentrations = dict(concentrations or {})
        if brine:
            merged_concentrations.update(brine)

        super().__init__(
            concentrations=merged_concentrations,
            temperature=temperature,
            final_time_days=final_time_days,
            database_path=database_path or _default_database_path(),
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
        default=0.5,
        help="Water activity below which the sandboxes inhibit methanogenesis.",
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
