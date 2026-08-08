"""
PFLOTRAN Input File Generator

Generates .in files for PFLOTRAN reactive transport simulations of
microbial redox networks in saline environments.

Usage:
    generator = PFLOTRANGenerator(
        concentrations={'Cl-': '2.68 T', 'Na+': '2.295 T'},
        aw_threshold=0.6,
        dimensions='1d',
    )
    generator.generate('my_simulation.in')

Sources:
    [1] O'Meara et al. (2024) JGR-Biogeosciences, doi:10.1029/2023JG007633
    [2] Furukawa et al. (2004) Limnol. Oceanogr., 49(6), 2058-2072
    [3] Boudreau (1997) Diagenetic Models and Their Implementation, Springer
    [4] Thompson et al. (1995) Estuaries, 18(3), 399-408
    See REFERENCES.md for full citations and parameter derivations.
"""

from datetime import datetime

from .pflotran_templates import (
    HEADER,
    PRIMARY_SPECIES,
    SECONDARY_SPECIES,
    MINERALS,
    MINERAL_KINETICS_AND_SORPTION,
    CHEMISTRY_OUTPUT,
    SOLVER,
    FLUID_PROPERTIES,
    MATERIAL_PROPERTIES,
    OUTPUT_OPTIONS,
    MICROBIAL_REACTIONS,
    GENERAL_REACTIONS,
)

# ═════════════════════════════════════════════════════════════════════
# Default parameter sets
# ═════════════════════════════════════════════════════════════════════

DEFAULT_RATE_CONSTANTS = {
    # Reaction rate constants [mol/(L·s)] for MICROBIAL_REACTION
    # [mol/(L·s)] without biomass tracking
    "fermentation": 6.00e-08,
    "dom_aerobic": 1.80e-07,  # [1]
    "fe_abiotic_oxidation": 1.00e-02,  # GENERAL_REACTION [1/(mol·L·s)]
    "fe_microbial_oxidation": 5.50e-05,
    "methylotrophic_methano": 9.10e-06,
    "hydrogenotrophic_methano": 7.20e-09,
    "acetate_aerobic": 3.00e-07,
    "hydrogen_oxidation": 1.50e-06,
    "fe_reduction": 2.25e-10,  # [1] calibrated, see REFERENCES.md
    "sulfate_reduction": 1.50e-09,
    "ebullition": 3.00e-08,
    "acetaclastic_methano": 1.50e-08,
    "methane_o2_oxidation": 1.50e-08,
    "methane_no3_oxidation": 1.50e-08,
    "methane_so4_oxidation": 1.50e-08,
    "methane_fe_oxidation": 1.50e-08,
}

DEFAULT_HALF_SATURATION = {
    # Half-saturation constants [mol/L] for Monod kinetics
    "dom1": {"fermentation": 5.00e-02, "aerobic": 1.00e-01},  # [1]
    "o2": {"standard": 1.00e-04, "fe_oxidation": 1.00e-08},
    "fe_plus2": 1.00e-04,
    "ch3oh": 1.00e-01,
    "h2": 1.00e-01,
    "hco3": 1.00e-01,
    "acetate": 4.00e-02,
    "fe_plus3": 1.00e-10,
    "so4": 1.00e-04,
    "ch4": 4.00e-02,
    "no3": 1.00e-04,
}

DEFAULT_THRESHOLDS = {
    # Threshold concentrations [mol/L] for inhibition switches
    "general": 1.10e-15,
    "very_low": 1.10e-16,
    "o2_inhibition": 1.00e-06,
    "acetate_inhibition": 8.00e-02,
    "cl_inhibition": 2.00e-01,  # seawater 1× Cl⁻ = 0.536 M
    "fe_inhibition": 1.00e-09,
    "h_plus_inhibition_1": 1.78e-06,
    "h_plus_inhibition_2": 2.88e-05,
    "h_plus_inhibition_3": 2.88e-07,
    "ch4_ebullition": 2.50e-03,
}

DEFAULT_INITIAL_CONCENTRATIONS = {
    # PFLOTRAN constraint strings: value + type code
    "DOM1": "5.00 T",
    "H+": "6.5 P",
    "O2(aq)": "2.d-4 T",
    "CO2(aq)": "2.d-4 T",
    "HCO3-": "4.00d-6 T",
    "Fe+++": "5.3d-5 M Fe(OH)3",
    "Fe++": "2.4d-3 T Fe(OH)2",
    "Mg++": "9.00d-04 T",
    "Ca++": "5.00d-04 T",
    "Na+": "2.00d-04 T",
    "K+": "2.00d-05 T",
    "CH3OH": "2.00d-03 T",
    "SO4--": "1.00d-03 T",
    "Cl-": "6.00d-04 Z",
}

# Atmospheric boundary condition
# CO2(aq): ~1.9e-5 M at 8°C, see REFERENCES.md (Henry's law derivation)
DEFAULT_ATMOSPHERIC_CONCENTRATIONS = {
    "H+": "5.0 P",
    "O2(aq)": "1.00d-04 T",
    "CO2(aq)": "1.906e-05 T",
    "HCO3-": "1.00d-05 T",
}

# Species present but set to trace levels in initial constraint
TRACE_SPECIES = [
    "NH4+",
    "Tracer",
    "Tracer2",
    "Tracer3",
    "CH4(aq)",
    "Acetate-",
    "H2(aq)",
    "HS-",
    "NO3-",
    "N2(aq)",
]

# Grid presets for different dimensionalities
#
# The 1D column is the default — it captures the vertical redox ladder
# (O₂ penetration → Fe reduction → sulfate reduction → methanogenesis)
# which is the dominant physics in a uniform sediment column.
#
# TODO: 2D and 3D presets should be expanded when lateral processes are added:
#   - 2D cross-section: tidal channel boundary on one side (time-varying
#     salinity/O₂ Dirichlet BC on WEST face), lateral salinity gradients
#     from evaporation ponds, or root O₂ injection at specific X positions.
#   - 3D: heterogeneous permeability fields, point-source root injection
#     at specific (X,Y) coordinates, or coupled surface water flow.
#   Each of these requires adding:
#     1. Additional BOUNDARY_CONDITION blocks in _build_regions_and_conditions()
#     2. Time-varying FLOW_CONDITION / TRANSPORT_CONDITION if tidal
#     3. Possibly multiple MATERIAL_PROPERTY zones (e.g. root zone vs bulk)
#     4. Region definitions for the lateral boundaries
# Rate keys of the reactions that produce methane. A salinity inhibition term
# has to be attached to these to have any effect on modelled methane; attaching
# it anywhere else -- as the AWINHIBIT sandboxes effectively do -- leaves the
# production pathways untouched.
METHANOGENESIS_RATE_KEYS = (
    "methylotrophic_methano",
    "hydrogenotrophic_methano",
    "acetaclastic_methano",
)

GRID_PRESETS = {
    "1d": {
        "grid_cells": "1 1 10",
        "cell_size_x": "1.0d0",
        "cell_size_y": "1.0d0",
        "cell_size_z": "0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0",
        "domain": (1.0, 1.0, 1.0),  # (Lx, Ly, Lz) in meters
    },
    "2d": {
        "grid_cells": "10 1 10",
        "cell_size_x": "0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0",
        "cell_size_y": "1.0d0",
        "cell_size_z": "0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0 0.1d0",
        "domain": (1.0, 1.0, 1.0),
    },
    "3d": {
        "grid_cells": "4 4 4",
        "cell_size_x": "0.25d0 0.25d0 0.25d0 0.25d0",
        "cell_size_y": "0.25d0 0.25d0 0.25d0 0.25d0",
        "cell_size_z": "0.25d0 0.25d0 0.25d0 0.25d0",
        "domain": (1.0, 1.0, 1.0),
    },
}


# ═════════════════════════════════════════════════════════════════════
# Generator class
# ═════════════════════════════════════════════════════════════════════


class PFLOTRANGenerator:
    """Generates PFLOTRAN .in files from configurable parameters.

    Parameters are organized into groups. Pass overrides as kwargs;
    anything not overridden uses the default from DEFAULT_* dicts above.
    """

    def __init__(
        self,
        # --- Concentrations (override individual species) ---
        concentrations=None,
        atmospheric_concentrations=None,
        # --- Kinetic parameters (override individual reactions) ---
        rate_constants=None,
        half_saturation=None,
        thresholds=None,
        # --- Inhibition mechanism toggles ---
        # Both ON by default (current behavior). For the double-counting
        # diagnostic (PR #50 item 4), run three variants:
        #   A: enable_cl_inhibition=True,  enable_aw_sandbox=False
        #   B: enable_cl_inhibition=False, enable_aw_sandbox=True
        #   C: enable_cl_inhibition=True,  enable_aw_sandbox=True   (default)
        enable_cl_inhibition=True,
        enable_aw_sandbox=True,
        # --- Reaction sandbox: water activity inhibition ---
        aw_threshold=0.5,
        aw_rate_constant=1.0e-10,
        aw_inhibition_type="THRESHOLD",
        # --- Domain geometry ---
        dimensions="1d",
        # --- Simulation control ---
        temperature=8.0,
        final_time_days=31,
        initial_timestep_hours=2.0,
        max_timestep_hours=12.0,
        # --- Chemistry configuration ---
        # Whether dissolved carbon dioxide is held in equilibrium with the
        # carbonate system, or carried as an independent primary species.
        #
        # Defaults to False, which is the historical behaviour of the sediment
        # column decks and is left alone here so that existing column results
        # stay reproducible.
        #
        # Setting it True moves CO2(aq) from the primary list to the secondary
        # list, so PFLOTRAN computes it from bicarbonate and pH through the
        # database reaction
        #
        #     CO2(aq) = HCO3- + H+ - H2O
        #
        # (hanford.dat, log K -6.3447 at 25 C). This matters because no reaction
        # in the network produces CO2(aq): every carbon-oxidising step yields
        # bicarbonate. Carried as a decoupled primary species, dissolved carbon
        # dioxide therefore never moves from its initial value, however much
        # carbon the organisms respire, and the model cannot predict a headspace
        # carbon dioxide concentration at all.
        #
        # Note that the same argument does NOT apply to CH4(aq), whose database
        # entry is a redox reaction rather than an acid-base one. Decoupling
        # methane is what allows the kinetic network to produce it, and must
        # stay.
        couple_carbonate=False,
        # --- Salinity inhibition on the reaction network itself ---
        # An inhibition term added directly to the network's methanogenesis
        # reactions, as opposed to the AWINHIBIT reaction sandboxes.
        #
        # This exists because the sandboxes do not inhibit the network. They add
        # their own parallel copies of the three methanogenesis pathways and
        # inhibit only those, at a rate constant of 1e-10 against the network's
        # 9.1e-6 for the methylotrophic route -- roughly ninety thousand times
        # smaller, before accounting for the sandbox's sixth-order rate law.
        # Raising the sandbox threshold until it is fully engaged in every
        # bottle changes the modelled methane by under one percent, because the
        # pathway it governs produces almost none of it.
        #
        # Pass a dict to switch this on, for example::
        #
        #     {"species": "Cl-", "threshold": 1.5, "interval": 1.0}
        #
        # ``threshold`` is in mol/L and ``interval`` is the width of the
        # transition in decades. TYPE SMOOTHSTEP is used rather than MONOD
        # deliberately: Monod inhibition is hyperbolic, so the most it can
        # deliver between the weakest and strongest brine here is roughly a
        # factor of twenty, whereas the measurements fall by four orders of
        # magnitude across the same range. A sigmoid can express a collapse; a
        # hyperbola cannot, at any half-saturation value.
        #
        # Defaults to None, leaving every existing deck unchanged.
        salinity_inhibition=None,
        # --- Paths ---
        database_path="/home/sshindad/miniconda/pflotran/md_test_files/hanford.dat",
    ):
        # Merge user overrides with defaults (user wins)
        self.rate_constants = {**DEFAULT_RATE_CONSTANTS, **(rate_constants or {})}
        self.half_saturation = {**DEFAULT_HALF_SATURATION, **(half_saturation or {})}
        self.thresholds = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
        self.concentrations = {
            **DEFAULT_INITIAL_CONCENTRATIONS,
            **(concentrations or {}),
        }
        self.atmospheric = {
            **DEFAULT_ATMOSPHERIC_CONCENTRATIONS,
            **(atmospheric_concentrations or {}),
        }

        # Water activity sandbox parameters
        self.aw_threshold = aw_threshold
        self.aw_rate_constant = aw_rate_constant
        self.aw_inhibition_type = aw_inhibition_type

        # Inhibition mechanism toggles
        self.enable_cl_inhibition = enable_cl_inhibition
        self.enable_aw_sandbox = enable_aw_sandbox

        # Chemistry configuration
        self.couple_carbonate = couple_carbonate
        self.salinity_inhibition = salinity_inhibition

        # Domain
        self.dimensions = dimensions.lower()
        grid = GRID_PRESETS[self.dimensions]
        self.grid_cells = grid["grid_cells"]
        self.cell_size_x = grid["cell_size_x"]
        self.cell_size_y = grid["cell_size_y"]
        self.cell_size_z = grid["cell_size_z"]
        self.domain = grid["domain"]

        # Simulation control
        self.temperature = temperature
        self.final_time_days = final_time_days
        self.initial_timestep_hours = initial_timestep_hours
        self.max_timestep_hours = max_timestep_hours
        self.database_path = database_path

    # ─────────────────────────────────────────────────────────────────
    # Helper: look up a half-saturation constant from flat or nested key
    # ─────────────────────────────────────────────────────────────────

    def _get_ks(self, key):
        """Resolve half-saturation constant from key or (key, subkey) tuple."""
        if isinstance(key, tuple):
            return self.half_saturation[key[0]][key[1]]
        return self.half_saturation[key]

    # ─────────────────────────────────────────────────────────────────
    # Section builders (each returns a string)
    # ─────────────────────────────────────────────────────────────────

    def _primary_species(self):
        """Species carried as independent primary unknowns in this deck."""
        primary = list(PRIMARY_SPECIES)
        if self.couple_carbonate and "CO2(aq)" in primary:
            primary.remove("CO2(aq)")
        return primary

    def _include_reaction(self, rxn):
        """Whether a reaction from the network belongs in this deck.

        A seam for subclasses that model a different physical system. The base
        generator includes every reaction.
        """
        return True

    def _build_mineral_kinetics_and_sorption(self):
        """Mineral kinetics, immobile species, gas species and sorption."""
        return MINERAL_KINETICS_AND_SORPTION

    def _build_chemistry_output(self):
        """Closing CHEMISTRY block: output requests and the database path."""
        return CHEMISTRY_OUTPUT.format(database_path=self.database_path)

    def _build_species_lists(self):
        """PRIMARY_SPECIES, DECOUPLED_EQUILIBRIUM_REACTIONS, SECONDARY_SPECIES, MINERALS

        When ``couple_carbonate`` is set, dissolved carbon dioxide moves out of
        the primary list and into the secondary list. See the constructor
        documentation for why.
        """
        primary = self._primary_species()
        secondary = list(SECONDARY_SPECIES)
        if self.couple_carbonate and "CO2(aq)" not in secondary:
            secondary.append("CO2(aq)")

        lines = ["\nPRIMARY_SPECIES"]
        for s in primary:
            lines.append(f"  {s}")
        lines.append("/")

        # The decoupled list is the primary list. Every primary species here is
        # a redox or acid-base species that the database would otherwise hold at
        # equilibrium; decoupling lets the kinetic reaction network drive them.
        lines.append("DECOUPLED_EQUILIBRIUM_REACTIONS")
        for s in primary:
            lines.append(f"  {s}")
        lines.append("/")

        lines.append("SECONDARY_SPECIES")
        for s in secondary:
            lines.append(f"  {s}")
        lines.append("/")

        lines.append("MINERALS")
        for m in MINERALS:
            lines.append(f"  {m}")
        lines.append("/")
        return "\n".join(lines)

    def _build_microbial_reaction(self, rxn):
        """Render one MICROBIAL_REACTION block from a reaction dict."""
        lines = [f'  # {rxn["comment"]}']
        lines.append("  MICROBIAL_REACTION")
        lines.append(f'    REACTION {rxn["reaction"]}')
        lines.append(
            f'    RATE_CONSTANT       {self.rate_constants[rxn["rate_key"]]:.2e}'
        )

        for m in rxn["monod"]:
            lines.append("    MONOD")
            lines.append(f'      SPECIES_NAME        {m["species"]}')
            lines.append(
                f'      HALF_SATURATION_CONSTANT {self._get_ks(m["ks_key"]):.2e}'
            )
            if "threshold_key" in m:
                lines.append(
                    f'      THRESHOLD_CONCENTRATION {self.thresholds[m["threshold_key"]]:.2e}'
                )
            else:
                lines.append(
                    f'      THRESHOLD_CONCENTRATION {m.get("threshold", 0.0):.2e}'
                )
            lines.append("    /")

        for inh in rxn["inhibition"]:
            # Skip Cl⁻ inhibition if disabled (for double-counting diagnostic)
            if inh["species"] == "Cl-" and not self.enable_cl_inhibition:
                continue
            lines.append("    INHIBITION")
            lines.append(f'      SPECIES_NAME        {inh["species"]}')
            lines.append("      TYPE MONOD")
            lines.append(
                f'      THRESHOLD_CONCENTRATION {self.thresholds[inh["threshold_key"]]:.2e}'
            )
            lines.append(f'      INHIBIT_{inh["direction"]}_THRESHOLD')
            lines.append("    /")

        lines.extend(self._build_salinity_inhibition(rxn))

        lines.append("  /")
        return "\n".join(lines)

    def _build_salinity_inhibition(self, rxn):
        """Extra inhibition lines for one reaction, or nothing.

        Applied only to the methanogenesis reactions. Inhibiting fermentation
        or the oxidation steps as well would suppress the whole carbon chain
        rather than the methanogens specifically, which is not what the salt
        stress in these incubations is understood to do.
        """
        spec = self.salinity_inhibition
        if not spec or rxn.get("rate_key") not in METHANOGENESIS_RATE_KEYS:
            return []

        lines = [
            "    INHIBITION",
            f'      SPECIES_NAME        {spec["species"]}',
            "      TYPE SMOOTHSTEP",
            f'      SMOOTHSTEP_INTERVAL {spec.get("interval", 1.0):.2f}',
            f'      THRESHOLD_CONCENTRATION {spec["threshold"]:.2e}',
            "      INHIBIT_ABOVE_THRESHOLD",
            "    /",
        ]
        return lines

    def _build_general_reaction(self, rxn):
        """Render one GENERAL_REACTION block."""
        lines = [f'  # {rxn["comment"]}']
        lines.append("  GENERAL_REACTION")
        lines.append(f'    REACTION {rxn["reaction"]}')
        lines.append(
            f'    FORWARD_RATE        {self.rate_constants[rxn["rate_key"]]:.2e}'
        )
        lines.append(f'    BACKWARD_RATE       {rxn["backward_rate"]:.2e}')
        lines.append("  /")
        return "\n".join(lines)

    def _build_all_reactions(self):
        """All MICROBIAL + GENERAL reactions."""
        blocks = []
        for rxn in MICROBIAL_REACTIONS:
            if not self._include_reaction(rxn):
                continue
            blocks.append(self._build_microbial_reaction(rxn))
        for rxn in GENERAL_REACTIONS:
            blocks.append(self._build_general_reaction(rxn))
        return "\n\n".join(blocks)

    def _build_reaction_sandbox(self):
        """REACTION_SANDBOX block for water activity inhibition.

        Three sandboxes targeting different methanogenesis pathways:
          AWINHIBIT        — hydrogenotrophic (4 H2 + HCO3- + H+ -> CH4 + 3 H2O)
          AWINHIBITACETATE — acetoclastic (Acetate- + H2O -> CH4 + HCO3-)
          AWINHIBITMETHYL  — methylotrophic (CH3OH + H2 -> CH4 + H2O)

        Parameters read by Fortran sandbox code:
          WATER_ACTIVITY_THRESHOLD — a_w below which reaction is fully inhibited
          RATE_CONSTANT — base rate [mol/(m³·s)] for the sandbox reaction
          INHIBITION_TYPE — THRESHOLD (binary) or SMOOTHSTEP (gradual)
        """
        sandbox_names = ["AWINHIBIT", "AWINHIBITACETATE", "AWINHIBITMETHYL"]
        lines = ["\nREACTION_SANDBOX"]
        for name in sandbox_names:
            lines.append(f"  {name}")
            lines.append(f"    WATER_ACTIVITY_THRESHOLD {self.aw_threshold:.1e}")
            lines.append(f"    RATE_CONSTANT {self.aw_rate_constant:.1e}")
            lines.append(f"    INHIBITION_TYPE {self.aw_inhibition_type}")
            lines.append("  /")
        lines.append("/")
        return "\n".join(lines)

    def _build_constraints(self):
        """Transport constraints: initial conditions + atmospheric boundary."""
        # Initial constraint
        lines = [
            "\n#=========================== transport constraints ============================",
            "CONSTRAINT initial",
            "  IMMOBILE",
            "    cellulose            8.00e+03",
            "    HRimm                1.00d-20",
            "  /",
            "  CONCENTRATIONS",
        ]

        # Parameterized species from self.concentrations. A species that has
        # been moved to the secondary list is computed by PFLOTRAN rather than
        # constrained, so it must not appear here.
        for species in self._primary_species():
            if species in self.concentrations:
                lines.append(f"    {species:20s}{self.concentrations[species]}")
            elif species in TRACE_SPECIES:
                lines.append(f"    {species:20s}1.00d-15 T")
            elif species == "H2O":
                lines.append(f'    {"H2O":20s}1.00d-03 T')

        lines.extend(
            [
                "  /",
                "  MINERALS",
                "    Fe(OH)3             9.6d-6  1.d2 m^2/m^3",
                "    Fe(OH)2             7.2d-1  1.d2 m^2/m^3",
                "    Rock(s)             0.5  5.0e3 m^2/m^3",
                "    MgCl2.H2O           1.0d-02  1.0e2 m^2/m^3",
                "  /",
                "END",
            ]
        )

        # Atmospheric constraint — only the species that differ from initial
        lines.extend(
            [
                "",
                "CONSTRAINT atmospheric",
                "  CONCENTRATIONS",
            ]
        )
        # Use atmospheric overrides; fill rest from initial or trace
        for species in PRIMARY_SPECIES:
            if species in self.atmospheric:
                lines.append(f"    {species:20s}{self.atmospheric[species]}")
            elif species in self.concentrations:
                lines.append(f"    {species:20s}{self.concentrations[species]}")
            elif species in TRACE_SPECIES:
                lines.append(f"    {species:20s}1.00d-15 T")
            elif species == "H2O":
                lines.append(f'    {"H2O":20s}1.00d-03 T')
        lines.extend(["  /", "END"])

        return "\n".join(lines)

    def _build_grid_and_time(self):
        """Discretization and time stepping."""
        return f"""\
#=========================== discretization ===================================
GRID
  TYPE structured
  ORIGIN 0.d0 0.d0 0.d0
  NXYZ {self.grid_cells}
  DXYZ
     {self.cell_size_x}
     {self.cell_size_y}
     {self.cell_size_z}
  /
END

PROC 1 1 1

#=========================== times ============================================
TIME
  FINAL_TIME {self.final_time_days} d
  INITIAL_TIMESTEP_SIZE {self.initial_timestep_hours:.1f}d0 h
  MAXIMUM_TIMESTEP_SIZE {self.max_timestep_hours:.1f}d0 h
END"""

    def _build_regions_and_conditions(self):
        """Regions, observations, flow/transport conditions, couplers."""
        lx, ly, lz = self.domain
        return f"""\
#=========================== regions ==========================================
REGION all
  COORDINATES
    0.d0 0.d0 0.d0
    {lx:.1f}d0 {ly:.1f}d0 {lz:.1f}d0
  /
END

REGION top_surface
  COORDINATES
    0.d0 0.d0 {lz:.1f}d0
    {lx:.1f}d0 {ly:.1f}d0 {lz:.1f}d0
  /
  FACE TOP
END

REGION bottom_surface
  COORDINATES
    0.d0 0.d0 0.d0
    {lx:.1f}d0 {ly:.1f}d0 0.d0
  /
  FACE BOTTOM
END

REGION center_obs
  COORDINATE {lx/2:.1f}d0 {ly/2:.1f}d0 {lz/2:.1f}d0 # noqa: E226
END

#=========================== observation points ===============================
OBSERVATION
  REGION center_obs
END

#=========================== transport conditions =============================
FLOW_CONDITION initial
  TYPE
    GAS_PRESSURE dirichlet
    GAS_SATURATION dirichlet
    TEMPERATURE dirichlet
  /
  GAS_PRESSURE 2.d5
  GAS_SATURATION 0.25
  TEMPERATURE {self.temperature:.1f}d0
/

FLOW_CONDITION atmospheric
  TYPE
    GAS_PRESSURE dirichlet
    GAS_SATURATION dirichlet
    TEMPERATURE dirichlet
  /
  GAS_PRESSURE 2.d5
  GAS_SATURATION 0.15
  TEMPERATURE {self.temperature:.1f}d0
/

TRANSPORT_CONDITION initial
  TYPE dirichlet
  CONSTRAINT_LIST
    0.d0 initial
  /
END

TRANSPORT_CONDITION atmospheric
  TYPE dirichlet
  CONSTRAINT_LIST
    0.d0 atmospheric
  /
END

#=========================== condition couplers ===============================
INITIAL_CONDITION
  TRANSPORT_CONDITION initial
  FLOW_CONDITION initial
  REGION all
END

BOUNDARY_CONDITION top_atm
  TRANSPORT_CONDITION atmospheric
  FLOW_CONDITION atmospheric
  REGION top_surface
END

#=========================== stratigraphy couplers ============================
STRATA
  REGION all
  MATERIAL soil1
END

END_SUBSURFACE"""

    # ─────────────────────────────────────────────────────────────────
    # Main assembly
    # ─────────────────────────────────────────────────────────────────

    def generate(self, filename="pflotran_input.in"):
        """Assemble and write the complete PFLOTRAN input file.

        Sections are built from templates (static text) and parameters
        (self.*). No inline PFLOTRAN strings in this method.
        """
        sections = [
            HEADER.format(timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            self._build_species_lists(),
            self._build_mineral_kinetics_and_sorption(),
            self._build_all_reactions(),
            self._build_reaction_sandbox() if self.enable_aw_sandbox else "",
            self._build_chemistry_output(),
            self._build_constraints(),
            SOLVER,
            self._build_grid_and_time(),
            FLUID_PROPERTIES,
            MATERIAL_PROPERTIES,
            OUTPUT_OPTIONS,
            self._build_regions_and_conditions(),
        ]

        content = "\n\n".join(sections)

        with open(filename, "w") as f:
            f.write(content)

        print(f"Generated PFLOTRAN input file: {filename}")
        print(f"  Dimensions: {self.dimensions}")
        print(f"  Grid: {self.grid_cells}")
        print(f"  Temperature: {self.temperature}°C")
        print(f"  a_w threshold: {self.aw_threshold}")
        print(f"  Cl⁻ inhibition: {'ON' if self.enable_cl_inhibition else 'OFF'}")
        print(f"  a_w sandbox: {'ON' if self.enable_aw_sandbox else 'OFF'}")
        return filename


# ═════════════════════════════════════════════════════════════════════
# CLI entry point
# ═════════════════════════════════════════════════════════════════════


def main():
    generator = PFLOTRANGenerator()
    generator.generate()


if __name__ == "__main__":
    main()
