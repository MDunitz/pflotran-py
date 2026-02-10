"""
PFLOTRAN input file template strings.

Static text blocks for sections that don't depend on configurable parameters.
Parameterized sections live in PFLOTRANGenerator methods.

PFLOTRAN constraint type codes:
    F = free ion concentration [mol/L]
    T = total aqueous component concentration [mol/L]
    P = pH
    PE = pe (for O2(aq) or H+ only)
    M = concentration from mineral equilibrium [mol/L]
    G = concentration from gas equilibrium [bars]
    L = log10 of concentration
    Z = charge balance

PFLOTRAN notation:
    d0 = e0 (Fortran double-precision scientific notation)
    e.g. 31.d0 = 31 days, 2.d-4 = 2e-4

Ref: https://documentation.pflotran.org/user_guide/cards/subsurface/constraint_card.html
"""

# ─────────────────────────────────────────────────────────────────────
# Simulation setup
# ─────────────────────────────────────────────────────────────────────

HEADER = """\
# Description: Microbial decomposition and redox network for saline environments
# Generated on: {timestamp}

SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE GENERAL
      OPTIONS
        ARITHMETIC_GAS_DIFFUSIVE_DENSITY
        ISOTHERMAL
        REPLACE_INIT_PARAMS_ON_RESTART
      /
    /
    SUBSURFACE_TRANSPORT
      MODE GIRT
    /
  /

END


SUBSURFACE

#=========================== useful transport parameters =======================

REFERENCE_LIQUID_DENSITY 1.d3

#=========================== chemistry ========================================
CHEMISTRY
ACTIVITY_WATER
ACTIVITY_COEFFICIENTS TIMESTEP
OUTPUT
  WATER_ACTIVITY_COEFFICIENT
/"""


# ─────────────────────────────────────────────────────────────────────
# Species definitions
# ─────────────────────────────────────────────────────────────────────

PRIMARY_SPECIES = [
    "DOM1",
    "H+",
    "O2(aq)",
    "HCO3-",
    "Fe+++",
    "Fe++",
    "NH4+",
    "Tracer",
    "Tracer2",
    "Tracer3",
    "CH4(aq)",
    "Acetate-",
    "H2(aq)",
    "Mg++",
    "Ca++",
    "Na+",
    "K+",
    "CH3OH",
    "SO4--",
    "HS-",
    "H2O",
    "Cl-",
    "CO2(aq)",
    "NO3-",
    "N2(aq)",
]

SECONDARY_SPECIES = [
    "OH-",
    "FeCO3+",
    "Fe(OH)4-",
    "Acetic_acid(aq)",
    "FeCH3COO+",
    "Fe(OH)2(aq)",
    "FeCO3(aq)",
    "CO3--",
    "CaHCO3+",
    "MgCl+",
]

MINERALS = ["Fe(OH)3", "Fe(OH)2", "Rock(s)", "MgCl2.H2O"]


# ─────────────────────────────────────────────────────────────────────
# Mineral kinetics, sorption, gas species
# ─────────────────────────────────────────────────────────────────────

MINERAL_KINETICS_AND_SORPTION = """\
MINERAL_KINETICS
      Fe(OH)3
        RATE_CONSTANT  1.d-6 mol/m^2-sec
      /
      Fe(OH)2
        RATE_CONSTANT  1.d-7 mol/m^2-sec
      /
      Rock(s)
        RATE_CONSTANT  0.0 mol/m^2-sec
      /
      MgCl2.H2O
        RATE_CONSTANT  1.d-6 mol/m^2-sec
      /
    /
    IMMOBILE_SPECIES
      cellulose
      HRimm
    /
    PASSIVE_GAS_SPECIES
      O2(g)
    /
    ACTIVE_GAS_SPECIES
      GAS_TRANSPORT_IS_UNVETTED
      CO2(g)
    /
    SORPTION
      ION_EXCHANGE_RXN
        CEC 2.00e+02
        CATIONS
          Fe++ 1.00d-01
          Fe+++ 3.00d-01
          Mg++ 1.10e+00
          Ca++ 4.10e+00
          Na+ 1.00e+00 REFERENCE
          K+ 9.00d-01
          H+ 1.10e+00
        /
      /
    /"""


# ─────────────────────────────────────────────────────────────────────
# Chemistry output + solver
# ─────────────────────────────────────────────────────────────────────

CHEMISTRY_OUTPUT = """\
  LOG_FORMULATION
  TRUNCATE_CONCENTRATION 1.00d-25
  DATABASE {database_path}

  OUTPUT
    PH
    FREE_ION
      O2(aq)
      CH4(aq)
      CO2(aq)
      HCO3-
      SO4--
      HS-
      Acetate-
      Cl-
      NO3-
      N2(aq)
    TOTAL
      Fe++
      Fe+++
      DOM1
      Tracer2
      Tracer3
      CH3OH
      Mg++
    MINERALS
      MgCl2.H2O
  /
END"""


SOLVER = """\
#=========================== solver options ===================================
NUMERICAL_METHODS TRANSPORT
LINEAR_SOLVER
  SOLVER DIRECT
END
NEWTON_SOLVER
  ATOL 1.d-20
  RTOL 1.d-15
  STOL 1.d-10
  ITOL 1.d-10
  NUMERICAL_JACOBIAN
END
END"""


# ─────────────────────────────────────────────────────────────────────
# Material and fluid properties
# ─────────────────────────────────────────────────────────────────────

FLUID_PROPERTIES = """\
#=========================== fluid properties =================================
FLUID_PROPERTY
  PHASE LIQUID
  DIFFUSION_COEFFICIENT 1.d-11
END

FLUID_PROPERTY
  PHASE gas
  DIFFUSION_COEFFICIENT 2.d-10
END"""


MATERIAL_PROPERTIES = """\
#=========================== material properties ==============================
MATERIAL_PROPERTY soil1
  ID 1
  POROSITY 0.97d0
  SOIL_COMPRESSIBILITY 1.d-07
  SOIL_REFERENCE_PRESSURE 201325.d0
  ROCK_DENSITY 2650.0d0
  SPECIFIC_HEAT 830.0d0
  THERMAL_CONDUCTIVITY_DRY 0.12037926674717922d0
  THERMAL_CONDUCTIVITY_WET 1.6082691464310437d0
  CHARACTERISTIC_CURVES sf01
  PERMEABILITY
    PERM_X 6.5870260083342112d-013
    PERM_Y 6.5870260083342112d-013
    PERM_Z 9.5870260083342112d-14
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


# ─────────────────────────────────────────────────────────────────────
# Output options
# ─────────────────────────────────────────────────────────────────────

OUTPUT_OPTIONS = """\
#=========================== output options ===================================
OUTPUT
  TIME_UNITS d

  SNAPSHOT_FILE
    FORMAT HDF5
    PERIODIC TIME 1 d
    VARIABLES
      NO_FLOW_VARIABLES
      NO_ENERGY_VARIABLES
      LIQUID_SATURATION
    /
  /

  OBSERVATION_FILE
    PERIODIC TIME 1 d
  /
/"""


# ─────────────────────────────────────────────────────────────────────
# Microbial reaction templates
# ─────────────────────────────────────────────────────────────────────
# Each reaction is a dict with: comment, reaction_str, rate_key,
# monod list, inhibition list, reaction_type (MICROBIAL or GENERAL)

MICROBIAL_REACTIONS = [
    {
        "comment": "fermentation",
        "reaction": "1.00e+00 DOM1 + 6.67e-01 H2O -> 3.33e-01 Acetate- + 3.33e-01 HCO3- + 6.67e-01 H+ + 6.67e-01 H2(aq) + 3.33e-01 Tracer",
        "rate_key": "fermentation",
        "monod": [
            {
                "species": "DOM1",
                "ks_key": ("dom1", "fermentation"),
                "threshold_key": "general",
            },
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
            {
                "species": "Acetate-",
                "threshold_key": "acetate_inhibition",
                "direction": "ABOVE",
            },
            {"species": "Cl-", "threshold_key": "cl_inhibition", "direction": "ABOVE"},
        ],
    },
    {
        "comment": "DOM aerobic respiration",
        "reaction": "1.00e+00 DOM1 + 1.00e+00 O2(aq) + 1.00e+00 H2O -> 1.00e+00 HCO3- + 1.00e+00 H+ + 1.00e+00 Tracer",
        "rate_key": "dom_aerobic",
        "monod": [
            {"species": "O2(aq)", "ks_key": ("o2", "standard"), "threshold": 0.0},
            {
                "species": "DOM1",
                "ks_key": ("dom1", "aerobic"),
                "threshold_key": "very_low",
            },
        ],
        "inhibition": [],
    },
    {
        "comment": "Fe(II) microbial oxidation",
        "reaction": "1.00e+00 Fe++ + 2.50e-01 O2(aq) + 1.00e+00 H+ -> 1.00e+00 Fe+++ + 5.00e-01 H2O",
        "rate_key": "fe_microbial_oxidation",
        "monod": [
            {"species": "O2(aq)", "ks_key": ("o2", "fe_oxidation"), "threshold": 0.0},
            {"species": "Fe++", "ks_key": "fe_plus2", "threshold_key": "general"},
        ],
        "inhibition": [],
    },
    {
        "comment": "methylotrophic methanogenesis",
        "reaction": "1.00e+00 CH3OH + 1.00e+00 H2(aq) -> 1.00e+00 CH4(aq) + 1.00e+00 H2O",
        "rate_key": "methylotrophic_methano",
        "monod": [
            {"species": "CH3OH", "ks_key": "ch3oh", "threshold_key": "general"},
            {"species": "H2(aq)", "ks_key": "h2", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
            {"species": "Cl-", "threshold_key": "cl_inhibition", "direction": "ABOVE"},
        ],
    },
    {
        "comment": "acetate aerobic respiration",
        "reaction": "1.00e+00 Acetate- + 2.00e+00 O2(aq) -> 2.00e+00 HCO3- + 2.00e+00 H+ + 2.00e+00 Tracer",
        "rate_key": "acetate_aerobic",
        "monod": [
            {"species": "Acetate-", "ks_key": "acetate", "threshold": 0.0},
            {
                "species": "O2(aq)",
                "ks_key": ("o2", "standard"),
                "threshold_key": "very_low",
            },
        ],
        "inhibition": [],
    },
    {
        "comment": "hydrogen oxidation",
        "reaction": "2.00e+00 H2(aq) + 1.00e+00 O2(aq) -> 2.00e+00 H2O",
        "rate_key": "hydrogen_oxidation",
        "monod": [
            {"species": "H2(aq)", "ks_key": "h2", "threshold_key": "general"},
            {"species": "O2(aq)", "ks_key": ("o2", "standard"), "threshold": 0.0},
        ],
        "inhibition": [],
    },
    {
        "comment": "Fe(III) reduction",
        "reaction": "1.00e+00 Acetate- + 8.00e+00 Fe+++ + 4.00e+00 H2O -> 2.00e+00 HCO3- + 8.00e+00 Fe++ + 9.00e+00 H+ + 2.00e+00 Tracer",
        "rate_key": "fe_reduction",
        "monod": [
            {"species": "Acetate-", "ks_key": "acetate", "threshold_key": "general"},
            {"species": "Fe+++", "ks_key": "fe_plus3", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
        ],
    },
    {
        "comment": "sulfate reduction",
        "reaction": "1.00e+00 Acetate- + 1.00e+00 SO4-- + 1.00e+00 H+ -> 2.00e+00 HCO3- + 1.00e+00 HS- + 1.00e+00 Tracer",
        "rate_key": "sulfate_reduction",
        "monod": [
            {"species": "Acetate-", "ks_key": "acetate", "threshold_key": "general"},
            {"species": "SO4--", "ks_key": "so4", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
            {
                "species": "Fe+++",
                "threshold_key": "fe_inhibition",
                "direction": "ABOVE",
            },
        ],
    },
    {
        "comment": "ebullition (CH4 degassing proxy)",
        "reaction": "1.00e+00 CH4(aq) -> 1.00e+00 Tracer2",
        "rate_key": "ebullition",
        "monod": [
            {"species": "CH4(aq)", "ks_key": "ch4", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "CH4(aq)",
                "threshold_key": "ch4_ebullition",
                "direction": "BELOW",
            },
        ],
    },
    {
        "comment": "methane oxidation (O2)",
        "reaction": "1.00e+00 CH4(aq) + 2.00e+00 O2(aq) -> 1.00e+00 HCO3- + 1.00e+00 H+ + 1.00e+00 H2O + 1.00e+00 Tracer",
        "rate_key": "methane_o2_oxidation",
        "monod": [
            {
                "species": "O2(aq)",
                "ks_key": ("o2", "standard"),
                "threshold_key": "general",
            },
            {"species": "CH4(aq)", "ks_key": "ch4", "threshold_key": "general"},
        ],
        "inhibition": [],
    },
    {
        "comment": "methane oxidation (NO3)",
        "reaction": "5.00e+00 CH4(aq) + 8.00e+00 NO3- + 3.00e+00 H+ -> 5.00e+00 HCO3- + 4.00e+00 N2(aq) + 9.00e+00 H2O + 1.00e+00 Tracer",
        "rate_key": "methane_no3_oxidation",
        "monod": [
            {"species": "NO3-", "ks_key": "no3", "threshold_key": "general"},
            {"species": "CH4(aq)", "ks_key": "ch4", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
        ],
    },
    {
        "comment": "methane oxidation (SO4)",
        "reaction": "1.00e+00 CH4(aq) + 1.00e+00 SO4-- -> 1.00e+00 HCO3- + 1.00e+00 HS- + 1.00e+00 H2O + 1.00e+00 Tracer",
        "rate_key": "methane_so4_oxidation",
        "monod": [
            {"species": "SO4--", "ks_key": "so4", "threshold_key": "general"},
            {"species": "CH4(aq)", "ks_key": "ch4", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
        ],
    },
    {
        "comment": "methane oxidation (Fe)",
        "reaction": "1.00e+00 CH4(aq) + 8.00e+00 Fe+++ + 3.00e+00 H2O -> 1.00e+00 HCO3- + 8.00e+00 Fe++ + 9.00e+00 H+ + 1.00e+00 Tracer",
        "rate_key": "methane_fe_oxidation",
        "monod": [
            {"species": "Fe+++", "ks_key": "fe_plus3", "threshold_key": "general"},
            {"species": "CH4(aq)", "ks_key": "ch4", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
        ],
    },
    {
        "comment": "hydrogenotrophic methanogenesis",
        "reaction": "4.00e+00 H2(aq) + 1.00e+00 HCO3- + 1.00e+00 H+ -> 1.00e+00 CH4(aq) + 3.00e+00 H2O",
        "rate_key": "hydrogenotrophic_methano",
        "monod": [
            {"species": "H2(aq)", "ks_key": "h2", "threshold_key": "general"},
            {"species": "HCO3-", "ks_key": "hco3", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
            {
                "species": "Fe+++",
                "threshold_key": "fe_inhibition",
                "direction": "ABOVE",
            },
            {
                "species": "H+",
                "threshold_key": "h_plus_inhibition_1",
                "direction": "ABOVE",
            },
            {"species": "Cl-", "threshold_key": "cl_inhibition", "direction": "ABOVE"},
        ],
    },
    {
        "comment": "acetoclastic methanogenesis",
        "reaction": "1.00e+00 Acetate- + 1.00e+00 H2O -> 1.00e+00 CH4(aq) + 1.00e+00 HCO3- + 1.00e+00 Tracer",
        "rate_key": "acetaclastic_methano",
        "monod": [
            {"species": "Acetate-", "ks_key": "acetate", "threshold_key": "general"},
        ],
        "inhibition": [
            {
                "species": "O2(aq)",
                "threshold_key": "o2_inhibition",
                "direction": "ABOVE",
            },
            {
                "species": "Fe+++",
                "threshold_key": "fe_inhibition",
                "direction": "ABOVE",
            },
            {
                "species": "H+",
                "threshold_key": "h_plus_inhibition_2",
                "direction": "ABOVE",
            },
            {
                "species": "H+",
                "threshold_key": "h_plus_inhibition_3",
                "direction": "BELOW",
            },
            {"species": "Cl-", "threshold_key": "cl_inhibition", "direction": "ABOVE"},
        ],
    },
]

# Fe(II) abiotic oxidation is a GENERAL_REACTION, not MICROBIAL
GENERAL_REACTIONS = [
    {
        "comment": "Fe(II) abiotic oxidation",
        "reaction": "1.00e+00 Fe++ + 2.50e-01 O2(aq) + 1.00e+00 H+ <-> 1.00e+00 Fe+++ + 5.00e-01 H2O + 1.00e+00 Tracer3",
        "rate_key": "fe_abiotic_oxidation",
        "backward_rate": 0.0,
    },
]
