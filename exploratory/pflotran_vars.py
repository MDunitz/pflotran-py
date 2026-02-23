############################
# Introduction
############################
"""
Anything that was not specified to be changed, is taken from the paper below:
O’Meara, T. A., Yuan, F., Sulman, B. N., Noyce, G. L., Rich, R., Thornton, P. E., & Megonigal, J. P. (2024). 
Developing a Redox Network for Coastal Saltmarsh Systems in the PFLOTRAN Reaction Model. 
Journal of Geophysical Research: Biogeosciences, 129(3), e2023JG007633. https://doi.org/10.1029/2023JG007633 

"""

############################
# List of Primary Species
############################
"""
  DOM1
  H+
  O2(aq)
  HCO3-
  Fe+++
  Fe++
  NH4+
  Tracer
  Tracer2
  Tracer3
  CH4(aq)
  Acetate-
  H2(aq)
  Mg++
  Ca++
  Na+
  K+
  CH3OH
"""

############################
# List of Secondary Species
############################
"""
  CO2(aq)
  OH-
  FeCO3+
  Fe(OH)4-
  Acetic_acid(aq)
  FeCH3COO+
  # FeIIIDOM1(aq)
  # FeIIDOM1(aq)
  # FeIIIAcetate(aq)
  Fe(OH)2(aq)
  FeCO3(aq)
  CO3--
  CaHCO3+
"""
############################
# The rate constants
# we would need to find the reaction rates, inhibition thresholds, etc. or find papers from which to find them.
############################

# 1) fermentation
"""
1.00e+00 DOM1  + 6.67e-01 H2O  -> 3.33e-01 Acetate-  + 3.33e-01 HCO3-  + 6.67e-01 H+  + 6.67e-01 H2(aq)  + 3.33e-01 Tracer 
"""

ferm_rxn_const = 6.00e-07

ferm_monod1_halfsat_conc = 5.00e-02
ferm_monod1_thresh_conc = 1.10e-15

ferm_o2_inhib_above = 1.00e-05
ferm_acetate_inhib_above = 4.00e-02

# 2) DOM aerobic respiration
"""
1.00e+00 DOM1  + 1.00e+00 O2(aq)  + 1.00e+00 H2O  -> 1.00e+00 HCO3-  + 1.00e+00 H+  + 1.00e+00 Tracer 
"""

dom_aerobic_rxn_const = 1.80e-07

dom_aerobic_o2_halfsat_conc = 1.00e-04
dom_aerobic_o2_thresh_conc = 0.00e+00

dom_aerobic_dom1_halfsat_conc = 1.00e-01
dom_aerobic_dom1_thresh_conc = 1.10e-16

# 3) Fe(II) abiotic oxidation
"""
1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  <-> 1.00e+00 Fe+++  + 5.00e-01 H2O + 1.00e+00 Tracer3
"""

fe2_abiotic_forward_rate = 1.00e-2
fe2_abiotic_backward_rate = 0.00e-06

# 4) Fe(II) microbial oxidation
"""
1.00e+00 Fe++  + 2.50e-01 O2(aq)  + 1.00e+00 H+  -> 1.00e+00 Fe+++  + 5.00e-01 H2O 
"""

fe2_microbial_rxn_const = 5.5e-05

fe2_microbial_o2_halfsat_conc = 1.00e-08
fe2_microbial_o2_thresh_conc = 0.00e+00

fe2_microbial_fe2_halfsat_conc = 1.00e-04
fe2_microbial_fe2_thresh_conc = 1.10e-15

# 5) Hydrogenotrophic methanogenesis
"""
4.00e+00 H2(aq)  + 1.00e+00 HCO3-  + 1.00e+00 H+  -> 1.00e+00 CH4(aq)  + 3.00e+00 H2O 
"""

hydro_methano_rxn_const = 7.2e-09

hydro_methano_h2_halfsat_conc = 1.00e-01
hydro_methano_h2_thresh_conc = 1.10e-15

hydro_methano_hco3_halfsat_conc = 1.00e-01
hydro_methano_hco3_thresh_conc = 1.10e-15

hydro_methano_o2_inhib_above = 1.00e-05
hydro_methano_fe3_inhib_above = 1.00e-10
hydro_methano_h_inhib_above = 1.78e-07

# 6) Acetate aerobic respiration
"""
1.00e+00 Acetate-  + 2.00e+00 O2(aq)  -> 2.00e+00 HCO3-  + 2.00e+00 H+  + 2.00e+00 Tracer 
"""

acetate_aerobic_rxn_const = 3.00e-07

acetate_aerobic_o2_halfsat_conc = 1.00e-04
acetate_aerobic_o2_thresh_conc = 0.00e+00

acetate_aerobic_acetate_halfsat_conc = 4.00e-02
acetate_aerobic_acetate_thresh_conc = 1.10e-16

# 7) Hydrogen oxidation
"""
2.00e+00 H2(aq)  + 1.00e+00 O2(aq)  -> 2.00e+00 H2O 
"""

h2_oxidation_rxn_const = 1.5e-06

h2_oxidation_h2_halfsat_conc = 1.00e-01
h2_oxidation_h2_thresh_conc = 1.10e-15

h2_oxidation_o2_halfsat_conc = 1.00e-04
h2_oxidation_o2_thresh_conc = 0.00e+00

# 8) Fe(III) reduction
"""
1.00e+00 Acetate-  + 8.00e+00 Fe+++  + 4.00e+00 H2O  -> 2.00e+00 HCO3-  + 8.00e+00 Fe++  + 9.00e+00 H+  + 2.00e+00 Tracer 
"""

fe3_reduction_rxn_const = 2.25e-08

fe3_reduction_acetate_halfsat_conc = 4.00e-2
fe3_reduction_acetate_thresh_conc = 1.10e-15

fe3_reduction_fe3_halfsat_conc = 1.00e-10
fe3_reduction_fe3_thresh_conc = 1.10e-15

fe3_reduction_o2_inhib_above = 1.00e-05

# 9) Ebullition
"""
1.00e+00 CH4(aq) -> 1.00e+00 Tracer2
"""

ebullition_rxn_const = 3.00e-08

ebullition_ch4_halfsat_conc = 4.00e-02
ebullition_ch4_thresh_conc = 1.10e-15

ebullition_ch4_inhib_below = 2.5e-3

# 10) Acetoclastic methanogenesis
"""
1.00e+00 Acetate-  + 1.00e+00 H2O  -> 1.00e+00 CH4(aq)  + 1.00e+00 HCO3-  + 1.00e+00 Tracer
"""

aceto_methano_rxn_const = 1.5e-8

aceto_methano_acetate_halfsat_conc = 4.00e-02
aceto_methano_acetate_thresh_conc = 1.10e-15

aceto_methano_o2_inhib_above = 1.00e-5
aceto_methano_fe3_inhib_above = 1.00e-10

aceto_methano_h_inhib_above = 2.88e-06
aceto_methano_h_inhib_below = 2.88e-06

# 11) Methylotrophic methanogenesis
"""
1.00e+00 CH3OH + 1.00e+00 H2(aq) -> 1.00e+00 CH4(aq) + 1.00e+00 H2O
"""

methyl_methano_rxn_const = 9.1e-06  # arbitrary number

methyl_methano_ch3oh_halfsat_conc = 1.00e-01
methyl_methano_ch3oh_thresh_conc = 1.10e-15

methyl_methano_h2_halfsat_conc = 1.00e-01
methyl_methano_h2_thresh_conc = 1.10e-15

methyl_methano_o2_inhib_above = 1.00e-05

############################
# Initial conditions
############################

# CONSTRAINT initial - Main initial conditions for the system
"""
Initial conditions for the sediment system
"""
# Immobile species concentrations
initial_cellulose = 8.00e+03
initial_hrimm = 1.00e-20

# Aqueous species concentrations (M - molarity, T - total, P - pH, G - gas equilibrium)
initial_dom1 = 5.00e-02  # T
initial_h_ph = 5.0  # P (pH)
initial_o2_aq = 2.0e-4  # T O2(g) - equilibrium with O2 gas
initial_hco3 = 4.00e-6  # T CO2(g) - equilibrium with CO2 gas
initial_fe3_conc = 5.3e-5  # M Fe(OH)3 - equilibrium with Fe(OH)3 mineral
initial_fe2_conc = 2.4e-3  # T Fe(OH)2
initial_nh4 = 1.00e-15  # T
initial_tracer = 1.00e-15  # T
initial_tracer2 = 1.00e-15  # T
initial_tracer3 = 1.00e-15  # T
initial_ch4_aq = 1.00e-15  # T
initial_acetate = 1.00e-15  # T
initial_h2_aq = 1.00e-15  # T
initial_mg = 5.00e-04  # T
initial_ca = 5.00e-04  # T
initial_na = 2.00e-03  # T
initial_k = 2.00e-05  # T
initial_ch3oh = 2.00e-15  # T

# Mineral volume fractions and surface areas
initial_fe_oh3_volfrac = 9.6e-6  # volume fraction
initial_fe_oh3_surface = 1.0e2   # m^2/m^3

initial_fe_oh2_volfrac = 7.2e-1  # volume fraction  
initial_fe_oh2_surface = 1.0e2   # m^2/m^3

initial_rock_volfrac = 0.5       # volume fraction
initial_rock_surface = 5.0e3     # m^2/m^3

# CONSTRAINT sed_air_interface - Sediment-air interface conditions
"""
Conditions at the sediment-air interface with atmospheric gas equilibrium
"""
sedair_dom1 = 5.00e-02  # T
sedair_h_ph = 5.0  # P (pH)
sedair_o2_aq = 0.2  # G O2(g) - gas equilibrium (atmospheric)
sedair_hco3 = 0.01  # G CO2(g) - gas equilibrium
sedair_fe3_conc = 0.37e-15  # M Fe(OH)3
sedair_fe2_conc = 0.37e-15  # M Fe(OH)2
sedair_nh4 = 1.00e-15  # T
sedair_tracer = 1.00e-15  # T
sedair_tracer2 = 1.00e-15  # T
sedair_tracer3 = 1.00e-15  # T
sedair_ch4_aq = 1.00e-15  # T
sedair_acetate = 1.00e-15  # T
sedair_h2_aq = 1.00e-15  # T
sedair_mg = 5.00e-15  # T
sedair_ca = 5.00e-15  # T
sedair_na = 2.00e-15  # T
sedair_k = 2.00e-15  # T
sedair_ch3oh = 2.00e-15  # T

# Mineral conditions at sediment-air interface
sedair_fe_oh3_volfrac = 0.875e-3  # volume fraction
sedair_fe_oh3_surface = 1.0e2     # m^2/m^3
sedair_fe_oh2_volfrac = 0.0e-20   # volume fraction
sedair_fe_oh2_surface = 1.0e2     # m^2/m^3
sedair_rock_volfrac = 0.5         # volume fraction
sedair_rock_surface = 5.0e3       # m^2/m^3

# CONSTRAINT drain - Drainage boundary conditions
"""
Low concentration conditions for drainage boundary
"""
drain_dom1 = 1.00e-15  # T
drain_h_ph = 7.0  # P (pH - neutral)
drain_o2_aq = 1.00e-15  # T O2(g)
drain_hco3 = 1.00e-15  # T CO2(g)
drain_fe3_conc = 1.00e-15  # T Fe(OH)3
drain_fe2_conc = 1.00e-15  # T
drain_nh4 = 1.00e-15  # T
drain_tracer = 1.00e-15  # T
drain_tracer2 = 1.00e-15  # T
drain_tracer3 = 1.00e-15  # T
drain_ch4_aq = 1.00e-15  # T
drain_acetate = 1.00e-15  # T
drain_h2_aq = 1.00e-15  # T
drain_mg = 1.00e-15  # T
drain_ca = 1.00e-15  # T
drain_na = 1.00e-15  # T
drain_k = 1.00e-15  # T
drain_ch3oh = 1.00e-15  # T

# CONSTRAINT recharge - Recharge/rainfall boundary conditions
"""
Very low concentration conditions for clean recharge water
"""
recharge_dom1 = 5.00e-20  # T
recharge_h_ph = 7.0  # P (pH - neutral)
recharge_o2_aq = 1.00e-15  # T O2(g)
recharge_hco3 = 1.00e-15  # T CO2(g)
recharge_fe3_conc = 1.00e-15  # T Fe(OH)3
recharge_fe2_conc = 0.37e-23  # T
recharge_nh4 = 1.00e-25  # T
recharge_tracer = 1.00e-25  # T
recharge_tracer2 = 1.00e-25  # T
recharge_tracer3 = 1.00e-25  # T
recharge_ch4_aq = 1.00e-25  # T
recharge_acetate = 1.00e-25  # T
recharge_h2_aq = 1.00e-25  # T
recharge_mg = 5.00e-25  # T
recharge_ca = 5.00e-25  # T
recharge_na = 2.00e-25  # T
recharge_k = 2.00e-25  # T
recharge_ch3oh = 2.00e-25  # T

############################
# Condition of the biomass
############################
"""
porosity = 0.97d0
SOIL_COMPRESSIBILITY 1.d-07              
SOIL_REFERENCE_PRESSURE 201325.d0
ROCK_DENSITY 2650.0d0
SPECIFIC_HEAT 830.0d0
THERMAL_CONDUCTIVITY_DRY 0.12037926674717922d0
THERMAL_CONDUCTIVITY_WET 1.6082691464310437d0


  PERMEABILITY
    PERM_X 6.5870260083342112d-013
    PERM_Y 6.5870260083342112d-013
    PERM_Z 9.5870260083342112d-14
    
"""

############################
# Time output
############################

# TIME block parameters
"""
Simulation time control parameters
"""
final_time = 100  # d (days)
initial_timestep_size = 1.0  # h (hours) 
maximum_timestep_size = 12.0  # h (hours)

# OUTPUT block parameters  
"""
Output frequency and format control
"""
screen_output_periodic = 10000  # output every 10000 timesteps to screen
file_output_periodic_time = 1.0  # d (days) - output every day
output_format = "TECPLOT POINT"  # output format

# Output variables control
"""
Variables to include/exclude from output
"""
include_flow_variables = False  # NO_FLOW_VARIABLES
include_energy_variables = False  # NO_ENERGY_VARIABLES
include_liquid_saturation = True  # LIQUID_SATURATION

# Observation file settings
"""
Observation point output settings  
"""
obs_output_periodic_time = 1.0  # d (days) - observation output frequency

# Flow condition time schedules
"""
Time schedules for different flow conditions
"""
# Drain condition - constant rate
drain_start_time = 0.0  # d
drain_end_time = 100.0  # d
drain_rate = -5.0e-09  # m^3/s (constant drainage rate)

# Recharge condition - single rain event
recharge_no_rain_start = 0.0   # d
recharge_no_rain_until = 20.0  # d
recharge_rain_start = 20.5     # d  
recharge_rain_end = 21.5       # d
recharge_rain_rate = 2.0e-5    # m^3/s
recharge_no_rain_after = 100.0 # d