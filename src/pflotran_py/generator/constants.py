"""Literature-sourced defaults for generated PFLOTRAN decks.

Runtime pipeline knobs (paths, species maps) live in ``pflotran_py.config``.
Post-processing physics (diffusion, viscosity) lives in
``pflotran_py.analysis.constants``. This module is the deck-generator
counterpart: kinetic and inhibition numbers that appear in ``.in`` files,
with the citation next to the value rather than inlined at every CLI flag.
"""

# ═════════════════════════════════════════════════════════════════════
# Sealed-bottle geometry and run setup
# ═════════════════════════════════════════════════════════════════════
#
# These match the constants the measurement pipeline uses to turn headspace
# concentrations into moles (saltyBiomass, experiments/analysis/incubations/
# constants.py). If those change, these must change with them, or the model and
# the measurement will be describing differently-sized bottles.

VIAL_VOLUME_L = 0.125  # total internal volume of the vial [L]
HEADSPACE_VOLUME_L = 0.100  # gas volume above the liquid [L]
LIQUID_VOLUME_L = VIAL_VOLUME_L - HEADSPACE_VOLUME_L  # brine + biomass [L]

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
# Water-activity inhibition (AWINHIBIT sandboxes)
# ═════════════════════════════════════════════════════════════════════
#
# ONE_MINUS_AW uses a_crit as the water activity where that pathway's rate
# hits zero: f = max(0, (a_w - a_crit) / (1 - a_crit)).
#
# Higher a_crit = more salt-sensitive (rate hits zero at wetter a_w).
#
# Ordering: acetoclastic > methylotrophic > hydrogenotrophic.
# Acetoclasts fail first under salt (Oren 1999, 2011). The other two are
# spaced 0.05 apart so the three pathways peel off across the bottle a_w
# range rather than sharing one cliff. Classic hypersaline surveys
# (McGenity 2010) find methylotrophic methanogens persist to still lower
# a_w than hydrogenotrophs; putting methyl in the middle is conservative
# relative to that literature, not a methane fit.
#
# The numbers themselves are round values on the Exp003/Exp004 meter
# span (0.824–1.000). The hydrogenotrophic floor sits just below the
# driest bottle so ONE_MINUS_AW stays nonzero across the incubations.
#
# Refs:
#   Oren, A. (1999). Bioenergetic aspects of halophilism.
#     Microbiol. Mol. Biol. Rev. 63:334–348.
#     https://doi.org/10.1128/MMBR.63.2.334-348.1999
#   Oren, A. (2011). Thermodynamic limits to microbial life at high salt
#     concentrations. Environ. Microbiol. 13:1908–1923.
#     https://doi.org/10.1111/j.1462-2920.2010.02365.x
#   McGenity, T.J. (2010). Methanogens and methanogenesis in hypersaline
#     environments. In Timmis (ed.), Handbook of Hydrocarbon and Lipid
#     Microbiology. Springer. https://doi.org/10.1007/978-3-540-77587-4_52

AW_INHIBITION_TYPE = "ONE_MINUS_AW"

AW_CRIT_HYDROGENOTROPHIC = 0.80
AW_CRIT_METHYLOTROPHIC = 0.85
AW_CRIT_ACETOCLASTIC = 0.90

# Shared CLI / constructor alias: hydrogenotrophic is the fallback when a
# pathway-specific a_crit is omitted.
AW_THRESHOLD = AW_CRIT_HYDROGENOTROPHIC
