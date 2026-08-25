"""Literature-sourced defaults for generated PFLOTRAN decks.

Runtime pipeline knobs (paths, species maps) live in ``pflotran_py.config``.
Post-processing physics (diffusion, viscosity) lives in
``pflotran_py.analysis.constants``. This module is the deck-generator
counterpart: kinetic and inhibition numbers that appear in ``.in`` files,
with the citation next to the value rather than inlined at every CLI flag.
"""

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
