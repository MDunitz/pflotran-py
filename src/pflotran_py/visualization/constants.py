"""Referenced physics/empirical constants for the visualization physics layer.

Separated from ``config.py`` (runtime pipeline settings: species map,
temperature, output paths) because these are fixed physical/empirical values
with literature provenance rather than tunable configuration. All values carry
``astropy.units`` so their dimensions travel with them.
"""

from astropy import units as u

# ═════════════════════════════════════════════════════════════════════
# Molecular diffusion coefficients
# ═════════════════════════════════════════════════════════════════════

# Molecular diffusion coefficients in water at 25 °C.
# Ref: Boudreau, B.P. (1997). Diagenetic Models and Their Implementation.
#      Springer. Table 4.3.
# NOTE: infinite-dilution freshwater values. In hypersaline brine the true
# tracer D is lower (higher viscosity, ion pairing); see project issue on
# brine-corrected diffusion. The Stokes-Einstein correction in physics.py
# corrects temperature only.
DIFFUSION_COEFFICIENTS_25C = {
    "CO2": 1.91e-9 * u.m**2 / u.s,
    "CH4": 1.49e-9 * u.m**2 / u.s,
    "SO4": 1.07e-9 * u.m**2 / u.s,
    "O2": 2.10e-9 * u.m**2 / u.s,
    "HCO3": 1.18e-9 * u.m**2 / u.s,
    "H2S": 2.12e-9 * u.m**2 / u.s,
    "Fe2": 0.72e-9 * u.m**2 / u.s,  # Fe²⁺
    "Cl": 2.03e-9 * u.m**2 / u.s,
    "Na": 1.33e-9 * u.m**2 / u.s,
}

# ═════════════════════════════════════════════════════════════════════
# Water viscosity — Vogel equation
# ═════════════════════════════════════════════════════════════════════

# Vogel-equation coefficients for the dynamic viscosity of pure liquid water,
# used by the Stokes-Einstein temperature correction in physics.py:
#     μ(T) = A · 10 ^ ( B / (T − C) )        [Pa·s], T absolute (K)
# Ref: standard Vogel water-viscosity correlation, e.g. Engineers Edge,
#      "Water - Density Viscosity Specific Weight"; also tabulated by the
#      Dortmund Data Bank (DDBST). Valid ~10–100 °C, deviation < ~0.8 %.
# The equivalent Celsius form is μ = A · 10 ^ (247.8 / (t_c + 133.15)), where
# 133.15 = 273.15 − C is the divergence temperature shifted into °C.
WATER_VISCOSITY_VOGEL_A = 2.414e-5 * u.Pa * u.s
WATER_VISCOSITY_VOGEL_B = 247.8 * u.K
WATER_VISCOSITY_VOGEL_C = 140.0 * u.K  # divergence temperature

# ═════════════════════════════════════════════════════════════════════
# Unit-conversion factors
# ═════════════════════════════════════════════════════════════════════

# Gradients are stored in M/m = mol/(L·m); Fick's law needs mol/(m³·m).
# Derive the factor with astropy instead of hard-coding 1000.
M_PER_M_TO_SI = (1.0 * u.mol / u.L / u.m).to(u.mol / u.m**3 / u.m).value  # 1000.0
