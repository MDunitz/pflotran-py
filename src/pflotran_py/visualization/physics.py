"""Concentration-gradient and diffusive-flux computation.

Physics constants use ``astropy.units`` so the values carry dimensions and the
M -> mol/m^3 conversion factor is derived rather than hard-coded. The
per-cell arithmetic runs on plain-float DataFrame columns (astropy Quantities
do not vectorize cleanly into pandas columns); units are re-attached at the
column-name level via ``columns.GRADIENT_UNITS`` / ``FLUX_UNITS``.
"""

import numpy as np
from astropy import units as u

from .columns import gradient_col, flux_col

# ═════════════════════════════════════════════════════════════════════
# Diffusion coefficients
# ═════════════════════════════════════════════════════════════════════

# Molecular diffusion coefficients in water at 25 °C.
# Ref: Boudreau, B.P. (1997). Diagenetic Models and Their Implementation.
#      Springer. Table 4.3.
# NOTE: infinite-dilution freshwater values. In hypersaline brine the true
# tracer D is lower (higher viscosity, ion pairing); see project issue on
# brine-corrected diffusion. Stokes-Einstein below corrects temperature only.
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

# Gradients are stored in M/m = mol/(L·m); Fick's law needs mol/(m³·m).
# Derive the factor with astropy instead of hard-coding 1000.
_M_PER_M_TO_SI = (1.0 * u.mol / u.L / u.m).to(u.mol / u.m**3 / u.m).value  # 1000.0


def stokes_einstein_correction(temperature_c, reference_c=25.0):
    """Temperature correction factor for diffusion coefficients.

    Stokes-Einstein scaling (D = k_B·T / (6·π·μ·r), so D ∝ T/μ):
        D(T) / D(T_ref) = (T / T_ref) × (μ_ref / μ_T)

    Water viscosity empirical fit (Vogel form, Sharqawy et al., 2010):
        μ(t) ≈ 2.414e-5 × 10^(247.8 / (t + 133.15))  [Pa·s], t in °C

    Parameters
    ----------
    temperature_c : float
        Temperature [°C].
    reference_c : float
        Reference temperature [°C] (default 25 °C for the Boudreau values).

    Returns
    -------
    float
        Dimensionless factor: D(T) = D(T_ref) × factor.
    """

    def viscosity(t_c):
        return 2.414e-5 * 10 ** (247.8 / (t_c + 133.15))

    t_k = (temperature_c * u.deg_C).to(u.K, equivalencies=u.temperature()).value
    t_ref_k = (reference_c * u.deg_C).to(u.K, equivalencies=u.temperature()).value
    return (t_k / t_ref_k) * (viscosity(reference_c) / viscosity(temperature_c))


# ═════════════════════════════════════════════════════════════════════
# Gradient calculation
# ═════════════════════════════════════════════════════════════════════


def _compute_point_gradient(
    coords, point, neighbor_values, center_value, neighbor_radius=0.5
):
    """Least-squares concentration gradient at a single point.

    Solve A·g = Δc for the local gradient over N neighbors within radius r:
        A  = [Δx_i, Δy_i, Δz_i, 1]   for i = 1..N
        Δc = [c_i − c_center]
        g  = [∂C/∂x, ∂C/∂y, ∂C/∂z, offset]

    Returns (∂C/∂x, ∂C/∂y, ∂C/∂z) or None if fewer than 3 neighbors.
    """
    distances = np.sqrt(np.sum((coords - point) ** 2, axis=1))
    neighbors_idx = np.where((distances > 0) & (distances <= neighbor_radius))[0]

    if len(neighbors_idx) < 3:
        return None

    A = np.column_stack([coords[neighbors_idx] - point, np.ones(len(neighbors_idx))])
    dc = neighbor_values[neighbors_idx] - center_value
    gradients, _, _, _ = np.linalg.lstsq(A, dc, rcond=None)
    return gradients[:3]


def calculate_gradients(df, species_map):
    """Compute concentration gradients for multiple species.

    Parameters
    ----------
    df : DataFrame
        Must contain 'X [m]', 'Y [m]', 'Z [m]', 'Time Index', and the
        concentration columns referenced in species_map.
    species_map : dict
        Maps short name to DataFrame column name, e.g.
        {'CO2': 'Free CO2(aq) [M]', 'CH4': 'Free CH4(aq) [M]'}.

    Returns
    -------
    DataFrame
        Copy of df with gradient columns added: {species}_grad_{x,y,z,magnitude}.
    """
    result = df.copy()

    for species in species_map:
        for component in ["x", "y", "z", "magnitude"]:
            result[gradient_col(species, component)] = 0.0

    for time_idx in df["Time Index"].unique():
        mask = result["Time Index"] == time_idx
        time_data = result.loc[mask]
        coords = time_data[["X [m]", "Y [m]", "Z [m]"]].values

        for idx, row in time_data.iterrows():
            point = row[["X [m]", "Y [m]", "Z [m]"]].values

            for species, col_name in species_map.items():
                grad = _compute_point_gradient(
                    coords,
                    point,
                    time_data[col_name].values,
                    row[col_name],
                )
                if grad is not None:
                    # Store the negated gradient (points toward lower
                    # concentration), so Fick's J = D·(stored) below is correct.
                    result.loc[idx, gradient_col(species, "x")] = -grad[0]
                    result.loc[idx, gradient_col(species, "y")] = -grad[1]
                    result.loc[idx, gradient_col(species, "z")] = -grad[2]

    for species in species_map:
        result[gradient_col(species, "magnitude")] = np.sqrt(
            result[gradient_col(species, "x")] ** 2
            + result[gradient_col(species, "y")] ** 2
            + result[gradient_col(species, "z")] ** 2
        )

    return result


def convert_to_flux(df, species_list, temperature_c=25.0):
    """Convert stored gradient columns to diffusive flux.

    Fick's First Law:
        J = −D(T)·∇C
    calculate_gradients() stores −∇C, so with that column the sign works out
    as J = +D(T)·(stored). The ×_M_PER_M_TO_SI factor converts the stored
    gradient from M/m = mol/(L·m) to mol/(m³·m) before applying D.

    Parameters
    ----------
    df : DataFrame
        Must already contain gradient columns from calculate_gradients().
    species_list : list of str
        Species short names (keys in DIFFUSION_COEFFICIENTS_25C).
    temperature_c : float
        Simulation temperature [°C] for the Stokes-Einstein correction.

    Returns
    -------
    DataFrame
        df with flux columns added: {species}_flux_{x,y,z,magnitude}.
    """
    temp_factor = stokes_einstein_correction(temperature_c)

    for species in species_list:
        d_t = (DIFFUSION_COEFFICIENTS_25C[species] * temp_factor).to(u.m**2 / u.s).value

        for component in ["x", "y", "z"]:
            g_col = gradient_col(species, component)
            f_col = flux_col(species, component)
            df[f_col] = d_t * df[g_col] * _M_PER_M_TO_SI

        df[flux_col(species, "magnitude")] = np.sqrt(
            df[flux_col(species, "x")] ** 2
            + df[flux_col(species, "y")] ** 2
            + df[flux_col(species, "z")] ** 2
        )

    return df
