"""
Shared utilities for PFLOTRAN visualization pipeline.

Functions:
    load_data             — Load processed DataFrame from pickle
    calculate_gradients   — Compute concentration gradients (∇C)
    convert_to_flux       — Convert gradients to diffusive flux (J = -D·∇C)
    column_name           — Generate consistent column names with units
    tooltip_entry         — Generate Bokeh tooltip entry matching column names

Units are defined once and propagated to column names and tooltips
so they cannot drift apart.
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


# ═════════════════════════════════════════════════════════════════════
# Unit definitions — single source of truth
# ═════════════════════════════════════════════════════════════════════

GRADIENT_UNITS = 'M/m'          # mol/L per meter
FLUX_UNITS = 'mol/(m²·s)'      # molar flux

# Molecular diffusion coefficients in water at 25°C [m²/s]
# Ref: Boudreau, B.P. (1997). Diagenetic Models and Their Implementation.
#      Springer. Table 4.3.
#
# Fick's First Law: J = -D · ∇C
#   J  = diffusive flux       [mol/(m²·s)]
#   D  = diffusion coefficient [m²/s]
#   ∇C = concentration gradient [mol/m³ / m]
#
# Note: ∇C from our calculation is in M/m = mol/(L·m).
#       Multiply by 1000 to convert to mol/(m³·m) before applying D.
DIFFUSION_COEFFICIENTS_25C = {
    'CO2':  1.91e-9,
    'CH4':  1.49e-9,
    'SO4':  1.07e-9,
    'O2':   2.10e-9,
    'HCO3': 1.18e-9,
    'H2S':  2.12e-9,
    'Fe2':  0.72e-9,   # Fe²⁺
    'Cl':   2.03e-9,
    'Na':   1.33e-9,
}


def stokes_einstein_correction(temperature_c, reference_c=25.0):
    """Temperature correction for diffusion coefficients.

    Stokes-Einstein equation scaling:
        D(T) / D(T_ref) ≈ (T / T_ref) × (μ_ref / μ_T)

    Simplified using water viscosity empirical fit (Sharqawy et al., 2010):
        μ(T) ≈ 2.414e-5 × 10^(247.8 / (T + 133.15))  [Pa·s]

    Parameters
    ----------
    temperature_c : float
        Temperature in °C.
    reference_c : float
        Reference temperature in °C (default 25°C for Boudreau values).

    Returns
    -------
    float
        Multiplicative correction factor. D(T) = D(25°C) × factor.
    """
    def viscosity(t_c):
        return 2.414e-5 * 10 ** (247.8 / (t_c + 133.15))

    t_k = temperature_c + 273.15
    t_ref_k = reference_c + 273.15
    return (t_k / t_ref_k) * (viscosity(reference_c) / viscosity(temperature_c))


# ═════════════════════════════════════════════════════════════════════
# Column naming — links computation to display
# ═════════════════════════════════════════════════════════════════════

def gradient_col(species, component):
    """Column name for a gradient component.

    Examples:
        gradient_col('CO2', 'x')         -> 'CO2_grad_x'
        gradient_col('CO2', 'magnitude') -> 'CO2_grad_magnitude'
    """
    return f'{species}_grad_{component}'


def flux_col(species, component):
    """Column name for a flux component.

    Examples:
        flux_col('CO2', 'x')         -> 'CO2_flux_x'
        flux_col('CO2', 'magnitude') -> 'CO2_flux_magnitude'
    """
    return f'{species}_flux_{component}'


def gradient_tooltips(species):
    """Bokeh tooltip entries for gradient columns of a species."""
    return [
        (f'{species} |∇C| [{GRADIENT_UNITS}]', f'@{gradient_col(species, "magnitude")}'),
        (f'{species} ∂C/∂x [{GRADIENT_UNITS}]', f'@{gradient_col(species, "x")}'),
        (f'{species} ∂C/∂y [{GRADIENT_UNITS}]', f'@{gradient_col(species, "y")}'),
        (f'{species} ∂C/∂z [{GRADIENT_UNITS}]', f'@{gradient_col(species, "z")}'),
    ]


def flux_tooltips(species):
    """Bokeh tooltip entries for flux columns of a species."""
    return [
        (f'{species} |J| [{FLUX_UNITS}]', f'@{flux_col(species, "magnitude")}'),
        (f'{species} Jx [{FLUX_UNITS}]', f'@{flux_col(species, "x")}'),
        (f'{species} Jy [{FLUX_UNITS}]', f'@{flux_col(species, "y")}'),
        (f'{species} Jz [{FLUX_UNITS}]', f'@{flux_col(species, "z")}'),
    ]


# ═════════════════════════════════════════════════════════════════════
# Data loading
# ═════════════════════════════════════════════════════════════════════

def load_data(filename='pflotran_data.pkl'):
    return pd.read_pickle(filename)


# ═════════════════════════════════════════════════════════════════════
# Gradient calculation
# ═════════════════════════════════════════════════════════════════════

def _compute_point_gradient(coords, point, neighbor_values, center_value, neighbor_radius=0.5):
    """Least-squares gradient at a single point.

    Gradient estimation equation:
        For N neighbors within radius r, solve A·g = Δc
        A = [Δx_i, Δy_i, Δz_i, 1]  for i = 1..N
        Δc = [c_i - c_center]
        g = [∂C/∂x, ∂C/∂y, ∂C/∂z, offset]

    Returns (∂C/∂x, ∂C/∂y, ∂C/∂z) or None if insufficient neighbors.
    """
    distances = np.sqrt(np.sum((coords - point)**2, axis=1))
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
        Must contain 'X [m]', 'Y [m]', 'Z [m]', 'Time Index', and
        the concentration columns referenced in species_map.
    species_map : dict
        Maps short name to DataFrame column name.
        Example: {'CO2': 'Free CO2(aq) [M]', 'CH4': 'Free CH4(aq) [M]'}

    Returns
    -------
    DataFrame
        Copy of df with gradient columns added:
        {species}_grad_x, _grad_y, _grad_z, _grad_magnitude
    """
    result = df.copy()

    # Initialize columns
    for species in species_map:
        for component in ['x', 'y', 'z', 'magnitude']:
            result[gradient_col(species, component)] = 0.0

    for time_idx in df['Time Index'].unique():
        mask = result['Time Index'] == time_idx
        time_data = result.loc[mask]
        coords = time_data[['X [m]', 'Y [m]', 'Z [m]']].values

        for idx, row in time_data.iterrows():
            point = row[['X [m]', 'Y [m]', 'Z [m]']].values

            for species, col_name in species_map.items():
                grad = _compute_point_gradient(
                    coords, point, time_data[col_name].values, row[col_name],
                )
                if grad is not None:
                    # Store negative gradient (convention: positive toward lower concentration)
                    result.loc[idx, gradient_col(species, 'x')] = -grad[0]
                    result.loc[idx, gradient_col(species, 'y')] = -grad[1]
                    result.loc[idx, gradient_col(species, 'z')] = -grad[2]

    # Magnitude
    for species in species_map:
        result[gradient_col(species, 'magnitude')] = np.sqrt(
            result[gradient_col(species, 'x')]**2 +
            result[gradient_col(species, 'y')]**2 +
            result[gradient_col(species, 'z')]**2
        )

    return result


def convert_to_flux(df, species_list, temperature_c=25.0):
    """Convert gradient columns to true diffusive flux.

    Fick's First Law:
        J = -D(T) × ∇C × 1000
        where the 1000 converts M (mol/L) to mol/m³

    Requires gradient columns from calculate_gradients() to already exist.

    Parameters
    ----------
    df : DataFrame
        Must contain gradient columns (species_grad_x, etc.)
    species_list : list of str
        Species short names (must be keys in DIFFUSION_COEFFICIENTS_25C)
    temperature_c : float
        Simulation temperature for Stokes-Einstein correction.

    Returns
    -------
    DataFrame
        df with additional flux columns:
        {species}_flux_x, _flux_y, _flux_z, _flux_magnitude
    """
    temp_factor = stokes_einstein_correction(temperature_c)

    for species in species_list:
        D_25 = DIFFUSION_COEFFICIENTS_25C[species]
        D_T = D_25 * temp_factor  # [m²/s] at simulation temperature

        # J = D × ∇C × 1000
        # gradient is in M/m, multiply by 1000 to get mol/m³/m, then by D
        for component in ['x', 'y', 'z']:
            g_col = gradient_col(species, component)
            f_col = flux_col(species, component)
            df[f_col] = D_T * df[g_col] * 1000.0

        df[flux_col(species, 'magnitude')] = np.sqrt(
            df[flux_col(species, 'x')]**2 +
            df[flux_col(species, 'y')]**2 +
            df[flux_col(species, 'z')]**2
        )

    return df


# ═════════════════════════════════════════════════════════════════════
# Backward compatibility aliases
# ═════════════════════════════════════════════════════════════════════

def calculate_single_species_gradient(coords, point, neighbor_data, current_val, flux_df, idx, species_prefix):
    """DEPRECATED: Use calculate_gradients() instead.

    Kept for backward compatibility with step3_flux.py / step4_plotflux.py.
    Note: writes to columns named {prefix}_flux_* but values are gradients.
    """
    grad = _compute_point_gradient(coords, point, neighbor_data, current_val)
    if grad is not None:
        flux_df.loc[idx, f'{species_prefix}_flux_x'] = -grad[0]
        flux_df.loc[idx, f'{species_prefix}_flux_y'] = -grad[1]
        flux_df.loc[idx, f'{species_prefix}_flux_z'] = -grad[2]


def calculate_concentration_gradients(df, species_names):
    """DEPRECATED: Use calculate_gradients() instead.

    Kept for backward compatibility. Writes to {prefix}_flux_* columns
    using the old naming convention.
    """
    flux_df = df.copy()

    for prefix in species_names:
        for dim in ['x', 'y', 'z', 'magnitude']:
            flux_df[f'{prefix}_flux_{dim}'] = 0.0

    for time_idx in df['Time Index'].unique():
        time_data = flux_df[flux_df['Time Index'] == time_idx].copy()
        coords = time_data[['X [m]', 'Y [m]', 'Z [m]']].values

        for idx, row in time_data.iterrows():
            point = row[['X [m]', 'Y [m]', 'Z [m]']].values
            for prefix, col in species_names.items():
                calculate_single_species_gradient(
                    coords, point, time_data[col].values, row[col],
                    flux_df, idx, prefix,
                )

    for prefix in species_names:
        flux_df[f'{prefix}_flux_magnitude'] = np.sqrt(
            flux_df[f'{prefix}_flux_x']**2 +
            flux_df[f'{prefix}_flux_y']**2 +
            flux_df[f'{prefix}_flux_z']**2
        )

    return flux_df
