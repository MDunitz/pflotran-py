### imports
import math
import streamlit as st
from astropy import units as u
from astropy.constants import R
from astropy import constants as const
import constants

### 1) title for streamlit app
st.title("Linear vs Exponential Decay Calculator (100 years)")

### 2) definining a few functions that we'll use for lin AND exp decay calculations
#   2a) linear release rates
def linear_release_rates(pct_total: float, years: float):
    frac_per_year  = pct_total / years / constants.percent_conversion
    frac_per_month = frac_per_year / constants.months_p_year
    frac_per_week  = frac_per_year / constants.weeks_p_year
    return frac_per_year, frac_per_month, frac_per_week

#   2b) exponential release rates
def exponential_release_rates(pct_total: float, years: float):
    # λ so that 1 – e^(–λ·years·12) = pct_total/100
    λ = -math.log(1 - pct_total/constants.percent_conversion) / (years * constants.months_p_year)
    # now project out:
    pct_year  = constants.percent_conversion * (1 - math.exp(-λ * constants.months_p_year))
    pct_month = constants.percent_conversion * (1 - math.exp(-λ * 1))
    pct_week  = constants.percent_conversion * (1 - math.exp(-λ * constants.weeks_p_year))
    return pct_year, pct_month, pct_week

#    2c) convert mg of carbon to moles
def moles_from_carbon_mass(mass_mg: float):
    mass = mass_mg * u.mg
    grams = mass.to(u.g)
    return (grams / (constants.molar_mass_C * u.g / u.mol)).value


#    2d) moles to co2 vol
def co2_volume_from_moles(moles_C: float):
    n = moles_C * u.mol
    T0 = constants.T * u.K
    P0 = 1 * const.atm
    vol_L = (n * R * T0 / P0).to(u.L).value
    vol_mL = (vol_L * u.L).to(u.mL).value
    return vol_L, vol_mL


#    2e) volume to concentration
def concentration_metrics(volume_mL: float, sample_mL: float = constants.sample_mL_default):
    fraction = (volume_mL / sample_mL)
    ppm = fraction * constants.frac_to_ppm
    return fraction * constants.percent_conversion, ppm


###  3) linear and exponential calculations

def calc_linear(pct_total: float, start_C: float, years: float):
    pr_y, pr_m, pr_w = linear_release_rates(pct_total, years)
    mg_month = start_C * pr_m
    moles    = moles_from_carbon_mass(mg_month)
    L, mL    = co2_volume_from_moles(moles)
    pct40, ppm = concentration_metrics(mL)
    return {
        'pct_per_year':  pr_y*constants.percent_conversion,
        'pct_per_month': pr_m*constants.percent_conversion,
        'pct_per_week':  pr_w*constants.percent_conversion,
        'carbon_mg':     mg_month,
        'moles_C':       moles,
        'liters_CO2':    L,
        'ml_CO2':        mL,
        'pct_40ml':      pct40,
        'ppm_in_1m3':    ppm,
    }

def calc_exponential(pct_total: float, start_C: float, years: float):
    p1y, p1m, p1w = exponential_release_rates(pct_total, years)
    mg_month = start_C * (p1m/constants.percent_conversion)
    mg_week  = start_C * (p1w/constants.percent_conversion)
    moles    = moles_from_carbon_mass(mg_month)
    L, mL    = co2_volume_from_moles(moles)
    pct40, ppm = concentration_metrics(mL)
    return {
        'pct_first_year':  p1y,
        'pct_first_month': p1m,
        'pct_first_week':  p1w,
        'carbon_mg_month': mg_month,
        'carbon_mg_week':  mg_week,
        'moles_C':         moles,
        'liters_CO2':      L,
        'ml_CO2':          mL,
        'pct_40ml':        pct40,
        'ppm_in_1m3':      ppm,
    }

# Streamlit inputs
num_years = st.number_input("Years to release:", min_value=1, max_value=10000, value=100, step=1)
pct_total_y = st.number_input("% released in your set years:", min_value=0.0, max_value=100.0, value=0.01, step=0.01)
biomass = st.number_input("Starting mg biomass:", min_value=0.0, value=1750.0, step=0.1)
start_carbon = st.number_input("Starting mg carbon:", min_value=0.0, value=875.0, step=0.1)
decay_model = st.selectbox("Decay model:", ["Linear", "Exponential"])

if st.button("Compute"):
    if decay_model == "Linear":
        res = calc_linear(pct_total_y, start_carbon, num_years)
        st.markdown(f"""
        **Assuming {pct_total_y:.2f}% released in your set years**

        - % lost in 1 yr: `{res['pct_per_year']:.6f} %`
        - % lost in 1 month: `{res['pct_per_month']:.6f} %`
        - % lost in 1 week: `{res['pct_per_week']:.6f} %`
        - Carbon used per month: `{res['carbon_mg']:.6f} mg C`
        - Moles C per month: `{res['moles_C']:.6e} mol`
        - CO₂ volume at STP: `{res['liters_CO2']:.6f} L` (`{res['ml_CO2']:.6f} mL`)
        - `{res['pct_40ml']:.6f} %` of a 40 mL headspace
        - `≈ {res['ppm_in_1m3']:.2f} ppm` in 1 m³ air
        """
        )
    else:
        res = calc_exponential(pct_total_y, start_carbon, num_years)
        st.markdown(f"""
        **Assuming {pct_total_y:.2f}% released in your set years**

        - % lost in first yr: `{res['pct_first_year']:.6f} %`
        - % lost in first month: `{res['pct_first_month']:.6f} %`
        - % lost in first week: `{res['pct_first_week']:.6f} %`
        - Carbon used in first month: `{res['carbon_mg_month']:.6f} mg C`
        - Carbon used in first week: `{res['carbon_mg_week']:.6f} mg C`
        - Moles C per month: `{res['moles_C']:.6e} mol`
        - CO₂ volume at STP: `{res['liters_CO2']:.6f} L` (`{res['ml_CO2']:.6f} mL`)
        - `{res['pct_40ml']:.6f} %` of a 40 mL headspace
        - `≈ {res['ppm_in_1m3']:.2f} ppm` in 1 m³ air
        """
        )
