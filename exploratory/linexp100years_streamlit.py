### imports
import math
import streamlit as st
from astropy import units as u
from astropy.constants import R
from astropy import constants as const

### title for streamlit app
st.title("Linear vs Exponential Decay Calculator (100 years)")

### function for pct calculations (linear)
def calc_linear(pct_100y, start_carbon):
    # Percentage fractions per time unit
    pct_per_year = pct_100y / 100
    pct_per_month = pct_per_year / 12
    pct_per_week = pct_per_year / 52

    # Mass of carbon released monthly
    carbon_mg_qty = start_carbon * u.mg
    carbon_mg_month = carbon_mg_qty * pct_per_month

    # Convert carbon mass to moles
    moles_C = carbon_mg_month.to(u.g) / (12.011 * u.g / u.mol)

    # Ideal gas law at STP: V = n R T / P
    T0 = 273.15 * u.K
    P0 = 1 * const.atm
    volume = (moles_C * R * T0 / P0).to(u.L)
    ml_CO2 = volume.to(u.mL)

    # Compute concentration
    fraction = (ml_CO2 / (40 * u.mL)).value
    ppm_in_1m3 = fraction * 1e6

    return {
        'pct_per_year': pct_per_year,
        'pct_per_month': pct_per_month,
        'pct_per_week': pct_per_week,
        'carbon_mg': carbon_mg_month.value,
        'moles_C': moles_C.value,
        'liters_CO2': volume.value,
        'ml_CO2': ml_CO2.value,
        'pct_40ml': fraction * 100,
        'ppm_in_1m3': ppm_in_1m3,
    }

### function for pct calculations (exponential)
def calc_exponential(pct_100y, start_carbon):
    # Decay constant based on 100 years = 1200 months
    decay_exp = math.log((100 - pct_100y) / 100) / -1200
    pct_first_year = 100 - (100 * math.exp(-decay_exp * 12))
    pct_first_month = 100 - (100 * math.exp(-decay_exp * 1))
    pct_first_week = 100 - (100 * math.exp(-decay_exp * 0.23))

    carbon_mg_qty = start_carbon * u.mg
    carbon_mg_month = carbon_mg_qty * (pct_first_month / 100)
    carbon_mg_week = carbon_mg_qty * (pct_first_week / 100)

    moles_C = carbon_mg_month.to(u.g) / (12.011 * u.g / u.mol)
    T0 = 273.15 * u.K
    P0 = 1 * const.atm
    volume = (moles_C * R * T0 / P0).to(u.L)
    ml_CO2 = volume.to(u.mL)

    fraction = (ml_CO2 / (40 * u.mL)).value
    ppm_in_1m3 = fraction * 1e6

    return {
        'pct_first_year': pct_first_year,
        'pct_first_month': pct_first_month,
        'pct_first_week': pct_first_week,
        'carbon_mg_month': carbon_mg_month.value,
        'carbon_mg_week': carbon_mg_week.value,
        'moles_C': moles_C.value,
        'liters_CO2': volume.value,
        'ml_CO2': ml_CO2.value,
        'pct_40ml': fraction * 100,
        'ppm_in_1m3': ppm_in_1m3,
    }

# Streamlit inputs
pct_100y = st.number_input("% released in 100 y:", min_value=0.0, max_value=100.0, value=0.01, step=0.01)
biomass = st.number_input("Starting mg biomass:", min_value=0.0, value=1750.0, step=0.1)
start_carbon = st.number_input("Starting mg carbon:", min_value=0.0, value=875.0, step=0.1)
decay_model = st.selectbox("Decay model:", ["Linear", "Exponential"])

if st.button("Compute"):
    if decay_model == "Linear":
        res = calc_linear(pct_100y, start_carbon)
        st.markdown(f"""
        **Assuming {pct_100y:.2f}% released in 100 years**

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
        res = calc_exponential(pct_100y, start_carbon)
        st.markdown(f"""
        **Assuming {pct_100y:.2f}% released in 100 years**

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
