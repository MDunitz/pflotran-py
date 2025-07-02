import math
import streamlit as st

st.title("Linear vs Exponential Decay Calculator (100 years)")

decay_model = st.selectbox("Decay model:", ["Linear", "Exponential"])
pct_100y = st.number_input("% released in 100 y:", min_value=0.0, max_value=100.0, value=0.01, step=0.01)
biomass = st.number_input("Starting mg biomass:", min_value=0.0, value=1750.0, step=0.1)
start_carbon = st.number_input("Starting mg carbon:", min_value=0.0, value=875.0, step=0.1)

if st.button("Compute"):
    if decay_model == "Linear":
        pct_per_year  = pct_100y / 100           
        pct_per_month = pct_per_year / 12
        pct_per_week = pct_per_year / 52           

        carbon_mg_month = start_carbon * pct_per_month
        moles_C   = carbon_mg_month / 1000 / 12.011
        liters_CO2 = moles_C * 22.414
        ml_CO2     = liters_CO2 * 1000
        pct_40ml   = ml_CO2 / 40 * 100
        ppm_in_1m3 = pct_40ml * 10000

        st.markdown(
            f"""
            **Working with the assumption that {pct_100y:.2f}% released in 100 years**

            - % lost in 1 yr: `{pct_per_year:.6f} %`
            - % lost in 1 month: `{pct_per_month:.6f} %`
            - % lost in 1 week: `{pct_per_week:.6f} %`
            - Carbon used per month: `{carbon_mg_month:.6f} mg C`
            - Moles C per month: `{moles_C:.6e} mol`
            - CO₂ volume at STP: `{liters_CO2:.6f} L`  (`{ml_CO2:.6f} mL`)
            - `{pct_40ml:.6f} %` of a 40 mL headspace
            - `≈ {ppm_in_1m3:.2f} ppm` in 1 m³ air
            """
        )
    else:
        decay_exp = math.log((100 - pct_100y) / 100) / -1200
        pct_first_year  = 100-(100*math.exp(-decay_exp*12))
        pct_first_month = 100-(100*math.exp(-decay_exp*1))
        pct_first_week = 100-(100*math.exp(-decay_exp*0.23))

        carbon_mg_firstmonth = start_carbon * pct_first_month
        carbon_mg_firstweek = start_carbon * pct_first_week
        moles_C   = carbon_mg_firstmonth / 1000 / 12.011
        liters_CO2 = moles_C * 22.414
        ml_CO2     = liters_CO2 * 1000
        pct_40ml   = ml_CO2 / 40 * 100
        ppm_in_1m3 = pct_40ml * 10000

        st.markdown(
            f"""
            **Working with the assumption that {pct_100y:.2f}% released in 100 years**

            - % lost in first yr: `{pct_first_year:.6f} %`
            - % lost in first month: `{pct_first_month:.6f} %`
            - % lost in first week: `{pct_first_week:.6f} %`
            - Carbon used first month: `{carbon_mg_firstmonth:.6f} mg C`
            - Carbon used first week: `{carbon_mg_firstweek:.6f} mg C`
            - Moles C per month: `{moles_C:.6e} mol`
            - CO₂ volume at STP: `{liters_CO2:.6f} L`  (`{ml_CO2:.6f} mL`)
            - `{pct_40ml:.6f} %` of a 40 mL headspace
            - `≈ {ppm_in_1m3:.2f} ppm` in 1 m³ air
            """
        )