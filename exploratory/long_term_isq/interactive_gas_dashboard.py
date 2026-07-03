"""
Interactive Gas Concentration Visualization Dashboard

This script creates an interactive Bokeh application that allows users to:
1. Select projection time period
2. Choose between methanogen, spirulina, or mix datasets
3. Toggle between treatment-organized or gas-organized plots
4. Update plots dynamically

Based on gas_concentration_visualization.ipynb
"""

import pandas as pd
import numpy as np
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    Select,
    Slider,
    Button,
    Div,
    HoverTool,
    RadioButtonGroup,
    CheckboxGroup,
)
from bokeh.palettes import Category10, Dark2
from scipy.optimize import curve_fit
import warnings
import os
from astropy import units as u
from astropy import constants as const

warnings.filterwarnings("ignore")

# Global variables to store data and plots
current_data = None
plots_container = None

# Define molar masses for unit conversion
molar_masses = {
    "CH4": 16.05 * u.g / u.mol,  # Methane
    "CO2": 44.01 * u.g / u.mol,  # Carbon dioxide
    "N2O": 46.01 * u.g / u.mol,  # Nitrous oxide
}

# Color palettes
gas_colors = {"CO2": "#1f77b4", "N2O": "#ff7f0e", "CH4": "#2ca02c"}
treatment_colors = Category10[10]  # Will use for treatments when grouped by gas


def parse_experiment_id(sample_id):
    """
    Parse experiment ID from sample ID format: exp_XX_YY_ZZ
    Returns experiment condition as exp_XX_YY
    """
    if pd.isna(sample_id) or not isinstance(sample_id, str):
        return "unknown"

    parts = sample_id.split("_")
    if len(parts) >= 4:
        return f"{parts[0]}_{parts[1]}_{parts[2]}"
    return sample_id


def create_experiment_label(experiment_id, water_activity):
    """
    Create experiment label using water activity if available, otherwise use experiment ID
    """
    if pd.isna(water_activity) or water_activity == "unknown":
        return experiment_id
    else:
        try:
            wa_value = float(water_activity)
            return f"water_activity = {wa_value:.3f}"
        except (ValueError, TypeError):
            return experiment_id


def ppm_to_mol_per_L(
    concentration_ppm, gas_type, temperature=298.15 * u.K, pressure=101325 * u.Pa
):
    """
    Convert ppm (parts per million by volume) to mol/L using astropy units
    """
    if pd.isna(concentration_ppm) or gas_type not in molar_masses:
        return np.nan

    # Convert ppm to fraction
    mole_fraction = concentration_ppm / 1e6

    # Use ideal gas law: n/V = P*χ/(R*T) where χ is mole fraction
    R = const.R  # Gas constant

    # Calculate molarity (mol/L)
    molarity = (pressure * mole_fraction / (R * temperature)).to(u.mol / u.L)

    return molarity.value


def load_and_process_data(filename):
    """Load and process CSV data"""
    try:
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            return None

        df = pd.read_csv(filename)
        print(f"Loaded {filename}: {len(df)} rows, {len(df.columns)} columns")
        print(f"Columns: {df.columns.tolist()}")
        print(
            f"Sample IDs: {df['Sample ID'].unique()[:10] if 'Sample ID' in df.columns else 'No Sample ID column'}"
        )

        # Convert numeric columns and handle non-numeric values
        df["Calc Concentration (ppm)"] = pd.to_numeric(
            df["Calc Concentration (ppm)"], errors="coerce"
        )
        df["days_since_start"] = pd.to_numeric(df["days_since_start"], errors="coerce")
        df["Water Activity"] = pd.to_numeric(df["Water Activity"], errors="coerce")

        # Drop rows with NaN values in critical columns
        df = df.dropna(subset=["Calc Concentration (ppm)", "days_since_start"])

        # Add experiment condition column
        df["experiment_condition"] = df["Sample ID"].apply(parse_experiment_id)

        # Calculate mean concentration for each condition, time point, and gas
        mean_data = (
            df.groupby(["experiment_condition", "days_since_start", "Gas"])
            .agg(
                {"Calc Concentration (ppm)": ["mean", "std"], "Water Activity": "mean"}
            )
            .reset_index()
        )

        # Flatten column names
        mean_data.columns = [
            "experiment_condition",
            "days_since_start",
            "Gas",
            "mean_conc_ppm",
            "std_conc_ppm",
            "mean_water_activity",
        ]

        # Create display labels using water activity when available
        mean_data["display_label"] = mean_data.apply(
            lambda row: create_experiment_label(
                row["experiment_condition"], row["mean_water_activity"]
            ),
            axis=1,
        )

        print(f"Processed data: {len(mean_data)} rows")
        print(
            f"Experiment conditions: {sorted(mean_data['experiment_condition'].unique())}"
        )
        print(f"Gases: {sorted(mean_data['Gas'].unique())}")

        # Convert to mol/L
        mean_data["mean_conc_molL"] = mean_data.apply(
            lambda row: ppm_to_mol_per_L(row["mean_conc_ppm"], row["Gas"]), axis=1
        )

        mean_data["std_conc_molL"] = mean_data.apply(
            lambda row: (
                ppm_to_mol_per_L(row["std_conc_ppm"], row["Gas"])
                if not pd.isna(row["std_conc_ppm"])
                else np.nan
            ),
            axis=1,
        )

        return mean_data

    except Exception as e:
        print(f"Error loading {filename}: {e}")
        import traceback

        traceback.print_exc()
        return None


# Curve fitting functions
def exponential_growth(x, a, b, c):
    """Exponential growth: a * exp(b * x) + c"""
    return a * np.exp(b * x) + c


def linear_model(x, a, b):
    """Linear model: a * x + b"""
    return a * x + b


def logistic_growth(x, a, b, c, d):
    """Logistic growth: a / (1 + exp(-b * (x - c))) + d"""
    return a / (1 + np.exp(-b * (x - c))) + d


def fit_best_curve(x_data, y_data):
    """Try different curve fitting models and return the best fit"""
    models = {
        "linear": (linear_model, [1, 1]),
        "exponential": (exponential_growth, [1, 0.1, 1]),
        "logistic": (logistic_growth, [1000, 0.1, 20, 0]),
    }

    best_model = None
    best_params = None
    best_r2 = -np.inf
    best_name = "linear"

    for name, (func, p0) in models.items():
        try:
            popt, _ = curve_fit(func, x_data, y_data, p0=p0, maxfev=1000)

            # Calculate R-squared
            y_pred = func(x_data, *popt)
            ss_res = np.sum((y_data - y_pred) ** 2)
            ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

            if r2 > best_r2:
                best_model = func
                best_params = popt
                best_r2 = r2
                best_name = name

        except Exception:
            continue

    return best_model, best_params, best_r2, best_name


def create_plots_by_treatment(data, projection_days):
    """Create plots organized by treatment (one plot per treatment)"""
    plots = []

    if data is None or len(data) == 0:
        return plots

    unique_conditions = sorted(data["experiment_condition"].unique())
    gases = sorted(data["Gas"].unique())

    # Take first 6 conditions for display
    selected_conditions = unique_conditions[:6]

    # Use large size as default
    plot_width = 1200
    plot_height = 700

    for condition in selected_conditions:
        # Get display label for this condition
        condition_data = data[data["experiment_condition"] == condition]
        display_label = (
            condition_data["display_label"].iloc[0]
            if len(condition_data) > 0
            else condition
        )

        # Create figure for this condition
        p = figure(
            title=f"Gas Concentrations - {display_label}",
            x_axis_label="Days since start",
            y_axis_label="Concentration (mol/L)",
            width=plot_width,
            height=plot_height,
            tools="pan,wheel_zoom,box_zoom,reset,save",
            toolbar_location="above",
        )

        # Increase axis label font sizes
        p.xaxis.axis_label_text_font_size = "16pt"
        p.yaxis.axis_label_text_font_size = "16pt"
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = "18pt"

        # Add hover tool
        hover = HoverTool(
            tooltips=[
                ("Treatment", f"{display_label}"),
                ("Gas", "@gas"),
                ("Day", "@x"),
                ("Concentration", "@y{0.000000} mol/L"),
                ("Model", "@model"),
            ]
        )
        p.add_tools(hover)

        # Filter data for this condition
        condition_data = data[data["experiment_condition"] == condition]

        # Plot each gas
        for gas in gases:
            if gas in gas_colors:
                gas_data = condition_data[condition_data["Gas"] == gas]

                if len(gas_data) > 0:
                    x_data = gas_data["days_since_start"].values
                    y_data = gas_data["mean_conc_molL"].values
                    y_err = gas_data["std_conc_molL"].values

                    # Remove NaN values
                    valid_mask = ~np.isnan(y_data)
                    x_data = x_data[valid_mask]
                    y_data = y_data[valid_mask]
                    y_err = y_err[valid_mask] if not np.isnan(y_err).all() else None

                    if len(x_data) > 0:
                        # Plot observed data points
                        p.circle(
                            x_data,
                            y_data,
                            size=8,
                            color=gas_colors[gas],
                            alpha=0.8,
                            legend_label=f"{gas}",
                        )

                        # Add error bars if available
                        if y_err is not None and not np.isnan(y_err).all():
                            p.segment(
                                x_data,
                                y_data - y_err,
                                x_data,
                                y_data + y_err,
                                line_color=gas_colors[gas],
                                line_alpha=0.5,
                            )

                        # Fit curve and project
                        if len(x_data) >= 3:
                            try:
                                best_model, best_params, r2, model_name = (
                                    fit_best_curve(x_data, y_data)
                                )

                                if best_model is not None:
                                    # Generate projection
                                    y_projection = best_model(
                                        projection_days, *best_params
                                    )
                                    y_projection = np.maximum(y_projection, 0)

                                    # Plot projection (without legend)
                                    p.line(
                                        projection_days,
                                        y_projection,
                                        line_color=gas_colors[gas],
                                        line_dash="dashed",
                                        line_width=2,
                                        alpha=0.7,
                                    )
                            except Exception:
                                pass

        # Customize legend
        p.legend.location = "top_left"
        p.legend.click_policy = "hide"
        p.legend.label_text_font_size = "14pt"

        plots.append(p)

    return plots


def create_plots_by_gas(data, projection_days):
    """Create plots organized by gas (one plot per gas)"""
    plots = []

    if data is None or len(data) == 0:
        return plots

    unique_conditions = sorted(data["experiment_condition"].unique())
    gases = sorted(data["Gas"].unique())

    # Use large size as default
    plot_width = 1400
    plot_height = 800

    for gas in gases:
        if gas in gas_colors:
            # Create figure for this gas
            p = figure(
                title=f"{gas} Concentrations - All Treatments",
                x_axis_label="Days since start",
                y_axis_label="Concentration (mol/L)",
                width=plot_width,
                height=plot_height,
                tools="pan,wheel_zoom,box_zoom,reset,save",
                toolbar_location="above",
            )

            # Increase axis label font sizes
            p.xaxis.axis_label_text_font_size = "16pt"
            p.yaxis.axis_label_text_font_size = "16pt"
            p.xaxis.major_label_text_font_size = "14pt"
            p.yaxis.major_label_text_font_size = "14pt"
            p.title.text_font_size = "18pt"

            # Add hover tool
            hover = HoverTool(
                tooltips=[
                    ("Treatment", "@treatment"),
                    ("Day", "@x"),
                    ("Concentration", "@y{0.000000} mol/L"),
                    ("Model", "@model"),
                ]
            )
            p.add_tools(hover)

            # Plot each treatment for this gas
            for i, condition in enumerate(
                unique_conditions[:10]
            ):  # Limit to 10 for visibility
                condition_data = data[
                    (data["experiment_condition"] == condition) & (data["Gas"] == gas)
                ]

                if len(condition_data) > 0:
                    # Get display label for this condition
                    display_label = condition_data["display_label"].iloc[0]

                    x_data = condition_data["days_since_start"].values
                    y_data = condition_data["mean_conc_molL"].values
                    y_err = condition_data["std_conc_molL"].values

                    # Remove NaN values
                    valid_mask = ~np.isnan(y_data)
                    x_data = x_data[valid_mask]
                    y_data = y_data[valid_mask]
                    y_err = y_err[valid_mask] if not np.isnan(y_err).all() else None

                    if len(x_data) > 0:
                        color = treatment_colors[i % len(treatment_colors)]

                        # Plot observed data points
                        p.circle(
                            x_data,
                            y_data,
                            size=8,
                            color=color,
                            alpha=0.8,
                            legend_label=f"{display_label} (observed)",
                        )

                        # Add error bars if available
                        if y_err is not None and not np.isnan(y_err).all():
                            p.segment(
                                x_data,
                                y_data - y_err,
                                x_data,
                                y_data + y_err,
                                line_color=color,
                                line_alpha=0.5,
                            )

                        # Fit curve and project
                        if len(x_data) >= 3:
                            try:
                                best_model, best_params, r2, model_name = (
                                    fit_best_curve(x_data, y_data)
                                )

                                if best_model is not None:
                                    # Generate projection
                                    y_projection = best_model(
                                        projection_days, *best_params
                                    )
                                    y_projection = np.maximum(y_projection, 0)

                                    # Plot projection (without legend)
                                    p.line(
                                        projection_days,
                                        y_projection,
                                        line_color=color,
                                        line_dash="dashed",
                                        line_width=2,
                                        alpha=0.7,
                                    )
                            except Exception:
                                pass

            # Customize legend
            p.legend.location = "top_left"
            p.legend.click_policy = "hide"
            p.legend.label_text_font_size = "14pt"

            plots.append(p)

    return plots


def update_plots():
    """Update plots based on current widget selections"""
    global current_data, plots_container

    # Get widget values
    dataset = dataset_select.value
    projection_days_max = projection_slider.value
    plot_organization = plot_org_radio.active  # 0 = by treatment, 1 = by gas

    # Load data
    filename = f"{dataset}_data.csv"
    print(f"Attempting to load: {filename}")
    current_data = load_and_process_data(filename)

    if current_data is None or len(current_data) == 0:
        # Show error message
        status_div.text = f"<div style='color: red;'>Error: Could not load {filename}. Please check if the file exists in the current directory.</div>"
        plots_container.children = []
        return

    # Create projection time points
    number_of_points = max(120, projection_days_max * 2)
    projection_days = np.linspace(0, projection_days_max, number_of_points)

    # Create plots based on organization choice
    if plot_organization == 0:  # By treatment
        plots = create_plots_by_treatment(current_data, projection_days)
        # Arrange all plots in a single column (6 plots stacked vertically)
        plots_container.children = plots
    else:  # By gas
        plots = create_plots_by_gas(current_data, projection_days)
        # Arrange vertically for gas plots
        plots_container.children = plots

    # Update status
    unique_conditions = len(current_data["experiment_condition"].unique())
    unique_gases = len(current_data["Gas"].unique())

    # Count how many experiments have water activity labels vs experiment IDs
    water_activity_labels = (
        current_data["display_label"].str.contains("water_activity").sum()
    )
    exp_id_labels = unique_conditions - (water_activity_labels > 0)

    status_div.text = f"""
    <div style='color: green;'>
    <b>Data loaded successfully!</b><br>
    Dataset: {dataset}<br>
    Treatments: {unique_conditions}<br>
    Gases: {unique_gases}<br>
    Projection period: {projection_days_max} days<br>
    Organization: {'By treatment' if plot_organization == 0 else 'By gas'}<br>
    <i>Note: Experiments with known water activity are labeled as "water_activity = [value]"</i>
    </div>
    """


# Create widgets
dataset_select = Select(
    title="Select Dataset:",
    value="methanogen",
    options=[
        ("methanogen", "Methanogen Data"),
        ("spirulina", "Spirulina Data"),
        ("mix", "Mix Data"),
    ],
)

projection_slider = Slider(
    title="Projection Days:", start=30, end=3650, value=60, step=10
)

plot_org_radio = RadioButtonGroup(
    labels=["Organize by Treatment", "Organize by Gas"], active=0
)

update_button = Button(label="Update Plots", button_type="primary")

status_div = Div(
    text="<div>Select parameters and click 'Update Plots' to begin.</div>",
    width=600,
    height=100,
)

# Container for plots
plots_container = column()

# Connect button to update function
update_button.on_click(update_plots)

# Create layout
header = Div(
    text="""
<h1>Interactive Gas Concentration Visualization Dashboard</h1>
<p>This dashboard allows you to explore gas concentration data with different visualization options.</p>
<p><strong>Features:</strong></p>
<ul>
<li><strong>Dataset Selection:</strong> Choose between methanogen, spirulina, or mix data</li>
<li><strong>Projection Control:</strong> Adjust the number of days for future projections (30-120 days)</li>
<li><strong>Organization Options:</strong> View data organized by treatment or by gas type</li>
<li><strong>Large Plots:</strong> All plots are displayed in large, easy-to-read format</li>
</ul>
""",
    width=800,
)

controls = column(
    dataset_select,
    projection_slider,
    plot_org_radio,
    update_button,
    status_div,
    width=350,
)

main_layout = column(header, row(controls, plots_container), sizing_mode="scale_width")

# Add to document
curdoc().add_root(main_layout)
curdoc().title = "Gas Concentration Dashboard"

# Initial load with default parameters
update_plots()
