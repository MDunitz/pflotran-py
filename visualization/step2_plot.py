#!/usr/bin/env python3
# PFLOTRAN Visualization Script
# Should run automatically after step1_extract.py
# Change  'variables_to_plot' list below to change which variables are visualized.
# After running, do: open multi_variable_concentration_3d.html

########################################################################
# 1) Load in the processed data (CSV) and other imports
########################################################################
from shared_utils import load_data

from utils_plotting import create_multi_variable_plot, create_single_variable_plot

# Re-export for callers that import plotting helpers from step2_plot.
__all__ = ["create_single_variable_plot", "create_multi_variable_plot", "main"]


########################################################################
# 0) Main function
########################################################################
def main(verbose=False):
    print("PFLOTRAN Visualization Starting...\n")

    df = load_data("pflotran_data.pkl")

    # CONFIGURATION: Change these variables to plot different data
    # Modify this list to visualize different variables from your simulation
    variables_to_plot = [
        "Total CH4(aq) [M]",
        "Total Acetate- [M]",
        "Total SO4-- [M]",
        "Gamma H2O",
    ]

    if verbose:
        print(f"Variables to plot: {variables_to_plot}")

    # Create multi-variable plot
    fig = create_multi_variable_plot(df, variables_to_plot)
    if fig is None:
        print("No requested variables found in data; skipping HTML export.")
        return

    # Export to HTML
    output_filename = "multi_variable_concentration_3d.html"
    fig.write_html(output_filename)
    print(f"Multi-variable 3D plot saved as: {output_filename}")
    print("to run, do open multi_variable_concentration_3d.html in terminal")

    # Also create individual plots for the first variable as an example, or if you want just one
    if variables_to_plot and variables_to_plot[0] in df.columns:
        single_var = variables_to_plot[0]
        single_fig = create_single_variable_plot(df, single_var)
        single_filename = f"single_variable_{single_var.replace(' ', '_').replace('[', '').replace(']', '')}.html"
        single_fig.write_html(single_filename)
        print(f"Single variable plot saved as: {single_filename}")

    print("\nVisualization complete! Open html file")


if __name__ == "__main__":
    main()
