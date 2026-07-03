#!/usr/bin/env python3
#### PFLOTRAN Visualization Script
# Should run automatically after step1_extract.py
# Change  'variables_to_plot' list below to change which variables are visualized.
# After running, do: open multi_variable_concentration_3d.html 

import os
import sys
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

########################################################################
###### 1) Load in the processed data (CSV)
########################################################################
def load_data(filename="pflotran_data.pkl"):
    if os.path.exists(filename):
        print(f"Loading data from: {filename}")
        df = pd.read_pickle(filename)
        print(f"Data loaded successfully! Shape: {df.shape}")
        return df
    elif os.path.exists(filename.replace('.pkl', '.csv')):
        csv_filename = filename.replace('.pkl', '.csv')
        print(f"Loading data from CSV: {csv_filename}")
        df = pd.read_csv(csv_filename)
        print(f"Data loaded successfully! Shape: {df.shape}")
        return df
    else:
        print(f"Data file not found: {filename}, run step1_extract.py first")
        return None

########################################################################
###### 2) Show which variables are available
########################################################################
def display_available_variables(df):
    if df is None:
        return
        
    print(f"\nAvailable variables for plotting:")
    variables = [col for col in df.columns if col not in ['Time Index', 'X [m]', 'Y [m]', 'Z [m]']]
    
    for i, var in enumerate(variables):
        # Get some basic statistics
        var_min = df[var].min()
        var_max = df[var].max()
        var_mean = df[var].mean()
        print(f"   {i:2d}: {var}")
        print(f"       Range: {var_min:.2e} to {var_max:.2e} (mean: {var_mean:.2e})")

########################################################################
###### 3) Create a single variable 3D plot with time animation, default for CH4(aq)
########################################################################
def create_single_variable_plot(df, variable):

    # Filter data for this variable
    plot_data = df[['X [m]', 'Y [m]', 'Z [m]', variable, 'Time Index']].copy()
    
    # Get consistent color scale
    var_min = plot_data[variable].min()
    var_max = plot_data[variable].max()
    
    # Handle case where all values are zero or very small
    if var_max <= 1e-15:
        var_max = 1e-15
        var_min = 0
    
    # Create the plot
    fig = px.scatter_3d(
        plot_data,
        x='X [m]',
        y='Y [m]',
        z='Z [m]',
        color=variable,
        animation_frame='Time Index',
        color_continuous_scale='Viridis',
        range_color=[var_min, var_max],
        labels={variable: f'{variable}'},
        title=f"3D Visualization of {variable} Over Time",
        opacity=0.7
    )
    
    # Improve layout
    fig.update_layout(
        scene=dict(
            xaxis_title='X [m]',
            yaxis_title='Y [m]',
            zaxis_title='Z [m]'
        ),
        coloraxis_colorbar=dict(
            tickformat='.2e' if var_max < 1e-6 else '.3f'
        ),
        updatemenus=[{
            'type': 'buttons',
            'showactive': False,
            'buttons': [
                {'label': 'Play', 'method': 'animate', 'args': [None, {'frame': {'duration': 100}}]},
                {'label': 'Pause', 'method': 'animate', 'args': [[None], {'frame': {'duration': 0}, 'mode': 'immediate'}]}
            ]
        }]
    )
    
    return fig

########################################################################
###### 4) Multi-variable 3D plot with time animation
########################################################################
def create_multi_variable_plot(df, variables_to_plot):
    # Filter available variables
    available_vars = [var for var in variables_to_plot if var in df.columns]
    
    if not available_vars:
        print(f"None of the requested variables found in data!")
        print(f"Requested: {variables_to_plot}")
        print(f"Available: {[col for col in df.columns if col not in ['Time Index', 'X [m]', 'Y [m]', 'Z [m]']]}")
        return None
    
    print(f"Creating multi-variable plot for: {available_vars}")
    
    # Create subplot grid
    n_plots = len(available_vars)
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols  # Ceiling division
    
    # Create subplots with 3D scenes
    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        specs=[[{'type': 'scatter3d'} for _ in range(n_cols)] for _ in range(n_rows)],
        subplot_titles=available_vars,
        vertical_spacing=0.15,
        horizontal_spacing=0.05
    )
    
    # Get unique time indices for frames
    time_indices = sorted(df['Time Index'].unique())
    
    # Create frames for animation
    frames = []
    for time_idx in time_indices:
        frame_data = []
        time_data = df[df['Time Index'] == time_idx]
        
        for i, var in enumerate(available_vars):
            row = (i // n_cols) + 1
            col = (i % n_cols) + 1
            
            # Calculate consistent color scale for this variable
            var_min = df[var].min()
            var_max = df[var].max()
            
            # Handle case where all values are zero or very small
            if var_max <= 1e-15:
                var_max = 1e-15
                var_min = 0
            
            # Better colorbar positioning
            colorbar_x = 1.02 if col == n_cols else 0.48
            colorbar_y = 0.75 if row == 1 else 0.25
            
            scatter = go.Scatter3d(
                x=time_data['X [m]'],
                y=time_data['Y [m]'],
                z=time_data['Z [m]'],
                mode='markers',
                marker=dict(
                    size=3,
                    color=time_data[var],
                    colorscale='Viridis',
                    cmin=var_min,
                    cmax=var_max,
                    colorbar=dict(
                        title=dict(
                            text=f"{var.split(' ')[1] if ' ' in var else var}",
                            font=dict(size=10)
                        ),
                        x=colorbar_x,
                        y=colorbar_y,
                        len=0.3,
                        thickness=15,
                        tickfont=dict(size=8),
                        tickformat='.2e' if var_max < 1e-6 else '.3f'
                    ),
                    opacity=0.7
                ),
                name=f"{var} (t={time_idx})"
            )
            frame_data.append(scatter)
        
        frames.append(go.Frame(data=frame_data, name=str(time_idx)))
    
    # Add initial data (first time step)
    initial_data = df[df['Time Index'] == time_indices[0]]
    for i, var in enumerate(available_vars):
        row = (i // n_cols) + 1
        col = (i % n_cols) + 1
        
        var_min = df[var].min()
        var_max = df[var].max()
        
        # Handle case where all values are zero or very small
        if var_max <= 1e-15:
            var_max = 1e-15
            var_min = 0
        
        # Better colorbar positioning
        colorbar_x = 1.02 if col == n_cols else 0.48
        colorbar_y = 0.75 if row == 1 else 0.25
        
        fig.add_trace(
            go.Scatter3d(
                x=initial_data['X [m]'],
                y=initial_data['Y [m]'],
                z=initial_data['Z [m]'],
                mode='markers',
                marker=dict(
                    size=3,
                    color=initial_data[var],
                    colorscale='Viridis',
                    cmin=var_min,
                    cmax=var_max,
                    colorbar=dict(
                        title=dict(
                            text=f"{var.split(' ')[1] if ' ' in var else var}",
                            font=dict(size=10)
                        ),
                        x=colorbar_x,
                        y=colorbar_y,
                        len=0.3,
                        thickness=15,
                        tickfont=dict(size=8),
                        tickformat='.2e' if var_max < 1e-6 else '.3f'
                    ),
                    opacity=0.7
                ),
                name=var
            ),
            row=row, col=col
        )
    
    # Update layout with animation controls
    fig.update_layout(
        title="Multi-Variable 3D Concentration Visualization Over Time",
        height=900,
        width=1200,
        showlegend=False,
        margin=dict(l=50, r=150, t=80, b=100),
        updatemenus=[{
            'type': 'buttons',
            'showactive': False,
            'x': 0.1,
            'y': 0.02,
            'buttons': [
                {'label': 'Play', 'method': 'animate', 'args': [None, {'frame': {'duration': 150}}]},
                {'label': 'Pause', 'method': 'animate', 'args': [[None], {'frame': {'duration': 0}, 'mode': 'immediate'}]}
            ]
        }],
        sliders=[{
            'steps': [
                {
                    'args': [[str(t)], {'frame': {'duration': 0}, 'mode': 'immediate'}],
                    'label': f't={t}',
                    'method': 'animate'
                } for t in time_indices
            ],
            'currentvalue': {'prefix': 'Time Index: '},
            'len': 0.8,
            'x': 0.1,
            'y': 0.02
        }]
    )
    
    # Update 3D scene properties for all subplots
    for i in range(1, n_plots + 1):
        fig.update_scenes(
            xaxis_title='X [m]',
            yaxis_title='Y [m]',
            zaxis_title='Z [m]',
            selector=dict(type='scene')
        )
    
    # Add frames to figure
    fig.frames = frames
    
    return fig

########################################################################
###### 0) Main function
########################################################################
def main():
    print("PFLOTRAN Visualization Starting...\n")
    
    # Load data
    print("Loading processed data...")
    df = load_data("pflotran_data.pkl")
    
    if df is None:
        return
    
    # Display available variables
    display_available_variables(df)
    
    # CONFIGURATION: Change these variables to plot different data
    # Modify this list to visualize different variables from your simulation
    variables_to_plot = [
        'Total CH4(aq) [M]',
        'Total Acetate- [M]',
        'Total SO4-- [M]',
        'Gamma H2O'
    ]
    
    print(f"Creating visualizations...")
    print(f"Variables to plot: {variables_to_plot}")
    print(f"\n💡 To change variables, edit the 'variables_to_plot' list in this script.")
    
    # Create multi-variable plot
    fig = create_multi_variable_plot(df, variables_to_plot)
    
    if fig is not None:
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
        
        print(f"\nVisualization complete! Open html file")

        
        # Show the plot
        try:
            fig.show()
        except Exception as e:
            print(f"  Could not display plot interactively: {e}")
            print(f"   Please open {output_filename} in your web browser.")
    
    else:
        print("Failed to create visualization. Please check your variable names.")

if __name__ == "__main__":
    main()
