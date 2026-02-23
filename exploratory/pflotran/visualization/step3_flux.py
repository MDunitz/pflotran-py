#!/usr/bin/env python3
"""
PFLOTRAN Flux Visualization Script
=================================
Visualizes CO2 and CH4 fluxes for the surface.

Should output an html with time fluxes.

Change these as you use:
*     for i, time_idx in enumerate(time_indices[:5]):
I put 5, change to the # of time steps you want to output. 

"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm

def install_dependencies():
    """Install required packages if not already installed."""
    packages = ['pandas', 'numpy', 'bokeh', 'tqdm', 'scipy']
    
    for package in packages:
        try:
            if package == 'bokeh':
                import bokeh
            elif package == 'scipy':
                import scipy
            else:
                __import__(package)
            print(f"{package} is already installed")
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
            print(f"{package} installed successfully")

def load_data(filename="pflotran_data.pkl"):
    """
    Load the processed PFLOTRAN data.
    
    Args:
        filename (str): Filename of the saved data
        
    Returns:
        pd.DataFrame: Combined dataframe, or None if file not found
    """
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
        print(f"Data file not found: {filename}")
        print("Please run step1_extract.py first to process the data.")
        return None

def calculate_concentration_gradients(df):
    """
    Calculate spatial gradients for CO2 and CH4 concentrations.
    
    Args:
        df (pd.DataFrame): The dataframe containing the data
        
    Returns:
        pd.DataFrame: DataFrame with added gradient columns
    """
    from scipy.spatial.distance import cdist
    
    # Create a copy of the dataframe
    flux_df = df.copy()
    
    # Species of interest for flux calculations
    co2_col = 'CO2(aq) [M]'
    ch4_total_col = 'Total CH4(aq) [M]'
    ch4_free_col = 'Free CH4(aq) [M]'
    
    # Initialize gradient columns
    flux_df['CO2_flux_x'] = 0.0
    flux_df['CO2_flux_y'] = 0.0  
    flux_df['CO2_flux_z'] = 0.0
    flux_df['CO2_flux_magnitude'] = 0.0
    flux_df['CH4_flux_x'] = 0.0
    flux_df['CH4_flux_y'] = 0.0
    flux_df['CH4_flux_z'] = 0.0
    flux_df['CH4_flux_magnitude'] = 0.0
    
    # Get unique time indices
    time_indices = flux_df['Time Index'].unique()
    
    print("Calculating concentration gradients and fluxes...")
    
    for time_idx in tqdm(time_indices, desc="Processing time steps"):
        # Get data for this time step
        time_data = flux_df[flux_df['Time Index'] == time_idx].copy()
        
        if len(time_data) < 8:  # Need enough points for gradient calculation
            continue
            
        # Create coordinate arrays
        coords = time_data[['X [m]', 'Y [m]', 'Z [m]']].values
        
        # Calculate gradients for each point
        for idx, row in time_data.iterrows():
            point = row[['X [m]', 'Y [m]', 'Z [m]']].values
            
            # Find neighboring points (within some distance threshold)
            distances = np.sqrt(np.sum((coords - point)**2, axis=1))
            neighbors_idx = np.where((distances > 0) & (distances <= 0.5))[0]  # Adjust threshold as needed
            
            if len(neighbors_idx) >= 3:  # Need at least 3 neighbors for gradient
                neighbor_coords = coords[neighbors_idx]
                neighbor_co2 = time_data.iloc[neighbors_idx][co2_col].values
                neighbor_ch4 = time_data.iloc[neighbors_idx][ch4_free_col].values
                
                # Simple finite difference approximation
                try:
                    # Calculate gradients using least squares fit
                    A = np.column_stack([neighbor_coords - point, np.ones(len(neighbor_coords))])
                    
                    # CO2 gradients
                    co2_gradients, _, _, _ = np.linalg.lstsq(A, neighbor_co2 - row[co2_col], rcond=None)
                    flux_df.loc[idx, 'CO2_flux_x'] = -co2_gradients[0] if len(co2_gradients) > 0 else 0
                    flux_df.loc[idx, 'CO2_flux_y'] = -co2_gradients[1] if len(co2_gradients) > 1 else 0
                    flux_df.loc[idx, 'CO2_flux_z'] = -co2_gradients[2] if len(co2_gradients) > 2 else 0
                    
                    # CH4 gradients
                    ch4_gradients, _, _, _ = np.linalg.lstsq(A, neighbor_ch4 - row[ch4_free_col], rcond=None)
                    flux_df.loc[idx, 'CH4_flux_x'] = -ch4_gradients[0] if len(ch4_gradients) > 0 else 0
                    flux_df.loc[idx, 'CH4_flux_y'] = -ch4_gradients[1] if len(ch4_gradients) > 1 else 0
                    flux_df.loc[idx, 'CH4_flux_z'] = -ch4_gradients[2] if len(ch4_gradients) > 2 else 0
                    
                except np.linalg.LinAlgError:
                    # If linear algebra fails, set fluxes to zero
                    pass
    
    # Calculate flux magnitudes
    flux_df['CO2_flux_magnitude'] = np.sqrt(
        flux_df['CO2_flux_x']**2 + flux_df['CO2_flux_y']**2 + flux_df['CO2_flux_z']**2
    )
    flux_df['CH4_flux_magnitude'] = np.sqrt(
        flux_df['CH4_flux_x']**2 + flux_df['CH4_flux_y']**2 + flux_df['CH4_flux_z']**2
    )
    
    # Print some statistics about the calculated fluxes
    print(f"CO2 flux statistics:")
    print(f"  Min: {flux_df['CO2_flux_magnitude'].min():.3e}")
    print(f"  Max: {flux_df['CO2_flux_magnitude'].max():.3e}")
    print(f"  Mean: {flux_df['CO2_flux_magnitude'].mean():.3e}")
    
    print(f"CH4 flux statistics:")
    print(f"  Min: {flux_df['CH4_flux_magnitude'].min():.3e}")
    print(f"  Max: {flux_df['CH4_flux_magnitude'].max():.3e}")
    print(f"  Mean: {flux_df['CH4_flux_magnitude'].mean():.3e}")
    
    return flux_df

def identify_surface_cells(df):
    """
    Identify surface cells (highest Z values for each X,Y position).
    
    Args:
        df (pd.DataFrame): The dataframe containing the data
        
    Returns:
        pd.DataFrame: DataFrame with surface indicators
    """
    surface_df = df.copy()
    surface_df['is_surface'] = False
    
    # For each X,Y position, find the maximum Z value
    for time_idx in df['Time Index'].unique():
        time_mask = df['Time Index'] == time_idx
        time_data = df[time_mask].copy()
        
        # Group by X,Y and find maximum Z for each group
        max_z_by_xy = time_data.groupby(['X [m]', 'Y [m]'])['Z [m]'].max().reset_index()
        max_z_by_xy.columns = ['X [m]', 'Y [m]', 'max_z']
        
        # Merge to identify surface cells
        time_data_with_max = time_data.merge(max_z_by_xy, on=['X [m]', 'Y [m]'])
        surface_mask = (time_data_with_max['Z [m]'] == time_data_with_max['max_z'])
        
        # Update the surface indicator using the original indices
        surface_indices = time_data_with_max[surface_mask].index
        surface_df.loc[surface_indices, 'is_surface'] = True
    
    return surface_df

def create_flux_visualization(flux_df):
    """
    Create interactive Bokeh visualization of CO2 and CH4 fluxes.
    
    Args:
        flux_df (pd.DataFrame): DataFrame with flux calculations
    """
    from bokeh.plotting import figure, save, output_file
    from bokeh.layouts import column, row
    from bokeh.models import HoverTool, ColorBar, ColumnDataSource, Select
    from bokeh.transform import linear_cmap
    from bokeh.palettes import Viridis256
    from bokeh.io import curdoc
    
    print("Creating flux visualization...")
    
    # Prepare data for different time steps
    time_indices = sorted(flux_df['Time Index'].unique())
    
    # Create the main flux visualization
    output_file("co2_ch4_flux_visualization.html")
    
    # Create plots for each time step
    plots = []
    
    for i, time_idx in enumerate(time_indices[:5]):  # Limit to first 5 time steps for testing
        time_data = flux_df[flux_df['Time Index'] == time_idx]
        
        # Create 2D projection plots (X-Y view, top-down)
        surface_data = time_data[time_data['is_surface'] == True] if 'is_surface' in time_data.columns else time_data
        
        if len(surface_data) == 0:
            surface_data = time_data
        
        # Check if there are any meaningful fluxes for this time step
        max_co2_flux = surface_data['CO2_flux_magnitude'].max()
        max_ch4_flux = surface_data['CH4_flux_magnitude'].max()
        
        print(f"Time step {time_idx}: Max CO2 flux = {max_co2_flux:.3e}, Max CH4 flux = {max_ch4_flux:.3e}")
        
        # Skip this time step if no meaningful fluxes (very relaxed thresholds to always show CH4)
        # Skip only if absolutely no fluxes (extremely relaxed threshold)
        if max_co2_flux <= 1e-25 and max_ch4_flux <= 1e-30:
            print(f"Skipping time step {time_idx} - no significant fluxes detected")
            continue
        
        # CO2 flux plot
        p1 = figure(
            title=f"CO2 Surface Flux - Time Step {time_idx}",
            x_axis_label="X [m]",
            y_axis_label="Y [m]",
            width=400,
            height=400,
            tools="pan,wheel_zoom,box_zoom,reset,save"
        )
        
        # Create color mapping for CO2 flux magnitude
        co2_flux_min = surface_data['CO2_flux_magnitude'].min()
        co2_flux_max = surface_data['CO2_flux_magnitude'].max()
        
        # Only plot if there are non-zero fluxes
        if co2_flux_max > co2_flux_min and co2_flux_max > 1e-15:
            color_mapper_co2 = linear_cmap('CO2_flux_magnitude', Viridis256, co2_flux_min, co2_flux_max)
            
            source_co2 = ColumnDataSource(surface_data)
            circles_co2 = p1.scatter(
                'X [m]', 'Y [m]', 
                size=15, 
                color=color_mapper_co2,
                source=source_co2,
                alpha=0.7
            )
            
            # Add arrows for flux direction
            arrow_scale = 1000  # Adjust as needed
            p1.segment(
                x0='X [m]', y0='Y [m]',
                x1=f'X [m]', y1=f'Y [m]',  # Will be computed in JavaScript
                source=source_co2,
                line_width=2,
                alpha=0.5,
                color='black'
            )
            
            # Add hover tool
            hover_co2 = HoverTool(
                tooltips=[
                    ("X [m]", "@{X [m]}{0.00}"),
                    ("Y [m]", "@{Y [m]}{0.00}"),
                    ("Z [m]", "@{Z [m]}{0.00}"),
                    ("CO2 Concentration [M]", "@{CO2(aq) [M]}"),
                    ("CO2 Flux Magnitude [M·m⁻¹]", "@CO2_flux_magnitude"),
                    ("CO2 Flux X [M·m⁻¹]", "@CO2_flux_x"),
                    ("CO2 Flux Y [M·m⁻¹]", "@CO2_flux_y"),
                    ("CO2 Flux Z [M·m⁻¹]", "@CO2_flux_z")
                ],
                formatters={
                    "@{CO2(aq) [M]}": "printf",
                    "@CO2_flux_magnitude": "printf",
                    "@CO2_flux_x": "printf",
                    "@CO2_flux_y": "printf",
                    "@CO2_flux_z": "printf"
                },
                renderers=[circles_co2]
            )
            p1.add_tools(hover_co2)
        
        # CH4 flux plot
        p2 = figure(
            title=f"CH4 Surface Flux - Time Step {time_idx}",
            x_axis_label="X [m]",
            y_axis_label="Y [m]",
            width=400,
            height=400,
            tools="pan,wheel_zoom,box_zoom,reset,save"
        )
        
        # Create color mapping for CH4 flux magnitude
        ch4_flux_min = surface_data['CH4_flux_magnitude'].min()
        ch4_flux_max = surface_data['CH4_flux_magnitude'].max()
        
        # Only plot if there are non-zero fluxes (very relaxed threshold for CH4 - always show if any detectable)
        if ch4_flux_max > ch4_flux_min and ch4_flux_max > 1e-25:
            color_mapper_ch4 = linear_cmap('CH4_flux_magnitude', Viridis256, ch4_flux_min, ch4_flux_max)
            
            source_ch4 = ColumnDataSource(surface_data)
            circles_ch4 = p2.scatter(
                'X [m]', 'Y [m]', 
                size=15, 
                color=color_mapper_ch4,
                source=source_ch4,
                alpha=0.7
            )
            
            # Add hover tool
            hover_ch4 = HoverTool(
                tooltips=[
                    ("X [m]", "@{X [m]}{0.00}"),
                    ("Y [m]", "@{Y [m]}{0.00}"),
                    ("Z [m]", "@{Z [m]}{0.00}"),
                    ("CH4 Concentration [M]", "@{Free CH4(aq) [M]}"),
                    ("CH4 Flux Magnitude [M·m⁻¹]", "@CH4_flux_magnitude"),
                    ("CH4 Flux X [M·m⁻¹]", "@CH4_flux_x"),
                    ("CH4 Flux Y [M·m⁻¹]", "@CH4_flux_y"),
                    ("CH4 Flux Z [M·m⁻¹]", "@CH4_flux_z")
                ],
                formatters={
                    "@{Free CH4(aq) [M]}": "printf",
                    "@CH4_flux_magnitude": "printf",
                    "@CH4_flux_x": "printf",
                    "@CH4_flux_y": "printf",
                    "@CH4_flux_z": "printf"
                },
                renderers=[circles_ch4]
            )
            p2.add_tools(hover_ch4)
        
        # Create row for this time step
        time_row = row(p1, p2)
        plots.append(time_row)
    
    # Create final layout
    layout = column(*plots)
    save(layout)
    
    print("CO2 and CH4 flux visualization saved as: co2_ch4_flux_visualization.html")


def main():
    """Main execution function."""
    print("Starting PFLOTRAN Flux Visualization...")
    
    # Install dependencies
    install_dependencies()
    
    # Load data
    df = load_data()
    if df is None:
        print("Cannot proceed without data. Exiting.")
        return
    
    # Check if we have the required gas species columns
    required_columns = ['CO2(aq) [M]', 'Free CH4(aq) [M]']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Missing required columns: {missing_columns}")
        print("Available columns:")
        for col in df.columns:
            if any(gas in col.lower() for gas in ['co2', 'ch4', 'methane']):
                print(f"   - {col}")
        return
    
    print(f"Dataset contains {len(df)} data points across {df['Time Index'].nunique()} time steps")
    
    # Print concentration statistics
    print(f"CO2 concentration statistics:")
    print(f"  Min: {df['CO2(aq) [M]'].min():.3e}")
    print(f"  Max: {df['CO2(aq) [M]'].max():.3e}")
    print(f"  Mean: {df['CO2(aq) [M]'].mean():.3e}")
    
    print(f"CH4 concentration statistics:")
    print(f"  Min: {df['Free CH4(aq) [M]'].min():.3e}")
    print(f"  Max: {df['Free CH4(aq) [M]'].max():.3e}")
    print(f"  Mean: {df['Free CH4(aq) [M]'].mean():.3e}")
    
    # Check if CH4 is at numerical floor
    if df['Free CH4(aq) [M]'].max() <= 1e-14:
        print("  Note: CH4 concentrations are at numerical floor - little/no methane production yet")
    
    # Calculate concentration gradients and fluxes
    print("Calculating flux fields...")
    flux_df = calculate_concentration_gradients(df)
    
    # Identify surface cells
    print("Identifying surface cells...")
    flux_df = identify_surface_cells(flux_df)
     # Create visualizations
    create_flux_visualization(flux_df)

    print("\nFlux visualization complete!")
    print("\nGenerated files:")
    print("   co2_ch4_flux_visualization.html - Interactive flux maps")
    print("\nTo view the results:")
    print("   Open the HTML file in your web browser")
    print("   Or use: open co2_ch4_flux_visualization.html")

if __name__ == "__main__":
    main()