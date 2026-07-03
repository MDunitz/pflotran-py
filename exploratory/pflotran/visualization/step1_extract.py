#!/usr/bin/env python3
### PFLOTRAN Data Extraction Script

# Should install + import dependencies
# Will show you what variables are available in the tecplot files
# Save the combined data for visualization in step2_plot.py, as a csv
# make sure to change these or else
"""
data_dir (str): Directory containing the .tec files
file_template (str): Template for file names with {} for numbering
n_files (int): Number of files to process
"""


import os
import sys
import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm


########################################################################
###### 1) Install and import dependencies
########################################################################
def install_dependencies():
    """Install required packages if not already installed."""
    packages = ['pandas', 'numpy', 'matplotlib', 'tqdm', 'plotly']
    
    for package in packages:
        try:
            __import__(package)
            print(f" {package} is already installed")
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
            print(f"{package} installed successfully")

########################################################################
###### 2) Get variable names from a sample tecplot file
########################################################################
# get the variable names from a sample tecplot file
# useful if you're working with a new PFLOTRAN output file
def get_variable_names(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("VARIABLES"):
                # Clean and split variable names
                variables = [v.strip().strip('"') for v in line.split("=")[1].split(",")]
                return variables
    return []

########################################################################
###### 3) Read tecplot files, return a dataframe to give structure
########################################################################
# It will also call function 2

def read_tec_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Extract variables
    var_line = next(l for l in lines if l.startswith("VARIABLES"))
    variables = [v.strip().strip('"') for v in var_line.split("=")[1].split(",")]
    
    # Extract data lines (after the ZONE line)
    data_start_idx = next(i for i, l in enumerate(lines) if l.startswith("ZONE"))
    data_lines = lines[data_start_idx + 1:]
    data = [list(map(float, l.strip().split())) for l in data_lines]
    
    df = pd.DataFrame(data, columns=variables)
    return df

########################################################################
###### 4) Extract pflotran data from all tecplot files
########################################################################
# Will also call functions 2 and 3
def extract_pflotran_data(data_dir=".", file_template="test29-{:03d}.tec", n_files=100):
    """
    Extract data from all PFLOTRAN tecplot files.
    
    Args:

        
    Returns:
        tuple: (combined_dataframe, variable_names)
    """
    print(f"Looking for files in: {os.path.abspath(data_dir)}")
    print(f"File template: {file_template}")
    print(f"Expected number of files: {n_files}")
    
    # Check if sample file exists to get variable names
    sample_filename = file_template.format(0)
    sample_filepath = os.path.join(data_dir, sample_filename)
    
    if not os.path.exists(sample_filepath):
        print(f"Sample file not found: {sample_filepath}")
        print("Available files in directory:")
        for f in os.listdir(data_dir):
            if f.endswith('.tec'):
                print(f"   {f}")
        return None, None
    
    # Get variable names from sample file
    variable_names = get_variable_names(sample_filepath)
    print(f"\n Variables found in tecplot files:")
    for i, var in enumerate(variable_names):
        print(f"   {i:2d}: {var}")
    
    # Storage for parsed data
    all_data = []
    
    print(f"\nReading {n_files} tecplot files...")
    
    # Read all files
    files_read = 0
    for i in tqdm(range(n_files), desc="Processing files"):
        filename = file_template.format(i)
        filepath = os.path.join(data_dir, filename)
        
        if os.path.exists(filepath):
            try:
                df = read_tec_file(filepath)
                df["Time Index"] = i
                all_data.append(df)
                files_read += 1
            except Exception as e:
                print(f"Error reading {filename}: {e}")
        else:
            print(f"File not found: {filename}")
    
    if not all_data:
        print("No data files were successfully read!")
        return None, None
    
    # Combine all dataframes
    print(f"\nCombining data from {files_read} files...")
    full_df = pd.concat(all_data, ignore_index=True)
    
    print(f"Data extraction complete!")
    print(f"   Combined dataset shape: {full_df.shape}")
    print(f"   Time indices: {sorted(full_df['Time Index'].unique())}")
    print(f"   Variables: {len(variable_names)}")
    
    return full_df, variable_names

########################################################################
###### 5) Save the combined dataframe
########################################################################
def save_data(df, filename="pflotran_data.pkl"):
    if df is not None:
        df.to_pickle(filename)
        print(f"💾 Data saved to: {filename}")
        
        # Also save as CSV for portability (though it will be larger)
        csv_filename = filename.replace('.pkl', '.csv')
        df.to_csv(csv_filename, index=False)
        print(f"💾 Data also saved as CSV: {csv_filename}")

########################################################################
###### 0) Main
########################################################################
def main():
    print("PFLOTRAN Data Extraction Starting...\n")
    
    # Install dependencies
    print("hecking dependencies...")
    install_dependencies()
    
    # Configuration - modify these as needed
    data_dir = "."
    file_template = "test29-{:03d}.tec"
    n_files = 100
    
    print(f"\nExtracting data...")
    
    # Extract data
    full_df, variable_names = extract_pflotran_data(
        data_dir=data_dir,
        file_template=file_template,
        n_files=n_files
    )
    
    if full_df is not None:
        print(f"\nSaving data...")
        save_data(full_df, "pflotran_data.pkl")
        
        print(f"\nNext steps:")
        print(f"   Run: python step2_plot.py")
        print(f"   This will create interactive 3D visualizations of your data.")
        
        # Import and run step2 automatically
        print(f"\nAutomatically running visualization script...")
        try:
            if os.path.exists("step2_plot.py"):
                import step2_plot
                step2_plot.main()
            else:
                print("step2_plot.py not found in current directory.")
                print("Please ensure step2_plot.py is in the same directory and run it manually.")
        except Exception as e:
            print(f"Error running visualization: {e}")
            print("Please run step2_plot.py manually.")
    else:
        print("Data extraction failed. Please check your file paths and try again.")

if __name__ == "__main__":
    main()
