#!/usr/bin/env python3
# PFLOTRAN Data Extraction Script

# Should install + import dependencies
# Will show you what variables are available in the tecplot files
# Save the combined data for visualization in step2_plot.py, as a csv

# ##################################################################
#  1) Install and import dependencies
# ##################################################################
import os
import pandas as pd


# ##################################################################
#  2) Get variable names from a sample tecplot file
# ##################################################################
# get the variable names from a sample tecplot file (there will only be one row)
# useful if you're working with a new PFLOTRAN output file
def extract_plfotran_tecplot_variable_names(filepath):
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("VARIABLES"):
                # Clean and split variable names
                variables = [
                    v.strip().strip('"') for v in line.split("=")[1].split(",")
                ]
                return variables
    return []


# ##################################################################
#  3) Read tecplot files, return a dataframe to give structure
# ##################################################################


def read_tec_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Extract variables
    var_line = next(ln for ln in lines if ln.startswith("VARIABLES"))
    variables = [v.strip().strip('"') for v in var_line.split("=")[1].split(",")]

    # Extract data lines (after the ZONE line)
    data_start_idx = next(i for i, ln in enumerate(lines) if ln.startswith("ZONE"))
    data_lines = lines[data_start_idx + 1 :]
    data = [list(map(float, ln.strip().split())) for ln in data_lines]

    df = pd.DataFrame(data, columns=variables)
    return df


# ##################################################################
#  4) Extract pflotran data from all tecplot files
# ##################################################################
# Will also call functions 2 and 3 (note: will only call 2 if you're doing troubleshooting, and you need to comment out stuff below)
def extract_pflotran_data_tec(
    data_dir=".", file_name_template="test29-{:03d}.tec", n_files=100, verbose=False
):

    # sample file exists to get variable names
    sample_filename = file_name_template.format(0)
    sample_filepath = os.path.join(data_dir, sample_filename)

    if verbose:
        # If you want to see the variable names, uncomment this block, and add it as something you return
        # Get variable names from sample file
        variable_names = extract_plfotran_tecplot_variable_names(sample_filepath)
        print("\n Variables found in tecplot files:")
        for i, var in enumerate(variable_names):
            print(f"   {i:2d}: {var}")

    # Storage for  data
    all_data = []

    # Read all files
    files_read = 0
    for i in range(n_files):
        filename = file_name_template.format(i)
        filepath = os.path.join(data_dir, filename)

        df = read_tec_file(filepath)
        df["Time Index"] = i
        all_data.append(df)
        files_read += 1

    # Combine all dataframes
    full_df = pd.concat(all_data, ignore_index=True)

    return full_df


# ##################################################################
#  5) Save the combined dataframe
# ##################################################################
def save_data(df, filename="pflotran_data.pkl"):
    df.to_pickle(filename)

    # Also save as CSV for convenience (though it will be larger)
    csv_filename = filename.replace(".pkl", ".csv")
    df.to_csv(csv_filename, index=False)


# ##################################################################
#  0) Main
# ##################################################################
def main(data_format, data_dir, file_name_template, n_files):
    data_dir = data_dir
    file_name_template = file_name_template
    n_files = n_files

    if data_format == "tec":
        print("\nExtracting data from tecplot files...")
        full_df = extract_pflotran_data_tec(
            data_dir=data_dir, file_name_template=file_name_template, n_files=n_files
        )
    else:
        print("\nExtracting data from hdf5 files is not yet implemented.")

    print("\nSaving data...")
    save_data(full_df, "pflotran_data.pkl")


if __name__ == "__main__":
    main()
