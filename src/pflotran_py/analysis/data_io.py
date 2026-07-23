"""Persistence for the extracted PFLOTRAN DataFrame (pickle + CSV)."""

import pandas as pd


def load_data(filename="pflotran_data.pkl"):
    """Load a processed DataFrame from pickle."""
    return pd.read_pickle(filename)


def save_data(df, filename="pflotran_data.pkl"):
    """Save the combined DataFrame as pickle, plus a CSV copy for convenience."""
    df.to_pickle(filename)
    csv_filename = filename.replace(".pkl", ".csv")
    df.to_csv(csv_filename, index=False)
