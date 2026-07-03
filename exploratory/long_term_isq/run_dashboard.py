#!/usr/bin/env python3
"""
Run the Interactive Gas Concentration Dashboard

This script starts a Bokeh server to run the interactive dashboard.
Usage: python run_dashboard.py
"""

import subprocess
import sys
import os


def run_dashboard():
    """Run the Bokeh dashboard server"""

    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dashboard_file = os.path.join(script_dir, "interactive_gas_dashboard.py")

    if not os.path.exists(dashboard_file):
        print(f"Error: {dashboard_file} not found!")
        sys.exit(1)

    print("Starting Interactive Gas Concentration Dashboard...")
    print("Dashboard will open in your default web browser.")
    print("Press Ctrl+C to stop the server.")
    print("-" * 50)

    try:
        # Run bokeh serve command
        cmd = [
            "bokeh",
            "serve",
            dashboard_file,
            "--show",  # Automatically open in browser
            "--port",
            "5006",  # Use port 5006
        ]

        subprocess.run(cmd, check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error running dashboard: {e}")
        print("\nMake sure you have bokeh installed:")
        print("pip install bokeh")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nDashboard stopped.")
    except FileNotFoundError:
        print("Error: 'bokeh' command not found.")
        print("Make sure you have bokeh installed:")
        print("pip install bokeh")
        sys.exit(1)


if __name__ == "__main__":
    run_dashboard()
