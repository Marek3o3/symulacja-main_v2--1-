# main.py
"""
Main script for satellite constellation simulation and visualization.
Run this script to generate both a static plot and an animation video.
"""
from orbit import propagate_orbits
from visualization import animate_coverage, plot_coverage_overlay
from config import simulation_mode, satellites

def filter_satellites(satellites, mode):
    return satellites

if __name__ == "__main__":
    import numpy as np
    # Filter satellites based on simulation_mode
    selected_sats = filter_satellites(satellites, simulation_mode)
    # Propagate orbits and calculate FOVs
    satellites_with_tracks, sky_datetime_objects = propagate_orbits(selected_sats)
    # Generate and save animation video (with coverage overlay)
    animate_coverage(satellites_with_tracks, sky_datetime_objects)
    # Plot static coverage overlay at the end
    plot_coverage_overlay(satellites_with_tracks)
    print("\nSimulation and visualization complete. Check the output video file and coverage plot in your folder.")