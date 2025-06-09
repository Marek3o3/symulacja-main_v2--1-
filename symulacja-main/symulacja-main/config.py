# config.py
"""
Configuration file for satellite constellation simulation.
Edit this file to change satellites, stations, map, and simulation settings.
"""

from datetime import datetime, timedelta
from skyfield.api import utc

# --- Simulation mode ---
# Choose 'MODIS', 'WorldView-3', or 'Both'
simulation_mode = 'WorldView-3'

# --- Earth parameters ---
R_earth_km = 6371
mu_km3_s2 = 398600.4418

# --- Satellite definitions (add more as needed) ---
base_satellite = {
    "altitude_km": 617,
    "inclination_deg": 98,
    "fov_deg": 1.2,  # For ~13.1 km swath at 617 km
    "color": "red",
    "tle_line1": "1 40115U 14045A   25141.50000000  .00000714  00000-0  39998-4 0  9992",
    "tle_line2": "2 40115  97.8743  47.2925 0001529  71.0323  289.1058 14.83390144560002"
}

satellites = {}
for i in range(200):
    sat = base_satellite.copy()
    # Modyfikacja TLE, aby każda satelita miała unikalną nazwę i lekko zmieniony parametr (np. RAAN)
    tle2_parts = sat["tle_line2"].split()
    # Zmieniamy RAAN (pole 3) o 1 stopień na satelitę
    raan = float(tle2_parts[2]) + i * (360/200)
    tle2_parts[2] = f"{raan:.4f}"
    sat["tle_line2"] = " ".join(tle2_parts)
    sat["color"] = "C" + str(i % 10)  # Kolory cyklicznie
    satellites[f"SAT_{i+1:03d}"] = sat

# --- Time configuration ---
start_datetime = datetime(2025, 5, 21, 0, 0, 0, tzinfo=utc)
simulation_duration_hours = 24  # skrócony czas symulacji
time_step_seconds = 60  # większy krok czasowy

# --- Reception Stations ---
reception_stations = {
    'Utqiagvik, USA': {'lat': 71.27, 'lon': -156.81, 'marker': 'o', 'color': 'green'},
    'Dundee, Scotland': {'lat': 56.40, 'lon': -3.18, 'marker': 's', 'color': 'blue'},
    'Chitose, Japan': {'lat': 42.77, 'lon': 141.62, 'marker': '^', 'color': 'red'},
    'Mojave, USA': {'lat': 35.05, 'lon': -118.15, 'marker': 'D', 'color': 'orange'},
    'Dubai, UAE': {'lat': 24.94, 'lon': 55.35, 'marker': 'v', 'color': 'purple'},
    'Paumalu, USA': {'lat': 21.67, 'lon': -158.03, 'marker': 'P', 'color': 'brown'},
    'Harmon, Guam': {'lat': 13.51, 'lon': 144.82, 'marker': '*', 'color': 'pink'},
    'Mwulire, Rwanda': {'lat': -1.96, 'lon': 30.39, 'marker': 'X', 'color': 'gray'},
    'Tahiti, French Polynesia': {'lat': -17.63, 'lon': -149.60, 'marker': 'o', 'color': 'olive'},
    'Awarua, New Zeland': {'lat': -46.52, 'lon': 168.38, 'marker': 's', 'color': 'cyan'},
    'Öjebyn, Sweden': {'lat': 65.33 , 'lon': 21.42, 'marker': '^', 'color': 'magenta'},
    'North Pole, USA': {'lat': 64.79, 'lon': -147.53, 'marker': 'D', 'color': 'teal'},
    'Guildford Pole, UK': {'lat': 51.24, 'lon': -0.61, 'marker': 'v', 'color': 'navy'},
    'Obihiro, Japan': {'lat': 42.59, 'lon': 143.45, 'marker': 'P', 'color': 'gold'},
    'Pendergrass, USA': {'lat': 34.17, 'lon': -83.67, 'marker': '*', 'color': 'lime'},
    'Accra, Ghana': {'lat': 5.74, 'lon': -0.30, 'marker': 'X', 'color': 'coral'},
    'Pretoria, South Africa': {'lat': -25.88, 'lon': 27.70, 'marker': 'o', 'color': 'maroon'},
    'Cordoba, Argentina': {'lat': -31.52, 'lon': -64.46, 'marker': 's', 'color': 'deepskyblue'},
    'Ushuaia, Argentina': {'lat': -54.51, 'lon': -67.11, 'marker': '^', 'color': 'crimson'},
    'Sodankylä, Finland': {'lat': 67.36, 'lon': 26.63, 'marker': 'D', 'color': 'orchid'},
    'Alice Springs, Australia': {'lat': -23.75, 'lon': 133.88, 'marker': 'v', 'color': 'slateblue'},
    'Mingenew, Australia': {'lat': -29.01, 'lon': 115.34, 'marker': 'P', 'color': 'darkorange'},
    'Stockholm, Sweden': {'lat': 59.33, 'lon': 18.06, 'marker': '*', 'color': 'turquoise'},
    'Dublin, Ireland': {'lat': 53.41, 'lon': 8.24, 'marker': 'X', 'color': 'firebrick'},
    'Portland, USA': {'lat': 45.52, 'lon': -122.67, 'marker': 'o', 'color': 'darkgreen'},
    'Columbus, USA': {'lat': 39.96, 'lon': -83.0, 'marker': 's', 'color': 'dodgerblue'},
    'Seoul, South Korea': {'lat': 37.55, 'lon': 126.99, 'marker': '^', 'color': 'indigo'},
    'Deadhorse, USA': {'lat': 70.20, 'lon': -148.45, 'marker': 'D', 'color': 'hotpink'},
    'Zallaq, Bahrain': {'lat': 26.04, 'lon': 50.48, 'marker': 'v', 'color': 'mediumseagreen'},
    'Kapolei, USA': {'lat': 21.33, 'lon': -158.08, 'marker': 'P', 'color': 'slategray'},
    'Singapore': {'lat': 1.35, 'lon': 103.81, 'marker': '*', 'color': 'peru'},
    'Cape Town, South Africa': {'lat': -33.92, 'lon': 18.42, 'marker': 'X', 'color': 'mediumvioletred'},
    'Dubbo, Australia': {'lat': -32.24, 'lon': 148.61, 'marker': 'o', 'color': 'darkred'},
    'Punta Arenas, Chile': {'lat': -53.16, 'lon': -70.90, 'marker': 's', 'color': 'steelblue'},
}

# --- Map configuration ---
use_custom_map = True  # Set to True to use your custom map
custom_map_path = 'custom_map.png'  # Place your PNG/GeoTIFF in the same folder
custom_map_extent = [-180, 180, -90, 90]  # [lon_min, lon_max, lat_min, lat_max]

# --- Animation parameters ---
animation_fps = 5  # mniej klatek na sekundę
animation_bitrate = 1800
fov_history_length = 5  # Number of previous FOVs to show with transparency
output_video_filename = "dual_satellite_animation.mp4"

# --- Coverage analysis ---
latitude_limit = 66.73  # Exclude coverage above/below +/- this latitude (e.g., 80 deg)

# --- Scenario Configuration ---
# Scenario A: Station outage simulation
disabled_stations = []  # Lista nazw stacji do wyłączenia, np. ['Utqiagvik, USA', 'Dundee, Scotland']

# Scenario B: Random transmission failure simulation  
enable_random_transmission_failure = False  # True to enable random transmission failures
transmission_failure_probability = 0.1  # Prawdopodobieństwo utraty transmisji (0.0 - 1.0)
random_seed = 42  # Seed dla powtarzalności wyników (None dla losowości)