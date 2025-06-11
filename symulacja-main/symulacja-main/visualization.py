# visualization.py
"""
Optimized calculation functions for satellite simulation: text output only with multiprocessing.
"""
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
from shapely.geometry import box, MultiPolygon
from shapely.ops import unary_union
from config import (
    reception_stations,
    latitude_limit,
    time_step_seconds,
    disabled_stations,
    enable_random_transmission_failure,
    transmission_failure_probability,
    random_seed,
)
import numpy as np
import pandas as pd
import math
import csv
from multiprocessing import Pool, cpu_count
import functools
import random


def calculate_elevation_angle(sat_lat, sat_lon, sat_alt_km, station_lat, station_lon):
    """Calculate elevation angle between satellite and ground station supporting full range 0°-180°."""
    # Convert to radians
    sat_lat_rad = np.radians(sat_lat)
    sat_lon_rad = np.radians(sat_lon)
    st_lat_rad = np.radians(station_lat)
    st_lon_rad = np.radians(station_lon)
    
    # Earth radius
    R_earth = 6371.0  # km
    
    # Calculate angular distance between satellite and station on Earth's surface (haversine formula)
    dlat = st_lat_rad - sat_lat_rad
    dlon = st_lon_rad - sat_lon_rad
    a = np.sin(dlat/2)**2 + np.cos(sat_lat_rad) * np.cos(st_lat_rad) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))  # Angular distance in radians
    
    # Calculate satellite distance from Earth center
    sat_distance_from_center = R_earth + sat_alt_km
    
    # Calculate distance from satellite to station using law of cosines
    satellite_to_station_distance = np.sqrt(
        sat_distance_from_center**2 + R_earth**2 - 
        2 * sat_distance_from_center * R_earth * np.cos(c)
    )
    
    if satellite_to_station_distance == 0:
        return 90.0  # Satellite directly overhead
    
    # Calculate elevation angle using law of sines in the triangle:
    # Earth center - Station - Satellite
    sin_elevation_angle = (sat_distance_from_center * np.sin(c)) / satellite_to_station_distance
    sin_elevation_angle = np.clip(sin_elevation_angle, -1, 1)  # Ensure valid range
    
    # Get the angle at station position in the triangle
    angle_at_station = np.arcsin(sin_elevation_angle)
    
    # Convert to elevation angle (angle above horizon)
    # When c < π/2: normal elevation (0° to 90°)
    # When c > π/2: satellite "behind" station (90° to 180°)
    if c <= np.pi/2:
        # Satellite approaching or at zenith
        elevation_deg = np.degrees(np.pi/2 - angle_at_station)
    else:
        # Satellite moving away (behind the station)
        elevation_deg = np.degrees(angle_at_station + np.pi/2)
    
    # Clamp to valid range and ensure non-negative
    elevation_deg = np.clip(elevation_deg, 0.0, 180.0)
    
    # Check if satellite is below horizon (very large angular distances)
    horizon_limit = np.arccos(R_earth / sat_distance_from_center)  # Horizon angle for this altitude
    if c > horizon_limit:
        return 0.0  # Below horizon
    
    return elevation_deg


def calculate_transmission_rate(elevation_deg):
    """Calculate transmission rate based on elevation angle using parabolic curve matching the chart."""
    # No transmission below 5 degrees or above 175 degrees
    if elevation_deg < 5 or elevation_deg > 175:
        return 0
    
    # Parameters matching the curve in the chart
    min_elevation = 5.0
    max_elevation = 175.0
    peak_elevation = 90.0  # Peak at 90 degrees
    min_rate = 480.0  # Mbps at 5° and 175°
    max_rate = 3000.0  # Mbps at 90°
    
    # Use quadratic formula for parabolic curve
    # The curve is symmetric around 90°, so we can use a parabola formula:
    # rate = a * (elevation - peak)² + max_rate
    # where 'a' is negative to create downward-opening parabola
    
    # Calculate 'a' using the constraint that at min_elevation (5°), rate = min_rate
    # min_rate = a * (5 - 90)² + max_rate
    # min_rate = a * 85² + max_rate
    # a = (min_rate - max_rate) / 85²
    a = (min_rate - max_rate) / ((min_elevation - peak_elevation) ** 2)
    
    # Calculate rate using parabolic formula
    rate = a * ((elevation_deg - peak_elevation) ** 2) + max_rate
    
    # Ensure rate doesn't go below minimum
    rate = max(rate, min_rate)
    
    return rate


def calculate_static_coverage(satellites):
    """Calculate total coverage statistics without visualization."""
    clip_box = box(-180, -latitude_limit, 180, latitude_limit)
    
    # Union of all FOVs clipped
    all_polys = []
    for params in satellites.values():
        for poly in params.get('fov_polygons_shapely', []):
            if poly and poly.is_valid and not poly.is_empty:
                clipped = poly.intersection(clip_box)
                if clipped.is_valid and not clipped.is_empty:
                    all_polys.append(clipped)

    if all_polys:
        coverage_union = unary_union(all_polys)
        if isinstance(coverage_union, MultiPolygon):
            total_area = sum(geom.area for geom in coverage_union.geoms)
        else:
            total_area = coverage_union.area
        
        # Calculate coverage as percentage of Earth's surface within latitude limits
        earth_surface_area = 360 * (2 * latitude_limit)  # rough approximation in degrees²
        coverage_percentage = (total_area / earth_surface_area) * 100
        
        print(f"Total Coverage Analysis:")
        print(f"Coverage area: {total_area:.2f} square degrees")
        print(f"Coverage percentage: {coverage_percentage:.2f}%")
        print(f"Latitude range: {latitude_limit}° to -{latitude_limit}°")
    else:
        print("No valid coverage polygons found.")


def process_frame(args):
    """Process a single frame for coverage calculations with realistic satellite-station pairing."""
    frame, satellites, sky_datetime_objects = args
    
    # Initialize random seed for reproducible results (Scenario B)
    if random_seed is not None:
        random.seed(random_seed + frame)  # Different seed for each frame but deterministic
    
    t = sky_datetime_objects[frame]
    ns = Nightshade(t, alpha=0)
    night_polys = list(ns.geometries())
    
    day_land = 0
    day_ocean = 0
    
    # Classify surface type
    def classify_surface(poly):
        for land_geom in cfeature.LAND.geometries():
            if poly.intersects(land_geom):
                return 'land'
        return 'ocean'
    
    # Check transmission status with elevation-dependent rates
    STATION_RADIUS_KM = 2935
    EARTH_RADIUS_KM = 6371
    all_potential_pairs = []
    
    # Dictionary to store best station for each satellite: sat_name -> (station_name, rate, elevation, distance)
    satellite_best_pairs = {}
    
    for sat_name, params in satellites.items():
        sat_lat = params['latitudes'][frame]
        sat_lon = params['longitudes'][frame]
        sat_alt = params.get('altitudes_actual_km', [617])[frame] if len(params.get('altitudes_actual_km', [617])) > frame else 617
        
        # Count FOV coverage
        curr_fov = params['fov_polygons_shapely'][frame]
        if curr_fov and curr_fov.is_valid and not curr_fov.is_empty:
            cen = curr_fov.centroid
            in_night = any(n.contains(cen) for n in night_polys)
            if not in_night:
                surf = classify_surface(curr_fov)
                if surf == 'land': 
                    day_land += 1
                else:           
                    day_ocean += 1
        
        # Find best station for this satellite
        best_station = None
        best_rate = 0
        best_elevation = 0
        best_distance = 0
        
        # Check station range with elevation-dependent transmission rates
        for st_name, st_info in reception_stations.items():
            # Scenario A: Skip disabled stations
            if st_name in disabled_stations:
                continue
                
            lat1, lon1 = np.radians(sat_lat), np.radians(sat_lon)
            lat2, lon2 = np.radians(st_info['lat']), np.radians(st_info['lon'])
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
            c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
            distance_km = EARTH_RADIUS_KM * c
            
            if distance_km <= STATION_RADIUS_KM:
                # Calculate elevation angle and transmission rate
                elevation = calculate_elevation_angle(
                    sat_lat, sat_lon, sat_alt,
                    st_info['lat'], st_info['lon']
                )
                transmission_rate = calculate_transmission_rate(elevation)
                
                if transmission_rate > 0:  # Only count if above 5° elevation
                    # Add to potential pairs (for original statistics)
                    all_potential_pairs.append({
                        'satellite': sat_name,
                        'station': st_name,
                        'elevation': elevation,
                        'rate_mbps': transmission_rate,
                        'distance_km': distance_km
                    })
                    
                    # Check if this is the best station for this satellite
                    if transmission_rate > best_rate:
                        best_rate = transmission_rate
                        best_station = st_name
                        best_elevation = elevation
                        best_distance = distance_km
        
        # Store best pair for this satellite
        if best_station is not None:
            satellite_best_pairs[sat_name] = (best_station, best_rate, best_elevation, best_distance)
    
    # Original logic: Select only ONE satellite for transmission (the one with highest rate)
    transmitting_pairs = []
    total_transmission_rate = 0
    
    if all_potential_pairs:
        # Find the satellite-station pair with the highest transmission rate
        best_pair = max(all_potential_pairs, key=lambda x: x['rate_mbps'])
        
        # Scenario B: Random transmission failure simulation
        transmission_successful = True
        if enable_random_transmission_failure:
            if random.random() < transmission_failure_probability:
                transmission_successful = False
        
        if transmission_successful:
            transmitting_pairs.append(best_pair)
            total_transmission_rate = best_pair['rate_mbps']
    
    # NEW: Realistic transmission logic - each satellite selects its best station
    # But each ground station can only handle ONE connection at a time
    realistic_transmitting_pairs = []
    aggregated_transmission_rate = 0
    simultaneous_transmissions = 0
    occupied_stations = set()  # Track which stations are already in use
    
    # Sort satellites by their best transmission rate (highest first) to prioritize better connections
    sorted_satellites = sorted(
        satellite_best_pairs.items(), 
        key=lambda x: x[1][1],  # Sort by rate (index 1 in the tuple)
        reverse=True
    )
    
    for sat_name, (station_name, rate, elevation, distance) in sorted_satellites:
        # Check if this station is already occupied
        if station_name in occupied_stations:
            continue  # Skip this satellite - its best station is already in use
            
        # Scenario B: Random transmission failure simulation
        transmission_successful = True
        if enable_random_transmission_failure:
            if random.random() < transmission_failure_probability:
                transmission_successful = False
        
        if transmission_successful:
            realistic_transmitting_pairs.append({
                'satellite': sat_name,
                'station': station_name,
                'elevation': elevation,
                'rate_mbps': rate,
                'distance_km': distance
            })
            aggregated_transmission_rate += rate
            simultaneous_transmissions += 1
            occupied_stations.add(station_name)  # Mark this station as occupied
    
    return {
        'frame': frame,
        'day_land': day_land,
        'day_ocean': day_ocean,
        'transmitting': len(transmitting_pairs) > 0,
        'transmitting_pairs': transmitting_pairs,
        'total_rate_mbps': total_transmission_rate,
        'potential_pairs_count': len(all_potential_pairs),
        # NEW: Realistic transmission metrics
        'realistic_transmitting_pairs': realistic_transmitting_pairs,
        'aggregated_transmission_rate': aggregated_transmission_rate,
        'simultaneous_transmissions': simultaneous_transmissions,
        'realistic_transmitting': len(realistic_transmitting_pairs) > 0
    }


def calculate_coverage_fast(satellites, sky_datetime_objects):
    """Fast calculation of coverage statistics using multiprocessing."""
    num_steps = len(sky_datetime_objects)
    
    print(f"Starting coverage calculations for {num_steps} time steps...")
    print(f"Using {cpu_count()} CPU cores for parallel processing")
    
    # Prepare arguments for multiprocessing
    args_list = [(frame, satellites, sky_datetime_objects) for frame in range(num_steps)]
    
    # Use multiprocessing to calculate coverage
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_frame, args_list)
    
    # Aggregate results
    total_day_land = sum(r['day_land'] for r in results)
    total_day_ocean = sum(r['day_ocean'] for r in results)
    total_transmit_steps = sum(1 for r in results if r['transmitting'])
    
    # Calculate data transmission with variable rates
    TOTAL_DATA_TB = 25.0
    total_transmitted_tb = 0
    
    # Calculate transmission statistics
    elevation_stats = []
    rate_stats = []
    potential_pairs_stats = []
    
    for result in results:
        potential_pairs_stats.append(result.get('potential_pairs_count', 0))
        if result['transmitting']:
            # Calculate data transmitted in this step with actual rates
            step_rate_mbps = result['total_rate_mbps']
            step_transmitted_tb = (step_rate_mbps * time_step_seconds) / 8 / 1e6  # Convert Mbps*s to TB
            total_transmitted_tb += step_transmitted_tb
            
            # Collect statistics
            for pair in result['transmitting_pairs']:
                elevation_stats.append(pair['elevation'])
                rate_stats.append(pair['rate_mbps'])
    
    data_left_tb = max(0, TOTAL_DATA_TB - total_transmitted_tb)
    total_transmit_seconds = total_transmit_steps * time_step_seconds
    
    # Calculate statistics
    avg_elevation = np.mean(elevation_stats) if elevation_stats else 0
    min_elevation = np.min(elevation_stats) if elevation_stats else 0
    max_elevation = np.max(elevation_stats) if elevation_stats else 0
    avg_rate = np.mean(rate_stats) if rate_stats else 0
    min_rate = np.min(rate_stats) if rate_stats else 0
    max_rate = np.max(rate_stats) if rate_stats else 0
    avg_potential_pairs = np.mean(potential_pairs_stats) if potential_pairs_stats else 0
    max_potential_pairs = np.max(potential_pairs_stats) if potential_pairs_stats else 0
    
    # NEW: Calculate realistic transmission statistics
    realistic_total_transmitted_tb = 0
    realistic_transmit_steps = 0
    realistic_elevation_stats = []
    realistic_rate_stats = []
    aggregated_rate_stats = []
    simultaneous_transmission_stats = []
    
    for result in results:
        if result.get('realistic_transmitting', False):
            realistic_transmit_steps += 1
            # Calculate data transmitted in this step with aggregated rates
            step_aggregated_rate_mbps = result['aggregated_transmission_rate']
            step_transmitted_tb = (step_aggregated_rate_mbps * time_step_seconds) / 8 / 1e6  # Convert Mbps*s to TB
            realistic_total_transmitted_tb += step_transmitted_tb
            
            # Collect aggregated statistics
            aggregated_rate_stats.append(step_aggregated_rate_mbps)
            simultaneous_transmission_stats.append(result['simultaneous_transmissions'])
            
            # Collect individual pair statistics
            for pair in result['realistic_transmitting_pairs']:
                realistic_elevation_stats.append(pair['elevation'])
                realistic_rate_stats.append(pair['rate_mbps'])
    
    # Calculate realistic statistics
    realistic_data_left_tb = max(0, TOTAL_DATA_TB - realistic_total_transmitted_tb)
    realistic_transmit_seconds = realistic_transmit_steps * time_step_seconds
    
    avg_realistic_elevation = np.mean(realistic_elevation_stats) if realistic_elevation_stats else 0
    min_realistic_elevation = np.min(realistic_elevation_stats) if realistic_elevation_stats else 0
    max_realistic_elevation = np.max(realistic_elevation_stats) if realistic_elevation_stats else 0
    avg_realistic_rate = np.mean(realistic_rate_stats) if realistic_rate_stats else 0
    min_realistic_rate = np.min(realistic_rate_stats) if realistic_rate_stats else 0
    max_realistic_rate = np.max(realistic_rate_stats) if realistic_rate_stats else 0
    
    avg_aggregated_rate = np.mean(aggregated_rate_stats) if aggregated_rate_stats else 0
    min_aggregated_rate = np.min(aggregated_rate_stats) if aggregated_rate_stats else 0
    max_aggregated_rate = np.max(aggregated_rate_stats) if aggregated_rate_stats else 0
    avg_simultaneous = np.mean(simultaneous_transmission_stats) if simultaneous_transmission_stats else 0
    max_simultaneous = np.max(simultaneous_transmission_stats) if simultaneous_transmission_stats else 0
    
    # Print results
    print(f"\n" + "="*70)
    print("PORÓWNANIE MODELI TRANSMISJI")
    print("="*70)
    print(f"Model obecny (najlepsza para):        {total_transmitted_tb:.2f} TB ({(total_transmitted_tb/TOTAL_DATA_TB)*100:.1f}%)")
    print(f"Model realistyczny (1 sat -> 1 st):  {realistic_total_transmitted_tb:.2f} TB ({(realistic_total_transmitted_tb/TOTAL_DATA_TB)*100:.1f}%)")
    print(f"Średnia zagregowana przepustowość:    {avg_aggregated_rate:.1f} Mbps")
    print(f"Maksymalna zagregowana przepustowość:  {max_aggregated_rate:.1f} Mbps")
    print(f"Średnia liczba równoczesnych transmisji: {avg_simultaneous:.1f}")
    print(f"Maksymalna liczba równoczesnych transmisji: {max_simultaneous}")
    print("="*70)

    print(f"\n=== COVERAGE ANALYSIS RESULTS ===")
    print(f"Total simulation steps: {num_steps}")
    print(f"Day-Land coverage events: {total_day_land}")
    print(f"Day-Ocean coverage events: {total_day_ocean}")
    print(f"Total coverage events: {total_day_land + total_day_ocean}")
    
    print(f"\n=== TRANSMISSION SELECTION STATISTICS ===")
    print(f"Average potential satellite-station pairs per step: {avg_potential_pairs:.1f}")
    print(f"Maximum potential pairs in single step: {max_potential_pairs}")
    print(f"Steps with transmission capability: {sum(1 for x in potential_pairs_stats if x > 0)}")
    print(f"Steps with actual transmission: {total_transmit_steps}")
    print(f"Transmission selection efficiency: {(total_transmit_steps / max(1, sum(1 for x in potential_pairs_stats if x > 0))) * 100:.1f}%")
    
    print(f"\n=== ELEVATION & TRANSMISSION STATISTICS ===")
    print(f"Average elevation angle: {avg_elevation:.1f}°")
    print(f"Elevation range: {min_elevation:.1f}° - {max_elevation:.1f}°")
    print(f"Average transmission rate: {avg_rate:.1f} Mbps")
    print(f"Transmission rate range: {min_rate:.1f} - {max_rate:.1f} Mbps")
    
    print(f"\n=== DATA TRANSMISSION RESULTS ===")
    print(f"Steps with active transmission: {total_transmit_steps}")
    print(f"Total transmission time: {total_transmit_seconds} seconds ({total_transmit_seconds/3600:.2f} hours)")
    print(f"Data transmitted: {total_transmitted_tb:.2f} TB")
    print(f"Data remaining: {data_left_tb:.2f} TB / {TOTAL_DATA_TB} TB")
    print(f"Transmission completion: {(total_transmitted_tb/TOTAL_DATA_TB)*100:.1f}%")
    
    # Save detailed results to CSV
    with open('coverage_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Metric', 'Value'])
        writer.writerow(['Total_simulation_steps', num_steps])
        writer.writerow(['Day_land_coverage_events', total_day_land])
        writer.writerow(['Day_ocean_coverage_events', total_day_ocean])
        writer.writerow(['Total_coverage_events', total_day_land + total_day_ocean])
        writer.writerow(['Transmission_steps', total_transmit_steps])
        writer.writerow(['Transmission_time_seconds', total_transmit_seconds])
        writer.writerow(['Transmission_time_hours', total_transmit_seconds/3600])
        writer.writerow(['Data_transmitted_TB', total_transmitted_tb])
        writer.writerow(['Data_remaining_TB', data_left_tb])
        writer.writerow(['Transmission_completion_percent', (total_transmitted_tb/TOTAL_DATA_TB)*100])
        writer.writerow(['Average_elevation_deg', avg_elevation])
        writer.writerow(['Min_elevation_deg', min_elevation])
        writer.writerow(['Max_elevation_deg', max_elevation])
        writer.writerow(['Average_transmission_rate_Mbps', avg_rate])
        writer.writerow(['Min_transmission_rate_Mbps', min_rate])
        writer.writerow(['Max_transmission_rate_Mbps', max_rate])
        writer.writerow(['Average_potential_pairs_per_step', avg_potential_pairs])
        writer.writerow(['Max_potential_pairs_single_step', max_potential_pairs])
        writer.writerow(['Steps_with_transmission_capability', sum(1 for x in potential_pairs_stats if x > 0)])
        writer.writerow(['Transmission_selection_efficiency_percent', (total_transmit_steps / max(1, sum(1 for x in potential_pairs_stats if x > 0))) * 100])
        
        # NEW: Realistic transmission metrics
        writer.writerow(['Realistic_data_transmitted_TB', realistic_total_transmitted_tb])
        writer.writerow(['Realistic_data_remaining_TB', realistic_data_left_tb])
        writer.writerow(['Realistic_transmission_completion_percent', (realistic_total_transmitted_tb/TOTAL_DATA_TB)*100])
        writer.writerow(['Realistic_transmission_steps', realistic_transmit_steps])
        writer.writerow(['Realistic_transmission_time_seconds', realistic_transmit_seconds])
        writer.writerow(['Realistic_transmission_time_hours', realistic_transmit_seconds/3600])
        writer.writerow(['Average_aggregated_rate_Mbps', avg_aggregated_rate])
        writer.writerow(['Min_aggregated_rate_Mbps', min_aggregated_rate])
        writer.writerow(['Max_aggregated_rate_Mbps', max_aggregated_rate])
        writer.writerow(['Average_simultaneous_transmissions', avg_simultaneous])
        writer.writerow(['Max_simultaneous_transmissions', max_simultaneous])
        writer.writerow(['Realistic_average_elevation_deg', avg_realistic_elevation])
        writer.writerow(['Realistic_min_elevation_deg', min_realistic_elevation])
        writer.writerow(['Realistic_max_elevation_deg', max_realistic_elevation])
        writer.writerow(['Realistic_average_transmission_rate_Mbps', avg_realistic_rate])
        writer.writerow(['Realistic_min_transmission_rate_Mbps', min_realistic_rate])
        writer.writerow(['Realistic_max_transmission_rate_Mbps', max_realistic_rate])
    
    # Save transmission time summary
    with open('czas_w_zasiegu.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sumaryczny_czas_transmisji_s'])
        writer.writerow([total_transmit_seconds])
    
    # Save detailed transmission log
    with open('transmission_details.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Frame', 'Satellite', 'Station', 'Elevation_deg', 'Rate_Mbps', 'Distance_km', 'Selected'])
        for result in results:
            if result['transmitting']:
                for pair in result['transmitting_pairs']:
                    writer.writerow([
                        result['frame'],
                        pair['satellite'],
                        pair['station'],
                        round(pair['elevation'], 2),
                        round(pair['rate_mbps'], 1),
                        round(pair['distance_km'], 1),
                        'YES'
                    ])
    
    # NEW: Save realistic detailed transmission log
    with open('realistic_transmission_details.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Frame', 'Satellite', 'Station', 'Elevation_deg', 'Rate_Mbps', 'Distance_km', 'Aggregated_Rate_Mbps', 'Simultaneous_Count'])
        for result in results:
            if result.get('realistic_transmitting', False):
                aggregated_rate = result['aggregated_transmission_rate']
                simultaneous_count = result['simultaneous_transmissions']
                for pair in result['realistic_transmitting_pairs']:
                    writer.writerow([
                        result['frame'],
                        pair['satellite'],
                        pair['station'],
                        round(pair['elevation'], 2),
                        round(pair['rate_mbps'], 1),
                        round(pair['distance_km'], 1),
                        round(aggregated_rate, 1),
                        simultaneous_count
                    ])
    
    print(f"\nResults saved to:")
    print(f"- 'coverage_results.csv' (summary)")
    print(f"- 'czas_w_zasiegu.csv' (transmission time)")
    print(f"- 'transmission_details.csv' (detailed transmission log)")
    print(f"- 'realistic_transmission_details.csv' (realistic transmission log)")


# Keep these functions for backward compatibility but make them no-ops
def plot_static(satellites, sky_datetime_objects):
    """Disabled - use calculate_static_coverage() instead."""
    print("Visual plotting disabled. Use calculate_static_coverage() for text results.")
    calculate_static_coverage(satellites)


def plot_coverage_overlay(satellites):
    """Disabled - use calculate_static_coverage() instead."""
    print("Visual plotting disabled. Use calculate_static_coverage() for text results.")
    calculate_static_coverage(satellites)


def animate_coverage(satellites, sky_datetime_objects):
    """Disabled - use calculate_coverage_fast() instead."""
    print("Animation disabled. Running fast calculations instead...")
    calculate_coverage_fast(satellites, sky_datetime_objects)


def animate_modis_only(satellites, sky_datetime_objects, output_filename="modis_animation.mp4"):
    """Disabled - use calculate_coverage_fast() instead."""
    print("Animation disabled. Running fast calculations instead...")
    calculate_coverage_fast(satellites, sky_datetime_objects)