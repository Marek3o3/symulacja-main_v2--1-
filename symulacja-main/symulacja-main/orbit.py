# orbit.py
"""
Orbit propagation and FOV calculation functions for satellite simulation.
"""
import numpy as np
import math
from shapely.geometry import box
from shapely.affinity import rotate
from skyfield.api import load, EarthSatellite
from config import R_earth_km, start_datetime, time_step_seconds, simulation_duration_hours
from datetime import timedelta

def propagate_orbits(satellites):
    ts = load.timescale()
    total_seconds = int(simulation_duration_hours * 3600)
    num_steps = total_seconds // time_step_seconds
    sky_datetime_objects = [start_datetime + timedelta(seconds=i * time_step_seconds) for i in range(num_steps)]
    sky_times = ts.utc(sky_datetime_objects)
    
    for name, params in satellites.items():
        satellite_obj = EarthSatellite(params["tle_line1"], params["tle_line2"], name, ts)
        geocentric = satellite_obj.at(sky_times)
        params['latitudes'] = geocentric.subpoint().latitude.degrees
        params['longitudes'] = geocentric.subpoint().longitude.degrees
        params['altitudes_actual_km'] = geocentric.subpoint().elevation.km
        lats = params['latitudes']
        lons = params['longitudes']
        alts_actual = params['altitudes_actual_km']
        fov_polygons_shapely = []
        h_actual_avg = np.mean(alts_actual)
        theta_ground = 2 * np.arcsin((R_earth_km / (R_earth_km + h_actual_avg)) * np.sin(np.deg2rad(params['fov_deg'] / 2)))
        delta_deg_fov_width = np.rad2deg(theta_ground)
        for i in range(len(lats)):
            cx, cy = lons[i], lats[i]
            angle_deg = 0
            if 0 < i < len(lats) - 1:
                dx = lons[i+1] - lons[i-1]; dy = lats[i+1] - lats[i-1]
                if abs(dx) > 180: dx = np.sign(dx) * (360 - abs(dx)) if dx !=0 else 0
                angle_deg = np.rad2deg(np.arctan2(dy, dx))
            elif i > 0 :
                dx = lons[i] - lons[i-1]; dy = lats[i] - lats[i-1]
                if abs(dx) > 180: dx = np.sign(dx) * (360 - abs(dx)) if dx !=0 else 0
                angle_deg = np.rad2deg(np.arctan2(dy, dx))
            fov_side_deg = delta_deg_fov_width
            current_half_lat_span = fov_side_deg / 2.0
            current_half_lon_span = (fov_side_deg / 2.0) / np.cos(np.deg2rad(cy)) if np.cos(np.deg2rad(cy)) > 0.05 else (fov_side_deg / 2.0)
            current_half_lon_span = min(current_half_lon_span, 60)
            try:
                rect_shapely = box(cx - current_half_lon_span, cy - current_half_lat_span,
                                   cx + current_half_lon_span, cy + current_half_lat_span)
                rotated_fov = rotate(rect_shapely, angle_deg, origin=(cx, cy), use_radians=False)
                fov_polygons_shapely.append(rotated_fov)
            except Exception:
                fov_polygons_shapely.append(None)
        params['fov_polygons_shapely'] = fov_polygons_shapely
        d_lon_diff = np.diff(params['longitudes'])
        params['jump_mask'] = np.concatenate(([False], np.abs(d_lon_diff) > 300))
    return satellites, sky_datetime_objects 