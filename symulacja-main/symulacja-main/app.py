import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box, Point, MultiPolygon
from shapely.affinity import rotate
# from shapely.ops import unary_union # Not used in animation part
from skyfield.api import load, EarthSatellite
from datetime import datetime, timedelta
from matplotlib.patheffects import withStroke

# --- Configuration & Parameters ---

# Earth parameters
R_earth_km = 6371
mu_km3_s2 = 398600.4418

# Satellite definitions (Using placeholder TLEs - replace with current ones for accuracy)
satellites = {
    "MODIS_Terra": {
        "altitude_km": 617,
        "fov_deg": 110, # Effective FOV for 2330km swath from 705km altitude
        "color": "blue",
        "tle_line1": "1 25994U 99068A   25141.50000000  .00000714  00000-0  39998-4 0  9991",
        "tle_line2": "2 25994  98.2041 100.6060 0001200  80.5935 279.5400 14.57107008270001"
    },
    "WorldView-3": {
        "altitude_km": 617,
        "fov_deg": 1.2, # For ~13.1 km swath at 617 km
        "color": "red",
        "tle_line1": "1 40115U 14045A   25141.50000000  .00000714  00000-0  39998-4 0  9992",
        "tle_line2": "2 40115  97.8743  47.2925 0001529  71.0323  289.1058 14.83390144560002"
    }
}

# Time configuration
start_datetime = datetime(2025, 5, 21, 0, 0, 0)
simulation_duration_hours = 3 # Shorter duration for faster animation generation
end_datetime = start_datetime + timedelta(hours=simulation_duration_hours)
time_step_seconds = 60 # Every 1 minute
total_seconds = int(simulation_duration_hours * 3600)
num_steps = total_seconds // time_step_seconds

# Reception Stations
reception_stations = {
    "Poznan_Station": {"lat": 52.4064, "lon": 16.9252, "marker": "o", "color": "green"},
    "Berlin_Station": {"lat": 52.5200, "lon": 13.4050, "marker": "s", "color": "purple"},
}

# Animation parameters
animation_fps = 10
animation_bitrate = 1800
fov_history_length = 5 # Number of previous FOVs to show with transparency

# --- Skyfield Setup ---
ts = load.timescale()
sky_datetime_objects = [start_datetime + timedelta(seconds=i * time_step_seconds) for i in range(num_steps)]
sky_times = ts.utc(sky_datetime_objects)


# --- Orbit Propagation & FOV Calculation ---
print(f"Simulating from {start_datetime} to {end_datetime} ({simulation_duration_hours} hours) with {num_steps} steps.")

for name, params in satellites.items():
    print(f"\nProcessing satellite: {name}")
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
    print(f"  {name}: Avg Alt={h_actual_avg:.2f} km, Ground FOV Width (approx)={delta_deg_fov_width:.2f} deg")

    for i in range(len(lats)):
        cx, cy = lons[i], lats[i]
        angle_deg = 0
        if 0 < i < len(lats) - 1:
            dx = lons[i+1] - lons[i-1]; dy = lats[i+1] - lats[i-1]
            if abs(dx) > 180: dx = np.sign(dx) * (360 - abs(dx)) if dx !=0 else 0
            angle_deg = np.rad2deg(np.arctan2(dy, dx))
        elif i > 0 : # For last point
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
        except Exception as e:
            fov_polygons_shapely.append(None)
            
    params['fov_polygons_shapely'] = fov_polygons_shapely
    # Pre-calculate jump masks for track plotting in animation
    d_lon_diff = np.diff(params['longitudes'])
    params['jump_mask'] = np.concatenate(([False], np.abs(d_lon_diff) > 300)) # Pad for consistent length with lons/lats

# --- Animation Setup ---
fig = plt.figure(figsize=(16, 9))
ax = plt.axes(projection=ccrs.PlateCarree())
map_boundary_box = box(-180, -90, 180, 90) # For clipping FOVs

print("\nStarting animation generation...")

def animate(frame_num):
    ax.clear() # Clear previous frame
    ax.set_global()
    ax.stock_img(zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='gray', zorder=1)
    ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False, zorder=2, color='lightgray', alpha=0.5)

    current_sim_time = sky_datetime_objects[frame_num]

    for name, params in satellites.items():
        lons_hist = params['longitudes'][:frame_num+1]
        lats_hist = params['latitudes'][:frame_num+1]
        jump_mask_hist = params['jump_mask'][:frame_num+1]
        
        # Plot ground track history up to current frame
        for i in range(len(lons_hist) - 1):
            if not jump_mask_hist[i+1]: # Check jump to next point
                 ax.plot([lons_hist[i], lons_hist[i+1]], [lats_hist[i], lats_hist[i+1]],
                         color=params['color'], linewidth=1.0, transform=ccrs.Geodetic(), zorder=3)
            # else: point already plotted by previous segment or will be by next valid segment

        # Plot FOV history (last few frames)
        for j in range(max(0, frame_num - fov_history_length), frame_num):
            poly = params['fov_polygons_shapely'][j]
            if poly and poly.is_valid and not poly.is_empty:
                visible_poly = poly.intersection(map_boundary_box)
                if visible_poly.is_empty: continue
                if isinstance(visible_poly, MultiPolygon):
                    for p_part in visible_poly.geoms:
                        if p_part.is_valid and not p_part.is_empty:
                            x, y = p_part.exterior.xy
                            ax.plot(x, y, color=params['color'], alpha=0.1, linewidth=0.5, transform=ccrs.Geodetic(), zorder=2)
                else:
                    x, y = visible_poly.exterior.xy
                    ax.plot(x, y, color=params['color'], alpha=0.1, linewidth=0.5, transform=ccrs.Geodetic(), zorder=2)

        # Plot current FOV prominently
        current_fov = params['fov_polygons_shapely'][frame_num]
        if current_fov and current_fov.is_valid and not current_fov.is_empty:
            visible_poly = current_fov.intersection(map_boundary_box)
            if not visible_poly.is_empty:
                if isinstance(visible_poly, MultiPolygon):
                    for p_part in visible_poly.geoms:
                         if p_part.is_valid and not p_part.is_empty:
                            x, y = p_part.exterior.xy
                            ax.plot(x, y, color=params['color'], alpha=0.7, linewidth=1.5, transform=ccrs.Geodetic(), zorder=4)
                else:
                    x, y = visible_poly.exterior.xy
                    ax.plot(x, y, color=params['color'], alpha=0.7, linewidth=1.5, transform=ccrs.Geodetic(), zorder=4)
    
    # Plot Reception Stations (static)
    for st_name, st_info in reception_stations.items():
        ax.plot(st_info["lon"], st_info["lat"],
                marker=st_info["marker"], color=st_info["color"], markersize=8,
                transform=ccrs.Geodetic(), linestyle='', zorder=5)
        ax.text(st_info["lon"] + 0.5, st_info["lat"] + 0.5, st_name,
                transform=ccrs.Geodetic(), color=st_info["color"], fontsize=8, zorder=6,
                path_effects=[withStroke(linewidth=1.5, foreground='white')])

    # Update title with current time
    ax.set_title(f"Dual Satellite Simulation: {current_sim_time.strftime('%Y-%m-%d %H:%M:%S')} UTC\nMODIS (Blue), WorldView-3 (Red)", fontsize=10)
    
    # Return list of artists to be redrawn if blit=True, though blit=False is simpler here
    return []


# Create and save animation
anim = FuncAnimation(fig, animate, frames=num_steps, interval=1000/animation_fps, blit=False)
writer = FFMpegWriter(fps=animation_fps, metadata=dict(artist='Satellite Simulator'), bitrate=animation_bitrate)

output_filename = "dual_satellite_animation.mp4"
try:
    anim.save(output_filename, writer=writer)
    print(f"\nAnimation saved as {output_filename}")
except Exception as e:
    print(f"\nError saving animation: {e}")
    print("Please ensure FFMpeg is installed and in your system's PATH.")

plt.close(fig) # Close the plot figure to free memory
print("Processing complete.")