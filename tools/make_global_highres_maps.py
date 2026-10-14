#!/usr/bin/env python3
import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Render original resolved topography, not a 59 m SLDEM2015 product.
with open("data/gravity.bin", "rb") as f:
    _degree = struct.unpack("<i", f.read(4))[0]
    _GM, R_MOON = struct.unpack("<dd", f.read(16))

# 1. Load ALL impact sites from summary
print("--> Loading impact data...")
s = pd.read_csv("results/summary.csv")
df_impact = s[s.impacted == 1]

# Convert longitudes from [-180, 180] to [0, 360] to match the topo file frame
imp_lon = (df_impact['impact_lon'].values) % 360
imp_lat = df_impact['impact_lat'].values
imp_time = df_impact['t_end_hr'].values

# 2. Load the binary topography grid WITHOUT downsampling [::4, ::4]
print("--> Loading raw topography binary...")
with open("data/topo.bin", "rb") as f:
    nlat, nlon = struct.unpack("<ii", f.read(8))
    topo = np.fromfile(f, "<f8", nlat*nlon).reshape(nlat, nlon)

# Calculate true height in km above mean radius
elev_km = (topo - R_MOON) / 1000.0

# 3. Render at high DPI; the underlying physical surface resolution is unchanged.
print("--> Generating crisp multi-panel figures...")
fig, ax = plt.subplots(figsize=(14, 8), dpi=200) # High physical sizing + high DPI

# Calculate percentile limits for terrain contrast balance
ev = np.percentile(np.abs(elev_km), 99)

# Show the background map at full resolution
im = ax.imshow(elev_km, extent=[0, 360, -90, 90], origin="upper",
               cmap="gist_earth", vmin=-ev, vmax=ev, aspect="auto")

# Overplot ALL impact points as crisp, clear cross markers color-coded by time
sc = ax.scatter(imp_lon, imp_lat, c=imp_time, cmap="plasma", 
                marker="x", s=60, linewidths=2.0, zorder=5, label="Impact Site")

# Formatting crisp labels and grid lines
ax.set_xlabel("Body-Fixed Longitude [deg]", fontsize=12, labelpad=10)
ax.set_ylabel("Latitude [deg]", fontsize=12)
ax.set_title("Global Satellite Impacts over Sampled LOLA Surface Grid", 
             fontsize=14, weight="bold", pad=15)
ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)
ax.grid(color='black', linestyle=':', alpha=0.4)

# Colorbars
cb_map = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb_map.set_label("Elevation Relative to Mean Radius [km]", fontsize=11)

cb_pts = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.12, shrink=0.6)
cb_pts.set_label("Satellite Orbital Lifespan Until Crash [hours]", fontsize=11)

# Save to disk as a lossless uncompressed PNG
output_path = "results/highres_global_topography.png"
fig.savefig(output_path, bbox_inches='tight')
plt.close(fig)

print(f"--> Rendered {output_path} from the original {nlat}x{nlon} shape grid.")
print("--> Higher image DPI does NOT constitute high-resolution impact modeling.")
