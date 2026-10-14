#!/usr/bin/env python3
"""Statistical analysis of a Monte-Carlo run, binned by INITIAL ALTITUDE.
Reads results/summary.csv (one row per satellite) and produces:
  results/montecarlo_stats.png   multi-panel statistics figure
  results/montecarlo_bins.csv    per-bin table
Quantities per altitude bin: impact fraction, satellite lifetime (time to
impact; survivors censored at max sim time), max eccentricity reached, orbits
completed, local interpolated terrain-normal impact angle, ground-relative terminal speed.
Usage: python tools/montecarlo_stats.py [summary_csv] [bin_width_km]"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv = sys.argv[1] if len(sys.argv) > 1 else "results/summary.csv"
bw = float(sys.argv[2]) if len(sys.argv) > 2 else 25.0
s = pd.read_csv(csv)
n = len(s)
tmax_d = s.t_end_hr.max() / 24.0
imp = s[s.impacted == 1]
surv = s[s.impacted == 0]

# altitude bins
lo, hi = np.floor(s.alt0_km.min()), np.ceil(s.alt0_km.max())
edges = np.arange(lo, hi + bw, bw)
s = s.assign(bin=pd.cut(s.alt0_km, edges, right=False))
rows = []
for b, g in s.groupby("bin", observed=True):
    gi = g[g.impacted == 1]
    rows.append(dict(
        alt_lo=b.left, alt_hi=b.right, n=len(g),
        impact_frac=g.impacted.mean(),
        life_med_d=(gi.t_end_hr.median() / 24.0) if len(gi) else np.nan,
        life_mean_d=(gi.t_end_hr.mean() / 24.0) if len(gi) else np.nan,
        ecc_med=g.ecc_max.median(), ecc_mean=g.ecc_max.mean(),
        orbits_med=g.orbit_count.median(),
        vterm_med=(gi.terminal_speed_ms.median() if len(gi) else np.nan)))
bins = pd.DataFrame(rows)
bins.to_csv("results/montecarlo_bins.csv", index=False)
bc = (bins.alt_lo + bins.alt_hi) / 2

print(f"N={n}  impacts={len(imp)} ({100*len(imp)/n:.1f}%)  "
      f"survivors={len(surv)}  max sim={tmax_d:.0f} d")
print(bins.to_string(index=False,
      float_format=lambda x: f"{x:.3f}"))

# ---------------- figure ----------------
fig, ax = plt.subplots(2, 3, figsize=(16, 9))

# 1. impact fraction per altitude bin
ax[0, 0].bar(bc, bins.impact_frac, width=bw*0.9, color="tab:red", alpha=0.75,
             edgecolor="k")
ax[0, 0].set(xlabel="initial altitude [km]", ylabel="impact fraction",
             ylim=(0, 1.05), title="Impact probability vs start altitude")
for x, nb in zip(bc, bins.n):
    ax[0, 0].text(x, 0.02, f"n={nb}", ha="center", va="bottom", fontsize=7, rotation=90)

# 2. lifetime vs initial altitude (scatter + binned median)
ax[0, 1].scatter(imp.alt0_km, imp.t_end_hr/24, s=10, c="tab:blue", label="impacted")
if len(surv):
    ax[0, 1].scatter(surv.alt0_km, surv.t_end_hr/24, s=12, marker="^",
                     c="tab:green", label=f"survived ({tmax_d:.0f} d)")
ax[0, 1].plot(bc, bins.life_med_d, "o-", c="k", lw=1.5, label="bin median")
ax[0, 1].set(xlabel="initial altitude [km]", ylabel="lifetime [days]",
             title="Satellite lifetime vs start altitude"); ax[0, 1].legend(fontsize=8)

# 3. lifetime distribution (impacted)
ax[0, 2].hist(imp.t_end_hr/24, bins=30, color="tab:purple", alpha=0.8)
ax[0, 2].set(xlabel="time to impact [days]", ylabel="count",
             title="Lifetime distribution (impacted)")

# 4. final / max eccentricity vs altitude
ax[1, 0].scatter(s.alt0_km, s.ecc_max, s=8, c=s.impacted, cmap="coolwarm")
ax[1, 0].plot(bc, bins.ecc_med, "o-", c="k", lw=1.5, label="bin median")
ax[1, 0].set(xlabel="initial altitude [km]", ylabel="max eccentricity",
             title="Max eccentricity vs start altitude"); ax[1, 0].legend(fontsize=8)

# 5. distribution of max eccentricity
ax[1, 1].hist(s.ecc_max, bins=30, color="tab:orange", alpha=0.8)
ax[1, 1].set(xlabel="max eccentricity", ylabel="count",
             title="Max-eccentricity distribution")

# 6. impact geometry: flight-path angle & terminal speed
if len(imp):
    ax[1, 2].hist(imp.impact_angle_deg, bins=25, color="tab:cyan", alpha=0.85)
    ax[1, 2].set(xlabel="local terrain-normal impact angle [deg]", ylabel="count",
                 title=f"Impact angle (v_term med "
                       f"{imp.terminal_speed_ms.median():.0f} m/s)")

fig.suptitle(f"Monte-Carlo statistics - N={n}, altitude [{lo:.0f},{hi:.0f}] km, "
             f"bin {bw:.0f} km", fontsize=14, weight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig("results/montecarlo_stats.png", dpi=130)
print("wrote results/montecarlo_stats.png and results/montecarlo_bins.csv")
