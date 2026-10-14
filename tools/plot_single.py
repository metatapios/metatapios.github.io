#!/usr/bin/env python3
"""Plots for a SINGLE trajectory from results/trajectories.csv.
  A results/single_run_panels.png    time series + ground track + phase (continuous)
  B results/single_run_3d_phase.png  3D orbit + 3D momentum hodograph
  C results/single_run_elements.png  osculating elements: precession & inclination
Impact marker added where the run crashed (from results/summary.csv).
Usage: python tools/plot_single.py [id] [traj_csv] [summary_csv]"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D            # noqa
from mpl_toolkits.mplot3d.art3d import Line3DCollection

sat_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
csv = sys.argv[2] if len(sys.argv) > 2 else "results/trajectories.csv"
summ = sys.argv[3] if len(sys.argv) > 3 else "results/summary.csv"
R_MOON_KM = 1737.4
GM = 4.902800e12

d = pd.read_csv(csv)
d = d[d.id == sat_id].reset_index(drop=True)
if d.empty:
    sys.exit(f"no rows for id={sat_id} in {csv}")
s = pd.read_csv(summ)
s = s[s.id == sat_id].iloc[0]
impacted = int(s.impacted) == 1

x, y, z = d.x.values, d.y.values, d.z.values
vx, vy, vz = d.vx.values, d.vy.values, d.vz.values
r = np.sqrt(x**2 + y**2 + z**2)
v_r = (x*vx + y*vy + z*vz) / r
r_km = r / 1000.0
t_d = d.t_hr.values / 24.0


def color_line(ax, xa, ya, t, cmap, lw=1.0, wrap=None):
    """Continuous line coloured by time via a LineCollection.
    wrap: if set, drop segments where |dx| > wrap (avoids longitude 0/360
    wraparound streaks across the map)."""
    xa = np.asarray(xa); ya = np.asarray(ya)
    pts = np.array([xa, ya]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
    carr = t[:-1]
    if wrap is not None:
        keep = np.abs(np.diff(xa)) <= wrap
        seg = seg[keep]; carr = np.asarray(carr)[keep]
    lc = LineCollection(seg, cmap=cmap, array=carr, linewidth=lw)
    ax.add_collection(lc)
    ax.set_xlim(xa.min(), xa.max()); ax.set_ylim(ya.min(), ya.max())
    return lc


def mark_impact(ax, xi, yi):
    if impacted:
        ax.plot(xi, yi, marker="*", color="red", ms=16, mec="k", zorder=5,
                label="impact")


# ---------- osculating orbital elements ----------
hx = y*vz - z*vy; hy = z*vx - x*vz; hz = x*vy - y*vx        # angular momentum
h = np.sqrt(hx**2 + hy**2 + hz**2)
inc = np.degrees(np.arccos(np.clip(hz/h, -1, 1)))          # inclination
raan = np.degrees(np.arctan2(hx, -hy))                     # ascending node
# eccentricity vector -> argument of periapsis
ex = (vy*hz - vz*hy)/GM - x/r
ey = (vz*hx - vx*hz)/GM - y/r
ez = (vx*hy - vy*hx)/GM - z/r
ecc = np.sqrt(ex**2 + ey**2 + ez**2)
nx, ny = -hy, hx                                           # node vector (z x h)
nmag = np.sqrt(nx**2 + ny**2)
argp = np.degrees(np.arccos(np.clip((nx*ex + ny*ey)/(nmag*ecc), -1, 1)))
argp = np.where(ez < 0, 360 - argp, argp)

# ---------- Figure A ----------
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
ax[0, 0].plot(t_d, d.alt_km, lw=0.7, color="tab:blue")
ax[0, 0].axhline(0, color="k", ls="--", lw=0.8); mark_impact(ax[0, 0], t_d[-1], d.alt_km.values[-1])
ax[0, 0].set(xlabel="time [days]", ylabel="altitude [km]",
             title="Altitude (perilune plunge)"); ax[0, 0].legend(loc="upper left")

ax[0, 1].plot(t_d, d.ecc, lw=0.7, color="tab:red")
ax[0, 1].set(xlabel="time [days]", ylabel="eccentricity", title="Eccentricity growth")

inv0 = d.invariant.iloc[0]
ax[0, 2].plot(t_d, (d.invariant - inv0)/abs(inv0), lw=0.7, color="tab:green")
ax[0, 2].set(xlabel="time [days]", ylabel="(I - I0)/|I0|",
             title="Invariant (Jacobi) conservation")
ax[0, 2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

lc = color_line(ax[1, 0], d.lon.values % 360, d.lat.values, t_d, "viridis", 0.8, wrap=180)
plt.colorbar(lc, ax=ax[1, 0], label="time [days]", fraction=0.04)
mark_impact(ax[1, 0], s.impact_lon % 360, s.impact_lat)
ax[1, 0].set(xlabel="body-fixed longitude [deg]", ylabel="latitude [deg]",
             xlim=(0, 360), ylim=(-90, 90), title="Ground track"); ax[1, 0].legend(loc="upper right")

lc = color_line(ax[1, 1], r_km, v_r, t_d, "plasma", 0.6)   # continuous phase curve
plt.colorbar(lc, ax=ax[1, 1], label="time [days]", fraction=0.04)
mark_impact(ax[1, 1], r_km[-1], v_r[-1])
ax[1, 1].set(xlabel="radius r [km]", ylabel="radial velocity  v_r [m/s]",
             title="Radial phase portrait (r, p_r)"); ax[1, 1].legend(loc="upper right")

lc = color_line(ax[1, 2], x/1000, vx, t_d, "cividis", 0.6)
plt.colorbar(lc, ax=ax[1, 2], label="time [days]", fraction=0.04)
mark_impact(ax[1, 2], x[-1]/1000, vx[-1])
ax[1, 2].set(xlabel="x [km]", ylabel="p_x = v_x [m/s]", title="Phase portrait (x, p_x)")

fig.suptitle(f"Lunar orbit - satellite id {sat_id}", fontsize=14, weight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.98]); fig.savefig("results/single_run_panels.png", dpi=130)
print("wrote results/single_run_panels.png")

# ---------- Figure B: 3D ----------
fig2 = plt.figure(figsize=(15, 7))
axo = fig2.add_subplot(1, 2, 1, projection="3d")
u, w = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
axo.plot_surface(R_MOON_KM*np.cos(u)*np.sin(w), R_MOON_KM*np.sin(u)*np.sin(w),
                 R_MOON_KM*np.cos(w), color="gray", alpha=0.25, linewidth=0)
pts = np.array([x/1000, y/1000, z/1000]).T.reshape(-1, 1, 3)
seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
lc3 = Line3DCollection(seg, cmap="viridis", array=t_d[:-1], linewidth=0.7)
axo.add_collection3d(lc3)
axo.scatter([x[0]/1000], [y[0]/1000], [z[0]/1000], color="green", s=40, label="start")
if impacted:
    axo.scatter([x[-1]/1000], [y[-1]/1000], [z[-1]/1000], color="red", marker="*", s=120,
                edgecolor="k", label="impact")
lim = R_MOON_KM*1.3
axo.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim, lim),
        xlabel="x [km]", ylabel="y [km]", zlabel="z [km]",
        title="3D orbit (inertial position space)")
axo.legend(); axo.set_box_aspect((1, 1, 1))

axm = fig2.add_subplot(1, 2, 2, projection="3d")
p = axm.scatter(vx, vy, vz, c=t_d, s=1.2, cmap="magma")
plt.colorbar(p, ax=axm, label="time [days]", fraction=0.03, pad=0.1)
axm.set(xlabel="p_x [m/s]", ylabel="p_y [m/s]", zlabel="p_z [m/s]",
        title="3D momentum space (hodograph)"); axm.set_box_aspect((1, 1, 1))
fig2.tight_layout(); fig2.savefig("results/single_run_3d_phase.png", dpi=130)
print("wrote results/single_run_3d_phase.png")

# ---------- Figure C: orbital elements (precession & inclination) ----------
fig3, ax3 = plt.subplots(1, 3, figsize=(16, 4.5))
ax3[0].plot(t_d, inc, lw=0.8, color="tab:purple")
ax3[0].set(xlabel="time [days]", ylabel="inclination [deg]",
           title=f"Inclination  (start {inc[0]:.2f}°, mean {inc.mean():.2f}°)")
ax3[1].plot(t_d, np.unwrap(np.radians(raan))*180/np.pi, lw=0.8, color="tab:brown")
ax3[1].set(xlabel="time [days]", ylabel="RAAN Omega [deg]",
           title="Nodal precession (ascending node)")
ax3[2].plot(t_d, np.unwrap(np.radians(argp))*180/np.pi, lw=0.8, color="tab:olive")
ax3[2].set(xlabel="time [days]", ylabel="arg. periapsis omega [deg]",
           title="Apsidal precession (periapsis)")
fig3.suptitle(f"Osculating elements - satellite id {sat_id}", fontsize=13, weight="bold")
fig3.tight_layout(rect=[0, 0, 1, 0.95]); fig3.savefig("results/single_run_elements.png", dpi=130)
print("wrote results/single_run_elements.png")

# rate estimates (linear fit over the run)
draan = np.polyfit(t_d, np.unwrap(np.radians(raan))*180/np.pi, 1)[0]
dargp = np.polyfit(t_d, np.unwrap(np.radians(argp))*180/np.pi, 1)[0]
print(f"inclination: {inc.min():.3f}..{inc.max():.3f} deg (swing {inc.max()-inc.min():.3f})")
print(f"nodal regression dOmega/dt  ~ {draan:+.3f} deg/day")
print(f"apsidal precession domega/dt ~ {dargp:+.3f} deg/day")
