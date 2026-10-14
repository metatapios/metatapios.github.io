#!/usr/bin/env python3
"""Build a Moon-map widget with a toggle between:
  (1) real topography (surface elevation vs mean radius), and
  (2) mascon map (free-air gravity anomaly, mGal) from GRGM1200B.
Ground track + impact of a trajectory are overlaid on both.
Outputs:
  results/moon_maps.html      standalone (open in browser)
  results/moon_maps_widget.html  inline fragment (base64 JPEGs + toggle)
Usage: python tools/make_map_widget.py [id]"""
import sys, struct, base64, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sat_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
R_MOON = 1737400.0
LMAX_ANOM = 300

# ---- overlay: ground track + impact ----
gt_lon = gt_lat = None; imp = None
try:
    import pandas as pd
    d = pd.read_csv("results/trajectories.csv"); d = d[d.id == sat_id]
    gt_lon = (d.lon.values % 360); gt_lat = d.lat.values
    s = pd.read_csv("results/summary.csv"); s = s[s.id == sat_id].iloc[0]
    if int(s.impacted) == 1:
        imp = (s.impact_lon % 360, s.impact_lat)
except Exception as e:
    print("overlay skipped:", e)

# ---- (1) topography from data/topo.bin ----
with open("data/topo.bin", "rb") as f:
    nlat, nlon = struct.unpack("<ii", f.read(8))
    topo = np.fromfile(f, "<f8", nlat*nlon).reshape(nlat, nlon)
elev_km = ((topo - R_MOON) / 1000.0)[::4, ::4]  # row0=+90 (DH), downsampled

# ---- (2) mascon free-air gravity anomaly from GRGM1200B ----
import pyshtools as pysh
clm = pysh.datasets.Moon.GRGM1200B(lmax=LMAX_ANOM)
GM = float(clm.gm); Rref = float(clm.r0)
gc = clm.copy()
for l in range(2, LMAX_ANOM + 1):
    gc.coeffs[:, l, :] *= (l + 1) * (GM / Rref**2) * 1e5   # -> mGal
gc.coeffs[:, 0, :] = 0.0; gc.coeffs[:, 1, :] = 0.0         # isolate anomalies
anom = pysh.SHCoeffs.from_array(gc.coeffs).expand(grid="DH").data[::4, ::4]  # downsampled


def render(data, cmap, label, title, vmin, vmax):
    fig, ax = plt.subplots(figsize=(4.0, 2.1))
    im = ax.imshow(data, extent=[0, 360, -90, 90], origin="upper",
                   cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    if gt_lon is not None:
        lo, la = gt_lon.copy(), gt_lat.copy()          # continuous line,
        brk = np.where(np.abs(np.diff(lo)) > 180)[0]   # broken at 0/360 seam
        lo = np.insert(lo, brk + 1, np.nan); la = np.insert(la, brk + 1, np.nan)
        ax.plot(lo, la, "-", lw=0.5, color="cyan", alpha=0.75)
    if imp is not None:
        ax.plot(imp[0], imp[1], "*", ms=13, color="red", mec="k", label="impact")
        ax.legend(loc="upper right", fontsize=7)
    ax.set(xlabel="body-fixed longitude [deg]", ylabel="latitude [deg]",
           xlim=(0, 360), ylim=(-90, 90), title=title)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label(label)
    fig.tight_layout()
    buf = io.BytesIO()
    try:
        fig.savefig(buf, format="jpg", dpi=82,
                    pil_kwargs={"quality": 20, "optimize": True})
        mime = "jpeg"
    except Exception:
        fig.savefig(buf, format="png", dpi=85); mime = "png"
    plt.close(fig)
    return mime, base64.b64encode(buf.getvalue()).decode()


ev = np.percentile(np.abs(elev_km), 99)
m_topo, b_topo = render(elev_km, "gist_earth", "elevation [km]",
                        "Real topography (LOLA shape)", -ev, ev)
m_anom, b_anom = render(anom, "seismic", "free-air anomaly [mGal]",
                        f"Mascon map (GRGM1200B, Lmax={LMAX_ANOM})", -400, 400)
print(f"topo img {len(b_topo)//1024} KB, anomaly img {len(b_anom)//1024} KB")

BODY = f"""
<div style="font-family:system-ui,sans-serif;max-width:920px;">
  <div style="margin:4px 0 8px;">
    <label style="margin-right:16px;cursor:pointer;">
      <input type="radio" name="moonmap" value="topo" checked> Real map (altitudes)</label>
    <label style="cursor:pointer;">
      <input type="radio" name="moonmap" value="mascon"> Mascons (gravity anomaly)</label>
  </div>
  <img id="mm_topo" src="data:image/{m_topo};base64,{b_topo}" style="width:100%;display:block;">
  <img id="mm_anom" src="data:image/{m_anom};base64,{b_anom}" style="width:100%;display:none;">
</div>
<script>
(function(){{
  var rs=document.getElementsByName('moonmap');
  function upd(){{
    var v=[...rs].find(r=>r.checked).value;
    document.getElementById('mm_topo').style.display=(v==='topo')?'block':'none';
    document.getElementById('mm_anom').style.display=(v==='mascon')?'block':'none';
  }}
  rs.forEach(r=>r.addEventListener('change',upd)); upd();
}})();
</script>
"""
open("results/moon_maps_widget.html", "w").write(BODY)
open("results/moon_maps.html", "w").write(
    "<!doctype html><meta charset=utf-8><title>Moon maps</title>" + BODY)
print("wrote results/moon_maps.html and results/moon_maps_widget.html")
