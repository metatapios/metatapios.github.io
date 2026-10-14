#!/usr/bin/env python3
"""Interactive 3D orbit with a selectable Moon skin:
  Sphere (plain, code default) | Topography | Mascons (gravity anomaly).
Sphere geometry is built in JS (small payload); only the two colour grids
(topography, anomaly) + the orbit are embedded.
Outputs results/orbit_sphere_3d.html + results/orbit_sphere_widget.html
Usage: python tools/make_orbit_sphere_widget.py [id]"""
import sys, json, struct
import numpy as np
import pandas as pd

sat_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
R_MOON = 1737400.0
NLAT = int(sys.argv[2]) if len(sys.argv) > 2 else 46      # sphere mesh rows
NLON = int(sys.argv[3]) if len(sys.argv) > 3 else 91      # sphere mesh cols
NPTS = int(sys.argv[4]) if len(sys.argv) > 4 else 700     # orbit sample cap
OUT = sys.argv[5] if len(sys.argv) > 5 else "both"        # widget|standalone|both
LMAX_ANOM = 300

# ---- orbit ----
d = pd.read_csv("results/trajectories.csv"); d = d[d.id == sat_id].reset_index(drop=True)
s = pd.read_csv("results/summary.csv"); s = s[s.id == sat_id].iloc[0]
step = max(1, len(d) // NPTS)
d = d.iloc[::step]
rnd = lambda a: [round(float(v), 1) for v in a]
ox, oy, oz = rnd(d.x/1000), rnd(d.y/1000), rnd(d.z/1000)
ot = [round(float(v), 3) for v in d.t_hr/24]
impacted = int(s.impacted) == 1

# ---- topography grid ----
with open("data/topo.bin", "rb") as f:
    tlat, tlon = struct.unpack("<ii", f.read(8))
    topo = np.fromfile(f, "<f8", tlat*tlon).reshape(tlat, tlon)
elev = (topo - R_MOON) / 1000.0     # km, row0=+90 (DH), lon 0..360

# ---- mascon anomaly grid ----
import pyshtools as pysh
clm = pysh.datasets.Moon.GRGM1200B(lmax=LMAX_ANOM)
GM = float(clm.gm); Rref = float(clm.r0); gc = clm.copy()
for l in range(2, LMAX_ANOM + 1):
    gc.coeffs[:, l, :] *= (l + 1) * (GM / Rref**2) * 1e5
gc.coeffs[:, 0, :] = 0.0; gc.coeffs[:, 1, :] = 0.0
anom = pysh.SHCoeffs.from_array(gc.coeffs).expand(grid="DH").data   # row0=+90, lon 0..360


def sample(grid, latd, lond):
    g_nlat, g_nlon = grid.shape
    j = int(round((90.0 - latd) / (180.0 / g_nlat))) % g_nlat
    i = int(round((lond % 360) / (360.0 / g_nlon))) % g_nlon
    return grid[j, i]


# sample both grids onto the sphere mesh (rows=lat, cols=lon)
COLT, COLM = [], []
for j in range(NLAT):
    latd = 90.0 - j * 180.0 / (NLAT - 1)
    rt, rm = [], []
    for i in range(NLON):
        lond = i * 360.0 / (NLON - 1)
        rt.append(int(round(float(sample(elev, latd, lond)))))
        rm.append(int(round(float(sample(anom, latd, lond)))))
    COLT.append(rt); COLM.append(rm)

data = dict(NLAT=NLAT, NLON=NLON, Rkm=R_MOON/1000.0,
            ox=ox, oy=oy, oz=oz, ot=ot, impacted=impacted,
            x0=ox[0], y0=oy[0], z0=oz[0], xi=ox[-1], yi=oy[-1], zi=oz[-1],
            colt=COLT, colm=COLM,
            tmin=float(np.floor(-np.percentile(np.abs(elev), 99))),
            tmax=float(np.ceil(np.percentile(np.abs(elev), 99))))
J = json.dumps(data, separators=(",", ":"))

BODY = """
<div style="font-family:system-ui,sans-serif;">
  <div style="margin:4px 0 6px;">
    <label style="margin-right:14px;cursor:pointer;"><input type="radio" name="skin" value="plain" checked> Sphere</label>
    <label style="margin-right:14px;cursor:pointer;"><input type="radio" name="skin" value="topo"> Topography</label>
    <label style="cursor:pointer;"><input type="radio" name="skin" value="mascon"> Mascons</label>
  </div>
  <div id="orb" style="width:100%;height:560px;"></div>
</div>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
(function(){
  var D=__DATA__;
  function build(){
    if(!window.Plotly){return setTimeout(build,60);}
    var SX=[],SY=[],SZ=[],R=D.Rkm;
    for(var j=0;j<D.NLAT;j++){
      var lat=(90-j*180/(D.NLAT-1))*Math.PI/180, rx=[],ry=[],rz=[];
      for(var i=0;i<D.NLON;i++){
        var lon=(i*360/(D.NLON-1))*Math.PI/180;
        rx.push(R*Math.cos(lat)*Math.cos(lon));
        ry.push(R*Math.cos(lat)*Math.sin(lon));
        rz.push(R*Math.sin(lat));
      }
      SX.push(rx);SY.push(ry);SZ.push(rz);
    }
    var plain=SX.map(r=>r.map(()=>0));
    var moon={type:'surface',x:SX,y:SY,z:SZ,surfacecolor:plain,
      colorscale:[[0,'#8a8a8a'],[1,'#8a8a8a']],showscale:false,
      cmin:0,cmax:1,hoverinfo:'skip',name:'Moon',lighting:{ambient:0.9}};
    var orbit={type:'scatter3d',mode:'lines',x:D.ox,y:D.oy,z:D.oz,
      line:{width:4,color:D.ot,colorscale:'Viridis',
        colorbar:{title:'days',thickness:12,x:1.02,len:0.45,y:0.78,yanchor:'middle'}},
      name:'orbit'};
    var start={type:'scatter3d',mode:'markers',x:[D.x0],y:[D.y0],z:[D.z0],
      marker:{size:5,color:'lime'},name:'start'};
    var traces=[moon,orbit,start];
    if(D.impacted)traces.push({type:'scatter3d',mode:'markers',x:[D.xi],y:[D.yi],z:[D.zi],
      marker:{size:7,color:'red',symbol:'diamond'},name:'impact'});
    var L={margin:{l:0,r:0,t:6,b:0},
      scene:{xaxis:{title:'x [km]'},yaxis:{title:'y [km]'},zaxis:{title:'z [km]'},aspectmode:'data'},
      legend:{x:0,y:1}};
    Plotly.newPlot('orb',traces,L,{responsive:true});

    function skin(v){
      if(v==='plain')
        Plotly.restyle('orb',{surfacecolor:[plain],
          colorscale:[[[0,'#8a8a8a'],[1,'#8a8a8a']]],showscale:[false],cmin:[0],cmax:[1]},[0]);
      else if(v==='topo')
        Plotly.restyle('orb',{surfacecolor:[D.colt],colorscale:['Earth'],
          reversescale:[false],showscale:[true],cmin:[D.tmin],cmax:[D.tmax],
          'colorbar.title':['km'],'colorbar.x':[1.02],'colorbar.len':[0.45],
          'colorbar.y':[0.22],'colorbar.yanchor':['middle'],'colorbar.thickness':[12]},[0]);
      else
        Plotly.restyle('orb',{surfacecolor:[D.colm],colorscale:['RdBu'],
          reversescale:[true],showscale:[true],cmin:[-400],cmax:[400],
          'colorbar.title':['mGal'],'colorbar.x':[1.02],'colorbar.len':[0.45],
          'colorbar.y':[0.22],'colorbar.yanchor':['middle'],'colorbar.thickness':[12]},[0]);
    }
    var rs=document.getElementsByName('skin');
    [].forEach.call(rs,function(r){r.addEventListener('change',function(){skin(this.value);});});
  }
  build();
})();
</script>
""".replace("__DATA__", J)

if OUT in ("widget", "both"):
    open("results/orbit_sphere_widget.html", "w").write(BODY)
if OUT in ("standalone", "both"):
    open("results/orbit_sphere_3d.html", "w").write(
        "<!doctype html><meta charset=utf-8><title>Orbit 3D</title>" + BODY)
print(f"mesh {NLAT}x{NLON}, orbit {len(ox)} pts, {len(BODY)//1024} KB, out={OUT}")
