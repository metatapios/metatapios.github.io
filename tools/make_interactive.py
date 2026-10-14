#!/usr/bin/env python3
"""Build an interactive (rotate/zoom/pan) 3D orbit with Plotly.
Writes results/orbit_3d.html (standalone) and results/orbit_widget.html (fragment).
Usage: python tools/make_interactive.py [id]"""
import sys, json
import numpy as np
import pandas as pd

sat_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
R = 1737.4
d = pd.read_csv("results/trajectories.csv")
d = d[d.id == sat_id].reset_index(drop=True)
s = pd.read_csv("results/summary.csv"); s = s[s.id == sat_id].iloc[0]
impacted = int(s.impacted) == 1

# downsample + round to keep the inline widget payload small
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 700
step = max(1, len(d) // npts)
d = d.iloc[::step].reset_index(drop=True)
rnd = lambda a: [round(float(v), 1) for v in a]
x, y, z = rnd(d.x/1000), rnd(d.y/1000), rnd(d.z/1000)
t = [round(float(v), 3) for v in d.t_hr/24]

# Moon sphere mesh (coarse, integer km)
u = np.linspace(0, 2*np.pi, 30); v = np.linspace(0, np.pi, 18)
ri = lambda M: [[int(round(c)) for c in row] for row in M]
sx = ri(R*np.outer(np.cos(u), np.sin(v)))
sy = ri(R*np.outer(np.sin(u), np.sin(v)))
sz = ri(R*np.outer(np.ones_like(u), np.cos(v)))

data = dict(x=x, y=y, z=z, t=t, sx=sx, sy=sy, sz=sz,
            x0=x[0], y0=y[0], z0=z[0],
            xi=x[-1], yi=y[-1], zi=z[-1], impacted=impacted)
J = json.dumps(data)

BODY = """
<div id="orbit3d" style="width:100%;height:560px;"></div>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
(function(){
  var D = __DATA__;
  function draw(){
    if(!window.Plotly){return setTimeout(draw,60);}
    var orbit={type:'scatter3d',mode:'lines',x:D.x,y:D.y,z:D.z,
      line:{width:4,color:D.t,colorscale:'Viridis',
        colorbar:{title:'time [days]',thickness:12}},name:'orbit'};
    var moon={type:'surface',x:D.sx,y:D.sy,z:D.sz,opacity:0.30,showscale:false,
      colorscale:[[0,'#888'],[1,'#aaa']],hoverinfo:'skip',name:'Moon'};
    var start={type:'scatter3d',mode:'markers',x:[D.x0],y:[D.y0],z:[D.z0],
      marker:{size:6,color:'green'},name:'start'};
    var traces=[moon,orbit,start];
    if(D.impacted){traces.push({type:'scatter3d',mode:'markers',
      x:[D.xi],y:[D.yi],z:[D.zi],
      marker:{size:8,color:'red',symbol:'diamond'},name:'impact'});}
    var L={margin:{l:0,r:0,t:10,b:0},
      scene:{xaxis:{title:'x [km]'},yaxis:{title:'y [km]'},zaxis:{title:'z [km]'},
        aspectmode:'data'},
      legend:{x:0,y:1}};
    Plotly.newPlot('orbit3d',traces,L,{responsive:true});
  }
  draw();
})();
</script>
""".replace("__DATA__", J)

open("results/orbit_widget.html", "w").write(BODY)
open("results/orbit_3d.html", "w").write(
    "<!doctype html><meta charset=utf-8><title>Lunar orbit</title>" + BODY)
print("wrote results/orbit_3d.html and results/orbit_widget.html "
      f"({len(x)} pts, impacted={impacted})")
