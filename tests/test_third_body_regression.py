#!/usr/bin/env python3
"""End-to-end third-body CSV, corrected diagnostic, and integration regression."""
import csv
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

executable, data_folder = sys.argv[1:3]
with tempfile.TemporaryDirectory() as d:
    d=Path(d)
    ephem=d/'external.csv'
    with ephem.open('w') as f:
        f.write('time_s,earth_x_m,earth_y_m,earth_z_m,sun_x_m,sun_y_m,sun_z_m,'
                'jupiter_x_m,jupiter_y_m,jupiter_z_m\n')
        # Constant external body positions in the common inertial frame.
        # Two days of padding covers Yoshida's non-monotonic stage times.
        for t in range(-172800, 172801, 1800):
            theta=2*math.pi*t/(27.321661*86400)
            earth_x=384400000*math.cos(theta)
            earth_y=384400000*math.sin(theta)
            f.write(f'{t},{earth_x:.8f},{earth_y:.8f},0,149597870700,0,0,'
                    f'778570000000,0,0\n')
    outputs=[]
    for i,extra in enumerate([[],['--kozai-lidov','on','--third-bodies','earth,sun,jupiter','--ephem-file',str(ephem)]]):
        out=d/f'out{i}.csv';trj=d/f'traj{i}.csv'
        cmd=[executable,'--grav',str(Path(data_folder)/'gravity.bin'),
             '--topo',str(Path(data_folder)/'topo.bin'),
             '--n','1','--lmax','10','--alt-min','100','--alt-max','100',
             '--days','0.25','--steps-per-orbit','200','--integrator','yoshida4',
             '--traj-count','1','--traj-dt-min','30','--seed','777',
             '--out',str(out),'--traj',str(trj),'--log',str(d/f'log{i}.txt')]+extra
        env=dict(os.environ,OMP_NUM_THREADS='1')
        subprocess.run(cmd,check=True,env=env,capture_output=True)
        with out.open() as f:
            row=next(csv.DictReader(f))
        with trj.open() as f:
            traj=list(csv.DictReader(f))
        outputs.append((row,traj))
    moon,moon_traj=outputs[0];tide,tide_traj=outputs[1]
    assert moon['impacted']=='0' and tide['impacted']=='0'
    assert len(moon_traj)>=8 and len(tide_traj)>=8
    assert math.isfinite(float(moon['invariant_drift_rel']))
    assert math.isnan(float(tide['invariant_drift_rel'])), 'old Jacobi must not be called conserved'
    assert math.isfinite(float(tide['lunar_jacobi_variation_rel']))
    assert float(tide['lunar_jacobi_variation_rel'])>0
    assert math.isfinite(float(tide['kl_incl_end_deg']))
    assert abs(float(tide['kl_jz_end']))<=1.0
    assert math.isnan(float(moon['kl_jz_end']))
    assert abs(float(moon['t_end_hr'])-6)<1e-8
    assert abs(float(tide['t_end_hr'])-6)<1e-8
    # The perturbations genuinely modify the integration, not merely output flags.
    for coord in ('x','y','z'):
        assert math.isfinite(float(tide_traj[-1][coord]))
    separation=math.sqrt(sum((float(tide_traj[-1][key])-float(moon_traj[-1][key]))**2
                             for key in ('x','y','z')))
    assert separation>0.01, f'external tidal forces did not modify orbit: {separation} m'
print('PASS end-to-end ephemeris, tidal orbit response and non-conserved Jacobi semantics')
