#!/usr/bin/env python3
"""Ensure the master KL switch truly controls the perturbations."""
import csv
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

exe, data_dir = sys.argv[1:3]
with tempfile.TemporaryDirectory() as td:
    td=Path(td)
    base=[exe,'--grav',str(Path(data_dir)/'gravity.bin'),
          '--topo',str(Path(data_dir)/'topo.bin'), '--n','1',
          '--lmax','10', '--days','0.1','--alt-min','100',
          '--alt-max','100', '--steps-per-orbit','120','--traj-count','1',
          '--seed','1234']
    env=dict(os.environ,OMP_NUM_THREADS='1')
    def simulate(label, args, ok=True):
        out=td/f'{label}.csv'
        cmd=base+['--out',str(out),'--traj',str(td/f'{label}_traj.csv'),
                 '--log',str(td/f'{label}.log')]+args
        run=subprocess.run(cmd,capture_output=True,text=True,env=env)
        assert (run.returncode==0)==ok, (label,run.returncode,run.stderr)
        if not ok: return run
        with out.open() as f: return next(csv.DictReader(f))
    baseline=simulate('default',[])
    disabled=simulate('explicit_off',['--kozai-lidov','off'])
    assert baseline==disabled, 'default must be precisely lunar-only'
    switched_on=simulate('enabled',['--kozai-lidov','on'])
    assert math.isnan(float(switched_on['invariant_drift_rel']))
    assert math.isfinite(float(switched_on['kl_jz0']))
    assert math.isnan(float(baseline['kl_jz0']))
    assert switched_on['lunar_jacobi_variation_rel']!=baseline['lunar_jacobi_variation_rel']
    assert 'requires --kozai-lidov on' in simulate('bad_body',[
        '--third-bodies','earth'],False).stderr
    assert 'while --kozai-lidov is off' in simulate('bad_ephem',[
        '--ephem-mode','analytic'],False).stderr
    assert 'on needs at least one' in simulate('bad_empty',[
        '--kozai-lidov','on','--third-bodies','none'],False).stderr
    assert 'accepts on or off' in simulate('bad_switch',[
        '--kozai-lidov','maybe'],False).stderr
    just_earth=simulate('just_earth',['--kozai-lidov','on',
                                      '--third-bodies','earth'])
    assert math.isfinite(float(just_earth['kl_jz0']))
    alias=simulate('alias',['--enable-kozai-lidov'])
    assert alias==switched_on, 'alias must match master switch exactly'
print('PASS: opt-in KL switch, default off, selected bodies, validation and aliases')
