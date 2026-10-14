#!/usr/bin/env python3
"""End-to-end deterministic event and terrain-angle regression test."""
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile

exe, data = sys.argv[1], sys.argv[2]
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    cmd = [exe, '--grav', str(Path(data) / 'gravity.bin'),
           '--topo', str(Path(data) / 'topo.bin'),
           '--n', '1', '--seed', '777', '--days', '12', '--lmax', '60',
           '--alt-min', '50', '--alt-max', '50', '--steps-per-orbit', '300',
           '--integrator', 'yoshida4', '--traj-count', '0',
           '--out', str(td / 'summary.csv'), '--traj', str(td / 'traj.csv'),
           '--log', str(td / 'run.log')]
    env = dict(os.environ, OMP_NUM_THREADS='1')
    subprocess.run(cmd, env=env, check=True, capture_output=True)
    with (td/'summary.csv').open() as f:
        r = next(csv.DictReader(f))
    assert int(r['impacted']) == 1, 'orbit did not impact'
    assert abs(float(r['t_end_hr']) - 161.205253464) < 0.02, 'impact time changed'
    assert 3.0 < float(r['impact_angle_deg']) < 5.0, 'invalid terrain normal angle'
    assert 0.0 < float(r['radial_angle_deg']) < 1.0, 'invalid radial angle'
    assert float(r['invariant_drift_rel']) < 5e-7, 'unexpected Jacobi drift'
print('PASS end-to-end lunar impact/time/terrain-normal regression')
