#!/usr/bin/env python3
"""End-to-end regression of append-only events, live statuses, and HDF5 output."""
import csv
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import h5py

exe = Path(sys.argv[1]).resolve()
data = Path(sys.argv[2]).resolve()

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

with tempfile.TemporaryDirectory(prefix='mascons_progress_test_') as temporary:
    root = Path(temporary)
    cp = root/'checkpoints'
    args = [str(exe),'--grav',str(data/'gravity.bin'),'--topo',str(data/'topo.bin'),
            '--n','6','--lmax','60','--seed','20261014','--alt-min','25',
            '--alt-max','500','--integrator','adaptive','--traj-count','2',
            '--checkpoint-dir',str(cp),'--checkpoint-every-hours','0.5',
            '--hdf5',str(root/'out.h5'),'--no-csv', '--log',str(root/'run.log'),
            '--progress-dir',str(root),'--progress-every-seconds','0.02',
            '--progress-sim-hours','0.02']
    env = dict(os.environ,OMP_NUM_THREADS='2')
    subprocess.run(args+['--days','0.05'],env=env,check=True,stdout=subprocess.DEVNULL)
    statuses = read(root/'particle_status.csv')
    events = read(root/'particle_events.csv')
    timeline = read(root/'progress.csv')
    assert len(statuses)==6 and len(events)==6
    assert all(x['status']=='survived' for x in statuses)
    assert all(float(x['sim_day'])==0.05 for x in statuses)
    assert all(float(x['energy_change_rel_signed']) == float(x['energy_change_rel_signed']) for x in statuses)
    assert all(float(x['numeric_invariant_drift_rel']) < 1.e-5 for x in statuses)
    assert all(float(x['last_accepted_local_error_ratio']) <= 1 for x in statuses)
    assert int(timeline[-1]['completed'])==6
    with h5py.File(root/'out.h5') as h:
        assert len(h['summary']) == 6
        assert len(h['trajectories/offsets'])==7
    # Resuming the same completed horizon should not repeat completed events.
    subprocess.run(args+['--days','0.05','--resume'],env=env,check=True,stdout=subprocess.DEVNULL)
    assert len(read(root/'particle_events.csv'))==6
    # Extending the horizon should record the new completion epoch once each.
    subprocess.run(args+['--days','0.08','--resume'],env=env,check=True,stdout=subprocess.DEVNULL)
    assert len(read(root/'particle_events.csv'))==12
    assert len(read(root/'progress.csv'))>=4
    lines=(root/'run.log').read_text().strip().splitlines()
    pattern=r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \| \d+/6 completed \| \d+\.\d{2}% \| ETA: (unknown|\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})$'
    assert lines and all(re.match(pattern, line) for line in lines),lines
    assert '6/6 completed' in lines[-1],lines[-1]
    assert '100.00%' in lines[-1],lines[-1]
    assert len(read(root/'particle_status.csv'))==6
    with h5py.File(root/'out.h5') as h:
        assert float(h['metadata'].attrs['days'])==0.08
    # With third bodies present, physical lunar Jacobi variations are NOT
    # falsely categorized as numerical conservation drift.
    forced = root/'forced'
    forced.mkdir()
    forced_args=[str(exe),'--grav',str(data/'gravity.bin'),'--topo',str(data/'topo.bin'),
                 '--n','2','--days','0.02','--lmax','60','--kozai-lidov','on',
                 '--ephem-mode','analytic','--progress-dir',str(forced),
                 '--log',str(forced/'run.log'),'--progress-every-seconds','0.02',
                 '--progress-sim-hours','0.02','--hdf5',str(forced/'out.h5'),'--no-csv']
    subprocess.run(forced_args,env=env,check=True,stdout=subprocess.DEVNULL)
    for row in read(forced/'particle_status.csv'):
        assert row['numeric_invariant_drift_rel']=='nan', row
        assert row['max_lunar_jacobi_variation_rel']!='nan'
    # Exact surface starts must appear as terminal lunar impacts in the
    # per-particle event ledger, not as unreported numerical failures.
    impacts = root/'impacts'
    impacts.mkdir()
    collision_args=[str(exe),'--grav',str(data/'gravity.bin'),'--topo',str(data/'topo.bin'),
        '--n','2','--days','0.02','--alt-min','0','--alt-max','0',
        '--lmax','60','--progress-dir',str(impacts),'--log',str(impacts/'run.log'),
        '--hdf5',str(impacts/'out.h5'),'--no-csv']
    subprocess.run(collision_args,env=env,check=True,stdout=subprocess.DEVNULL)
    assert [r['status'] for r in read(impacts/'particle_status.csv')]==['impacted','impacted']
    impact_events=read(impacts/'particle_events.csv')
    assert len(impact_events)==2
    assert all(r['outcome']=='impacted' for r in impact_events)
    assert all(float(r['sim_day'])==0 for r in impact_events)
print('live progress, restart, HDF5 archive and physical energy labels: PASS')
