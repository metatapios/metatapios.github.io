#!/usr/bin/env python3
"""Validate the native, compressed HDF5 archive and checkpoint extension semantics."""
import csv
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np

exe, data = sys.argv[1:3]
env = dict(os.environ, OMP_NUM_THREADS='1')

def invoke(folder, label, days, cp=None, resume=False, trajectories=3, csv_on=True):
    args=[exe,'--grav',str(Path(data)/'gravity.bin'),
          '--topo',str(Path(data)/'topo.bin'),'--n','3','--seed','20261014',
          '--lmax','25','--alt-min','25','--alt-max','50',
          '--steps-per-orbit','200','--days',str(days),
          '--traj-count',str(trajectories),'--traj-dt-min','6',
          '--kozai-lidov','off', '--hdf5',str(folder/(label+'.h5')),
          '--out',str(folder/(label+'_sum.csv')),
          '--traj',str(folder/(label+'_traj.csv')),
          '--log',str(folder/(label+'.log'))]
    if not csv_on: args.append('--no-csv')
    if cp:args+=['--checkpoint-dir',str(cp),'--checkpoint-every-hours','0.2']
    if resume:args.append('--resume')
    subprocess.run(args,env=env,check=True,capture_output=True,text=True)

with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp); cp=root/'checkpoints'
    invoke(root,'stage1',.10,cp)
    invoke(root,'continued',.20,cp,True,csv_on=False)
    invoke(root,'direct',.20)
    assert not (root/'continued_sum.csv').exists(), '--no-csv produced CSV'
    for label in ['stage1','continued','direct']:
        with h5py.File(root/(label+'.h5'),'r') as f:
            assert f.attrs['schema'].decode()=='lunar_orbits_hdf5_v1'
            assert f['metadata'].attrs['kozai_lidov']==b'off'
            assert f['metadata'].attrs['seed']==20261014
            assert f['summary'].compression=='gzip'
            assert f['trajectories/samples'].compression=='gzip'
            assert f['summary'].shape==(3,)
            ids=f['summary']['id']; assert np.array_equal(ids,np.arange(3))
            offsets=f['trajectories/offsets'][:]
            assert len(offsets)==4 and offsets[0]==0
            assert offsets[-1]==len(f['trajectories/samples'])
            assert np.all(np.diff(offsets)>0)
            for idx in range(3):
                samples=f['trajectories/samples'][int(offsets[idx]):int(offsets[idx+1])]
                assert samples[0]['t_hr']==0
                assert np.all(np.diff(samples['t_hr'])>0)
                assert samples[-1]['t_hr']<=.20*24+1e-6
            assert not (root/(label+'.h5.tmp')).exists()
    with h5py.File(root/'stage1.h5') as a, h5py.File(root/'continued.h5') as b, h5py.File(root/'direct.h5') as d:
        assert a['metadata'].attrs['run_signature']==b['metadata'].attrs['run_signature']
        assert b['metadata'].attrs['restart_mode']==b'resumed'
        assert np.allclose(b['summary']['inv0_si'],d['summary']['inv0_si'],rtol=1e-14)
        assert np.allclose(b['summary']['t_end_hr'],4.8)
        off_a=a['trajectories/offsets'][:];off_b=b['trajectories/offsets'][:]
        for i in range(3):
            x=a['trajectories/samples'][int(off_a[i]):int(off_a[i+1])]
            y=b['trajectories/samples'][int(off_b[i]):int(off_b[i+1])]
            assert len(y)>len(x)
            # Exact checkpoint history restoration of the prior stage.
            for field in x.dtype.names:
                np.testing.assert_array_equal(x[field],y[field][:len(x)])
    # Missing sampling (`traj_count=0`) is a valid empty ragged dataset.
    invoke(root,'zero_samples',.02,trajectories=0,csv_on=False)
    with h5py.File(root/'zero_samples.h5','r') as f:
        assert len(f['trajectories/samples'])==0
        assert np.array_equal(f['trajectories/offsets'][:],[0,0,0,0])
    with open(root/'direct_sum.csv') as fh:
        direct_csv=list(csv.DictReader(fh))
    with h5py.File(root/'direct.h5') as f:
        for idx,row in enumerate(direct_csv):
            rr=f['summary'][idx]
            assert int(row['impacted'])==rr['impacted']
            assert int(row['n_steps'])==rr['n_steps']
            assert math.isclose(float(row['ecc_max']),rr['ecc_max'],abs_tol=1e-8)
            assert math.isclose(float(row['invariant_drift_rel']),rr['invariant_drift_rel'],rel_tol=2e-6,abs_tol=1e-14)
            assert math.isclose(float(row['t_end_hr']),rr['t_end_hr'],abs_tol=1e-8)
print('PASS HDF5 compression, typed fields, schema, restart history, no-csv, zero samples, CSV agreement')
