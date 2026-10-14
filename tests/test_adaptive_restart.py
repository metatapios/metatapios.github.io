#!/usr/bin/env python3
"""Regression for variable-step DOPRI HDF5 restart over the 25--500 km domain."""
import pathlib,sys,tempfile,subprocess
import numpy as np
import h5py
exe=pathlib.Path(sys.argv[1]); data=pathlib.Path(sys.argv[2])
with tempfile.TemporaryDirectory(prefix='mascons_dopri_restart_') as d:
 root=pathlib.Path(d)
 base=[str(exe),'--grav',str(data/'gravity.bin'),'--topo',str(data/'topo.bin'),
       '--n','4','--lmax','40','--seed','20261014','--alt-min','25',
       '--alt-max','500','--integrator','adaptive','--rtol','1e-10',
       '--atol-pos-m','0.01','--atol-vel-mps','1e-5',
       '--adaptive-max-dt-s','90','--event-search','adaptive',
       '--event-max-depth','32','--traj-count','0','--kozai-lidov','off',
       '--checkpoint-every-hours','1','--no-csv']
 def run(tag,days,resume=False):
  out=root/tag;out.mkdir(exist_ok=True)
  cmd=base+['--days',str(days),'--checkpoint-dir',str(out/'snap'),
            '--hdf5',str(out/'results.h5'),'--log',str(out/'log.txt')]
  if resume:cmd+=['--resume']
  result=subprocess.run(cmd,capture_output=True,text=True,timeout=60)
  if result.returncode: raise AssertionError(result.stdout[-1200:]+result.stderr[-1200:])
  states=[]
  for i in range(4):
   p=out/'snap'/f'particle_{i:08d}.h5'
   assert p.exists(),p
   with h5py.File(p) as f:
    assert int(f.attrs['finished'])==1
    states.append(np.asarray(f['state_t_q_v_a']))
  with h5py.File(out/'results.h5') as f:
   assert len(f['summary'])==4
  return np.asarray(states)
 partial=run('split',0.075)
 assert all(abs(s[0]-6480)<1e-6 for s in partial)
 completed=run('split',0.15,True)
 direct=run('direct',0.15)
 assert np.max(np.abs(completed[:,0]-direct[:,0]))<1e-7
 delta=np.linalg.norm(completed[:,1:4]-direct[:,1:4],axis=1)
 assert max(delta)<5.,delta
 print('PASS adaptive HDF5 restart 25--500 km; max split/full separation [m]:',max(delta))
