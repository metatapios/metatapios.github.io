#!/usr/bin/env python3
"""Short reproducible degree-300 adaptive vs. fixed-step validation."""
import argparse, subprocess, time, csv, json, pathlib
import h5py, numpy as np
ROOT=pathlib.Path(__file__).resolve().parent.parent
B=ROOT/'build/lunarmc'
O=ROOT/'validation_adaptive'
O.mkdir(exist_ok=True)
rows=[]
for alt in (25, 150, 500):
 for typ in ('adaptive','yoshida4'):
  name=f'alt{alt}_{typ}'
  folder=O/name
  folder.mkdir(parents=True,exist_ok=True)
  args=[str(B),'--grav',str(ROOT/'data/gravity.bin'),'--topo',str(ROOT/'data/topo.bin'),
    '--lmax','max','--n','1','--seed','20261014','--alt-min',str(alt),'--alt-max',str(alt),
    '--days','0.5','--integrator',typ,'--steps-per-orbit','960',
    '--event-search','adaptive','--event-max-depth','32',
    '--kozai-lidov','off','--traj-count','1','--traj-dt-min','30',
    '--hdf5',str(folder/'out.h5'),'--no-csv','--log',str(folder/'log.txt'),
    '--checkpoint-dir',str(folder/'checkpoints'),'--checkpoint-every-hours','4']
  if typ=='adaptive':
    args+=['--rtol','1e-10','--atol-pos-m','0.01','--atol-vel-mps','1e-5',
        '--adaptive-max-dt-s','90','--adaptive-low-steps','480',
        '--adaptive-high-steps','120']
  start=time.perf_counter()
  r=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
  runtime=time.perf_counter()-start
  if r.returncode: raise RuntimeError(name+':\n'+r.stdout[-2500:])
  with h5py.File(folder/'checkpoints/particle_00000000.h5') as h:
    print(name,'datasets',list(h.keys()))
    print('attributes',list(h.attrs.keys()))
    x=np.asarray(h['state']).flatten() if 'state' in h else []
  # Summary CSV not written; read HDF5's output summary via its declared schema
  with h5py.File(folder/'out.h5') as h:
    print('RESULT keys',list(h.keys()))
  rows.append(dict(altitude_km=alt,integrator=typ,elapsed_s=runtime))
  print('RESULT',name,runtime,'sec',flush=True)
with open(O/'benchmark.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
