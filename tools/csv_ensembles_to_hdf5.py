#!/usr/bin/env python3
"""Archive existing September/October comparison CSV ensembles as compressed HDF5.

A one-off *conversion* of previously rounded CSV values; not a numerical
rerun, and these converted files are NOT valid restart checkpoints.
The output schema deliberately matches native writer field names.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SUMMARY_FIELDS = [
    ('id','i4','id'),('lat0_deg','f8','lat0'),('lon0_deg','f8','lon0'),
    ('alt0_km','f8','alt0_km'),('heading_deg','f8','heading'),
    ('inc0_deg','f8','inc0'),('impacted','i4','impacted'),
    ('t_end_hr','f8','t_end_hr'),('impact_lat_deg','f8','impact_lat'),
    ('impact_lon_deg','f8','impact_lon'),('impact_angle_deg','f8','impact_angle_deg'),
    ('radial_angle_deg','f8','radial_angle_deg'),
    ('terminal_speed_ms','f8','terminal_speed_ms'),
    ('terrain_elev_m','f8','terrain_elev_m'),
    ('vx_term_ms','f8','vx_term'),('vy_term_ms','f8','vy_term'),
    ('vz_term_ms','f8','vz_term'),('alt_min_km','f8','alt_min_km'),
    ('alt_max_km','f8','alt_max_km'),('ecc_max','f8','ecc_max'),
    ('orbit_count','i4','orbit_count'),('n_steps','i8','n_steps'),
    ('inv0_si','f8','inv0'),('invariant_drift_rel','f8','invariant_drift_rel'),
    ('lunar_jacobi_variation_rel','f8','lunar_jacobi_variation_rel'),
    ('kl_incl0_deg','f8','kl_incl0_deg'),('kl_incl_end_deg','f8','kl_incl_end_deg'),
    ('kl_jz0','f8','kl_jz0'),('kl_jz_end','f8','kl_jz_end')
]
TRAJ_FIELDS = [
    ('t_hr','f8','t_hr'),('lat_deg','f8','lat'),('lon_deg','f8','lon'),
    ('alt_km','f8','alt_km'),('jacobi_si','f8','invariant'),
    ('ecc','f8','ecc'),('x_m','f8','x'),('y_m','f8','y'),
    ('z_m','f8','z'),('vx_ms','f8','vx'),('vy_ms','f8','vy'),
    ('vz_ms','f8','vz'),('kl_incl_deg','f8','kl_incl_deg'),
    ('kl_jz','f8','kl_jz')
]

def to_structured(df, spec):
    dest = np.empty(len(df),dtype=[(name,typ) for name,typ,_ in spec])
    for name,typ,source in spec:
        if source in df.columns:
            dest[name] = df[source].to_numpy(dtype=dest.dtype[name])
        else:
            dest[name] = -1 if typ.startswith('i') else np.nan
    return dest

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for chunk in iter(lambda:stream.read(1<<20),b''):
            h.update(chunk)
    return h.hexdigest()

def add_run(root, name, summary_csv, trajectory_csv):
    ss=pd.read_csv(summary_csv)
    ts=pd.read_csv(trajectory_csv)
    ss=ss.sort_values('id',kind='stable').reset_index(drop=True)
    ts=ts.sort_values('id',kind='stable').reset_index(drop=True)
    n=len(ss)
    if not np.array_equal(ss['id'].to_numpy(),np.arange(n)):
        raise ValueError(f'{name}: id must enumerate 0..N-1')
    if len(ts) and ((ts['id']<0).any() or (ts['id']>=n).any()):
        raise ValueError(f'{name}: invalid trajectory particle id')
    counts=np.bincount(ts['id'],minlength=n)
    offsets=np.concatenate(([0],np.cumsum(counts,dtype=np.uint64)))
    if len(ts):
        for i in range(n):
            t=ts.iloc[int(offsets[i]):int(offsets[i+1])]['t_hr'].to_numpy()
            if np.any(np.diff(t)<0):
                raise ValueError(f'{name}: nonmonotonic times for particle {i}')
    grp=root.create_group(name)
    grp.attrs['conversion']='legacy_csv_to_hdf5_v1'
    grp.attrs['native_checkpoints_compatible']=False
    grp.attrs['conversion_warning']='Data reconstructed from rounded CSV; no additional numerical precision is available.'
    grp.attrs['summary_source_file']=str(summary_csv.name)
    grp.attrs['trajectory_source_file']=str(trajectory_csv.name)
    grp.attrs['summary_source_sha256']=sha256(summary_csv)
    grp.attrs['trajectory_source_sha256']=sha256(trajectory_csv)
    grp.attrs['n_particles']=n
    grp.create_dataset('summary',data=to_structured(ss,SUMMARY_FIELDS),compression='gzip',shuffle=True,compression_opts=4)
    trajgrp=grp.create_group('trajectories')
    trajgrp.create_dataset('samples',data=to_structured(ts,TRAJ_FIELDS),compression='gzip',shuffle=True,compression_opts=4)
    trajgrp.create_dataset('offsets',data=offsets,compression='gzip',shuffle=True,compression_opts=4)
    trajgrp.attrs['layout']='samples[offsets[i]:offsets[i+1]] are particle i; missing values in older versions are NaN'
    print(f'{name}: {n} particles, {len(ts)} samples, {int(ss.impacted.sum())} impacts')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('comparison_dir',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.parent.mkdir(exist_ok=True,parents=True)
    tmp=args.output.with_suffix('.h5.tmp')
    with h5py.File(tmp,'w') as output:
        output.attrs['schema']='lunar_orbits_converted_ensembles_v1'
        output.attrs['description']='Previously computed 100-particle lunar orbit comparisons, CSV fidelity (not new integrations)'
        output.attrs['units']='Time hours, distances km or m, velocity m/s; see field names'
        runs=output.create_group('runs')
        for label,sumfile,trajfile in [
            ('older_degree40','old_summary.csv','old_trajectories.csv'),
            ('corrected_degree40','new_summary.csv','new_trajectories.csv'),
            ('older_degree60','old_l60_summary.csv','old_l60_trajectories.csv'),
            ('corrected_degree60','new_l60_summary.csv','new_l60_trajectories.csv')]:
            add_run(runs,label,args.comparison_dir/sumfile,args.comparison_dir/trajfile)
    tmp.replace(args.output)
    print(f'Wrote {args.output}: {args.output.stat().st_size / 1024**2:.2f} MiB')
