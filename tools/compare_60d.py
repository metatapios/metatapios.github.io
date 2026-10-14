#!/usr/bin/env python3
"""Reproducible old-CSV versus corrected-HDF5 lunar orbit comparison.

Requires numpy, pandas, h5py, matplotlib. Original September outputs are read
unchanged; the corrected archive stores every orbit in compressed HDF5.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GM = 4.902800122e12


def load(old_s, old_t, new_h):
    osum = pd.read_csv(old_s)
    od = pd.read_csv(old_t)
    old = {}
    for i, sub in od.groupby('id', sort=False):
        old[int(i)] = sub.copy()
    with h5py.File(new_h) as h:
        ns = pd.DataFrame.from_records(h['summary'][:])
        a = h['trajectories/samples'][:]
        offsets = h['trajectories/offsets'][:]
        new = {i:pd.DataFrame.from_records(a[offsets[i]:offsets[i+1]])
               for i in range(len(offsets)-1)}
    return osum, old, ns, new


def paired_outcomes(o, n, dest):
    s = o[['id','impacted','t_end_hr','impact_lat','impact_lon']].rename(
        columns={'impacted':'old_impacted','t_end_hr':'old_t_end_hr',
                 'impact_lat':'old_impact_lat_deg','impact_lon':'old_impact_lon_deg'})
    t = n[['id','impacted','t_end_hr','impact_lat_deg','impact_lon_deg']].rename(
        columns={'impacted':'new_impacted','t_end_hr':'new_t_end_hr',
                 'impact_lat_deg':'new_impact_lat_deg','impact_lon_deg':'new_impact_lon_deg'})
    match = s.merge(t, on='id', validate='one_to_one')
    match.to_csv(dest, index=False)
    return match


def finish(fig, out):
    fig.tight_layout()
    fig.savefig(out, dpi=165, bbox_inches='tight')
    plt.close(fig)


def plot_survival(o, n, folder, horizon=60):
    fig, ax = plt.subplots(figsize=(9.6,5.4))
    for tag, s in [('Original September',o), ('Corrected October',n)]:
        events = np.sort(np.asarray(s.loc[s.impacted==1,'t_end_hr'])/24)
        times = np.r_[0, events, horizon]
        remaining = np.r_[100,100-np.arange(1,len(events)+1),100-len(events)]
        ax.step(times, remaining, where='post',label=f'{tag}: {len(events)} impacts')
    ax.set(xlim=(0,horizon),ylim=(0,102),xlabel='Elapsed time [days]',ylabel='Satellites still orbiting',
           title='100 identical starting orbits | 25–50 km | degree-60 lunar gravity')
    ax.grid(alpha=.3); ax.legend()
    finish(fig,folder/'01_survival.png')


def plot_sampled(old,new,folder,variable,filename,ylab):
    ids=(27,56,73,22)
    fig,ax=plt.subplots(figsize=(9.6,5.4))
    for i in ids:
        so=old[i]; sn=new[i]
        x=so['t_hr'].to_numpy()/24
        ax.plot(x,so[variable].to_numpy(),alpha=.75,lw=1.3,label=f'#{i} original')
        ax.plot(sn.t_hr.to_numpy()/24,sn[variable].to_numpy(),ls='--',alpha=.8,lw=1.5,
                label=f'#{i} corrected')
    ax.set(xlim=(0,60),xlabel='Elapsed time [days]',ylabel=ylab,title='Matched individual satellite orbits (solid: original; dashed: corrected)')
    ax.grid(alpha=.28);ax.legend(ncol=2,fontsize=8)
    finish(fig,folder/filename)


def jacobi_drift_ensemble(old,new,folder):
    fig,ax=plt.subplots(figsize=(9.6,5.4))
    days=np.linspace(0,60,241)
    for label,traj,col in [('Original September',old,'invariant'),('Corrected October',new,'jacobi_si')]:
        # Running *global maximum* of the relative deviation from each particle's initial J.
        levels=np.zeros_like(days)
        for df in traj.values():
            if len(df)==0:continue
            t=np.asarray(df.t_hr)/24
            J=np.asarray(df[col]); denom=max(abs(J[0]),1.0)
            err=np.maximum.accumulate(abs(J-J[0])/denom)
            # Once a particle impacts its previous maximum remains part of the
            # ensemble maximum; do not extrapolate its ongoing trajectory.
            val=np.interp(days,t,err,left=0,right=err[-1])
            levels=np.maximum(levels,val)
        ax.plot(days,levels,label=label)
    ax.set(xlim=(0,60),yscale='log',xlabel='Elapsed time [days]',
           ylabel='Maximum accrued |J(t)−J(0)|/|J(0)|',
           title='Jacobi conservation: historical maximum across all 100 orbits')
    ax.grid(alpha=.3,which='both');ax.legend()
    finish(fig,folder/'04_jacobi.png')


def node_and_apse(df):
    q=np.stack([df[k].to_numpy() for k in ('x','y','z')],axis=1) if 'x' in df else \
      np.stack([df[k].to_numpy() for k in ('x_m','y_m','z_m')],axis=1)
    v=np.stack([df[k].to_numpy() for k in ('vx','vy','vz')],axis=1) if 'vx' in df else \
      np.stack([df[k].to_numpy() for k in ('vx_ms','vy_ms','vz_ms')],axis=1)
    h=np.cross(q,v)
    node=np.stack([-h[:,1],h[:,0]],axis=1)
    omega=np.degrees(np.unwrap(np.arctan2(node[:,1],node[:,0])))
    r=np.linalg.norm(q,axis=1)
    ev=np.cross(v,h)/GM-q/r[:,None]
    ep=np.linalg.norm(ev[:,:2],axis=1)
    peri=np.degrees(np.unwrap(np.arctan2(ev[:,1],ev[:,0])))
    peri[ep<0.0015]=np.nan  # projected ellipse is ill-conditioned near a circle
    return omega,peri


def precession(old,new,folder):
    ids=(27,56,73)
    for kind,col,filename,title in [('node',0,'05_node_precession.png','Osculating ascending-node precession'),
                                     ('apsidal',1,'06_apsidal_precession.png','Projected eccentricity-vector orientation')]:
        fig,ax=plt.subplots(figsize=(9.6,5.4))
        for i in ids:
            for tag, data, ls in [('original',old,'-'),('corrected',new,'--')]:
                d=data[i]
                z=node_and_apse(d)[col]
                if np.count_nonzero(np.isfinite(z))<5:continue
                first=z[np.isfinite(z)][0]
                ax.plot(np.asarray(d.t_hr)/24,z-first,ls=ls,lw=1.25,label=f'#{i} {tag}')
        ax.set(xlim=(0,60),xlabel='Elapsed time [days]',ylabel='Unwrapped orientation change [degrees]',title=title)
        ax.grid(alpha=.3);ax.legend(ncol=2,fontsize=8)
        if col==1:
            fig.text(.12,.01,'Ill-conditioned close to circular or projected-circular orbits; masked when projected eccentricity <0.0015.',fontsize=8)
        finish(fig,folder/filename)


def impact_map(o,n,folder):
    fig=plt.figure(figsize=(10.4,5.9))
    ax=fig.add_subplot(111,projection='mollweide')
    for label,s,lat,lon,mark in [('Original September',o,'impact_lat','impact_lon','o'),
                                 ('Corrected October',n,'impact_lat_deg','impact_lon_deg','x')]:
        d=s.loc[s.impacted==1]
        L=((np.asarray(d[lon],dtype=float)+180)%360)-180
        B=np.asarray(d[lat],dtype=float)
        ax.scatter(np.deg2rad(L),np.deg2rad(B),s=16,marker=mark,
                   alpha=.75,label=f'{label} ({len(d)})')
    ax.grid(alpha=.32);ax.legend(loc='lower left',bbox_to_anchor=(.0,-.18),ncol=2,fontsize=9)
    ax.set_title('Lunar impact locations | east-positive body-fixed longitude')
    finish(fig,folder/'07_impact_locations.png')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--old-summary',required=True)
    p.add_argument('--old-trajectories',required=True)
    p.add_argument('--new-hdf5',required=True)
    p.add_argument('--output-dir',required=True)
    c=p.parse_args(); folder=Path(c.output_dir); folder.mkdir(parents=True,exist_ok=True)
    o,old,n,new=load(c.old_summary,c.old_trajectories,c.new_hdf5)
    match=paired_outcomes(o,n,folder/'paired_particle_outcomes.csv')
    plot_survival(o,n,folder)
    plot_sampled(old,new,folder,'ecc','02_eccentricity.png','Osculating eccentricity (monopole convention)')
    plot_sampled(old,new,folder,'alt_km','03_altitude.png','Altitude above local terrain [km]')
    jacobi_drift_ensemble(old,new,folder)
    precession(old,new,folder)
    impact_map(o,n,folder)
    print('Old impacts:',int(o.impacted.sum()),'new impacts:',int(n.impacted.sum()))
    print('Old/new paired impacts:',sum((match.old_impacted==1)&(match.new_impacted==1)))
    print('Median impact day old/new:',o.loc[o.impacted==1,'t_end_hr'].median()/24,
          n.loc[n.impacted==1,'t_end_hr'].median()/24)
    print('Figures:',*sorted(x.name for x in folder.glob('*.png')))

if __name__=='__main__':main()
