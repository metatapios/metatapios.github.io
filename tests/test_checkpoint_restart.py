#!/usr/bin/env python3
"""Executable restart/extension, validation, and interrupted-process regressions."""
import csv
import h5py
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

exe, data = sys.argv[1:3]
env = dict(os.environ, OMP_NUM_THREADS='1')

def rows(file):
    with open(file) as fh:
        return list(csv.DictReader(fh))

def run(root, label, days, cp=None, resume=False, check=True, extra=()):
    command = [exe, '--grav', str(Path(data)/'gravity.bin'),
               '--topo', str(Path(data)/'topo.bin'),
               '--n', '2', '--seed', '20261014', '--lmax', '25',
               '--alt-min', '25', '--alt-max', '50', '--steps-per-orbit', '180',
               '--traj-count', '2', '--traj-dt-min', '10',
               '--kozai-lidov', 'off', '--days', str(days),
               '--out', str(root/(label+'.csv')),
               '--traj', str(root/(label+'_traj.csv')),
               '--log', str(root/(label+'.log'))]
    if cp:
        command += ['--checkpoint-dir', str(cp), '--checkpoint-every-hours', '0.15']
    if resume:
        command += ['--resume']
    command += list(extra)
    return subprocess.run(command, env=env, capture_output=True, text=True, check=check)

with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp); cp=root/'snapshots'
    run(root, 'stage1', .13, cp)
    short=rows(root/'stage1.csv')
    assert sorted(p.name for p in cp.glob('particle_*.h5') if '.prev.h5' not in p.name) == [
        'particle_00000000.h5', 'particle_00000001.h5']
    before=(cp/'particle_00000000.h5').read_bytes()
    run(root, 'same', .13, cp, resume=True)
    assert (root/'same.csv').read_bytes()==(root/'stage1.csv').read_bytes()
    assert (root/'same_traj.csv').read_bytes()==(root/'stage1_traj.csv').read_bytes()
    assert (cp/'particle_00000000.h5').read_bytes()==before

    # Resuming a finished survivor must extend, not regenerate from t=0.
    run(root, 'extended', .27, cp, resume=True)
    extended=rows(root/'extended.csv')
    run(root, 'reference', .27)
    reference=rows(root/'reference.csv')
    for old, new, ref in zip(short, extended, reference):
        assert old['id']==new['id']==ref['id']
        assert float(new['t_end_hr'])==.27*24.
        assert float(new['inv0'])==float(old['inv0'])
        assert int(new['n_steps'])>int(old['n_steps'])
        assert abs(float(new['ecc_max'])-float(ref['ecc_max']))<5e-5
        assert abs(float(new['alt_min_km'])-float(ref['alt_min_km']))<.3
    for ident in ('0','1'):
        a=[x for x in rows(root/'extended_traj.csv') if x['id']==ident]
        t=[float(x['t_hr']) for x in a]
        assert t[0]==0 and all(y>x for x,y in zip(t,t[1:])),t
        assert max(t)>3.12  # original stop 0.13 days

    # All physics-input and numerical settings are part of the run identity.
    for changed in [('--steps-per-orbit','181'),('--seed','999'),('--omega','0')]:
        bad=run(root, 'wrong_settings', .30, cp, resume=True, check=False, extra=changed)
        assert bad.returncode!=0 and 'manifest' in bad.stderr,bad.stderr
    bad=run(root,'shortened',.1,cp,resume=True,check=False)
    assert bad.returncode!=0 and 'shorten' in bad.stderr,bad.stderr
    # Detect damaged checkpoint snapshots without silently restarting the orbit.
    snapshot=cp/'particle_00000000.h5';original=snapshot.read_bytes()
    previous_file=cp/'particle_00000000.prev.h5'
    older_original=previous_file.read_bytes()
    snapshot.write_bytes(b"corrupt-hdf5")
    previous_file.write_bytes(b"corrupt-hdf5")
    bad=run(root,'corrupt',.30,cp,resume=True,check=False)
    assert bad.returncode!=0 and 'both checkpoint generations' in bad.stderr,bad.stderr
    snapshot.write_bytes(original)
    previous_file.write_bytes(older_original)
    # The latest snapshot can fail without losing the previous generation.
    with h5py.File(cp/'particle_00000000.prev.h5','r') as f:
        assert f.attrs['schema']==b'lunar_particle_restart_hdf5_v2'
    # A failed latest must automatically fall back to the older generation.
    snapshot.write_bytes(b'corrupt-hdf5')
    recovered=run(root,'recover',.30,cp,resume=True)
    assert recovered.returncode==0
    with h5py.File(snapshot,'r') as f:
        assert f.attrs['schema']==b'lunar_particle_restart_hdf5_v2'
    # Keep the earlier uninterrupted checkpoint for the adaptive test.

    # Adaptive DOPRI restart state: next proposed step and diagnostics.
    cp2=root/'adaptive'
    run(root,'adapt1',.08,cp2,extra=('--integrator','adaptive'))
    run(root,'adapt2',.16,cp2,resume=True,extra=('--integrator','adaptive'))
    assert all(float(x['t_end_hr'])==.16*24 for x in rows(root/'adapt2.csv'))

    # NaN (KL fields) must survive serialization with third-body forces enabled.
    cp3=root/'tides'
    run(root,'tides1',.03,cp3,extra=('--kozai-lidov','on','--ephem-mode','analytic'))
    run(root,'tides2',.06,cp3,resume=True,
        extra=('--kozai-lidov','on','--ephem-mode','analytic'))
    assert all(float(x['t_end_hr'])==.06*24 for x in rows(root/'tides2.csv'))

    # Persist a completed impact event, without resurrecting it on extension.
    cp4=root/'impact'
    special=('--n','1','--lmax','60','--seed','777','--alt-min','50',
             '--alt-max','50','--steps-per-orbit','300','--traj-count','0')
    run(root,'hit1',7,cp4,extra=special)
    hit1=rows(root/'hit1.csv')[0]
    assert hit1['impacted']=='1' and float(hit1['t_end_hr'])<7*24
    snap=(cp4/'particle_00000000.h5').read_bytes()
    run(root,'hit2',12,cp4,resume=True,extra=special)
    assert rows(root/'hit2.csv')[0]==hit1
    assert (cp4/'particle_00000000.h5').read_bytes()==snap

    # SIGTERM an active job after a nonzero-time snapshot, then resume it.
    cp5=root/'interrupted'
    cmd=[exe,'--grav',str(Path(data)/'gravity.bin'),'--topo',str(Path(data)/'topo.bin'),
         '--n','1','--lmax','40','--seed','12','--days','1.5',
         '--alt-min','45','--alt-max','45','--steps-per-orbit','1500',
         '--traj-count','1','--checkpoint-dir',str(cp5),
         '--checkpoint-every-hours','0.03','--out',str(root/'aborted.csv'),
         '--traj',str(root/'aborted_traj.csv'),'--log',str(root/'aborted.log')]
    proc=subprocess.Popen(cmd,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    progress=False
    deadline=time.monotonic()+20
    while time.monotonic()<deadline:
        chk=cp5/'particle_00000000.h5'
        if chk.exists():
            try:
                # The checkpoint is atomically replaced by a complete file.
                with h5py.File(chk,'r') as f:
                    t=float(f['state_t_q_v_a'][0]);progress=t>3600*.04
            except (OSError,ValueError,IndexError):
                pass
            if progress:break
        if proc.poll() is not None:break
        time.sleep(.006)
    proc.terminate()
    try:proc.wait(timeout=5)
    except subprocess.TimeoutExpired:proc.kill();proc.wait()
    assert progress,'no noninitial checkpoint was captured before terminating'
    command=cmd.copy();command[command.index('--days')+1]='.23'
    # In-flight stored time is earlier than 0.23 days.
    command.append('--resume')
    end=subprocess.run(command,env=env,capture_output=True,text=True)
    assert end.returncode==0,end.stderr
    assert abs(float(rows(root/'aborted.csv')[0]['t_end_hr'])-.23*24)<1e-8
print('PASS checkpoint restart, extension, impact persistence, corruption, adaptive, KL, and SIGTERM')
