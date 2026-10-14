#!/usr/bin/env python3
"""CLI regression for --lmax max and invalid requested degree."""
import subprocess, sys, tempfile
from pathlib import Path
exe,datadir=sys.argv[1:3]
with tempfile.TemporaryDirectory() as temp:
    d=Path(temp)
    common=[exe,'--grav',str(Path(datadir)/'gravity.bin'),
            '--topo',str(Path(datadir)/'topo.bin'),
            '--n','1','--days','0.002','--alt-min','25',
            '--alt-max','50','--steps-per-orbit','80',
            '--traj-count','0','--out',str(d/'summary.csv'),
            '--traj',str(d/'traj.csv'),'--log',str(d/'run.log')]
    run=subprocess.run(common+['--lmax','max'],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    assert 'gravity: lmax=300' in run.stdout,run.stdout
    bad=subprocess.run(common+['--lmax','600'],capture_output=True,text=True)
    assert bad.returncode!=0
    assert 'requested harmonic degree exceeds data file' in bad.stderr
print('PASS --lmax max resolves available 300 degrees; excessive 600 rejected')
