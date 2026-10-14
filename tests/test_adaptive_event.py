#!/usr/bin/env python3
"""Compare new clearance-aware event search with dense fixed probes on a known impact."""
import csv,os,subprocess,sys,tempfile
from pathlib import Path
exe,data=sys.argv[1:3]
with tempfile.TemporaryDirectory() as tmp:
    rows=[]
    for mode in ('sampled','adaptive'):
        d=Path(tmp);outfile=d/(mode+'.csv')
        cmd=[exe,'--grav',str(Path(data)/'gravity.bin'),'--topo',str(Path(data)/'topo.bin'),
             '--n','1','--seed','777','--days','12','--lmax','60','--alt-min','50','--alt-max','50',
             '--steps-per-orbit','300','--integrator','yoshida4','--traj-count','0',
             '--event-search',mode,'--event-substeps','16',
             '--out',str(outfile),'--traj',str(d/(mode+'_tr.csv')),'--log',str(d/(mode+'.log'))]
        proc=subprocess.run(cmd,env=dict(os.environ,OMP_NUM_THREADS='1'),capture_output=True,text=True)
        if proc.returncode!=0:raise RuntimeError(proc.stderr+'\n'+proc.stdout)
        with outfile.open() as f:rows.append(next(csv.DictReader(f)))
    assert all(int(r['impacted'])==1 for r in rows),rows
    assert abs(float(rows[0]['t_end_hr'])-float(rows[1]['t_end_hr']))*3600<0.01, rows
print('adaptive clearance search matches 16-subdivision impact within 0.01s')
