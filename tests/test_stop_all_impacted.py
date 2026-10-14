"""Integration regression for optional explicit all-impacted stop logging.

Particle propagation already stops individually on surface impact, regardless
of this flag; the option records the ensemble termination criterion explicitly.
"""
import csv
import pathlib
import subprocess
import sys
import tempfile

binary, data = sys.argv[1:3]
with tempfile.TemporaryDirectory() as td:
    d = pathlib.Path(td)
    summaries = []
    for name, extra in [('opt_in', ['--stop-when-all-impacted']), ('default', [])]:
        summary = d / (name + '.csv')
        log = d / (name + '.log')
        cmd = [binary, '--grav', str(pathlib.Path(data) / 'gravity.bin'),
               '--topo', str(pathlib.Path(data) / 'topo.bin'),
               '--n', '3', '--lmax', '20', '--days', '0.08',
               '--steps-per-orbit', '100', '--alt-min', '0', '--alt-max', '0',
               '--seed', '99', '--traj-count', '0', '--kozai-lidov', 'off',
               '--out', str(summary), '--traj', str(d / (name + '_t.csv')),
               '--log', str(log), *extra]
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        with summary.open() as f:
            summaries.append(list(csv.DictReader(f)))
        logtxt = proc.stdout  # Verbose details are intentionally console-only.
        if name == 'opt_in':
            assert 'ensemble stop: all 3 particles impacted' in logtxt, logtxt
        else:
            assert 'ensemble stop:' not in logtxt, logtxt
    assert len(summaries[0]) == len(summaries[1]) == 3
    assert all(int(row['impacted']) == 1 for row in summaries[0])
    assert summaries[0] == summaries[1], 'stop policy should not change orbits'
print('optional all-impacted logging regression PASS')
