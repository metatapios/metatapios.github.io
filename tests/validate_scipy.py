#!/usr/bin/env python3
"""Independent potential/acceleration test against scipy.special.lpmv.
Usage: python tests/validate_scipy.py build/gravity_probe data/gravity.bin
No pyshtools required. SciPy uses Condon-Shortley sign, removed here.
"""
import math
from pathlib import Path
import struct
import subprocess
import sys

import numpy as np
from scipy.special import lpmv, gammaln

probe, binary = sys.argv[1:3]
L = 40
with open(binary, 'rb') as f:
    nfile, = struct.unpack('<i', f.read(4))
    gm, radius = struct.unpack('<dd', f.read(16))
    if nfile < L:
        raise ValueError('gravity dataset has fewer than 40 degrees')
    cs = np.zeros((L + 1, L + 1))
    ss = np.zeros_like(cs)
    for n in range(nfile + 1):
        for m in range(n + 1):
            c, s = struct.unpack('<dd', f.read(16))
            if n <= L:
                cs[n, m], ss[n, m] = c, s


def ref_potential(q):
    x, y, z = q
    r = np.linalg.norm(q)
    sintheta = z / r
    lon = math.atan2(y, x)
    total = 0.0
    for n in range(L + 1):
        inner = 0.0
        for m in range(n + 1):
            norm = math.sqrt((2 if m else 1) * (2 * n + 1) *
                             math.exp(float(gammaln(n - m + 1) - gammaln(n + m + 1))))
            plm = (-1)**m * norm * float(lpmv(m, n, sintheta))
            inner += plm * (cs[n, m] * math.cos(m * lon) + ss[n, m] * math.sin(m * lon))
        total += (radius / r)**n * inner
    return gm / r * total


def cpp(q):
    cmd = [probe, binary, str(L), *map(str, q)]
    return np.array([float(v) for v in subprocess.check_output(cmd, text=True).split()])


points = [
    np.array([radius + 50000., 0, 0]),
    np.array([50000., radius + 50000., 190000.]),
    np.array([radius*.75, radius*.43, radius*.65]),
    np.array([-radius*.4, radius*.7, -radius*.65]),
    np.array([-200000., 350000., radius + 30000.]),
]
max_pot = max_acc = 0.0
for q in points:
    pot_cpp, *a_cpp = cpp(q)
    pot_sci = ref_potential(q)
    fd = np.array([(ref_potential(q + np.eye(3)[i]*10) -
                    ref_potential(q - np.eye(3)[i]*10))/20.0 for i in range(3)])
    max_pot = max(max_pot, abs(pot_cpp - pot_sci)/abs(pot_sci))
    max_acc = max(max_acc, float(np.linalg.norm(np.asarray(a_cpp)-fd)/np.linalg.norm(fd)))
print(f'Independent SciPy spherical-harmonic potential relative error: {max_pot:.3e}')
print(f'Independent SciPy Cartesian acceleration relative error: {max_acc:.3e}')
if max_pot > 1e-11 or max_acc > 3e-7:
    raise SystemExit('FAILED independent harmonic reference comparison')
print('PASS independent SciPy degree-40 gravity comparison')
