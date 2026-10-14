#!/usr/bin/env python3
"""Build a synthetic (J2-only, spherical surface) dataset without overwriting the
real GRGM1200B/LOLA lunar data. For integration smoke tests only.
Example:
  python tools/make_synth_data.py --out-dir data_synthetic
  ./build/lunarmc --grav data_synthetic/gravity.bin --topo data_synthetic/topo.bin --lmax 2
"""
import argparse
from pathlib import Path
import struct
import math

parser = argparse.ArgumentParser()
parser.add_argument('--out-dir', default='data_synthetic')
a = parser.parse_args()
out = Path(a.out_dir)
out.mkdir(parents=True, exist_ok=True)
GM=4.902794e12
R=1737400.0
lmax=2
Cbar20=-2.033e-4/math.sqrt(5.0)
with (out/'gravity.bin').open('wb') as f:
    f.write(struct.pack('<idd', lmax, GM, R))
    for n in range(lmax+1):
        for m in range(n+1):
            C=1.0 if (n,m)==(0,0) else Cbar20 if (n,m)==(2,0) else 0.0
            f.write(struct.pack('<dd',C,0.0))
ny,nx=180,360
with (out/'topo.bin').open('wb') as f:
    f.write(struct.pack('<ii',ny,nx))
    f.write(struct.pack('<d',R)*(ny*nx))
print('Synthetic test-only data written to',out)
