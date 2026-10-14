#!/usr/bin/env python3
"""Export GRGM1200B gravity coefficients and lunar topography to binary blobs
that the C++ propagator reads. Run once (needs pyshtools + network the first
time to fetch the datasets).

gravity.bin : int32 lmax, double GM, double R,
              then n=0..lmax, m=0..n : double C, double S   (4pi-normalized)
topo.bin    : int32 nlat, int32 nlon, then nlat*nlon doubles = absolute radius [m]
              pyshtools DH grid: row j -> lat = 90 - j*180/nlat (N to S),
              col i -> lon = i*360/nlon.  Matches src/topography.hpp.
"""
import argparse
import struct
import numpy as np
import pyshtools as pysh


def export_gravity(path, lmax):
    clm = pysh.datasets.Moon.GRGM1200B(lmax=lmax)
    C = clm.coeffs[0]  # [n, m]
    S = clm.coeffs[1]
    GM = float(clm.gm)
    R = float(clm.r0)
    with open(path, "wb") as f:
        f.write(struct.pack("<i", lmax))
        f.write(struct.pack("<dd", GM, R))
        for n in range(lmax + 1):
            for m in range(n + 1):
                f.write(struct.pack("<dd", float(C[n, m]), float(S[n, m])))
    print(f"gravity.bin: lmax={lmax} GM={GM:.6e} R={R:.1f} -> {path}")


def export_topo(path, ltopo):
    # LOLA_shape = absolute lunar radius; matches the validation scripts.
    topo = pysh.datasets.Moon.LOLA_shape(lmax=ltopo)
    coeffs = pysh.SHCoeffs.from_array(topo.coeffs)
    grid = coeffs.expand(grid="DH")         # absolute radius [m], DH sampling
    data = np.asarray(grid.data, dtype="<f8")   # shape (nlat, nlon), row0 = +90 lat
    nlat, nlon = data.shape
    with open(path, "wb") as f:
        f.write(struct.pack("<ii", nlat, nlon))
        data.tofile(f)
    print(f"topo.bin: {nlat}x{nlon} R in [{data.min()/1e3:.2f},"
          f" {data.max()/1e3:.2f}] km -> {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grav", default="data/gravity.bin")
    ap.add_argument("--topo", default="data/topo.bin")
    ap.add_argument("--lmax", type=int, default=300)
    ap.add_argument("--ltopo", type=int, default=180)
    args = ap.parse_args()
    import os
    os.makedirs(os.path.dirname(args.grav) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.topo) or ".", exist_ok=True)
    export_gravity(args.grav, args.lmax)
    export_topo(args.topo, args.ltopo)
