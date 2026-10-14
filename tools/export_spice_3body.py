#!/usr/bin/env python3
"""Export dated Moon->Earth/Sun/Jupiter ephemerides for lunarmc.

Needs user-supplied NASA NAIF/SPICE kernels (never bundled) and spiceypy:
  python -m pip install spiceypy
  python tools/export_spice_3body.py --kernels naif0012.tls pck00011.tpc de440.bsp \
     --start-utc '2026-10-13T00:00:00 UTC' --days 366 \
     --step-min 60 --out data/third_bodies_2026.csv

JPL positions are transformed from J2000 to a FIXED Cartesian basis equal to
IAU_MOON at the start epoch. The C++ lunar gravity uses a constant spin rate
omega about the z-axis thereafter, so lunar libration/pole variation is NOT
captured. This is an explicit limitation for high-precision/year-long analysis.
"""
import argparse
import csv
import math
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kernels', nargs='+', required=True, help='local NAIF SPK, PCK and LSK kernel files')
    p.add_argument('--start-utc', required=True, help='e.g. 2026-10-13T00:00:00 UTC')
    p.add_argument('--days', type=float, required=True, help='length of simulation [days]')
    p.add_argument('--step-min', type=float, default=60.0, help='sampling cadence [minutes]')
    p.add_argument('--padding-hours', type=float, default=48.0,
                   help='time padding on both sides for Yoshida stage times')
    p.add_argument('--jupiter-target', default='JUPITER BARYCENTER',
                   help='SPICE target name; default planetary system barycenter')
    p.add_argument('--out', required=True)
    args = p.parse_args()
    if not (args.days > 0 and args.step_min > 0 and args.padding_hours >= 0):
        p.error('days/step-min must be positive and padding-hours nonnegative')
    try:
        import spiceypy as spice
    except ImportError as exc:
        raise SystemExit('spiceypy is required: python -m pip install spiceypy') from exc

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        for kernel in args.kernels:
            spice.furnsh(kernel)
        epoch = spice.str2et(args.start_utc)
        # Once and only once: map inertial J2000 axes to the Moon's body frame
        # at the start epoch; keep this transformation FIXED for all samples.
        inertial_to_start_moon = spice.pxform('J2000', 'IAU_MOON', epoch)
        step = 60.0 * args.step_min
        pad = args.padding_hours * 3600.0
        duration = args.days * 86400.0
        first = -pad
        count = int(math.ceil((duration + 2 * pad) / step))
        with out.open('w', newline='') as f:
            f.write('# source=NAIF SPICE geometric states (ABCORR=NONE), observer=MOON\n')
            f.write(f'# initial_UTC={args.start_utc} (et={epoch:.8f} s since J2000 TDB)\n')
            f.write('# frame=J2000 rotated into fixed IAU_MOON axes at initial epoch\n')
            f.write('# body_center=MOON; lengths=m; t=simulation seconds\n')
            f.write('# Jupiter target=' + args.jupiter_target + '\n')
            w = csv.writer(f)
            w.writerow(['time_s','earth_x_m','earth_y_m','earth_z_m',
                        'sun_x_m','sun_y_m','sun_z_m',
                        'jupiter_x_m','jupiter_y_m','jupiter_z_m'])
            targets = ('EARTH','SUN',args.jupiter_target)
            for i in range(count+1):
                t = first + i*step
                values=[f'{t:.8f}']
                for body in targets:
                    position, _ = spice.spkpos(body, epoch+t, 'J2000', 'NONE', 'MOON')
                    # spice uses kilometers; input C++ uses meters.
                    for row in inertial_to_start_moon:
                        values.append(f'{1000.0*sum(row[k]*position[k] for k in range(3)):.16e}')
                w.writerow(values)
        print(f'Exported {count+1} rows to {out}.')
        print('NOTE: the C++ lunar field retains its constant-spin approximation, '
              'not a time-varying JPL lunar orientation.')
    finally:
        spice.kclear()


if __name__ == '__main__':
    main()
