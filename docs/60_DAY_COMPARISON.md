# Controlled 60-day comparison: unchanged September code vs corrected 14 October solver

**Status:** Both actual numerical integrations completed successfully on 14 October 2026. These are independent direct 60-day integrations, not extensions from the old 30-day analysis CSV files (which cannot serve as integrator restart state). A separate 60.00 -> 60.01 day checkpoint-resume smoke test passed using the corrected solver.

## Identical experiment

- 100 low lunar orbits from identical pseudorandom seed `20261014`, with initial local-terrain altitude uniform on [25,50] km, independent random surface positions and tangential directions.
- Lunar gravity through spherical-harmonic degree 60 with the original supplied `gravity.bin`; identical `topo.bin` at 603 x 603 cells; constant spin `2.661699e-6 rad/s`.
- Nominally 240 steps per initial monopole orbit, near-surface dense factor 1; up to 60 simulated days or collision, whichever comes first.
- No third-body forces. Kozai--Lidov off. Trajectories sampled approximately every 30 minutes for all 100 particles.
- **Original September code unedited**: the original executable's option `--integrator yoshida4` incorrectly runs its Verlet kernel. The original gravity normalization also differs from the corrected code. Both results must be labeled by code revision, not interpreted as convergence of two identical force models.
- **Corrected October code**: actual Yoshida4, corrected spherical harmonics and integration-stage timing, HDF5 results, two-generation HDF5 checkpoints for each particle, interval 6 simulated hours. Flag `--stop-when-all-impacted` enabled; three particles survived the 60-day horizon.

## Measured results

| Metric | Original September | Corrected October |
|---|---:|---:|
| Impacts by 30 days | 58 | 78 |
| Impacts by 60 days | 83 | **97** |
| Still orbiting at 60 days | 17 | **3** |
| Median time among impacting satellites [day] | 19.626 | 12.756 |
| First impact [day] | 2.141 | 1.596 |
| Latest impact [day] | 58.276 | 57.948 |
| Median maximum osculating eccentricity | 0.01968 | 0.02004 |
| Median maximum relative Jacobi drift | 9.61e-6 | 3.78e-7 |
| Maximum relative Jacobi drift | 1.38e-5 | 6.52e-7 |

Every one of the original version's 83 impacts also impacted in the corrected simulation by day 60. Another 14 particles impacted only in the corrected simulation; no particles impacted exclusively in the older run. For the 83 matching impact IDs, the **median separation of body-fixed collision sites** was approximately **20.55 degrees** and the **median absolute impact-time difference** was **1.17 days**.

**Cautions:** (1) the two versions use different effective gravity/integration and cannot isolate which correction causes the difference; (2) this is one deterministic N=100 ensemble, not a statistically representative sampling of all lunar orbital designs; (3) harmonic-degree and terrain-grid convergence have not been established; (4) only the lunar-only Jacobi invariant is expected to be conserved in these tests; (5) apsidal phase angles become numerically ill-conditioned as eccentricity approaches zero; (6) 240 steps per orbit is not a guaranteed mission-accuracy configuration.

## Reproduce the runs

Run the original binary from the *unchanged* 2026-09-10 source directory, where the original `data` files are present:

```sh
OMP_NUM_THREADS=2 ./build/lunarmc --n 100 --lmax 60 --days 60 \
  --steps-per-orbit 240 --alt-min 25 --alt-max 50 --omega 2.661699e-6 \
  --integrator yoshida4 --dense-factor 1 --traj-count 100 \
  --traj-dt-min 30 --seed 20261014 \
  --out old_60d_summary.csv --traj old_60d_trajectories.csv \
  --log old_60d_run.log
```

Run the **corrected** executable from the root of this package:

```sh
OMP_NUM_THREADS=4 ./build/lunarmc --n 100 --lmax 60 --days 60 \
  --steps-per-orbit 240 --alt-min 25 --alt-max 50 --omega 2.661699e-6 \
  --integrator yoshida4 --dense-factor 1 --traj-count 100 \
  --traj-dt-min 30 --seed 20261014 --kozai-lidov off \
  --stop-when-all-impacted --checkpoint-dir results/checkpoints \
  --checkpoint-every-hours 6 --hdf5 results/new_60d.h5 \
  --no-csv --log results/new_60d_run.log
```

To extend **the corrected 60-day checkpoint set** to, e.g., 90 days, preserve all options and inputs, change `--days 60` to `--days 90`, add `--resume` and set a new `--hdf5`/`--log` path. Never reuse old 30-day CSVs as checkpoints; they lack exact integrator state.

## Generate figures

```sh
python tools/compare_60d.py \
    --old-summary old_60d_summary.csv \
    --old-trajectories old_60d_trajectories.csv \
    --new-hdf5 results/new_60d.h5 \
    --output-dir results/figures_60d
```

This produces a paired-particle outcome table and seven figures: survival, eccentricity, altitude, running maximum Jacobi drift, nodal precession, projected eccentricity orientation, and lunar impact locations. Plot coordinates are lunar body-fixed east-positive longitude; osculating angles are computed from the inertial Cartesian sampled states.

## GitHub delivery

The source package is prepared for publication at `https://github.com/metatapios/mascons`. In the execution environment used to generate this release, attempts to connect to `github.com` failed at DNS resolution. Consequently, **no push was performed**. The ZIP is the local release artifact, not proof of a remote commit. Copy its contents into a local clone, inspect `git diff`, commit, and push from an authenticated machine. See `PUBLISHING.md` for commands.
