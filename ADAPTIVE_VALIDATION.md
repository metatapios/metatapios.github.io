# Adaptive controller: short, actual numerical measurements

Date: 2026-10-14. Corrected degree-300 gravity and rotating Moon, no third bodies.
All comparison starts use the same seed 20261014 and constant initial altitude
(25, 150 or 500 km above terrain), one particle evolved for 0.5 day.

| Initial altitude | DOPRI RK5(4) elapsed | Yoshida4 960 elapsed | DOPRI accepted steps | Yoshida4 steps |
|:--|--:|--:|--:|--:|
| 25 km | 4.070 s | 5.306 s | 2,534 | 6,257 |
| 150 km | 1.496 s | 4.346 s | 1,006 | 5,646 |
| 500 km | 0.946 s | 3.448 s | 551 | 4,373 |

**Endpoint accuracy check (12 h):** At 25 km, default DOPRI differed by
0.246 m in position from an independent Yoshida4 3840-step reference;
with tighter controls (`rtol=1e-12`, abs position=0.0001 m,
abs velocity=1e-7 m/s, cap=45 s) this difference was 0.003 m.
At 500 km default DOPRI versus the tightened DOPRI control differed by
0.221 m in position and 0.000148 m/s in velocity.

There was no impact in these three short tests. These results establish that
adaptive time selection and HDF5 restart operate, not one-year scientific
convergence. Actual timestep/impact convergence must be repeated with
third-body forces and dated, validated ephemerides across an orbit ensemble.
The impact detector's bounds remain heuristic near grid interpolation limits.

All 12 CTest checks pass, including corrected gravity, fixed Yoshida4,
DOPRI continuous-extension endpoint and tolerance response, impact regression,
HDF5 output, and split/restart agreement.

## HPC limitations

The current checkpoint implementation produces two HDF5 files per particle;
25,000 trajectories would create 50,000 files. With `--traj-count 0`, this
retains all particle states, diagnostics and impacts, but *not* dense time
series for all orbits. For fully resolved one-year analysis, adopt sharded
append-only trajectory output and sharded checkpoints to control filesystem
metadata and I/O. This is not yet implemented in this release.
