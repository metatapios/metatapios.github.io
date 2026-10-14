# HPC production note: 25,000 lunar orbits on 387 CPU cores

**Target**: 25,000 deterministic randomized starts; altitude 25–500 km above interpolated lunar terrain; propagate until first impact or 365.25 simulated days; full available spherical harmonics (degree/order 300); optional dated Earth/Sun/Jupiter forces.

## Deployment status

The present code is **C++17 + OpenMP**, not MPI. Its inner `#pragma omp parallel for schedule(dynamic)` distributes particles among threads **in one shared-memory process**. If 387 cores are spread over several machines/nodes, one process cannot automatically use them all. A Slurm/PBS job-array or MPI-work distributor would need global particle IDs, consistent initialization, controlled writer ownership, and tested restart semantics; none is supplied as a production-ready distributed scheduler.

**Do not** simply launch several invocations with smaller `--n` and expect them to reproduce the 25,000-particle seeded sample: the current master RNG draws initial conditions sequentially, and `--n` and seed conventions are part of reproducibility. A global-ID partitioning feature is a prerequisite for a fully equivalent partitioned run.

## Recommended preflight

1. Obtain and verify a *dated* Moon-centered ephemeris, including reference frame, Moon rotation/libration, and coverage of all possible integration stage times. The synthetic analytic mode is not a substitute.
2. Build on compute nodes with the HPC compiler, OpenMP and thread-safe/compatible HDF5 libraries. Run `ctest`; confirm CMake found HDF5 and OpenMP.
3. Run a one-year **pilot** (suggested: 1,000 starts) to measure impact lifetimes, force-evaluation throughput, checkpoint overhead and memory on the actual nodes. No end-to-end benchmark on 387 cores is included in this release.
4. Compare adaptive tolerance `1e-10` against tightened `1e-11` and `1e-12`, and run selected fixed Yoshida-4 controls. Compare impact/survival classification, event times/sites, long-lived orbit phase and statistical distributions. Degree-300 full-year convergence is **not established** by earlier degree-60/60-day tests.
5. Only after pilot I/O and convergence tests choose the final launch configuration and wall-time request.

## Existing code defaults and caveats

- `run_adaptive_25000_1year.sh`: `N=25000`, `DAYS=365.25`, `LMAX=max`, `ALT_MIN=25`, `ALT_MAX=500`, `INTEGRATOR=adaptive`, `THIRD_BODIES=on`, `TRAJ_COUNT=0`.
- `run_converged_ensemble.sh` selects seed `20261014`, `rtol=1e-10`, `atol-pos-m=0.01`, `atol-vel-mps=1e-5`, event search `adaptive`, `event-max-depth=32`, `impact-tol-s=0.001`.
- `MASCONS_EPHEM_FILE=data/validated_spice.csv` is **required** when third bodies are on unless a synthetic demonstration is explicitly authorized with `MASCONS_ALLOW_TOY=1`.
- Checkpoint interval `MASCONS_CHECKPOINT_SIM_HOURS=12` is in **simulated hours**; it is NOT a 12-hour wall-clock cadence. Two generations are written **for every particle**, potentially 50,000 HDF5 files.
- The checkpoint writer is serialized for ordinary non-thread-safe HDF5; this can dominate HPC runtime with frequent snapshots.
- `TRAJ_COUNT=0` suppresses sampled trajectory arrays to avoid enormous checkpoints. Particle summaries and impacts still go to final `orbits.h5`.
- The final consolidated `orbits.h5` is written after **all particle integrations finish successfully**, not as an incrementally appended live archive. Interrupted runs rely on their checkpoint files for recovery.

## Concise live log

`run.log` is designed for people monitoring an HPC allocation:

```text
2026-10-14 16:10:00 | 135/25000 completed | 0.54% | ETA: 2026-10-15 22:40:00
```

Counts and timestamps shown are illustrative. Poll every 600 **wall-clock** seconds and write on a changed count; otherwise write a line every 1800 seconds. The ETA depends on the rate of new completions and can be misleading when the slower long-lived particles dominate the remaining job. Do not equate successful completion count with accepted integration steps.

## Storage architecture change recommended before full allocation

Implement and test single-writer HDF5 **shards**, two rolling checkpoint generations per shard, and append-only trajectory output separate from restart state. This architecture is **not implemented in the release**. Preserve all completed impact records and never mix checkpoints created with different numerical configurations.

Runtime projections from short benchmark runs are *planning estimates only*; they do not include validated 387-core scaling or the cost of checkpoint contention. Never use them as a guaranteed HPC reservation estimate.
