# LunarMasconsOrbit

A C++17/OpenMP Monte Carlo propagator for lunar satellites in a **rotating, non-axisymmetric lunar gravitational field**, with spherical harmonics, a terrain-intersection model, and optional Earth/Sun/Jupiter third-body accelerations.

This repository provides the **14 October 2026 adaptive degree-300 release**, including native HDF5 results, two-generation restart checkpoints, and a minimal live progress log.

> **Research status:** The corrected gravity evaluator, Yoshida-4 integrator, adaptive DOPRI5(4) solver, impact detection, and checkpoint/restart system have automated tests. A production study of **25,000 randomized orbits, initially 25–500 km above the Moon, evolved for up to 365.25 days or impact**, is *planned*, **not completed or certified converged**. Validate results with timestep, gravitational-degree, topographic, and ephemeris tests before scientific use.

## Documentation

- **[Complete technical guide (PDF)](LunarMasconsOrbit_Guide.pdf)** — equations, numerical methods, command-line options, inputs, validation and limitations. Editable source: [`docs/LunarMasconsOrbit_Guide.tex`](docs/LunarMasconsOrbit_Guide.tex).
- [HPC production deployment: 387 cores](docs/HPC_PRODUCTION.md) — supported parallelism, checkpoints, output and job scheduling caveats.
- [Live progress log and telemetry](docs/LIVE_PROGRESS_SCHEMA.md) — exact log semantics and cadence.
- [Third-body models and ephemerides](docs/THIRD_BODY_MODEL.md).
- [Adaptive integration validation](ADAPTIVE_VALIDATION.md).
- [Historical 60-day comparison](docs/60_DAY_COMPARISON.md) and [KL trial](docs/KL_60DAY_EXPERIMENT.md): demonstrations, **not** validation of the proposed one-year 25,000-orbit study.
- [Changelog](CHANGELOG.md) · [Publishing instructions](PUBLISHING.md) · [Data provenance](DATA_PROVENANCE.txt).

## Included numerical methods

| Component | Implemented behaviour |
| --- | --- |
| Lunar gravity | Non-axisymmetric spherical harmonics through the maximum **degree/order 300** in the bundled `data/gravity.bin`; `--lmax max` reads that limit. |
| Time stepping | Adaptive Dormand–Prince 5(4) with scaled position/velocity local-error control; fixed-step Yoshida-4 remains for independent convergence comparisons. |
| Surface impacts | Adaptive terrain-clearance event search with fail-closed handling of unresolved near-grazing cases; impact terminates the particle. |
| Third bodies | Optional differential Earth/Sun/Jupiter acceleration; physical runs need a validated dated ephemeris CSV. The `--kozai-lidov` switch enables the third-body model; it is **not** a separate phenomenological KL force. |
| Final output | **One consolidated HDF5 analysis file** written only on successful run completion. |
| Checkpoints | **Two HDF5 snapshots per particle** (current and previous), with restart integrity checks. Sharded checkpoints have **not** been implemented. |
| Live progress | A four-field `run.log` line whenever the 10-minute poll sees new completions; otherwise a 30-minute heartbeat. |

The bundled terrain model is a coarse LOLA-derived interpolation grid; degree-300 gravity does not imply kilometre-scale surface realism or mission-grade impact-location accuracy.

## Build and test

Requirements: C++17 compiler, CMake ≥ 3.15 and an HDF5 C development library to enable HDF5 output and restarts. OpenMP is recommended. Python, `h5py` and SciPy are used by selected tests and analysis scripts.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

CMake can build a CSV-only binary without HDF5, but the launch scripts below require the native HDF5-enabled build. OpenMP is optional at compile time, but without it the propagator runs serially.

## Simple cluster launcher

Edit [`cluster_setup.sh`](cluster_setup.sh) to choose cores, particles, duration, Slurm resources and numerical settings. Submit with `bash run_cluster.sh --submit`, or use `bash run_cluster.sh` inside an existing single-node allocation. Build and startup errors go to `console.log`; live progress goes to `run.log` in the selected output directory. See [cluster quick start](docs/CLUSTER_QUICKSTART.md).

## Smoke test (no dated ephemeris required)

```bash
OMP_NUM_THREADS=4 \
MASCONS_N=8 MASCONS_DAYS=1 \
MASCONS_THIRD_BODIES=off \
MASCONS_OUTPUT_DIR=results/smoke \
./run_converged_ensemble.sh
```

The launcher defaults to degree 300 (`MASCONS_LMAX=max`), altitude 25–500 km, adaptive integration and seed `20261014`. It writes `results/smoke/orbits.h5` after successful completion, plus rolling restart files and the simple progress log.

## Planned 25,000-particle / one-year experiment

The launcher below provides the *intended experiment configuration*, **not a claim that the current parallel I/O architecture is production-scaled**:

```bash
OMP_NUM_THREADS=8 \
MASCONS_EPHEM_FILE=data/validated_spice.csv \
MASCONS_OUTPUT_DIR=results/one_year \
./run_adaptive_25000_1year.sh
```

You must first supply a **validated Moon-centered ephemeris** with matching time/frame conventions; the release does not ship one. For a synthetic *demonstration only*, explicitly use `MASCONS_ALLOW_TOY=1` instead. Do not publish toy-ephemeris results as dated lunar predictions.

For an existing run, use the **same physics and numerical settings** and set `MASCONS_RESUME=1`; the checkpoint signature rejects incompatible changes. Use separate output/checkpoint directories when comparing tolerances or integrators. See the PDF guide for the full interface.

Adaptive defaults are `rtol=1e-10`, position absolute tolerance `0.01 m`, velocity absolute tolerance `1e-5 m/s`, maximum step `90 s`, minimum step `0.02 s`, and altitude-dependent maximum-step protection. **These are initial numerical settings, not demonstrated year-long error bounds.** A degree-300 convergence campaign is required before scaling to 25,000 orbits.

## Simple live log

`results/<run>/run.log` contains only:

```text
YYYY-MM-DD HH:MM:SS | N/25000 completed | XX.XX% | ETA: YYYY-MM-DD HH:MM:SS
```

“Completed” means **impacted** or **reached the requested simulation horizon**; failed trajectories do not count. Every **600 wall-clock seconds** the monitor checks the completion total, writes on change, and otherwise writes once every **1800 seconds**. The ETA is provisional and can be `unknown`. For example:

```bash
tail -f results/one_year/run.log
```

Additional structured telemetry is saved separately as `progress.csv`, `particle_status.csv` and `particle_events.csv`. Those files do not alter the concise `run.log` format.

## Data, reproducibility and storage

`data/gravity.bin` and `data/topo.bin` are included for running the code; the data provenance note is at [`DATA_PROVENANCE.txt`](DATA_PROVENANCE.txt). Initialization uses a fixed master seed, uniform initial altitude in the specified range, isotropic initial position and randomized in-plane heading; the initial orbit speed is circular in the lunar monopole approximation.

The current executable is **OpenMP-only; it has no MPI distribution layer**. A 387-core allocation spanning multiple nodes is not equivalent to `OMP_NUM_THREADS=387`. Splitting a run into independent jobs also requires careful global particle-ID and RNG handling to preserve a consistent 25,000-particle ensemble. See [HPC production deployment](docs/HPC_PRODUCTION.md).

The current final HDF5 file is emitted **at the end**, not incrementally. Restart files are **per particle**, with two generations; at 25,000 particles that can mean up to **50,000 checkpoint files**. For a large multi-node campaign, use a tested partition/orchestration workflow or first implement sharded HDF5 output/checkpointing. Do not push runtime files, checkpoints, or private ephemeris kernels into the source repository.

## Project status and license

This archive is prepared for [`metatapios/mascons`](https://github.com/metatapios/mascons), but **no remote commit or push has been made from this package**. The upstream repository may contain changes not present here: review the diff before publishing. The archive does **not** establish a software license; choose and include a `LICENSE` file before inviting unrestricted reuse.
