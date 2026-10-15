# Cluster quick start

This launcher runs one shared-memory OpenMP process on one node. It does not distribute work across nodes.

## 1. Prepare the software

Clone this repository on the cluster. Load the site-provided C++17 compiler, CMake, OpenMP and HDF5 development modules. Module names differ by cluster; ask the administrators for their names. Load these modules before submission so the job inherits their environment.

## 2. Edit cluster_setup.sh

Set `CORES` to the number of allocated CPU cores on one node. Set `PARTITION`, `ACCOUNT`, `WALLTIME` and `MEMORY` for Slurm. Set particle count, duration, output directory and physics options in the same file. Start with the default 20-particle pilot.

Every listed setting can be overridden with an environment variable. `EXTRA_ARGS` holds remaining executable options; run `build/lunarmc --help` for the complete interface. Configuration files are shell scripts and must be trusted.

## 3. Submit on Slurm

From the repository:

```bash
bash run_cluster.sh --submit
```

For example, override the setup without editing it:

```bash
CORES=64 MASCONS_N=25000 MASCONS_DAYS=365.25 \
MASCONS_OUTPUT_DIR=results/year_run bash run_cluster.sh --submit
```

This requests exactly one node, one task and 64 CPUs per task. Choose a core count available on one node. Do not use `mpirun` or start multiple copies against the same output directory.

For another scheduler, request one node and the desired cores using that scheduler, then run inside the allocation:

```bash
CORES=32 bash run_cluster.sh
```

The launcher configures and builds the program, requires OpenMP and HDF5, sets thread placement, then starts propagation. Compilation uses at most `min(BUILD_JOBS, CORES)` processes. Do not launch production computation on a login node.

## 4. Monitor

With the default output directory:

```bash
tail -f results/cluster_run/console.log
tail -f results/cluster_run/run.log
```

`console.log` captures compilation and startup errors. `run.log` appears once initialization succeeds. The default poll is 10 seconds; a changed completion count creates a line, otherwise a heartbeat appears about every 30 seconds. These are wall-clock intervals. Slurm also writes `mascons-JOBID.console.log` in the repository.

`command.txt` records the resolved executable command. `setup-used.sh` copies the setup file; environment overrides are reflected in `command.txt`. Results are written to `orbits.h5` after successful completion. Checkpoints are written per particle, with two generations. With 25,000 particles this can create 50,000 files and serialize HDF5 I/O; validate storage and scaling on a pilot first.

## 5. Resume

After the previous job has stopped, use the same output directory and numerical/physics settings:

```bash
MASCONS_RESUME=1 bash run_cluster.sh --submit
```

Increase `MASCONS_DAYS` to extend the horizon. Never run two jobs against the same output/checkpoint directories. Use a different output directory for independent experiments.

Third bodies default to off. For physical external tides, set `MASCONS_THIRD_BODIES=on` and `MASCONS_EPHEM_FILE` to a validated Moon-centered ephemeris CSV. `MASCONS_ALLOW_TOY=1` permits demonstration ephemerides only.

Validation on macOS: 20 particles, degree 300, four threads, 0.01 days completed with HDF5, checkpoints and progress logs. Slurm submission arguments were checked with a stub; actual scheduling requires testing on your cluster.
