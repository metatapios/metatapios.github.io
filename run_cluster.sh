#!/usr/bin/env bash
# Usage: bash run_cluster.sh [--submit] [path/to/setup.sh]
# --submit submits one Slurm task; otherwise run inside an existing allocation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
SUBMIT=0
if [[ "${1:-}" == --submit ]]; then SUBMIT=1; shift; fi
SETUP="${1:-$ROOT/cluster_setup.sh}"
[[ $# -le 1 ]] || { echo 'Usage: run_cluster.sh [--submit] [setup.sh]' >&2; exit 2; }
SETUP="$(cd "$(dirname "$SETUP")" && pwd)/$(basename "$SETUP")"
set -a
source "$SETUP"
set +a
[[ "$CORES" =~ ^[1-9][0-9]*$ && "$BUILD_JOBS" =~ ^[1-9][0-9]*$ ]] || {
  echo 'CORES and BUILD_JOBS must be positive integers.' >&2; exit 2;
}
if (( SUBMIT )); then
  OPTIONS=(--nodes=1 --ntasks=1 --cpus-per-task="$CORES" --time="$WALLTIME"
           --mem="$MEMORY" --job-name=mascons --chdir="$ROOT"
           --output="$ROOT/mascons-%j.console.log" --export=ALL)
  [[ -z "$PARTITION" ]] || OPTIONS+=(--partition="$PARTITION")
  [[ -z "$ACCOUNT" ]] || OPTIONS+=(--account="$ACCOUNT")
  exec sbatch "${OPTIONS[@]}" "$ROOT/run_cluster.sh" "$SETUP"
fi
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  [[ "${SLURM_NNODES:-1}" == 1 && "${SLURM_NTASKS:-1}" == 1 ]] || {
    echo 'Use exactly one node and one task. This solver uses OpenMP, not MPI.' >&2; exit 2;
  }
  (( CORES <= ${SLURM_CPUS_PER_TASK:-1} )) || {
    echo 'CORES exceeds SLURM_CPUS_PER_TASK.' >&2; exit 2;
  }
fi
mkdir -p "$MASCONS_OUTPUT_DIR"
# Capture configuration, compilation and startup errors even before run.log exists.
printf 'Console log: %s/console.log\n' "$MASCONS_OUTPUT_DIR"
exec >> "$MASCONS_OUTPUT_DIR/console.log" 2>&1
printf 'Started %s; cores=%s; setup=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$CORES" "$SETUP"
cp "$SETUP" "$MASCONS_OUTPUT_DIR/setup-used.sh"
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release \
  -DMASCONS_REQUIRE_PARALLEL_HDF5=ON
COMPILE_JOBS=$BUILD_JOBS
(( COMPILE_JOBS <= CORES )) || COMPILE_JOBS=$CORES
cmake --build "$ROOT/build" --parallel "$COMPILE_JOBS"
export OMP_NUM_THREADS="$CORES" OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=spread
exec bash "$ROOT/run_converged_ensemble.sh" "${EXTRA_ARGS[@]}"
