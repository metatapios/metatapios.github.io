#!/usr/bin/env bash
# Production-scale launcher; explicitly invoke only after convergence validation.
# This starts 25,000 independent particles, each evolving to impact or day 365.25.
# WARNING: current per-particle HDF5 checkpoint approach creates 50,000 snapshot
# files and requires a suitable parallel filesystem; shard across jobs for HPC.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MASCONS_N="${MASCONS_N:-25000}"
export MASCONS_DAYS="${MASCONS_DAYS:-365.25}"
export MASCONS_LMAX="${MASCONS_LMAX:-max}"
export MASCONS_ALT_MIN="${MASCONS_ALT_MIN:-25}"
export MASCONS_ALT_MAX="${MASCONS_ALT_MAX:-500}"
export MASCONS_INTEGRATOR="${MASCONS_INTEGRATOR:-adaptive}"
export MASCONS_TRAJ_COUNT="${MASCONS_TRAJ_COUNT:-0}"
export MASCONS_THIRD_BODIES="${MASCONS_THIRD_BODIES:-on}"
export MASCONS_OUTPUT_DIR="${MASCONS_OUTPUT_DIR:-results/adaptive_25000_1year_l300}"
exec "$ROOT/run_converged_ensemble.sh"
