#!/usr/bin/env bash
# Physically dated Earth/Sun/Jupiter tides require a SPICE-generated ephemeris CSV.
# Example: MASCONS_EPHEM_FILE=data/ephemeris_2026.csv ./run_kl_100_60d.sh
# Toy alternative: MASCONS_ALLOW_TOY=1 ./run_kl_100_60d.sh
# Resume: MASCONS_RESUME=1 MASCONS_ALLOW_TOY=1 ./run_kl_100_60d.sh
set -euo pipefail
cd "$(dirname "$0")"
: "${MASCONS_OUTDIR:=results/kl300_100_60d}"
: "${MASCONS_DAYS:=60}"
: "${MASCONS_N:=100}"
: "${MASCONS_SEED:=20261014}"
: "${MASCONS_STEPS:=960}"
: "${MASCONS_CP_HOURS:=24}"
: "${MASCONS_EPHEM_FILE:=}"
: "${MASCONS_ALLOW_TOY:=0}"
: "${MASCONS_RESUME:=0}"
: "${OMP_NUM_THREADS:=4}"
export OMP_NUM_THREADS
if [[ -n "$MASCONS_EPHEM_FILE" ]]; then
  EPHEM_ARGS=(--ephem-file "$MASCONS_EPHEM_FILE")
elif [[ "$MASCONS_ALLOW_TOY" == "1" ]]; then
  echo 'WARNING: Running undated synthetic Earth/Sun/Jupiter ephemerides; NOT a scientific/date-specific forecast.' >&2
  EPHEM_ARGS=(--ephem-mode analytic)
else
  echo 'FATAL: supply MASCONS_EPHEM_FILE (SPICE) or explicitly select MASCONS_ALLOW_TOY=1.' >&2
  exit 2
fi
mkdir -p "$MASCONS_OUTDIR"
RESUME_ARGS=()
if [[ "$MASCONS_RESUME" == "1" ]]; then RESUME_ARGS=(--resume); fi
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/lunarmc \
  --n "$MASCONS_N" --seed "$MASCONS_SEED" \
  --alt-min 25 --alt-max 50 --days "$MASCONS_DAYS" --lmax max \
  --steps-per-orbit "$MASCONS_STEPS" --integrator yoshida4 \
  --event-search adaptive --event-max-depth 32 --impact-tol-s 0.001 \
  --kozai-lidov on --third-bodies earth,sun,jupiter \
  "${EPHEM_ARGS[@]}" \
  --traj-count "$MASCONS_N" --traj-dt-min 30 \
  --checkpoint-dir "$MASCONS_OUTDIR/checkpoints" \
  --checkpoint-every-hours "$MASCONS_CP_HOURS" \
  --stop-when-all-impacted \
  --hdf5 "$MASCONS_OUTDIR/orbits.h5" --no-csv \
  --log "$MASCONS_OUTDIR/run.log" \
  "${RESUME_ARGS[@]}"
