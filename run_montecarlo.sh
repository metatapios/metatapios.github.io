#!/usr/bin/env bash
# Usage: ./run_montecarlo.sh [N] [DAYS] [LMAX] [ALT_MIN] [ALT_MAX] [BIN_KM]
set -euo pipefail
cd "$(dirname "$0")"
: "${OMP_NUM_THREADS:=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
export OMP_NUM_THREADS
N=${1:-25000}; DAYS=${2:-366}; LMAX=${3:-300}
ALT_MIN=${4:-25}; ALT_MAX=${5:-250}; BIN=${6:-10}
# The third-body force switch is OFF by default for reproducible lunar-only runs.
# Scientific example: MASCONS_KOZAI_LIDOV=on MASCONS_EPHEM_FILE=data/my_spice.csv ./run_montecarlo.sh
# Toy ephemeris: MASCONS_KOZAI_LIDOV=on MASCONS_EPHEM_MODE=analytic ./run_montecarlo.sh
KOZAI_LIDOV=${MASCONS_KOZAI_LIDOV:-off}
THIRD_BODIES=${MASCONS_THIRD_BODIES:-all}
EPHEM_FILE=${MASCONS_EPHEM_FILE:-}
EPHEM_MODE=${MASCONS_EPHEM_MODE:-}
EXTRA=("${@:7}")
THIRD_ARGS=(--kozai-lidov "$KOZAI_LIDOV")
if [[ "$KOZAI_LIDOV" == on ]]; then
  THIRD_ARGS+=(--third-bodies "$THIRD_BODIES")
  if [[ -n "$EPHEM_FILE" ]]; then THIRD_ARGS+=(--ephem-file "$EPHEM_FILE")
  elif [[ "$EPHEM_MODE" == analytic ]]; then THIRD_ARGS+=(--ephem-mode analytic)
  else
    echo 'Production external tides require MASCONS_EPHEM_FILE with dated NAIF SPICE positions.' >&2
    echo 'For demonstration only, explicitly set MASCONS_EPHEM_MODE=analytic.' >&2
    exit 2
  fi
elif [[ "$KOZAI_LIDOV" != off ]]; then
  echo 'MASCONS_KOZAI_LIDOV must be on or off' >&2; exit 2
elif [[ -n "$EPHEM_FILE" || -n "$EPHEM_MODE" || "${MASCONS_THIRD_BODIES:-none}" != none ]]; then
  echo 'External force settings supplied but MASCONS_KOZAI_LIDOV=off' >&2; exit 2
fi
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
if [[ ! -f data/gravity.bin || ! -f data/topo.bin ]]; then
  echo 'Real data not present; export using tools/export_data.py (requires pyshtools).' >&2
  exit 2
fi
mkdir -p results
./build/lunarmc --n "$N" --days "$DAYS" --lmax "$LMAX" \
  --alt-min "$ALT_MIN" --alt-max "$ALT_MAX" \
  --integrator yoshida4 --traj-count 12 --steps-per-orbit 1000 \
  --dense-factor 1 --out results/summary.csv "${THIRD_ARGS[@]}" "${EXTRA[@]}"
python tools/montecarlo_stats.py results/summary.csv "$BIN"
python tools/make_global_highres_maps.py
echo 'DONE: New outputs in results/; the archived 2026 results are NOT valid for this corrected solver.'
