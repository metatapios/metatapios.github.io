#!/usr/bin/env bash
# Reproducible experiment launcher. A production run is scientifically converged
# only after comparing different timestep, gravity-degree and terrain settings.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
N="${MASCONS_N:-100}"
DAYS="${MASCONS_DAYS:-60}"
LMAX="${MASCONS_LMAX:-max}"
STEPS="${MASCONS_STEPS:-480}"
SEED="${MASCONS_SEED:-20261014}"
THIRD="${MASCONS_THIRD_BODIES:-off}"
METHOD="${MASCONS_INTEGRATOR:-adaptive}"
TRAJ_COUNT="${MASCONS_TRAJ_COUNT:-0}"
OUT="${MASCONS_OUTPUT_DIR:-results/${METHOD}_lmax${LMAX}_seed${SEED}}"
CHECKPOINTS="${MASCONS_CHECKPOINT_DIR:-$OUT/checkpoints}"
FLAGS=(--n "$N" --days "$DAYS" --seed "$SEED" --lmax "$LMAX"
  --alt-min "${MASCONS_ALT_MIN:-25}" --alt-max "${MASCONS_ALT_MAX:-500}"
  --integrator "$METHOD" --steps-per-orbit "$STEPS"
  --event-search adaptive --event-max-depth 32 --impact-tol-s 0.001
  --traj-count "$TRAJ_COUNT" --traj-dt-min "${MASCONS_TRAJ_DT_MIN:-30}" --checkpoint-dir "$CHECKPOINTS"
  --checkpoint-every-hours "${MASCONS_CHECKPOINT_SIM_HOURS:-12}" --hdf5 "$OUT/orbits.h5" --no-csv
  --log "$OUT/run.log" --progress-dir "$OUT"
  --progress-every-seconds "${MASCONS_PROGRESS_EVERY_SECONDS:-600}"
  --progress-heartbeat-seconds "${MASCONS_PROGRESS_HEARTBEAT_SECONDS:-1800}"
  --progress-sim-hours "${MASCONS_PROGRESS_SIM_HOURS:-6}"
  --stop-when-all-impacted)
if [[ "$METHOD" == adaptive ]]; then
  FLAGS+=(--rtol "${MASCONS_RTOL:-1e-10}"
    --atol-pos-m "${MASCONS_ATOL_POS_M:-0.01}"
    --atol-vel-mps "${MASCONS_ATOL_VEL_MPS:-1e-5}"
    --adaptive-max-dt-s "${MASCONS_ADAPTIVE_MAX_DT_S:-90}"
    --adaptive-min-dt-s "${MASCONS_ADAPTIVE_MIN_DT_S:-0.02}"
    --adaptive-low-steps "${MASCONS_ADAPTIVE_LOW_STEPS:-480}"
    --adaptive-high-steps "${MASCONS_ADAPTIVE_HIGH_STEPS:-120}")
fi
if [[ "$THIRD" == on ]]; then
  FLAGS+=(--kozai-lidov on --third-bodies "${MASCONS_BODIES:-earth,sun,jupiter}")
  if [[ -n "${MASCONS_EPHEM_FILE:-}" ]]; then
    FLAGS+=(--ephem-file "$MASCONS_EPHEM_FILE")
  elif [[ "${MASCONS_ALLOW_TOY:-0}" == 1 ]]; then
    FLAGS+=(--ephem-mode analytic)
    echo "WARNING: arbitrary-phase analytic third-body ephemerides are only a demonstration." >&2
  else
    echo "External tides requested: set MASCONS_EPHEM_FILE or MASCONS_ALLOW_TOY=1 for demo." >&2
    exit 2
  fi
else
  FLAGS+=(--kozai-lidov off)
fi
if [[ "${MASCONS_RESUME:-0}" == 1 ]]; then FLAGS+=(--resume);fi
mkdir -p "$OUT"
# Preserve the resolved exact command in a reproducibility text file.
printf '%q ' "$ROOT/build/lunarmc" "${FLAGS[@]}" > "$OUT/command.txt"
printf '\n' >> "$OUT/command.txt"
exec "$ROOT/build/lunarmc" "${FLAGS[@]}"
