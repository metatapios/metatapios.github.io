#!/usr/bin/env bash
# Single-orbit run + all plots (panels, elements, maps).
# Usage: ./run_single.sh [ALT_KM] [DAYS] [LMAX] [SEED] [EXTRA lunarmc args...]
set -e
cd "$(dirname "$0")"

ALT=${1:-100}
DAYS=${2:-90}
LMAX=${3:-120}
SEED=${4:-777}
shift 4 2>/dev/null || true

cmake --build build -j >/dev/null 2>&1 || {
  cmake -B build -S . -DCMAKE_BUILD_TYPE=Release >/dev/null; cmake --build build -j >/dev/null; }

if [ ! -f data/gravity.bin ]; then
  echo "== export gravity + topography =="; python tools/export_data.py --lmax 300 --ltopo 180
fi

echo "== single run: alt=$ALT km days=$DAYS lmax=$LMAX seed=$SEED =="
mkdir -p results
./build/lunarmc --n 1 --alt-min "$ALT" --alt-max "$ALT" --days "$DAYS" \
  --lmax "$LMAX" --seed "$SEED" --traj-count 1 --traj-dt-min 1 "$@"

echo "== plots =="
python tools/plot_single.py 0
python tools/make_map_widget.py 0
echo "== done: results/single_run_*.png, results/moon_maps.html =="
