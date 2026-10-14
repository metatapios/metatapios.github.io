#!/usr/bin/env bash
# Regenerate LunarMasconsOrbit_Guide.pdf from its LaTeX source.
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
src="$root/docs/LunarMasconsOrbit_Guide.tex"
target="$root/LunarMasconsOrbit_Guide.pdf"
if ! command -v pdflatex >/dev/null 2>&1; then
  echo 'pdflatex is required to regenerate the technical PDF guide.' >&2
  exit 2
fi
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
for run in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$scratch" "$src" >/dev/null
done
cp "$scratch/LunarMasconsOrbit_Guide.pdf" "$target"
echo "Updated $target"
