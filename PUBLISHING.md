# Publishing this release to GitHub

**Destination:** <https://github.com/metatapios/mascons>

This archive is **ready for review and publication**, but has **not** been pushed, committed or merged to the remote repository by this process. Review the upstream version first, as it may have newer commits or an existing project layout.

## Safe workflow (authenticated workstation)

```bash
git clone https://github.com/metatapios/mascons.git
cd mascons
# Extract the ZIP into a separate staging directory, then copy reviewed files
# into the clone. Do not overwrite or delete the .git directory.
git status
git diff --stat
git diff
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
git add README.md CHANGELOG.md PUBLISHING.md ADAPTIVE_VALIDATION.md   .gitignore CMakeLists.txt DATA_PROVENANCE.txt LunarMasconsOrbit_Guide.pdf   docs/ src/ tests/ tools/ notebooks/ data/ lunar_orbit.yml run_*.sh
git commit -m "Update degree-300 adaptive lunar-orbit code and HPC docs"
git push origin HEAD
```

Use a feature branch/pull request when branch protection or review policy requires it. The commit example assumes you have inspected every diff and have the authority to publish the input data.

## Contents and exclusions

The release includes the corrected C++ solver, launchers, regression tests, bounded lunar-gravity input to degree 300, coarse topography, PDF guide and editable LaTeX source, and Markdown documentation.

Do **not** commit compiled binaries, checkpoints, HDF5 simulation results, scheduler logs, or date-dependent/proprietary SPICE kernels without reviewing their licensing and provenance. These belong in HPC storage, a GitHub Release, or a research data archive as appropriate.

The code currently runs one OpenMP process; use caution before advertising 387-core multinode scaling. Sharded checkpoints and MPI are not yet implemented. The full 25,000-orbit, one-year scientific result has not yet been run and convergence-validated.

**Licensing:** No upstream license has been selected in these files. Add an explicit LICENSE with the copyright holders' approval rather than assuming an open-source license.
