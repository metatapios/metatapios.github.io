## 2026-10-14 — Concise HPC progress log

- `run.log` contains only `YYYY-MM-DD HH:MM:SS | N/Total completed | P% | ETA: timestamp/unknown`.
- Check completion counts every 10 wall-clock minutes; write on change or after 30 idle minutes.
- Successful completion is impact or reaching simulated end time; exclude numerical failures.
- ETA uses only new completions in the current process, avoiding false throughput from restored checkpoints.
- Verbose physics diagnostics are console-only; detailed CSV/HDF5 remain separately available.
- New `--progress-heartbeat-seconds` option; default poll 600 s and heartbeat 1800 s.
- Added deterministic cadence-policy test and log-format restart checks. No orbit dynamics changes.

## 2026-10-14 — live per-particle progress / consolidated final HDF5

- Confirmed that one compressed HDF5 archive is written at successful job completion and contains every particle summary (trajectory arrays only for selected IDs).
- Added append-only run progress and particle event CSV logs and an atomic current-status CSV with one row for every particle, updated while OpenMP workers integrate.
- Logs status, simulated days, integrated particle-days, steps, altitude, cumulative maximum eccentricity, energy and Jacobi changes, numerical invariant drift (lunar-only), adaptive local-error ratio, impact sites and failures.
- Progress wall-clock cadence and per-particle simulated-time reporting interval are independently configurable; restart preserves previous logs and avoids re-emitting already completed terminal events when continuing the same horizon.
- Existing corrected gravity, adaptive DOPRI5(4), Yoshida4 and first-terrain-contact algorithm are unchanged apart from exposing the DOPRI accepted local-error ratio.
- Added `live_progress_and_resume` end-to-end regression test.
- Sharded checkpoints and streaming final HDF5 trajectory output are NOT implemented here.

## 2026-10-14 — degree-300 adaptive (25–500 km) release

- Corrected the planned production altitude distribution from 25–50 to **25–500 km**;
  dedicated launch script targets **25,000 particles, 365.25 days**.
- Added DOPRI5(4) component-scaled absolute/relative error control, FSAL
  endpoint reuse and continuous fourth-degree event interpolation.
- Removed the old broad low-altitude terrain cap from the adaptive controller;
  replaced it with a smooth altitude-dependent maximum timestep and a tighter
  ceiling within 10 km of the surface. Kept independent adaptive event search.
- Kept Yoshida4 and corrected gravity evaluator unchanged for controls.
- Added solver test for dense interpolation endpoints and tolerance response;
  performed short degree-300 benchmark at 25, 150, and 500 km.
- Checkpoint run fingerprint version changed; do not mix old checkpoint files.
- No long-run ensemble convergence or dated SPICE run is claimed.

## 2026-10-14: Adaptive event-search / convergence release

- Retained verified corrected time-aware Yoshida4 and harmonic gravity evaluator.
- Added optional clearance-aware adaptive search (`--event-search adaptive`, now default); earliest-candidate-first recursive subinterval examination with conservative terrain/speed bound and fail-closed unresolved grazing events.
- Kept original fixed-probe event finder with `--event-search sampled` for historical bitwise-comparison tests.
- Added event-max-depth (default 32), new config fingerprint fields and native HDF5 output metadata attributes.
- Corrected initial clearance=0 immediate impact.
- Added `run_converged_ensemble.sh`, explicit SPICE-vs-undated forcing guard, and high-resolution numerical prescription in updated PDF guide.
- All 11 automated tests pass, including new adaptive-vs-fixed-probe event regression.
- Verified 20-particle, 60-day step/event test; full 100-particle convergence comparisons reported separately in validation results.
- Does not claim the 603x603 terrain or 300-degree gravity expansion represent converged mission-level physical data.

## 2026-10-14 — maximum-degree option and KL campaign

- Added `--lmax max`: read the exact highest degree from `gravity.bin` before
  generating any restart signature. The supplied real-data binary has degree 300,
  equivalent to 45,451 fully normalized harmonic coefficient pairs, **not 600**.
- Kept corrected gravity normalization and Yoshida-4 time staging unchanged.
- Reused the existing physically correct Earth/Sun/Jupiter differential third-body
  forces (there is no independent artificial Kozai--Lidov acceleration).
- Documented 100-particle 60-day experiment and SPICE limitations.
- Added explicit maximum-degree and overspecified-degree regression tests.

## 2026-10-14 — Explicit all-impacted termination and 60-day campaign

- Added `--stop-when-all-impacted`: logs the actual termination reason.
  The independent-particle solver already terminates each trajectory at
  impact, so this switch does not change results or computational cost.
- Added an automated regression confirming that opt-in and default trajectories
  coincide and that a three-particle all-impact ensemble terminates early.
- Updated README, LaTeX source and PDF technical guide for this option.
- Added reproducible 100-satellite, 60-day, degree-60 comparison instructions.
- Preserved corrected gravity, Yoshida-4 and adaptive integrator source bytes.

## 2026-10-14 — HDF5 checkpoint V2 (two generations per particle)

- Preserved the corrected rotating non-axisymmetric gravity and the true
  Yoshida-4 / Verlet / DOPRI5(4) integrator implementations unchanged.
- Added native typed, compressed per-particle HDF5 restart checkpoints.
- Each particle retains the current and preceding snapshots, atomically
  rotated after staging and checksum verification.
- Read-time schema, run signature, identity and full-state checksum validation.
- Automatically falls back to the previous valid generation if current is
  corrupt or missing, without resetting impacted particles.
- Snapshots include all restart controls, result diagnostics, permanent impact
  records, and previously recorded trajectory samples if enabled.
- Serialized HDF5 checkpoint calls across OpenMP workers for non-thread-safe
  HDF5 builds; CSV-only use remains possible without HDF5.
- Extended restart tests for both generations, deliberate corruption, process
  interruption, old result comparison, adaptive steps and impact persistence.
- V2 format is not backward compatible with legacy V1 text `.chk` snapshots.

# 2026-10-13 — Rewritten technical PDF guide

- Restored `LunarMasconsOrbit_Guide.pdf` at project root as a 12-page reference manual for the corrected source and optional Earth/Sun/Jupiter module.
- Added maintainable `docs/LunarMasconsOrbit_Guide.tex` and `docs/rebuild_guide.sh`.
- Documented the exact code defaults and output fields, normalized lunar harmonic potential, non-averaged external tides and Jacobi variation, SPICE conventions, solver stages, impact geometry, and required convergence studies.
- Preserved clear caveats: supplied coarse LOLA topography is not a 59 m SLDEM2015 DEM, dated SPICE-frame fidelity has not been benchmarked, and the original Monte Carlo predictions have not been rerun.
- This documentation-only release does not change the C++ physics implementation or scientific input datasets.

---

# 2026-10-13 — Optional KL master switch and updated manuscript

- Added explicit `--kozai-lidov on|off` switch (`off` by default), plus
  `--enable-kozai-lidov` / `--disable-kozai-lidov` aliases.
- A bare `--kozai-lidov on` selects all three perturbers and the undated
  analytic example ephemeris (NOT for scientific runs). `--ephem-file` enables
  dated external forcing with an independently supplied ephemeris.
- External `--third-bodies` and ephemeris settings cannot silently enable
  the tide while the master switch is off; contradictions fail fast.
- Production script now uses `MASCONS_KOZAI_LIDOV=on|off` and requires
  explicit ephemeris settings to activate tides.
- Added `tests/test_kl_switch.py`, including default/off identity, on/off
  force semantics, aliases, input validation and body selection.
- Updated accompanying Elsevier CAS manuscript equations for inertial-frame
  Hamiltonian, indirect term, Jacobi variation, KL quadrupole/diagnostic,
  correct integration and impact-angle equations. Marked archival results
  and unsupported 59-m claims as not yet reproduced by corrected solver.

# Additional third-body dynamics release — 2026-10-13

- Added opt-in Earth, Sun and Jupiter differential gravitational forces,
  including the essential indirect lunar-origin acceleration.
- Implemented explicit approximate analytic orbit mode and a Moon-centered
  ephemeris CSV reader with cubic interpolation, bounds checks and SI units.
- Added optional `tools/export_spice_3body.py` to generate CSV positions from
  locally provided NASA/NAIF SPICE ephemeris, PCK and time kernels.
- All integrators sample the same external forces at the correct time stages.
- Added osculating KL indicators and corrected conservation reporting: the
  lunar-only Jacobi integral ceases to be conserved when external tides act.
- Fixed an additional termination bug in which survivors at an exact stop
  epoch could loop indefinitely on a floating-point remainder.
- Added `third_body_tests` and `lunar_third_body_end_to_end` CTest coverage.
- Extended run scripts and README to expose source-of-ephemeris choices.
- SPICE position *export* was not executed here (spiceypy/kernels unavailable);
  synthetic ephemeris loading, interpolation and complete simulated orbits
  **were** tested. Accurate lunar libration/frame matching remains pending.
- Historical 25,000-orbit statistics remain uncorrected and must be rerun.

---

# Corrective scientific/software release — 2026-10-13

## Verified software fixes

- Fixed 4pi-normalized geodetic associated Legendre sectoral seed:
  `Pbar_11 = sqrt(3)*cos(latitude)`, not `sqrt(3/2)*cos(latitude)`.
  All m>0 gravitational contributions are affected. Added input checks,
  near-pole Cartesian gradient fallback, independent SciPy reference tests.
- Fixed integrator dispatch: `--integrator yoshida4` now invokes the actual
  fourth-order Yoshida composition. Historic `--integrator yoshida` is
  explicitly mapped to `yoshida4`; unknown names cause a clear error.
- Correctly evaluates rotation of the gravity field at each force stage;
  the old method froze the lunar rotation angle during each complete step.
- Replaced the non-terminating adaptive Verlet step-doubling algorithm with
  explicit-failure, error-controlled Dormand–Prince RK5(4). Default rtol 1e-10.
- Reintegrates/brackets forward-time topographic surface impacts. Position,
  velocity and time are now evaluated consistently at the same event time.
- Calculates a local bilinear surface normal and ground-relative impact angle;
  additionally reports `radial_angle_deg` for comparison.
- Implements formerly unused `--dense-factor`/`--dense-envelope-km` near-surface
  time-step controls and explicitly disclaims symplecticity for variable dt.
- Validates CLI options, data headers and grid sizes; propagates thread-local
  orbit failures instead of generating partial apparently valid output.
- Prevents synthetic test-data tool from overwriting actual lunar datasets.
- Removes stale notebook outputs and the outdated PDF guide generator.
- Corrects mislabeled surface image (higher rendering DPI does NOT improve
  physical topographic resolution), and clarifies plotting angle semantics.
- Removes archived production CSVs from the NEW code release; these are
  **numerically superseded** and remain accessible in the original tarball.

## Validation

CMake/CTest: unit + independent SciPy reference + one impact end-to-end
regression, all passing. Yoshida 300/600/1200 steps-per-orbit refinement
shows the expected fourth-order drop in Jacobi drift in a degree-60, 50-km
impacting case. At rtol=1e-10, adaptive RK5(4) agrees on impact timing.
See `validation/validation_report.md` for the actual figures.

## Not implemented (requires external science inputs or expanded modeling)

- No separate 59 m SLDEM2015 tiles/dataset or true high-resolution local
  terrain refinement was present; these cannot be reconstructed from a coarse
  603x603 grid.
- Earth/Sun/Jupiter gravity is now optional and tested as described above;
  no solar radiation pressure or physical models of cratering, ejecta,
  or optical flashes are present.
- No recomputation of the original 25,000-orbit, 366-day science tables.
- No proof of harmonic-degree convergence at N=300 or comprehensive ensemble
  uncertainty quantification. Further work is needed before paper submission.
