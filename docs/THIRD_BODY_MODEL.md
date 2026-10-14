# Earth, Sun and Jupiter tides / Kozai–Lidov dynamics

**Release: 2026-10-13.** This is an extension of the corrected lunar-mascon
C++ solver. The archived 25,000-particle statistics were **not** regenerated.

## Master switch: opt-in external perturbations

The controlling CLI option is **`--kozai-lidov on|off`** (default: `off`).
With `off` all external accelerations are zero and KL diagnostics are `nan`.
With `on`, all three perturbers are enabled unless `--third-bodies` explicitly
chooses a nonempty subset. Enabling the switch without an ephemeris file
selects the undated **analytic demo ephemeris** and emits a warning; this is
not suitable for physical results. To use dated SPICE positions supply
`--kozai-lidov on --ephem-file FILE`.

To prevent mistakes, external body/ephemeris options are rejected when the
switch is `off`. Equivalent aliases are `--enable-kozai-lidov` and
`--disable-kozai-lidov`.

This is a **switch for external differential gravity**, not a request to add
an independently averaged KL term; the latter would double-count the same
quadrupole contribution. Actual KL oscillations are not guaranteed.

## Dynamics

The satellite position **r** and all perturbing-body positions **R_b(t)** are
expressed in the **same Moon-centered, nonrotating Cartesian axes**, in meters.
The total inertial acceleration is

    r'' = R_z(omega*t) grad U_Moon[R_z(-omega*t) r]
          + sum_b GM_b [(R_b-r)/|R_b-r|^3 - R_b/|R_b|^3],
          b in {Earth, Sun, Jupiter}.

The last term is the **indirect** acceleration accounting for acceleration of
the Moon's origin by the same body. Simply adding GM_b(R_b-r)/|R_b-r|^3
would be incorrect by many orders of magnitude, particularly for the Sun.
The expansion for |r| << |R_b| is

    a_tide,b = GM_b / R_b^3 [3 (n_b · r) n_b - r] + O((r/R_b)^2),
    n_b = R_b/R_b.

This is the familiar quadrupole perturbation underlying Kozai–Lidov dynamics.
We **do not** add an independent Kozai–Lidov acceleration, or use a secular
quadrupole average, because doing so on top of full three-body gravity would
double-count the effect. The code resolves both the short-period variations and
whatever genuine secular KL evolution is present in the chosen configuration.
The lunar mascon precession and collisions can suppress classical KL cycles;
**a Kozai oscillation is not guaranteed** over a 366-day orbit.

The code implements source-term forces in all three integrators:
`verlet`, `yoshida4`, and error-controlled `adaptive` RK5(4).

## Run examples

Build/test:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Approximate toy ephemerides: 

```sh
./build/lunarmc --n 20 --lmax 60 --days 30 \
  --kozai-lidov on --third-bodies earth,sun,jupiter --ephem-mode analytic \
  --integrator yoshida4 --out results/toy_three_body.csv
```

**This toy model is for debugging only.** Earth is placed in an eccentric
Keplerian orbit about the Moon with an inclined orbital plane; the Earth-Sun
and Sun-Jupiter vectors come from independent approximate eccentric Kepler
orbits with **arbitrary phases**. The planets' correct actual positions on any
date are *not* claimed.

### Dated JPL/SPICE positions (recommended)

Provide appropriate NAIF leapseconds, Moon-orientation PCK and planetary SPK
kernels (not distributed in this archive):

```sh
python -m pip install spiceypy
python tools/export_spice_3body.py \
  --kernels /path/to/naif0012.tls /path/to/pck00011.tpc /path/to/de440.bsp \
  --start-utc '2026-10-13T00:00:00 UTC' --days 366 \
  --step-min 60 --padding-hours 48 --out data/third_bodies_spice.csv

./build/lunarmc --n 10 --lmax 60 --days 30 \
  --kozai-lidov on --third-bodies earth,sun,jupiter \
  --ephem-file data/third_bodies_spice.csv \
  --integrator yoshida4 --out results/spice_three_body.csv
```

Or invoke the production script with `MASCONS_KOZAI_LIDOV=on`,
`MASCONS_THIRD_BODIES=earth,sun,jupiter`, and
`MASCONS_EPHEM_FILE=data/third_bodies_spice.csv` in the environment.
For this workflow **the user must supply and validate the SPICE kernels**.
`spiceypy` and the SPICE kernels were unavailable in the release environment,
so export from *real* SPICE kernels was not tested here. The C++ table reader
and third-order interpolation were tested using synthetic ephemerides.

The exported file contains 10 comma-separated columns:

    time_s,earth_x_m,earth_y_m,earth_z_m,sun_x_m,sun_y_m,sun_z_m,jupiter_x_m,jupiter_y_m,jupiter_z_m

`time_s=0` is the chosen start epoch. All coordinates are **geometric**
Moon-to-body positions (SPICE light-time/aberration flag `NONE`, not apparent
positions) expressed in the fixed inertial Cartesian basis obtained by
applying the `J2000 -> IAU_MOON` transformation **at the initial epoch only**.
SPICE reports km; the export converts to m. Jupiter's planetary barycenter
is used unless `--jupiter-target` is changed. The Moon-centered ephemerides
are evaluated by 4-point Lagrange interpolation; all requested integration
stage times must lie inside the file's coverage. 48 h padding avoids
out-of-coverage negative Yoshida internal stages at run boundaries.

**Important precision caveat:** The lunar spherical-harmonic field and
collision grid still rotate rigidly about the fixed z-axis with constant
`--omega`, rather than following the high-precision lunar libration and pole
history from JPL. At a given epoch the inertial axes should also be checked
against the exact GRGM1200B/LOLA coefficient reference frame (mean Earth/polar
axis versus principal axes). The supplied code does **not** correct these
convention differences. Consequently even with JPL third-body ephemerides,
year-long mission-grade/peer-review quantitative claims require an additional
lunar-orientation/frame validation and comparison with an established orbit
propagator (e.g. SPICE/GMAT/Orekit).

## Output interpretation

The lunar-only Jacobi expression

    J_M = |v|^2/2 - U_Moon(r,t) - omega (r x v)_z

is conserved only when the lunar gravity field rotates uniformly and there
are **no third-body forces**. The third-body terms introduce explicit
additional time dependence, hence:

- `invariant_drift_rel` is `nan` with nonzero third bodies: a *conservation*
  diagnostic would be physically wrong.
- `lunar_jacobi_variation_rel` tracks the normalized maximum variation in
  this reference quantity; it includes **physical** perturbations and numerical
  error, and is never an accuracy estimator when third bodies are enabled.
- `kl_incl0_deg`, `kl_incl_end_deg` are osculating orbit inclinations relative
  to the **instantaneous outer orbital normal** of the first enabled body in
  order Earth, Sun, Jupiter.
- `kl_jz0`, `kl_jz_end`, and trajectory CSV `kl_jz` provide the classical
  normalized test-particle diagnostic `sqrt(1-e^2) cos(i_relative)` for bound
  osculating orbits (otherwise `nan`). This quantity is **not conserved**
  with mascons, multiple bodies, or moving orbital planes; secular averaging
  is needed before interpreting it as a KL constant.
- `invariant` in trajectory output is also only the **lunar-only reference**
  Jacobi quantity when third bodies are on.

A useful hierarchy at low lunar altitude is Earth ≫ Sun ≫ Jupiter. A simple
quadrupole-order estimate at radius ~1800 km gives external tidal forces of
order `1e-5`, `1e-7`, `1e-12` m/s² respectively (body distances vary). The
strong lunar gravity harmonics and collision constraints usually dominate
short-term spacecraft behavior. Please compare `--kozai-lidov off` with `--kozai-lidov on --third-bodies
earth`, `earth,sun`, and `all` at identical initial conditions.

## Verification

- Exact direct-minus-indirect differential force against independent formula.
- Vanishing perturbation at the Moon's center; quadrupole far-field limit.
- Three independent body contributions and selection switches.
- Synthetic polynomial table interpolation, time bounds and invalid data.
- End-to-end patched orbit compared to the same unperturbed orbit, including
  nonconservation semantics, trajectory fields, and end-time termination.
- All original tests (gravity normalization/gradient, numerical integrations,
  surface collision, independent SciPy) continue to run unchanged.

All verification is **code-level**, not a high-precision NAIF/mission benchmark,
and not a recomputation of the published Monte Carlo statistics.

Reference manuals:
- NAIF SPICE `SPKEZR`/`SPKPOS`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html
- NAIF frame conventions: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
- JPL planetary satellite GM: https://ssd.jpl.nasa.gov/sats/phys_par/
