# Earth/Sun/Jupiter 100-particle low lunar orbit experiment (October 2026)

The full 100-particle, 60-day degree-60 ensemble with the undated analytic
Earth/Sun/Jupiter orbit demo is complete (97 impacted, three survivors).
A maximum-degree 300 run was attempted on the same 100 initial conditions;
31 trajectories completed before the compute-limited run was suspended,
with the remaining 69 restartable using the included two-generation HDF5
checkpoints (if transferred from the separately provided results package).

For an *actual dated* third-body experiment, regenerate a Moon-centered
SPICE ephemeris file and change to `--ephem-file`; **do not reuse** the toy
checkpoints when the physical ephemeris is changed.

See `README.md`, `docs/THIRD_BODY_MODEL.md`, and the updated PDF guide.
The 12-page illustrated report and the status/resume instructions accompany
the separate simulation results package.
