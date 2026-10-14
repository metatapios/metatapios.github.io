# Release notes — 2026-10-14

This package is based on `mascons_degree300_adaptive_simple_progress_20261014.zip`; the **numerical C++ source, tests, scripts, and PDF guide are unchanged** from that release. Markdown documentation has been reorganized for GitHub.

- Degree-300 lunar gravity + 25–500 km initial altitude; adaptive RK5(4), fixed Yoshida-4 control.
- Third-body Earth/Sun/Jupiter forcing requires a scientifically validated dated ephemeris for physical results.
- One final HDF5 output at completion; dual per-particle HDF5 restart checkpoints, no distributed sharding yet.
- `run.log` checks completion every ten wall-clock minutes, writes when the count changes, and otherwise records a heartbeat after thirty minutes.
- Planned large run: 25,000 randomized starts, 365.25 simulated days, up to 387 allocated CPU cores. **Not yet carried out or validated at production scale.**

For detailed implementation changes and prior demonstrations, see [`../CHANGELOG.md`](../CHANGELOG.md), [`../ADAPTIVE_VALIDATION.md`](../ADAPTIVE_VALIDATION.md), and the technical PDF guide.
