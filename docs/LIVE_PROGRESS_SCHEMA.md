# Simple log and detailed telemetry

`run.log` is intentionally **only** the concise completion-progress line:

```text
YYYY-MM-DD HH:MM:SS | N/Total completed | XX.XX% | ETA: YYYY-MM-DD HH:MM:SS
```

Before enough completions are observed in the current process, `ETA: unknown`.
The ETA is extrapolated from newly completed trajectories per wall-clock second
since this process started, so it can be biased by early short-lived orbits.
Completed includes impact and successful arrival at the requested cutoff;
numerical failure is excluded. All 25,000 trajectories remain independent,
so a single simulated clock time is not defined.

The monitor *polls* every 600 wall-clock seconds; on any changed completed
count it writes immediately at the next poll, otherwise once every 1800
seconds it writes a heartbeat. An initial fresh-run line and a final line are
also emitted. Resume appends to the existing log rather than truncating it.

CLI: `--progress-every-seconds 600 --progress-heartbeat-seconds 1800`.
Launcher environment: `MASCONS_PROGRESS_EVERY_SECONDS=600` and
`MASCONS_PROGRESS_HEARTBEAT_SECONDS=1800`.

`progress.csv` remains an optional detailed aggregate timeline, updated when
the simple log writes. `particle_status.csv` is a per-particle atomically
replaced diagnostic snapshot updated on the same cadence, and
`particle_events.csv` remains a separate append-only terminal-event ledger.
These are **not** written to `run.log`.

Verbose configuration and physics diagnostics now go to stdout, not `run.log`.
Scientific results remain in the single final `orbits.h5`, while restart
checkpoints keep their existing two-generation-per-particle HDF5 scheme.
