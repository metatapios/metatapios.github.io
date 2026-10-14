#pragma once
// Scheduler for terse run.log updates, independent of orbital/physics code.
// Call at fixed wall-clock polling intervals. Emit on any completion-count
// change, or after one heartbeat with no change. A final status is mandatory.
#include <cmath>
inline bool should_write_progress(int completed, int last_completed,
                                  double seconds_since_last_write,
                                  double heartbeat_seconds, bool final=false) {
    return final || last_completed < 0 || completed != last_completed ||
           seconds_since_last_write >= heartbeat_seconds;
}
