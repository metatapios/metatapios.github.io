# Edit this file before launching. Values can also be overridden in the environment.
# Load your cluster's compiler, CMake, OpenMP and HDF5 modules before launching.
CORES=${CORES:-8}                         # OpenMP threads on ONE node
BUILD_JOBS=${BUILD_JOBS:-4}               # Concurrent compiler processes
PARTITION=${PARTITION:-}                 # Slurm partition; empty uses cluster default
ACCOUNT=${ACCOUNT:-}                     # Slurm account; empty uses cluster default
WALLTIME=${WALLTIME:-24:00:00}
MEMORY=${MEMORY:-16G}                    # Slurm memory per node
MASCONS_N=${MASCONS_N:-20}
MASCONS_DAYS=${MASCONS_DAYS:-1}
MASCONS_LMAX=${MASCONS_LMAX:-max}
MASCONS_ALT_MIN=${MASCONS_ALT_MIN:-25}    # km
MASCONS_ALT_MAX=${MASCONS_ALT_MAX:-500}   # km
MASCONS_SEED=${MASCONS_SEED:-20261014}
MASCONS_OUTPUT_DIR=${MASCONS_OUTPUT_DIR:-results/cluster_run}
MASCONS_CHECKPOINT_DIR=${MASCONS_CHECKPOINT_DIR:-$MASCONS_OUTPUT_DIR/checkpoints}
MASCONS_RESUME=${MASCONS_RESUME:-0}       # 1 resumes the same output/checkpoints
MASCONS_CHECKPOINT_SIM_HOURS=${MASCONS_CHECKPOINT_SIM_HOURS:-12}
MASCONS_INTEGRATOR=${MASCONS_INTEGRATOR:-adaptive}
MASCONS_STEPS=${MASCONS_STEPS:-480}       # Fixed-step steps per orbit
MASCONS_RTOL=${MASCONS_RTOL:-1e-10}
MASCONS_ATOL_POS_M=${MASCONS_ATOL_POS_M:-0.01}
MASCONS_ATOL_VEL_MPS=${MASCONS_ATOL_VEL_MPS:-1e-5}
MASCONS_ADAPTIVE_MIN_DT_S=${MASCONS_ADAPTIVE_MIN_DT_S:-0.02}
MASCONS_ADAPTIVE_MAX_DT_S=${MASCONS_ADAPTIVE_MAX_DT_S:-90}
MASCONS_ADAPTIVE_LOW_STEPS=${MASCONS_ADAPTIVE_LOW_STEPS:-480}
MASCONS_ADAPTIVE_HIGH_STEPS=${MASCONS_ADAPTIVE_HIGH_STEPS:-120}
MASCONS_TRAJ_COUNT=${MASCONS_TRAJ_COUNT:-0}
MASCONS_TRAJ_DT_MIN=${MASCONS_TRAJ_DT_MIN:-30}
MASCONS_THIRD_BODIES=${MASCONS_THIRD_BODIES:-off}
MASCONS_BODIES=${MASCONS_BODIES:-earth,sun,jupiter}
MASCONS_EPHEM_FILE=${MASCONS_EPHEM_FILE:-} # Validated Moon-centered ephemeris CSV
MASCONS_ALLOW_TOY=${MASCONS_ALLOW_TOY:-0} # Demonstration only
MASCONS_PROGRESS_EVERY_SECONDS=${MASCONS_PROGRESS_EVERY_SECONDS:-10}
MASCONS_PROGRESS_HEARTBEAT_SECONDS=${MASCONS_PROGRESS_HEARTBEAT_SECONDS:-30}
MASCONS_PROGRESS_SIM_HOURS=${MASCONS_PROGRESS_SIM_HOURS:-6}
# Any remaining lunarmc options. These override earlier options when repeated.
# All supported options: ./build/lunarmc --help
EXTRA_ARGS=(--event-search adaptive --event-max-depth 32 --impact-tol-s 0.001
            --event-substeps 4 --adaptive-scale-height-km 75
            --dense-factor 1 --dense-envelope-km 11 --omega 2.661699e-6)
