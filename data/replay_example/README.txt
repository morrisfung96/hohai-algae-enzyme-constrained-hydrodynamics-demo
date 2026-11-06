REAL-STATE SINGLE-STEP REGRESSION FIXTURE

real_source_pool_steps.npz and its JSON metadata contain 36 small cases from
2019--2022: three source formulations and three retained state timestamps per
year. The cases jointly cover 1704 grid-cell steps, including the repaired
ecological-loss to dry-weight fallback.

State values come from retained full-grid ecological exports. Initial pools
come from the retained corrected 30-cycle initialization for each source
formulation. For the independent arithmetic tests, those initial pool values
are paired with selected state timestamps without advancing between them.
These cases test numerical equations; they are not station predictions,
continuous trajectories, or independent field validation.

Expected arrays were produced by executing only the named original function
definitions and frozen defaults. The metadata records source paths and SHA256
identities. The public tests compare all three next pools and all six diagnostic
vectors exactly, then check inventory conservation using distinct real source
and receiving control-volume arrays.

For a continuous replay, use the separate yearly full input archives with
scripts/replay_lake.py. They contain all ecological source columns, consecutive
timestamps, the 4723-cell wet geometry and velocities, normalized station
weights, and the corrected full-grid initial pools. No FACHB905 network or LP
producer is needed by this downstream replay.
