Frozen ecological arithmetic
============================

state.py retains the original full per-grid R18 ecological generator and its
transitive arithmetic helpers. production.py retains the conservative
capacity-constrained production projection. config.py retains V16Config.
The single default scenario is the adopted OE1 allocation with the formal
R18 dynamic-capacity defaults; candidate ranking and development controllers
are excluded. frozen_parameters.json preserves the source configuration
values; omitted fields use the original dataclass defaults.

The response object only needs lookup(**environment_and_rho) returning the
precomputed response schema. No FACHB905 model, SBML, cobra, optimization
solver, or LP-generation fallback is used by these modules.

The generator yields each pre-update source row and Bm_next/Bo_next. A lake
controller may send back transported Bm/Bo inventories for the next step,
exactly as in the archived source. Consistent control-volume conversion,
wet-domain handling, initial-field analysis, time support and station
aggregation remain responsibilities of that controller.

Regression command, from the repository root:
  python -m unittest discover -s src/lake_mc/ecology -p test_fixture.py

The fixture uses 16 retained environmental hours and a retained initial
profile from one 2020 grid, together with response lookups recorded from the
released 7555-row table under exact R18 query keys. It covers night/day
production, losses, capacity projection and state balance. The original
producer and extracted update matched with zero numerical difference.
This arithmetic regression does not establish a fresh full-lake trajectory
reproduction or recover historical solver-native water fluxes.

See extraction_provenance.json for original source hashes and line segments.
