REPRESENTATIVE METABOLIC CALCULATION RESULTS

The three CSV files are unchanged retained numerical summaries. manifest.json
gives their source hashes, row counts and column names. No network, solver,
stoichiometric matrix, GPR or enzyme-capacity table is included.

metabolic_period_overview.csv: five representative environmental states,
growth/MC-LR capacity, exchange usage and ratios to Typical-daylight.
Growth rates are h^-1. MC-LR capacity and absolute pathway/exchange fluxes use
mmol gDW^-1 h^-1 on the retained metabolic flux basis. Usage fractions and
relative_to_typical values are dimensionless. Light is W m^-2; DIN and DIP
are mg N L^-1 and mg P L^-1. Timestamps retain their explicit UTC offset.

metabolic_selected_pathways.csv: eight pathway summaries in each of five
states (40 rows). edge_id is a retained pathway-map label, not a network model.

metabolic_severe_light_diagnostic.csv: two representative states and two
explicit diagnostic_constraint_challenge states. These four rows are not an
expanded observational data set or a new response-table interpolation.

These representative calculations used rho=0.9. Formal lake simulations use
state-dependent growth commitment. The Low-N-daylight DIN=0.007305 mg N/L
lies below the formal response-table minimum 0.0532; the direct-solver summary
must not be described as an in-domain lookup result. The retained dark zeros
follow the adopted diagnostic/night conventions and do not establish a new
independent dark-physiology calculation. No new simulations were performed
for this release.
