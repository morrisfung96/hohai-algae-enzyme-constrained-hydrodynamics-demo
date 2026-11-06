ARCHIVED METABOLIC RESPONSE INPUT

metabolic_response.csv.gz is the unchanged formal 7555-row, 16-column table.
SHA256: c68f04181b4825ce6dae1781d5c7b88e426d2945f1566fdd50468ede1942178d
data_dictionary.json records units, support, source hashes and field meanings.
The public calculation chain begins with this table. It does not reconstruct
or solve the nonpublic FACHB905 source network.

Query interface
  from lake_mc.response import ResponseConfig, PrecomputedMetabolicResponse
  response = PrecomputedMetabolicResponse(
      "data/response/metabolic_response.csv.gz",
      ResponseConfig(require_native_v16_response=True))
  result = response.lookup(temp_C=24.5, light_Wm2=100,
      DIN_mgN_L=0.3, DIP_mgP_L=0.003, rho_growth_commitment=0.7)

The archived implementation first canonicalizes the table and aggregates exact
coordinate duplicates using maximum numeric capacities. It sorts the five
coordinates through the pandas groupby operation. It standardizes coordinates
by population standard deviations (ddof=0), with zero scales replaced by 1 and
the rho scale bounded below by 0.10. The default query selects four neighbors
using scipy.spatial.cKDTree. If selected neighbors contain non-night rows,
selected exact-night rows are discarded. Outputs use normalized 1/(d+1e-9)
weights. The original NumPy fallback is retained; its tie selection/order may
differ from cKDTree. Both paths have distinct archived regression cases.

Only light is clipped below at zero and rho to [0,1]. Temperature, DIN and DIP
are not clipped to the table's marginal ranges. Outside-support queries still
use the same nearest-neighbor rule; this is neither a new direct LP solution
nor evidence that the archived response covers a full Cartesian domain.
Light <= 1e-9 W m^-2 triggers a strict zero response before querying or caching.
This is the adopted lake-query night gate, not an independent dark-growth LP.

Capacity is mmol MC-LR gDW^-1 h^-1. The toxin column is NOT multiplied by 1000.
Raw biomass reaction flux times 1000 gives the archived growth rate in h^-1.
Every source row is labelled as an equality-constrained solve. Numerical audit
finds one retained exception with a nonzero raw equality residual; see
source_constraint_audit.json. The source bytes are preserved for historical
reproduction, and this known source-data issue is not silently corrected.

Historical cache identity
Rounded keys use temperature/light to 3 decimals, DIN to 5, DIP to 6, rho to 3,
then k_env. The first unrounded query populating a rounded key is returned by
later queries sharing it. The cache holds up to 300000 keys without eviction.
Changing query order can therefore change results. Cold-cache results must
not be labelled exact historical replay merely because the table is unchanged.

The 2022 R0 recovery is archived as recovered_worker_cache_2022.csv.gz, with
its retained recovery-query records, original forward audit and a manifest.
Use RecoveredCacheArchive(...).for_grid(surface, grid_id) for a strict context.
It is only for the 123 listed official station-neighborhood grids in the R0
reproduction audit. Missing non-night keys raise an error. Do not apply it to
full-lake runs, other years, R1 or counterfactuals, or merge worker contexts.
The original forward state/process audit passed at <1e-10 on all 123 grids.
The public cache export repeats only response interpolation; it does not
repeat the upstream LP or the ecology integration.

query_regression_cases.json contains expectations from the isolated original
query class: 63 cases per backend, including support points, irregular-support
interpolation, night gates, clipping and queries outside marginal support.
