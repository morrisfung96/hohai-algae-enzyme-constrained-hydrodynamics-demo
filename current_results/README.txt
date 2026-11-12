Current archived result assets

This directory accompanies the revised English submission. It provides a 67,452-row precomputed environment/commitment response, the current 164 station-date evaluation, selected computed capacity diagnostics and source data. It does not contain the original FACHB905 metabolic network, full enzyme/proteomics tables, model-construction notebooks, upstream LP/response generators or raw workbooks.

Run python scripts/verify_current_results.py to verify this manifest and rescore the saved current predictions. Run python scripts/query_current_response.py --help for queries of the saved response table. Neither command refits observations or regenerates ecological states.

The current dynamic OLS MAE is 0.3580401138314635 on log10(C+0.01); its 164-row identity is separate from supplementary water-column configurations and the legacy release. The original 7,555-row table and four-year frozen replay remain archived under their existing paths and v0.1.0 input archives. Their replay verification does not certify a complete rerun of this current configuration.

Initial main-configuration toxin pools are declared cycle-30 states from 72-hour forcing replay. None of 12 runs met the strict convergence criterion; these are not equilibrium states. See observation_initialization_summary.json and manuscript/current_methods.txt.

Numerical CSV tokens were preserved. Machine-local path cells were replaced by opaque retained-reference identifiers; Chinese comparison labels were given English aliases where specified. Original package/member SHA256 values preserve provenance without distributing local path lists or private model payloads.

The named partial helper is scripts/numerical_capacity_interface.py. Check its 13 frozen query witnesses with:
python scripts/verify_capacity_query.py current_results/data/response/environment_response.csv.gz
The helper exposes saved-response queries and selected ecological target equations, not a full ecological controller. Its equation specification is manuscript/numerical_capacity_interface_specification.txt; parameter values are current_results/numerical_capacity_interface_parameters.json.
