# Hohai algae enzyme-constrained hydrodynamics demo

This repository accompanies **Linking Metabolic Capacity and Multiphase
Inventories to Algal Biomass and Microcystin Distributions in Lake Erie**.
It connects archived environment-dependent metabolic capacities to lake-scale
ecology, offline hydrodynamic exchange, toxin inventories, and surface observations.

This repository provides reusable ecological arithmetic, a four-year
**frozen-ecological-state replay**, and paper-result reproduction. The frozen
replay has been compared over all four years and three source variants against
79,296 station-hour records, 3,336 station days, and 9,900 mass-ledger steps;
numeric differences were zero. It does not regenerate the full-lake ecological
trajectory from environmental forcing. See `docs/reproduction_scope.md`.

The public computational starting point is the **precomputed metabolic response
table**. The underlying FACHB905 network, its reconstruction notebooks, complete
enzyme parameter tables, and upstream LP generation programs are not included.
Key metabolic formulation details belong to the article's Supporting Information.

## Installation

Use Python 3.11 or newer. From the repository directory:

```sh
python -m venv .venv
# Activate the virtual environment for your operating system.
python -m pip install -r requirements-verified.txt
python -m pip install --no-deps -e .
```

## Reproduction entry points

```sh
python scripts/verify_assets.py
python scripts/validate_response.py
python scripts/reproduce_metrics.py
python scripts/reproduce_figures.py
python -m unittest discover -s tests -v
python -m unittest discover -s src/lake_mc/ecology -p test_fixture.py -v
```

The first commands verify archived inputs, validate the response interface,
rescore frozen predictions, and assemble approved manuscript artwork. Rescoring
does not refit coefficients or rerun lake simulations. The figure command preserves
the four approved main figures; its manifest distinguishes artwork assembly from
new numerical plotting. See `docs/reproduction_scope.md` for the lake-replay status.

The response source contains **one documented numerical constraint exception**.
Archive identity, original acceptance tolerances, and query reproduction pass;
the separate strict growth-equality test remains an expected failure. This is
reported in `docs/known_issues.md` and is not a certification of strict LP
feasibility for every archived row. The response values have not been modified.

## Four-year frozen-state replay

Download a yearly `lake_inputs_YYYY.zip` attachment from the corresponding release,
then use:

```sh
python scripts/unpack_inputs.py lake_inputs_2022.zip
python scripts/replay_lake.py --inputs external_data/2022 --output results/replay_2022 --variant dynamic
```

Available source variants are `dynamic`, `fixed`, and `proportional`. Omit
`--max-hours` for the complete archived window. The yearly archives contain the
complete numerical inputs and reference outputs; example fixtures are additional
unit-regression assets and do not replace the full archives.

## Reader-facing model names

The biomass comparison uses **Metabolic capacity** and **Empirical capacity**
in the article and Figure 3. Historical identifiers remain in original data
columns, configuration keys, and provenance records. The depth-retention surface
expression is a separate observation setting. See `docs/model_names.md`.

Release v0.1.1 updates these names and Figure 3 labels. Its numerical code,
frozen data, configurations, and four yearly input archives are unchanged from
v0.1.0. Download the same annual archives using `manifests/external_inputs.json`.

## Data and model identity

* `data/response/`: the adopted 7,555-row response and its data dictionary.
* `data/metabolic_diagnostics/`: computed representative-state pathway summaries.
* `data/frozen_evaluation/`: the adopted prediction and observation tables.
* `data/paper_diagnostics/`: frozen stage, phase, and mapping diagnostics.
* `data/grid_contracts/`: wet-domain geometry and station mappings.
* `configs/`: frozen ecological configuration.
* `artwork/approved/`: approved main and supplementary artwork.
* `manifests/`: file checksums, input identities, and reproduction status.

Each validation retains its own sample identity: total-MC evaluation uses 164
station-date rows, the 2022 surface-biomass comparison uses 64 rows, and the
four-year biomass table uses 181 rows. They are not pooled into one score.
Legacy internal names in archived columns identify source lineage; reader-facing
outputs describe their scientific meaning.

The formal depth–retention observation expression changes sampling-layer
representation, not water-column biomass or toxin mass. The archived ordinary
LOYO point predictions, quantile mapping, and two-stage risk analysis remain
separate outputs.

## Hydrodynamic forcing

The lake model is driven by LEOFS/FVCOM fields through the documented offline
exchange implementation. This repository does not rebuild the historical
operational LEOFS system. Native historical face-volume flux recovery and causal
transport attribution are separate provenance questions, documented in
`docs/hydrodynamic_lineage.md`.

## Citation and reuse

Please cite the associated article and identify the release/commit used.
The article DOI and author metadata will be linked when finalized; no DOI has
been minted by preparing this directory. See `NOTICE` for reuse and third-party
rights. No project-wide open-source license has yet been selected by the authors.
