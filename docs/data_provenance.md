# Data provenance

The release uses the manuscript and Supporting Information revised on
2026-11-04. Their source checksums are recorded in `manifests/manuscript_identity.json`.
They are reference identities, not a claim that manuscript files are distributed.

The metabolic response is a frozen computed result. The growth-rate columns are
in h^-1; the MC-LR capacity is in mmol gDW^-1 h^-1. The original biomass reaction
requires a factor of 1,000 to express growth in h^-1. That conversion is not
applied to the toxin-capacity column. The response inventory and query data
dictionary define the supported coordinates and archived solver diagnostics.

The metabolic table is supported at irregular environment/commitment nodes.
Marginal coordinate ranges must not be interpreted as a fully populated grid.
Representative direct-LP pathway diagnostics are computed outputs with a
different role from the runtime response table and the approved network artwork.

The observational evaluation originates from the field observations and frozen
joins described in the article: surface particulate plus dissolved total MC,
LOD/2 censoring, 0.75-m sampling representation, and complete local-day support
for the 164-row MC comparison. Internal process concentrations are MC-LR model
quantities; the mapping targets observed total MC. They are distinct quantities.

Observation sources are attributed in the article to NOAA NCEI/GLERL and the
supporting field data providers; light inputs are derived from NASA POWER.
The repository does not assign a new blanket license to those source datasets.
Publication-specific source citations remain those of the manuscript and SI.

`manifests/asset_sources.json` pairs public file checksums with source checksums.
Metadata path normalization is recorded separately from byte-for-byte copies.
No personal confirmation emails, database exports, private model files, or
bundled third-party source trees are required for the public response interface.
