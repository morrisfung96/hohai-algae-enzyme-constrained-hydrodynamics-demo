# R29 hydrodynamic lineage boundary

## Frozen project identity

The formal research line remains the Delft3D/lake-scale environmental
ecological--metabolic--hydrodynamic coupled simulation.  The conservative ecological core, all
ecological parameters, the decline gate, and the held-out 2019 validation year
remain frozen.

R29 is water-only.  It may audit hydrodynamic input lineage and water-volume
continuity.  It may not perform biomass attribution until conservative water
continuity has passed.

## Three distinct identities

1. **Formal coupled-model lineage** -- the Delft3D/lake-scale research model
   and its ecological/metabolic coupling.
2. **Archived hydrodynamic input identity** -- the retained 2020, 2021 and 2022
   LEOFS NetCDF state files.  Their own global metadata consistently record
   `title=LEOFS` and `source=FVCOM_4.3`.
3. **Public reference source** -- the downloaded public FVCOM 4.3 source used
   during R28 to interpret possible meanings of `omega`, `WTS`, edge fluxes and
   boundary variables.  It is not proven to be the exact historical source
   tree, executable, compiler build, namelist, forcing set or runtime used to
   generate the archived files.

Identity 2 does not imply identity 3.  Consequently the project must not claim
that the metabolic-capacity ecological model itself is an FVCOM 4.3 model, or that compiling the public FVCOM 4.3
tree reproduces the historical hydrodynamics.

## Status of work already completed

- R29C's closed-domain test validates only generic sign conventions,
  double-precision accumulation and cancellation bookkeeping.
- R29D demonstrates only that part of a public reference implementation can be
  compiled locally in serial mode.
- Neither result is a historical run, a Delft3D run, a LEOFS reproduction, a
  water-continuity acceptance result, or evidence supporting biomass
  attribution.
- The public-source compile extension is stopped after the front-module
  preflight; it must not be advanced as the formal R29 path.

## Active R29 decision

The retained six-hour state snapshots do not contain verified, time-integrated
solver-native horizontal face, vertical-interface, open-boundary and
source/sink volume transports.  R28B already established that those transports
cannot be uniquely reconstructed with acceptable continuity from the retained
snapshots.

Therefore the active status is:

`R29_HISTORICAL_NATIVE_WATER_CONTINUITY = INPUT_LINEAGE_BLOCKED`

and:

`R28C_BIOMASS_ATTRIBUTION = BLOCKED`

No ecological parameter tuning, decline-gate candidate, biomass attribution,
or 2019 validation run is authorized by the reference-source work.

## Publication wording

If the R28/R29 diagnostic is mentioned, use wording equivalent to:

> The archived LEOFS hydrodynamic state files identify FVCOM 4.3 as their data
> source. Public FVCOM v4.3 source code was consulted only as a reference for
> interpreting archived variable semantics; it was not treated as the verified
> historical executable or as the lineage of the formal coupled model.

