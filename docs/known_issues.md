# Known numerical and reproduction issues

## One archived growth-equality exception

The original 7,555-row response has one row that does not satisfy a strict
1e-10 raw growth-residual criterion. At environment
`ENR3000089_24p9357_1p55_0p358134_0p00287857`, growth commitment is 0.1,
the growth target is approximately 4.60592e-7 raw, and the archived solution
is approximately zero. Its residual is approximately -4.60592e-7 raw,
or -0.000460592 h^-1 after the documented growth-unit conversion.
The deviation is approximately 100% of this row's small growth target.

This row is present in the original formal solver sidecar under the same model
and configuration hashes; it was not introduced by release packaging. The
producer's previously declared absolute acceptance threshold was 2e-6 raw,
which the row passes. The actual optimizer feasibility tolerance was not
recorded and is **UNKNOWN**. Passing the archival guard does not establish
strict satisfaction of the mathematical equality.

The row can be selected by an exact single-neighbour query. It was selected
zero times in the archived 2022 cache of 88,191 grid-specific query keys across
123 grids. This is a cache-subset result, not a full-lake or four-year influence
test. See `data/response/source_constraint_audit.json` for the exact record.

The release preserves the response checksum and values. The strict scientific
test remains an explicit expected failure; the archive identity/original-guard
checks and downstream numerical regressions are reported separately. Fixing
the response would require a new upstream source version and renewed downstream
comparison, not silently editing this frozen table.

## Full upstream trajectory scope

The package reproduces all archived downstream windows from frozen ecological
states. Its per-grid ecology regression is not a verification of a fresh
full-lake ecological trajectory. See `reproduction_scope.md` for the interfaces
and modules actually included.
