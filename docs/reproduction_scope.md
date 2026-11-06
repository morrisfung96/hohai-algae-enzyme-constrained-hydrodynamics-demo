# Reproduction scope

The computational modules form a demonstration and frozen-result reproduction
release. The archived input and result identities are retained throughout.

| Layer | Included material | Verification |
| --- | --- | --- |
| Metabolic response | 7,555 computed rows, lookup implementation, units, diagnostics | Original-class numerical query regression; one separately documented strict constraint exception |
| Ecological arithmetic | Adopted production projection, full per-grid metabolic-capacity ecological update, configuration | Retained 16-hour day/night driver regression with zero numerical difference |
| Ecological trajectory | Four complete frozen ecological-state bundles | Source checksums preserved; no new ecological simulation or selection |
| Toxin source, inventories and offline exchange | Portable numerical code, full annual geometry/velocities and retained initial pools | Four years × three variants, complete archived windows, zero numerical difference from retained reference outputs |
| Observation evaluation | Frozen 164-MC, 64-surface-biomass and 181-biomass tables, coefficients and metrics code | Strict common-sample identity and archived metric checks |
| Figures | Four approved main figures, supplementary network artwork, numerical tables | Checksummed artwork assembly; not a new renderer for every panel |

The full ecological lake controller, including the original initialization
analysis and time/transport scheduling, is not repackaged as a new standalone
end-to-end executable here. The extracted per-grid generator exposes the same
transported-state feedback interface as the original arithmetic. The annual
frozen-state inputs support exact downstream replay without changing this scope.

The underlying FACHB905 network and its reconstruction/LP-generation workflow
are not distributed. SI Sections S2.1–S2.8 describe the model formulation; the
published response table is the public metabolic interface.

The archived 2022 recovered response cache covers 123 retained grids. It is
context-specific and is not a substitute for an unarchived full-lake cache.
The strict cache interface stops on unknown non-night keys.

The LEOFS/FVCOM history fields drive the offline exchange implementation.
Historical solver-native horizontal/vertical/open-boundary flux recovery is not
claimed. No new causal water-transport attribution is introduced by this release.
