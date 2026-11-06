"""Explicit one-to-one inventory transfer between named control volumes.

Concentration (or concentration per unit time) is first converted to inventory
on the source volume. The receiver retains that inventory and derives its own
concentration. This defines an inventory interface, not an overlap remap or a
reconstruction of historical hydrodynamic fluxes. Missing source volumes are
errors; callers must never silently substitute receiver volumes.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


class InventoryContractError(ValueError):
    """The requested inventory transfer has no complete volume/identity contract."""


@dataclass(frozen=True)
class InventoryHandoff:
    source_ids: tuple[str, ...]
    receiver_ids: tuple[str, ...]
    source_volume_m3: np.ndarray
    receiver_volume_m3: np.ndarray

    def __post_init__(self):
        src = tuple(map(str, self.source_ids))
        dst = tuple(map(str, self.receiver_ids))
        if not src or len(set(src)) != len(src) or len(set(dst)) != len(dst):
            raise InventoryContractError('Empty or duplicate grid identities')
        if set(src) != set(dst):
            raise InventoryContractError('Source and receiver identities differ; an explicit remap is required')
        object.__setattr__(self, 'source_ids', src)
        object.__setattr__(self, 'receiver_ids', dst)
        for name, ids in [('source_volume_m3', src), ('receiver_volume_m3', dst)]:
            vol = np.array(getattr(self, name), dtype=np.float64, copy=True)
            if vol.shape != (len(ids),):
                raise InventoryContractError(f'{name}: one volume per grid is required')
            bad = ~np.isfinite(vol) | (vol <= 0)
            if bad.any():
                examples = [ids[i] for i in np.flatnonzero(bad)[:8]]
                raise InventoryContractError(f'{name}: {int(bad.sum())} nonpositive or undefined volumes; grids={examples}')
            vol.setflags(write=False)
            object.__setattr__(self, name, vol)
        positions = {gid: i for i, gid in enumerate(src)}
        index = np.array([positions[gid] for gid in dst], dtype=np.int64)
        index.setflags(write=False)
        object.__setattr__(self, '_receiver_order', index)

    def to_inventory(self, source_concentration: np.ndarray) -> np.ndarray:
        """Return receiver-ordered inventory; a concentration rate yields inventory/time."""
        value = np.asarray(source_concentration, dtype=np.float64)
        if value.ndim < 1 or value.shape[-1] != len(self.source_ids):
            raise InventoryContractError('Last dimension must match the source grid order')
        if not np.isfinite(value).all() or (value < 0).any():
            raise InventoryContractError('Source concentration must be finite and nonnegative')
        inventory = value * self.source_volume_m3
        if not np.isfinite(inventory).all():
            raise InventoryContractError('Inventory conversion overflow')
        return inventory[..., self._receiver_order]

    def to_concentration(self, source_concentration: np.ndarray) -> np.ndarray:
        """Return receiver concentration without changing the transferred inventory."""
        return self.to_inventory(source_concentration) / self.receiver_volume_m3
