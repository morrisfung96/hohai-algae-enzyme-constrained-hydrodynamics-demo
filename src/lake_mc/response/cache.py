"""Strict import of the archived 2022 R0 grid-specific response context."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .surface import ResponseSurfaceV16

KEY_COLUMNS = ('key_temp_C', 'key_light_Wm2', 'key_DIN_mgN_L',
               'key_DIP_mgP_L', 'key_rho', 'key_k_env')
TEXT_FIELDS = ('matched_env_ids', 'growth_constraint_mode',
               'toxin_solution_status', 'response_mode')
NUMERIC_FIELDS = ('matched_env_distance', 'mu_max_h', 'mu_committed_h',
                  'toxin_capacity_upper_per_gDW_h', 'mu_toxin_solution_h',
                  'growth_constraint_target_raw', 'biomass_flux_toxin_solution_raw',
                  'growth_constraint_residual_raw')


class RecoveredCacheArchive:
    """Load once, then select each original worker/grid context explicitly.

    This retained asset is valid only for the 2022 R0 station-neighborhood
    reproduction audit. It is not an observation mapping or fitted correction.
    """

    def __init__(self, path: str | Path, manifest: str | Path | None = None):
        self.path = Path(path)
        manifest_path = Path(manifest) if manifest else self.path.with_name(
            'recovered_worker_cache_2022_manifest.json')
        self.manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        digest = hashlib.sha256(self.path.read_bytes()).hexdigest()
        if digest != self.manifest['cache_sha256']:
            raise ValueError('Recovered cache SHA256 does not match its manifest')
        self.df = pd.read_csv(self.path, float_precision='round_trip', keep_default_na=False)
        required = {'grid_id', *KEY_COLUMNS, *TEXT_FIELDS, *NUMERIC_FIELDS}
        missing = required.difference(self.df.columns)
        if missing:
            raise ValueError(f'Recovered cache missing columns: {sorted(missing)}')
        if len(self.df) != self.manifest['cache_rows']:
            raise ValueError('Recovered cache row count does not match its manifest')
        self.grids = tuple(sorted(self.df.grid_id.unique()))
        if list(self.grids) != self.manifest['grids']:
            raise ValueError('Recovered cache grid coverage does not match its manifest')
        if self.df.duplicated(['grid_id', *KEY_COLUMNS]).any():
            raise ValueError('Recovered cache contains ambiguous grid/key entries')

    def for_grid(self, surface: ResponseSurfaceV16, grid_id: str) -> 'RecoveredGridResponse':
        """Select a strict grid context without merging it into another cache."""
        if grid_id not in self.grids:
            raise ValueError(f'Grid {grid_id!r} is outside the retained 2022 R0 cache coverage')
        table_sha = hashlib.sha256(surface.path.read_bytes()).hexdigest()
        if table_sha != self.manifest['response_table_sha256']:
            raise ValueError('Response-table lineage does not match the recovered cache')
        if surface.cfg.exact_night_light_threshold_Wm2 != 1e-9:
            raise ValueError('Night gate differs from the retained recovered-cache configuration')
        return RecoveredGridResponse(surface, grid_id, self.df[self.df.grid_id.eq(grid_id)])


class RecoveredGridResponse:
    """Response-like object that fails closed outside one archived grid context."""

    def __init__(self, surface: ResponseSurfaceV16, grid_id: str, rows: pd.DataFrame):
        self.surface = surface
        self.grid_id = grid_id
        self.lookup_cache = {}
        for record in rows.to_dict(orient='records'):
            key = tuple(float(record[c]) for c in KEY_COLUMNS[:-1]) + (int(record[KEY_COLUMNS[-1]]),)
            self.lookup_cache[key] = {
                **{c: str(record[c]) for c in TEXT_FIELDS},
                **{c: float(record[c]) for c in NUMERIC_FIELDS},
            }

    def __getattr__(self, name):
        return getattr(self.surface, name)

    def lookup(self, *, temp_C: float, light_Wm2: float, DIN_mgN_L: float,
               DIP_mgP_L: float, rho_growth_commitment: float,
               k_env: int = 4) -> dict[str, float | str]:
        light = max(float(light_Wm2), 0.0)
        if light <= self.surface.cfg.exact_night_light_threshold_Wm2:
            return self.surface.lookup(temp_C=temp_C, light_Wm2=light,
                DIN_mgN_L=DIN_mgN_L, DIP_mgP_L=DIP_mgP_L,
                rho_growth_commitment=rho_growth_commitment, k_env=k_env)
        rho = float(np.clip(rho_growth_commitment, 0.0, 1.0))
        key = ResponseSurfaceV16._round_key(temp_C, light, DIN_mgN_L, DIP_mgP_L, rho, k_env)
        if key not in self.lookup_cache:
            raise ValueError(f'Query left the retained 2022 R0 context for {self.grid_id}: {key}')
        return self.lookup_cache[key]
